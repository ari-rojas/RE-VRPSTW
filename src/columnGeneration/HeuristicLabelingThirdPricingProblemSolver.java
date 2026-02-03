package columnGeneration;

import model.EVRPTW;
import model.EVRPTW.Arc;
import model.EVRPTW.Vertex;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Set;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;

import branchAndPrice.ChargingTimeInequality;
import branchAndPrice.FixArc;
import branchAndPrice.RemoveArc;


/**
 * This class provides a heuristic solver for the ng-SPPRC pricing problem
 * It uses a relaxed dominance rule (but an exact dominance rule)
 */
public final class HeuristicLabelingThirdPricingProblemSolver extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem> {

	public Vertex[] vertices = dataModel.vertices; 			//vertices of the instance
	public PriorityQueue<Vertex> nodesToProcess; 			//labels that need be processed
	public final int numCols = 400; 						//maximum number of routes (columns) allowed
	public int[] infeasibleArcs; 						//arcs that cannot be used by branching
	public final int similarityThreshold = 5; 				//diversification of columns


	/**
	 * Labeling algorithm to solve the ng-SPPRC
	 */
	public HeuristicLabelingThirdPricingProblemSolver(EVRPTW dataModel, PricingProblem pricingProblem) {
		super(dataModel, pricingProblem);
		this.name="ExactLabelingSolver"; //Set a name for the solver
		this.infeasibleArcs = new int[dataModel.numArcs];
		this.nodesToProcess = new PriorityQueue<Vertex>(dataModel.V, new SortVertices());
	}

	/** Runs the labeling algorithm. */
	public void runLabeling() {
		//initialization
		int[] remain_energy = new int[dataModel.gamma + 1]; Arrays.fill( remain_energy, dataModel.E);
		Label initialLabel = new Label(dataModel.C+1, -1, 0, -pricingProblem.dualCost, dataModel.Q, vertices[dataModel.C+1].closing_tw, remain_energy, 0, new boolean[dataModel.C], new boolean[dataModel.C], new boolean[pricingProblem.subsetRowCuts.size()], new HashSet<Integer>(pricingProblem.subsetRowCuts.size()));
		this.nodesToProcess.add(vertices[dataModel.C+1]);
		initialLabel.index = 0;
		vertices[dataModel.C+1].unprocessedLabels.add(initialLabel);
		dataModel.infeasibleArcs = this.infeasibleArcs;

		//Labeling algorithm 
		
		while (!nodesToProcess.isEmpty() && System.currentTimeMillis()<timeLimit) {
			ArrayList<Label> labelsToProcessNext = labelsToProcessNext();
			for(Label currentLabel: labelsToProcessNext) {
				
				boolean isDominated;
				if (currentLabel.vertex == 0) isDominated = checkDominanceAtDepot(currentLabel);
				else isDominated = checkDominance(currentLabel);
				if(isDominated) continue;
				
				else if (currentLabel.vertex > 0) { // The non-dominated label is marked as processed and extended along its incoming arcs.
					currentLabel.index = vertices[currentLabel.vertex].processedLabels.size();
					vertices[currentLabel.vertex].processedLabels.add(currentLabel);
				
					for(Arc a: dataModel.graph.incomingEdgesOf(currentLabel.vertex)) {
						if(a.tail > 0 && a.head <= dataModel.C && a.modifiedCost >= -dataModel.precision) continue;
						if(infeasibleArcs[a.id] > 0) continue;

						Label extendedLabel = extendLabel(currentLabel, a);
						if (extendedLabel!=null) { //verifies if the extension is feasible
							updateNodesToProcess(extendedLabel);
						}
					}
				} else { // The currentVertex is 0, the label corresponds to a complete route.
					
					List<Label> labels_to_remove = new ArrayList<>();
					for (Label processedLabel: vertices[0].processedLabels){ if (isDominatedAtDepot(processedLabel, currentLabel)) labels_to_remove.add(processedLabel); }
					vertices[0].processedLabels.removeAll(labels_to_remove);

					vertices[0].processedLabels.add(currentLabel);

				}
			}
		}
	}


	/**
	 * Selects a set of labels to process (the one with the most remaining load)
	 */
	public ArrayList<Label> labelsToProcessNext(){

		ArrayList<Label> labelsToProcessNext = new ArrayList<Label>();
		Vertex currentVertex = nodesToProcess.poll();
		while(true) {
			Label currentLabel = currentVertex.unprocessedLabels.poll();
			if(labelsToProcessNext.isEmpty()) labelsToProcessNext.add(currentLabel);
			else {
				boolean isDominated = false;
				for(Label L2: labelsToProcessNext) {
					if (currentVertex.id == 0) isDominated = isDominatedAtDepot(currentLabel, L2);
					else isDominated = isDominated(currentLabel, L2);
					if(isDominated) break;
				}
				if(!isDominated) labelsToProcessNext.add(currentLabel);
			}
			if(currentVertex.unprocessedLabels.isEmpty() || (currentVertex.id<=dataModel.C && currentVertex.unprocessedLabels.peek().remainingLoad<currentLabel.remainingLoad)) break;
		}

		if(!currentVertex.unprocessedLabels.isEmpty()) nodesToProcess.add(currentVertex);
		return labelsToProcessNext;
	}


	/** Given a new (non-dominated) label, updates the nodes to be processed. */
	public void updateNodesToProcess(Label extendedLabel) {
		
		Vertex currentVertex = vertices[extendedLabel.vertex];
		currentVertex.unprocessedLabels.add(extendedLabel);
		if (currentVertex.unprocessedLabels.size() == 1) nodesToProcess.add(currentVertex);
	}

	/**
	 * Label extension procedure
	 */
	public Label extendLabel(Label currentLabel, Arc arc) {

		int source = arc.tail;
		if (source>=1 && source<=dataModel.C)
			if (currentLabel.unreachable[source-1]) return null;

		double reducedCost = currentLabel.reducedCost+arc.modifiedCost;

		boolean[] eta = currentLabel.eta.clone();
		HashSet<Integer> srcIndices = new HashSet<Integer>(currentLabel.srcIndices);
		for(int srcIndex: vertices[source].SRCIndices) {
			if(currentLabel.eta[srcIndex]) {
				eta[srcIndex] = false;
				int dualIndex = dataModel.C+dataModel.last_charging_period+srcIndex;
				reducedCost-=pricingProblem.dualCosts[dualIndex];
				srcIndices.remove(srcIndex);
			}
			else {eta[srcIndex]=true; srcIndices.add(srcIndex);}
		}
		reducedCost = Math.floor(reducedCost*10000)/10000;

		//only negative reduced cost labels
		if (source==0 && reducedCost>= pricingProblem.reducedCostThreshold-dataModel.precision) return null;

		int remainingLoad = currentLabel.remainingLoad-vertices[source].load;
		int remainingTime = currentLabel.remainingTime-arc.time;
		if(remainingTime>vertices[source].closing_tw) remainingTime = vertices[source].closing_tw;

		int[] remainingEnergy = new int[dataModel.gamma + 1];
		remainingEnergy[0] = currentLabel.remainingEnergy[0]-arc.energy; if (remainingEnergy[0] < 0) return null;
		for (int gam = 1; gam <= dataModel.gamma; gam++){
			if (currentLabel.remainingEnergy[gam-1] - arc.energy_deviation < currentLabel.remainingEnergy[gam]){ remainingEnergy[gam] = currentLabel.remainingEnergy[gam-1] - arc.energy - arc.energy_deviation; }
			else { remainingEnergy[gam] = currentLabel.remainingEnergy[gam] - arc.energy; }
			if (remainingEnergy[gam] < 0) return null;
		}
		
		int chargingTime = dataModel.f_inverse[dataModel.E-remainingEnergy[dataModel.gamma]];

		//Quick check
		if(source>0 && remainingTime-dataModel.graph.getEdge(0, source).time<vertices[0].opening_tw) return null;

		//Check whether the extension is actually feasible
		if(remainingTime<vertices[source].opening_tw || chargingTime>= (int) (remainingTime/10)) return null;

		boolean[] unreachable = Arrays.copyOf(currentLabel.unreachable.clone(), currentLabel.unreachable.length);
		if(source>0) unreachable[source-1] = true;

		//Mark unreachable customers and ng-path cycling restrictions
		if(source>0) {
			
			for (Arc c: dataModel.graph.incomingEdgesOf(source)) {
				if(c.tail==0 || unreachable[c.tail-1]) continue;
				//unreachable
				if (remainingLoad-vertices[c.tail].load<0 || remainingTime-c.time<vertices[c.tail].opening_tw || remainingEnergy[dataModel.gamma]-c.min_energy < 0 || 
					Math.min(remainingTime-c.time, vertices[c.tail].closing_tw)-dataModel.graph.getEdge(0, c.tail).time<vertices[0].opening_tw ||
					remainingEnergy[dataModel.gamma]-c.min_energy - dataModel.graph.getEdge(0, c.tail).min_energy<0) {
					unreachable[c.tail-1] = true;
				}

			}
		}
		Label extendedLabel = new Label(source, arc.id, currentLabel.index, reducedCost, remainingLoad, remainingTime, remainingEnergy, chargingTime,unreachable, currentLabel.ng_path, eta, srcIndices);
		return extendedLabel;

	}


	/**
	 * When the CG procedure terminates, the close function is invoked. 
	 */
	@Override
	public void close() {
		if(this.pricingProblemInfeasible) {
			for (int i = 0; i < vertices.length; i++) {
				vertices[i].processedLabels = new ArrayList<Label>(dataModel.numArcs);
				vertices[i].unprocessedLabels =  new PriorityQueue<Label>(dataModel.numArcs, new Label.SortLabels());
			}
		}else {
			for (int i = 0; i < vertices.length; i++) {
				vertices[i].processedLabels = new ArrayList<Label>(dataModel.numArcs);
				vertices[i].unprocessedLabels =  new PriorityQueue<Label>(dataModel.numArcs, new Label.SortLabels());
				vertices[i].SRCIndices = new ArrayList<>();
			}
		}
		this.nodesToProcess = new PriorityQueue<Vertex>(new SortVertices());
	}

	/**
	 * This method produces zero or more columns. 
	 */
	@Override
	protected List<Route> generateNewColumns() {

		//Solve the problem and check the solution

		long startTime = System.currentTimeMillis();

		this.runLabeling(); 									//runs the labeling algorithm
		ArrayList<Label> filtered_labels = pricingProblem.charging_pricing_filtering(vertices[0].processedLabels);
		ArrayList<Route> newRoutes = this.charging_pricing(filtered_labels);

		long totalTime = System.currentTimeMillis()-startTime;
		dataModel.heuristicPricingTime+=totalTime;
		if (dataModel.print_log) logger.debug("Time solving (heuristically) the pricing problem (s): " + getTimeInSeconds(totalTime)); 

		if (newRoutes.isEmpty()) { pricingProblemInfeasible=true; this.objective=Double.MAX_VALUE; }
		else {pricingProblemInfeasible=false;}

		if (dataModel.print_log) {
				logger.debug("Finished heuristic pricing: "+vertices[0].processedLabels.size()+" processed, "+vertices[0].unprocessedLabels.size()+" unprocessed.");
				logger.debug("Found " + newRoutes.size() + " columns");
		}

		close(); //restart
		return newRoutes;
	}

	public ArrayList<Route> charging_pricing(ArrayList<Label> labels){

		ArrayList<Route> newRoutes = new ArrayList<>();

		for (Label label: labels){

			int departureTime = (int) (label.remainingTime/10);
			int chargingTime = label.chargingTime;
			double reducedCost = label.reducedCost;

			if (chargingTime<departureTime){ // Only performs pricing for the routes that will have at least one column with reduced cost

				int load = dataModel.Q - label.remainingLoad;
				int energy = dataModel.E-label.remainingEnergy[dataModel.gamma];
				
				// Retrieve route information
				boolean isElementary = true;
				HashMap<Integer, Integer> route=new HashMap<Integer, Integer>(dataModel.C);
				ArrayList<Integer> arcs = new ArrayList<Integer>(dataModel.C);
				
				int currentVertex = 0; Label currentLabel = label.clone();
				int cost = 0;

				while(currentVertex!=dataModel.C+1) {
					Arc currentArc = dataModel.arcs[currentLabel.nextArc];
					cost+=currentArc.cost; int nextVertex = currentArc.head;
					if (currentVertex>=1 && currentVertex<=dataModel.C) {
						if(route.containsKey(currentVertex)) { route.replace(currentVertex, route.get(currentVertex)+1); isElementary = false;}
						else route.put(currentVertex, 1);
					}

					currentLabel = vertices[nextVertex].processedLabels.get(currentLabel.nextLabelIndex);
					currentVertex = nextVertex; arcs.add(currentArc.id);
				}

				//Gets the route sequence (of customers)
				int[] routeSequence = new int[arcs.size()-1];
				int counter = 0;
				for(Integer arc: arcs) {
					if(counter>=routeSequence.length) break;
					routeSequence[counter] = dataModel.arcs[arc].head; counter++;
				}
				
				BitSet last_charging_ts = pricingProblem.last_charging_periods.get(label.index);
				for (int t = last_charging_ts.nextSetBit(0); t >= 0; t = last_charging_ts.nextSetBit(t + 1)){

					double rc = reducedCost + pricingProblem.charging_reducedCosts.get(chargingTime).get(t);
					if (rc < -dataModel.precision){
						int initial = t-chargingTime+1;

						Route column = new Route("heuristicLabeling", false, route, routeSequence, pricingProblem, cost, departureTime, energy, load, rc, arcs, initial, chargingTime);
						newRoutes.add(column);
					}

				}
			
			}
			
		}

		return newRoutes;

	}

	/**
	 * Finds disjoint block of routes (to diversify)
	 */
	public List<Route> disjointBlocks(List<Route> newRoutes){

		if(newRoutes.isEmpty()) return newRoutes;
		Collections.sort(newRoutes, new Comparator<Route>() {
			public int compare(Route a, Route b){
				if(a.reducedCost>b.reducedCost) return 1;
				if(a.reducedCost<b.reducedCost) return -1;
				return 0;
			}
		});
		this.objective = newRoutes.get(0).reducedCost;

		//Diversify routes
		int blocks = 10;
		List<Route> disjointRoutes = new ArrayList<Route>(this.numCols);
		int[][] blocksWithCustomer = new int[dataModel.C][blocks];
		for(Route route: newRoutes) {
			for (int j = 0; j < blocks; j++) {
				int similarity = 0;
				for (int i: route.route.keySet()) {similarity+=blocksWithCustomer[i-1][j];}
				if (similarity<= similarityThreshold) {
					for (int i: route.route.keySet()) {blocksWithCustomer[i-1][j]+=1;}
					disjointRoutes.add(route);
					break;
				}
			}
		}
		return disjointRoutes;
	}

	/**
	 * When the Pricing Problem is solved, the set objective function gets invoked first. 
	 */
	@Override
	protected void setObjective() {
		// Este sí
		pricingProblem.reducedCostThreshold = 0.0;
		pricingProblem.bestReducedCost = -Double.MAX_VALUE;
		//Update the objective function with the new dual values
		for (int a = 0; a < dataModel.numArcs; a++) {
			Arc arc = dataModel.arcs[a];
			if (arc.tail>=1 && arc.tail<=dataModel.C) //routing arcs
				arc.modifiedCost = arc.cost-pricingProblem.dualCosts[arc.tail-1];
			else if(arc.tail== 0) arc.modifiedCost = arc.cost; //arcs from the depot source
			else if(arc.tail>dataModel.V) arc.modifiedCost = -pricingProblem.dualCosts[arc.tail-3];
			else arc.modifiedCost = 0;
		}

		//Check charging time branching decisions
		int i=0;
		for(ChargingTimeInequality branching: pricingProblem.branchesOnChargingTimes) {
			if(!branching.lessThanOrEqual) pricingProblem.reducedCostThreshold+= pricingProblem.dualCosts[dataModel.C+dataModel.last_charging_period+pricingProblem.subsetRowCuts.size()+i];
			i++;
		}
	}

	/**
	 * Verifies if a label is dominated. Returns true if it is, false otherwise.
	 * If the label is dominated it is discarded
	 * If the label is not dominated, the existing labels dominated by the label is discarded
	 * @param label to which check dominance
	 */
	public boolean checkDominance(Label newLabel) {

		/* // DELETE BLOCK LATER
		int[] lookup_route = new int[]{0,12,9,3,20,10,1}; // DELETE LATER
		int[] nl_sequence = get_route_sequence(newLabel); // DELETE LATER
		boolean is_nl_subset = false; // DELETE LATER
		if (nl_sequence.length <= lookup_route.length){
			is_nl_subset = sequence_is_subset(nl_sequence, lookup_route);
		} */

		Vertex currentVertex = vertices[newLabel.vertex];
		ArrayList<Label> labelsToDelete = new ArrayList<Label>();
		for(Label existingLabel: currentVertex.unprocessedLabels) {

			/* // DELETE BLOCK LATER
			boolean existing_is_discarded = false; 
			int[] el_sequence = get_route_sequence(existingLabel); 
			boolean is_el_subset = false;
			if (el_sequence.length <= lookup_route.length){ 
				is_el_subset = sequence_is_subset(el_sequence, lookup_route);
			} */

			if(isDominated(existingLabel, newLabel)) {
				//existing_is_discarded = true; // DELETE LATER
				labelsToDelete.add(existingLabel);
			}
			
		}
		currentVertex.unprocessedLabels.removeAll(labelsToDelete);
		if(currentVertex.unprocessedLabels.isEmpty()) nodesToProcess.remove(currentVertex);

		//boolean new_is_discarded = false; // DELETE LATER
		for(Label existingLabel: currentVertex.processedLabels) {
			//int[] el_sequence = get_route_sequence(existingLabel); // DELETE LATER
			if(isDominated(newLabel, existingLabel)) return true;
			
		}

		return false;
	}

	/**
	 * Verifies if a label is dominated. Returns true if it is, false otherwise.
	 * If the label is dominated it is discarded
	 * If the label is not dominated, the existing labels dominated by the label is discarded
	 * @param label to which check dominance
	 */
	public boolean checkDominanceAtDepot(Label newLabel) {
		Vertex currentVertex = vertices[newLabel.vertex];

		ArrayList<Label> labelsToDelete = new ArrayList<Label>();
		for(Label existingLabel: currentVertex.unprocessedLabels) {
			if(isDominatedAtDepot(existingLabel, newLabel)) labelsToDelete.add(existingLabel);
		}

		currentVertex.unprocessedLabels.removeAll(labelsToDelete);
		if(currentVertex.unprocessedLabels.isEmpty()) nodesToProcess.remove(currentVertex);

		for(Label existingLabel: currentVertex.processedLabels) {
			if(isDominatedAtDepot(newLabel, existingLabel))
				return true;
		}
		return false;
	}


	/**
	 * Verifies if L1 is (strongly) dominated by L2
	 * @param L1, L2 labels
	 */
	public boolean isDominated(Label L1, Label L2) {
		
		/* int[] nl_sequence = get_route_sequence(L1); // DELETE LATER
		int[] el_sequence = get_route_sequence(L2); // DELETE LATER */

		if(L1.vertex>dataModel.C) { //charging time vertices
			if (L2.chargingTime>L1.chargingTime) return false;
			if (L2.reducedCost-L1.reducedCost>dataModel.precision) return false;
			return true;

		}else { 					//customer vertices

			if (L1.vertex>0 && L2.remainingLoad<L1.remainingLoad) return false; 	//load
			if (L2.reducedCost-L1.reducedCost>dataModel.precision) return false; 	//reduced cost
			if (L2.remainingTime<L1.remainingTime) return false; 					//time
			
			for (int gam = 0; gam <= dataModel.gamma; gam ++){
				if (L2.remainingEnergy[gam]<L1.remainingEnergy[gam]) return false;				 //energy
			}
			
			//reducedCost
			double reducedCostL2 = 0;
			if (L1.vertex>0) {
				for(int i: L2.srcIndices) {
					if(!L1.eta[i]) {
						SubsetRowInequality src = pricingProblem.subsetRowCuts.get(i);
						if(!L2.unreachable[src.cutSet[0]-1] || !L2.unreachable[src.cutSet[1]-1] || !L2.unreachable[src.cutSet[2]-1]) {
							int dualIndex = dataModel.C+dataModel.last_charging_period+i;
							reducedCostL2+=pricingProblem.dualCosts[dualIndex];
						}
					}
					if (L2.reducedCost-reducedCostL2-L1.reducedCost>dataModel.precision) return false;
				}
			}

			if (L2.reducedCost-reducedCostL2-L1.reducedCost>dataModel.precision) return false;

			/* // Ng-paths and unreachable resources
			Vertex currentVertex = vertices[L1.vertex];
			if (currentVertex.id > 0) {
				for(int i=1; i<=dataModel.C; i++) {
					
					//boolean check_binaries = (L2.ng_path[i-1] || L2.unreachable[i-1]) && !(L1.ng_path[i-1] || L1.unreachable[i-1]);
					boolean other_way = L2.ng_path[i-1] && (!L1.unreachable[i-1] && !L1.ng_path[i-1]); // Dani's way
					if (other_way) {
						return false;
					}
				}
			} */

			return true;
		}
	}

	/**
	 * Verifies if L1 is (strongly) dominated by L2 (Heuristic rule)
	 * @param L1, L2 labels
	 */
	public boolean isDominatedAtDepot(Label L1, Label L2) {

		if (L1.vertex>0 && L1.remainingLoad > L2.remainingLoad) return false; 		// Load
		if ((int)(L1.remainingTime/10) > (int)(L2.remainingTime/10)) return false; 	// Departure Time
		if (L1.chargingTime < L2.chargingTime) return false;						// Charging Time
		if (L2.reducedCost > L1.reducedCost + dataModel.precision) return false;	// Reduced Cost
		
		return true;
	}


	public int[] get_route_sequence(Label label) {

		Label new_label = label.clone();

		ArrayList<Integer> arcs = new ArrayList<Integer>(dataModel.C);
		int currentVertex = new_label.vertex;
		while(currentVertex!=dataModel.C+1) {
			Arc currentArc = dataModel.arcs[new_label.nextArc];
			int nextVertex = currentArc.head;

			new_label = vertices[nextVertex].processedLabels.get(new_label.nextLabelIndex);
			if(currentArc.tail>=0 && currentArc.tail<=dataModel.C) arcs.add(currentArc.id);
			currentVertex = nextVertex;
		}

		//Gets the route sequence (of customers)
		if (arcs.size() > 0){

			int[] routeSequence = new int[arcs.size()];
			routeSequence[0] = label.vertex;

			int counter = 1;
			for(Integer arc: arcs) {
				if(counter>=routeSequence.length) break;
				routeSequence[counter] = dataModel.arcs[arc].head;
				counter++;
			}

			return routeSequence;
		} else {

			return new int[0];
		}
	}

	public boolean sequence_is_subset(int[] sequence, int[] lookup_sequence){

		int m = sequence.length;
        int n = lookup_sequence.length;

        // Compare the last m elements of n2 with n1
        for (int i = 0; i < m; i++) {
            if (lookup_sequence[n - m + i] != sequence[i]) {
                return false;
            }
        }
        return true;
	}

	/**
	 * Listen to branching decisions. The pricing problem is changed by the branching decisions.
	 * @param bd BranchingDecision
	 */
	@Override
	public void branchingDecisionPerformed(BranchingDecision bd) {
		if(bd instanceof FixArc) { 			//Fixing one arc
			FixArc fixArcDecision = (FixArc) bd;
			for(int infeasibleArc: fixArcDecision.infeasibleArcs) this.infeasibleArcs[infeasibleArc] ++;
		}else if(bd instanceof RemoveArc) {//Removing one arc
			RemoveArc removeArcDecision= (RemoveArc) bd;
			infeasibleArcs[removeArcDecision.arc] ++;
		}
	}

	/**
	 * When the Branch-and-Price algorithm backtracks, branching decisions are reversed.
	 * @param bd BranchingDecision
	 */
	@Override
	public void branchingDecisionReversed(BranchingDecision bd) {
		if(bd instanceof FixArc) { 			//Fixing one arc
			FixArc fixArcDecision = (FixArc) bd;
			for(int infeasibleArc: fixArcDecision.infeasibleArcs) this.infeasibleArcs[infeasibleArc] --;
		}else if(bd instanceof RemoveArc) {//Removing one arc
			RemoveArc removeArcDecision= (RemoveArc) bd;
			infeasibleArcs[removeArcDecision.arc] --;
		}
	}

	/**
	 * Returns the time in seconds (and considering two decimals)
	 */
	public double getTimeInSeconds(double time) {
		double realTime = time*0.001;
		realTime = Math.floor(realTime*100)/100; //two decimals
		return realTime;
	}

	/**
	 * @return a negative integer, zero, or a positive integer as this object is less than, equal to, or greater than the specified object.
	 */
	public class SortVertices implements Comparator<Vertex> {

		@Override
		public int compare(Vertex vertex1, Vertex vertex2) {

			if(vertex2.id==0 && (vertex1.id>0 && vertex1.id<=dataModel.C)) return -1;
			if(vertex1.id==0 && (vertex2.id>0 && vertex2.id<=dataModel.C)) return 1;

			if(vertex1.id<dataModel.V && vertex2.id>=dataModel.V) return -1;
			if(vertex1.id>=dataModel.V && vertex2.id<dataModel.V) return 1;
			if(vertex1.id>=dataModel.V && vertex2.id>=dataModel.V) {
				if(vertex1.id>vertex2.id) return -1;
				else return 1;
			}

			Label L1 = vertex1.unprocessedLabels.peek();
			Label L2 = vertex2.unprocessedLabels.peek();
			if(L1.remainingLoad>L2.remainingLoad) return -1;
			if(L1.remainingLoad<L2.remainingLoad) return 1;
			if(L1.remainingEnergy[dataModel.gamma]>L2.remainingEnergy[dataModel.gamma]) return -1;
			if(L1.remainingEnergy[dataModel.gamma]<L2.remainingEnergy[dataModel.gamma]) return 1;
			if(L1.remainingTime>L2.remainingTime) return -1;
			if(L1.remainingTime<L2.remainingTime) return 1;
			if(L1.reducedCost<L2.reducedCost) return -1;
			if(L1.reducedCost>L2.reducedCost) return 1;
			return 0;
		}
	}
}
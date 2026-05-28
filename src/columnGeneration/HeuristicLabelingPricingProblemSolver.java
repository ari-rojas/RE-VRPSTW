package columnGeneration;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.function.BiPredicate;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;
import branchAndPrice.ChargingTimeInequality;
import branchAndPrice.FixArc;
import branchAndPrice.RemoveArc;
import model.EVRPTW;
import model.EVRPTW.Arc;
import model.EVRPTW.Vertex;

/**
 * This class provides a heuristic solver for the ng-SPPRC pricing problem
 * It considers only the min-cost arcs and uses a relaxed dominance rule
 */
public final class HeuristicLabelingPricingProblemSolver extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem> {

	public Vertex[] vertices = dataModel.vertices; 						//vertices of the instance
	public PriorityQueue<Vertex> nodesToProcess; 						//labels that need be processed
	public final int numCols = 400; 									//maximum number of routes (columns) allowed
	public int[] infeasibleArcs; 										//arcs that cannot be used by branching
	public final int similarityThreshold = 5; 							//for the disjoint columns diversification strategy

	private int Gamma;
	private int depotID;

	/** Heuristic Labeling algorithm to solve the ng-SPPRC. */
	public HeuristicLabelingPricingProblemSolver(EVRPTW dataModel, PricingProblem pricingProblem) {
		super(dataModel, pricingProblem);
		this.name="HeuristicLabelingSolver"; //Set a name for the solver
		this.infeasibleArcs = new int[dataModel.numArcs];
		this.nodesToProcess = new PriorityQueue<Vertex>(dataModel.numVertices, new SortVertices());
		this.Gamma = dataModel.gamma;
		this.depotID = dataModel.C+1;
	}

	/**
	 * Runs the labeling algorithm
	 */
	public void runLabeling() {

		// Initialization
		int[] remain_energy = new int[Gamma + 1]; Arrays.fill( remain_energy, dataModel.E);
		Label initialLabel = new Label(0, -pricingProblem.dualCost, dataModel.Q, vertices[dataModel.C+1].closing_tw, remain_energy, 0,new boolean[dataModel.C], new boolean[dataModel.C], new boolean[pricingProblem.subsetRowCuts.size()], new HashSet<Integer>(pricingProblem.subsetRowCuts.size()));
		initialLabel.index = 0; initialLabel.vertex = depotID; initialLabel.vertex = depotID; initialLabel.nextArc = depotID;
		this.nodesToProcess.add(vertices[depotID]);
		vertices[depotID].unprocessedLabels.add(initialLabel);
		
		////////////////////////////////////////////
		/// Routing Labeling
		////////////////////////////////////////////
		
		long startTime = System.currentTimeMillis();
		while (!nodesToProcess.isEmpty() && System.currentTimeMillis()<timeLimit) {
			ArrayList<Label> labelsToProcessNext = routingLabelsToProcessNext();
			Set<Arc> incomingArcs = new HashSet<Arc>(dataModel.graph.incomingEdgesOf(labelsToProcessNext.get(0).vertex));
			incomingArcs.removeIf(arc -> infeasibleArcs[arc.id] > 0);
			
			for (Label currentLabel: labelsToProcessNext) {
				
				boolean isDominated = checkDominance(currentLabel);
				if(isDominated) continue;
				
				currentLabel.index = vertices[currentLabel.vertex].processedLabels.size();
				vertices[currentLabel.vertex].processedLabels.add(currentLabel);
				
				for(Arc a: incomingArcs) {
					Label extendedLabel = extendLabel(currentLabel, a);
					if (extendedLabel!=null) { // Verifies if the extension is feasible
						extendedLabel.vertex = a.tail;
						extendedLabel.nextArc = a.id;

						Vertex extendedVertex = vertices[extendedLabel.vertex];
						extendedVertex.unprocessedLabels.add(extendedLabel);
						if (extendedVertex.id > 0 && extendedVertex.unprocessedLabels.size() == 1) nodesToProcess.add(extendedVertex);
						
					}
				}
			}
		}

		////////////////////////////////////////////
		/// Middle point
		////////////////////////////////////////////
		
		if (System.currentTimeMillis()>=timeLimit) vertices[0].unprocessedLabels.clear();
		if (!vertices[0].unprocessedLabels.isEmpty()) nodesToProcess.add(vertices[0]);

		/////////////////////////////////////
		/// Charging Scheduling Labeling
		////////////////////////////////////
		
		while (!nodesToProcess.isEmpty() && System.currentTimeMillis()<timeLimit) {
			ArrayList<Label> labelsToProcessNext = chargingLabelsToProcessNext();
			
			for (Label currentLabel: labelsToProcessNext) {
				
				boolean isDominated = checkDominance(currentLabel);
				if(isDominated) continue;
				
				currentLabel.index = vertices[currentLabel.vertex].processedLabels.size();
				vertices[currentLabel.vertex].processedLabels.add(currentLabel);
				
				for (Arc a: dataModel.graph.incomingEdgesOf(currentLabel.vertex)) {
					Label extendedLabel = extendLabelChargingTime(currentLabel, a);
					if (extendedLabel!=null) { //verifies if the extension is feasible
						
						extendedLabel.vertex = a.tail;
						extendedLabel.nextArc = a.id;

						Vertex extendedVertex = vertices[extendedLabel.vertex];
						extendedVertex.unprocessedLabels.add(extendedLabel);
						if (extendedVertex.id != dataModel.V && extendedVertex.unprocessedLabels.size() == 1) nodesToProcess.add(extendedVertex);
					}
				}
			}
		}

		long totalTime = System.currentTimeMillis()-startTime;
		dataModel.heuristicPricingTime+=totalTime;
		if (dataModel.print_log) logger.debug("Time solving (exactly) the pricing problem (s): " + getTimeInSeconds(totalTime)); 
	}

	/**
	 * Selects a set of labels to process (the one with the most remaining load)
	 */
	public ArrayList<Label> routingLabelsToProcessNext(){

		ArrayList<Label> labelsToProcessNext = new ArrayList<Label>();
		Vertex currentVertex = nodesToProcess.poll();
		
		while(true) {
			Label currentLabel = currentVertex.unprocessedLabels.poll();
			if(labelsToProcessNext.isEmpty()) labelsToProcessNext.add(currentLabel);
			else {
				boolean isDominated = false;
				for(Label L2: labelsToProcessNext) {
					
					isDominatedRouting(currentLabel, L2);
					if(isDominated) break;
				}
				if (!isDominated) labelsToProcessNext.add(currentLabel);
			}
			if(currentVertex.unprocessedLabels.isEmpty() || (currentVertex.unprocessedLabels.peek().remainingLoad<currentLabel.remainingLoad)) break;
		}

		if(!currentVertex.unprocessedLabels.isEmpty()) nodesToProcess.add(currentVertex);
		return labelsToProcessNext;
	}

	public ArrayList<Label> chargingLabelsToProcessNext(){

		ArrayList<Label> labelsToProcessNext = new ArrayList<Label>();
		Vertex currentVertex = nodesToProcess.poll();
		
		while (!currentVertex.unprocessedLabels.isEmpty()) labelsToProcessNext.add(currentVertex.unprocessedLabels.poll());
		return labelsToProcessNext;
		
	}

	/**
	 * This class is invoked when only nonelementary routes are found. 
	 * @return true if the maximum size per neighborhood has been reached
	 */
	public boolean enlargeNeighborhoods(List<Route> nonElementaryRoutes) {

		boolean enlarged = false;
		ArrayList<Integer> cyclingVertics = new ArrayList<Integer>();
		for(Route route:nonElementaryRoutes) {
			ArrayList<Integer> visitedCustomers = new ArrayList<Integer>(dataModel.C);
			boolean[] visited = new boolean[dataModel.C];
			for(int arc: route.arcs) {
				int head = dataModel.arcs[arc].head;
				if(head<dataModel.C) {
					visitedCustomers.add(head);
					if(visited[head-1]) { //there is a cycle
						for (int i = visitedCustomers.size()-2; i >=0; i--) {
							int node = visitedCustomers.get(i);
							if(node == head) {cyclingVertics.add(head); break;}
							else if(!dataModel.vertices[node].neighbors.contains(head) && dataModel.vertices[node].neighbors.size()<=dataModel.DeltaMax) { 
								dataModel.vertices[node].neighbors.add(head); 
								enlarged = true;
								if (dataModel.print_log) logger.debug("Adding: " + head + " to the neighborhood of: "+node + " (size=" + dataModel.vertices[node].neighbors.size()+")");
							}
						}
					}else visited[head-1] = true;
				}
			}
		}
		return enlarged;
	}

	/**
	 * Label extension procedure
	 */
	public Label extendLabel(Label currentLabel, Arc arc) {

		int source = arc.tail;
		if (source>=1 && source<=dataModel.C)
			if (currentLabel.unreachable[source-1]|| currentLabel.ng_path[source-1]) return null;

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

		////////////////////////////////////////////
		/// Bounding Procedure
		////////////////////////////////////////////

		if (source == 0){
			double min_rc = reducedCost + pricingProblem.charging_bounds.get(chargingTime).get((int)(remainingTime/10)); 
			if (min_rc >= -dataModel.precision)  return null;
		}

		boolean[] unreachable = Arrays.copyOf(currentLabel.unreachable.clone(), currentLabel.unreachable.length);

		//Mark unreachable customers and ng-path cycling restrictions
		if(source>0) {
			unreachable[source-1] = true;
			for (int i: vertices[source].unreachable) unreachable[i-1] = true;
			
			for (Arc c: dataModel.graph.incomingEdgesOf(source)) {
				if(c.tail==0 || unreachable[c.tail-1]) continue;
				//unreachable
				if (remainingLoad-vertices[c.tail].load<0 || remainingTime-c.min_time<vertices[c.tail].opening_tw || remainingEnergy[dataModel.gamma]-c.min_energy < 0 || 
					Math.min(remainingTime-c.min_time, vertices[c.tail].closing_tw)-dataModel.graph.getEdge(0, c.tail).min_time<vertices[0].opening_tw ||
					remainingEnergy[dataModel.gamma]-c.min_energy - dataModel.graph.getEdge(0, c.tail).min_energy<0) {
					unreachable[c.tail-1] = true;
				}
			}
		} else {
			remainingTime = (int) (remainingTime/10);
		}

		Label extendedLabel = new Label(currentLabel.index, reducedCost, remainingLoad, remainingTime, remainingEnergy, chargingTime,unreachable, currentLabel.ng_path, eta, srcIndices);
		return extendedLabel;

	}

	/**
	 * Label extension procedure
	 */
	public Label extendLabelChargingTime(Label currentLabel, Arc arc) {

		int source = arc.tail;

		if(arc.head==0 && (source-dataModel.V<currentLabel.chargingTime || source-dataModel.V>=currentLabel.remainingTime)) return null;
		if(source == dataModel.V && currentLabel.chargingTime>0 ) return null;

		double reducedCost = currentLabel.reducedCost+arc.modifiedCost;
		reducedCost = Math.floor(reducedCost*10000)/10000;
		int chargingTime = currentLabel.chargingTime;
		if(source!=dataModel.V) {
			chargingTime-=1;
			if(chargingTime<0) return null;
		}

		if (source == dataModel.V){
			double rc = currentLabel.reducedCost;
			if (rc >= -dataModel.precision) return null; // Only negative reduced costs labels will get to the source node
		}

		Label extendedLabel = new Label(currentLabel.index, reducedCost, currentLabel.remainingLoad, currentLabel.remainingTime, currentLabel.remainingEnergy, chargingTime , currentLabel.unreachable, currentLabel.ng_path, currentLabel.eta, currentLabel.srcIndices);
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
	 * When the CG procedure terminates, the close function is invoked. 
	 */
	public void restart() {
		for (int i = 0; i < vertices.length; i++) {
			vertices[i].processedLabels = new ArrayList<Label>(dataModel.numArcs);
			vertices[i].unprocessedLabels =  new PriorityQueue<Label>(dataModel.numArcs, new Label.SortLabels());
		}
		this.nodesToProcess = new PriorityQueue<Vertex>(new SortVertices());
	}

	/**
	 * This method produces zero or more columns. 
	 */
	@Override
	protected List<Route> generateNewColumns() {

		List<Route> newRoutes=new ArrayList<>(this.numCols);  			//list of routes
		
		this.runLabeling();

		if(vertices[dataModel.V].unprocessedLabels.isEmpty()) {
			pricingProblemInfeasible=true; this.objective=Double.MAX_VALUE;
			
		} else {
			this.pricingProblemInfeasible=false;
			for (Label label: vertices[dataModel.V].unprocessedLabels) {
				int departureTime = label.remainingTime;
				int load = dataModel.Q - label.remainingLoad;
				if (label.reducedCost<=-dataModel.precision) {		//generate new column if it has negative reduced cost
					
					HashMap<Integer, Integer> route=new HashMap<Integer, Integer>(dataModel.C); int cost = 0; int energy = dataModel.E-label.remainingEnergy[dataModel.gamma]; double reducedCost = label.reducedCost;
					ArrayList<Integer> arcs = new ArrayList<Integer>(dataModel.C);
					int initialChargingTime = dataModel.arcs[label.nextArc].head-dataModel.V; int chargingTime = 0;
					int currentVertex = label.vertex;
					while(currentVertex!=dataModel.C+1) {
						Arc currentArc = dataModel.arcs[label.nextArc];
						cost+=currentArc.cost;
						int nextVertex = currentArc.head;
						if (currentVertex>=1 && currentVertex<=dataModel.C) {
							route.put(currentVertex, 1);
						}else if(currentVertex!=dataModel.V && currentVertex!=0) chargingTime++;

						label = vertices[nextVertex].processedLabels.get(label.nextLabelIndex);
						if(currentArc.tail>=0 && currentArc.tail<=dataModel.C) arcs.add(currentArc.id);
						currentVertex = nextVertex;
					}

					//Gets the route sequence (of customers)
					int[] routeSequence = new int[arcs.size()-1];
					int counter = 0;
					for(Integer arc: arcs) {
						if(counter>=routeSequence.length) break;
						routeSequence[counter] = dataModel.arcs[arc].head;
						counter++;
					}
					Route column = new Route("exactLabeling", false, route, routeSequence, pricingProblem, cost, departureTime, energy, load, reducedCost, arcs, initialChargingTime, chargingTime);
					newRoutes.add(column);
				}
			}

		}

		if (dataModel.print_log) {
				logger.debug("Finished exact pricing: "+vertices[0].processedLabels.size()+" processed, "+vertices[0].unprocessedLabels.size()+" unprocessed.");
				logger.debug("Found " + newRoutes.size() + " columns");
		}
		
		close();
		return newRoutes;
	}

	/**
	 * When the Pricing Problem is solved, the set objective function gets invoked first. 
	 */
	@Override
	protected void setObjective() {

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
			if(branching.startCharging) dataModel.graph.getEdge(dataModel.V, dataModel.V+branching.timestep).modifiedCost-=pricingProblem.dualCosts[dataModel.C+dataModel.last_charging_period+pricingProblem.subsetRowCuts.size()+i];
			else dataModel.graph.getEdge(dataModel.V+branching.timestep,0).modifiedCost-=pricingProblem.dualCosts[dataModel.C+dataModel.last_charging_period+pricingProblem.subsetRowCuts.size()+i];
			if(!branching.lessThanOrEqual) pricingProblem.reducedCostThreshold+= pricingProblem.dualCosts[dataModel.C+dataModel.last_charging_period+pricingProblem.subsetRowCuts.size()+i];
			i++;
		}
	}

	public BiPredicate<Label, Label> getDominanceChecker(int vx_id){

		if (vx_id == 0) return this::isDominatedDepot;
		if (vx_id <= depotID) return this::isDominatedRouting;
		
		return this::isDominatedCharging;
	}

	public boolean checkDominance(Label newLabel) {
		
		/* // DELETE BLOCK LATER
		int[] lookup_route = new int[]{0,12,9,3,20,10,1}; // DELETE LATER
		int[] nl_sequence = get_route_sequence(newLabel); // DELETE LATER
		boolean is_nl_subset = false; // DELETE LATER
		if (nl_sequence.length <= lookup_route.length){
			is_nl_subset = sequence_is_subset(nl_sequence, lookup_route);
			} */
		
		Vertex currentVertex = vertices[newLabel.vertex];
		BiPredicate<Label, Label> isDominatedChecker = getDominanceChecker(currentVertex.id);

		ArrayList<Label> labelsToDelete = new ArrayList<Label>();
		for(Label existingLabel: currentVertex.unprocessedLabels) {

			/* // DELETE BLOCK LATER
			boolean existing_is_discarded = false; 
			int[] el_sequence = get_route_sequence(existingLabel); 
			boolean is_el_subset = false;
			if (el_sequence.length <= lookup_route.length){ 
				is_el_subset = sequence_is_subset(el_sequence, lookup_route);
			} */

			if(isDominatedChecker.test(existingLabel, newLabel)) {
				//existing_is_discarded = true; // DELETE LATER
				labelsToDelete.add(existingLabel);
			}
			
		}
		currentVertex.unprocessedLabels.removeAll(labelsToDelete);
		if(currentVertex.unprocessedLabels.isEmpty()) nodesToProcess.remove(currentVertex);

		//boolean new_is_discarded = false; // DELETE LATER
		for(Label existingLabel: currentVertex.processedLabels) {
			//int[] el_sequence = get_route_sequence(existingLabel); // DELETE LATER
			if(isDominatedChecker.test(newLabel, existingLabel)) return true;
		}

		return false;
	}

	/**
	 * Verifies if L1 is (strongly) dominated by L2
	 * @param L1, L2 labels
	 */
	public boolean isDominatedDepot(Label L1, Label L2) {
		
		/* int[] nl_sequence = get_route_sequence(L1); // DELETE LATER
		int[] el_sequence = get_route_sequence(L2); // DELETE LATER */

		if (L2.reducedCost-L1.reducedCost>dataModel.precision) return false; 	//reduced cost
		if (L2.remainingTime<L1.remainingTime) return false; 					//departure time
		if (L2.chargingTime>L1.chargingTime) return false;						//charging time

		return true;
	}

	/**
	 * Verifies if L1 is (strongly) dominated by L2
	 * @param L1, L2 labels
	 */
	public boolean isDominatedRouting(Label L1, Label L2) {
		
		/* int[] nl_sequence = get_route_sequence(L1); // DELETE LATER
		int[] el_sequence = get_route_sequence(L2); // DELETE LATER */

		if (L2.remainingLoad<L1.remainingLoad) return false; 	//load
		if (L2.reducedCost-L1.reducedCost>dataModel.precision) return false; 	//reduced cost
		if (L2.remainingTime<L1.remainingTime) return false; 					//time
		
		for (int gam=0; gam<=Gamma; gam++){
			if (L2.remainingEnergy[gam]<L1.remainingEnergy[gam]) return false; //energy
		}
		
		//reducedCost
		double reducedCostL2 = 0;
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

		if (L2.reducedCost-reducedCostL2-L1.reducedCost>dataModel.precision) return false;

		return true;
	}

	public boolean isDominatedCharging(Label L1, Label L2) {
		
		if (L2.chargingTime>L1.chargingTime) return false;
		if (L2.reducedCost-L1.reducedCost>dataModel.precision) return false;
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
		}else if(bd instanceof RemoveArc) {	//Removing one arc
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

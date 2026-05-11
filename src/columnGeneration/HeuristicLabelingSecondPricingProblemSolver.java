package columnGeneration;

import model.EVRPTW;
import model.EVRPTW.Arc;
import model.EVRPTW.PPArc;
import model.EVRPTW.Vertex;
import model.EVRPTW.PPVertex;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;

import org.jgrapht.graph.DirectedWeightedMultigraph;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;

import branchAndPrice.FixArc;
import branchAndPrice.RemoveArc;


/**
 * This class provides a heuristic solver for the ng-SPPRC pricing problem
 * It uses a relaxed dominance rule
 */
public final class HeuristicLabelingSecondPricingProblemSolver extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem> {

	public PPVertex[] PPvertices = dataModel.PPvertices; 			//vertices of the instance
	public Vertex[] vertices = dataModel.vertices;
	public PriorityQueue<PPVertex> nodesToProcess; 			//labels that need be processed
	public final int numCols = 400; 						//maximum number of routes (columns) allowed
	public int[] infeasibleArcs; 						//arcs that cannot be used by branching
	public final int similarityThreshold = 5; 				//diversification of columns

	public double bestReducedCost;

	// Identifiers for the differnt types of vertices in the Pricing Problem Graph
	public static final byte C0 = EVRPTW.C0; 	 	// Customer depot nodes
	public static final byte C1 = EVRPTW.C1; 		// Non-first customer nodes
	public static final byte Tt = EVRPTW.Tt; 		// Charging time period nodes
	public static final byte Source = EVRPTW.Source;	// Dummy source node

	// Identifiers for the different types of arcs in the Pricing Problem Graph
	public static final byte AR0 = EVRPTW.AR0;	// Routing arcs between customer depot nodes and non-first customer nodes
	public static final byte AR1 = EVRPTW.AR1;	// Routing arcs between non-first customer nodes
	public static final byte AC1 = EVRPTW.AC1;	// Charging Scheduling arcs to finishing charging times
	public static final byte AC2 = EVRPTW.AC2;	// Charging Scheduling arcs between consecutive charging time periods
	public static final byte AC3 = EVRPTW.AC3;	// Charging Scheduling arcs to starting charging times

	public final int depotID;

	/**
	 * Labeling algorithm to solve the ng-SPPRC
	 */
	public HeuristicLabelingSecondPricingProblemSolver(EVRPTW dataModel, PricingProblem pricingProblem) {
		super(dataModel, pricingProblem);
		this.name="HeuristicLabelingSolver"; //Set a name for the solver
		this.infeasibleArcs = new int[dataModel.numArcs];
		this.nodesToProcess = new PriorityQueue<PPVertex>(dataModel.V, new SortVertices());
		this.depotID = dataModel.T_startID;
	}

	/**
	 * Runs the labeling algorithm
	 */
	public void runLabeling() {

		this.bestReducedCost = Double.MAX_VALUE;
		//Initialization
		int[] remain_energy = new int[dataModel.gamma + 1]; Arrays.fill( remain_energy, dataModel.E);
		Label initialLabel = new Label(0, -pricingProblem.dualCost, dataModel.Q, vertices[dataModel.C+1].closing_tw, remain_energy, 0,new boolean[dataModel.C], new boolean[dataModel.C], new boolean[pricingProblem.subsetRowCuts.size()], new HashSet<Integer>(pricingProblem.subsetRowCuts.size()));
		initialLabel.index = 0; initialLabel.vertex = depotID; initialLabel.nextArc = depotID;
		this.nodesToProcess.add(PPvertices[depotID]);
		PPvertices[depotID].unprocessedLabels.add(initialLabel);

		//Labeling algorithm
		long startTime = System.currentTimeMillis();
		while (!nodesToProcess.isEmpty() && System.currentTimeMillis()<timeLimit) {
			ArrayList<Label> labelsToProcessNext = labelsToProcessNext();
			for(Label currentLabel: labelsToProcessNext) {
				boolean isDominated = checkDominance(currentLabel);
				if(isDominated) continue;
				else {currentLabel.index = PPvertices[currentLabel.vertex].processedLabels.size(); PPvertices[currentLabel.vertex].processedLabels.add(currentLabel);}
				
				for(PPArc a: dataModel.PPgraph.incomingEdgesOf(currentLabel.vertex)) {
					if(infeasibleArcs[a.id] > 0) continue;
					Label extendedLabel;
					if(a.arc_type <= AR1) extendedLabel = extendLabel(currentLabel, a.routing_arc, a.arc_type, a.modifiedCost);
					else extendedLabel = extendLabelChargingTime(currentLabel, PPvertices[a.tail_vertex_id].node_number, a.arc_type, a.modifiedCost);
					if (extendedLabel!=null) { //verifies if the extension is feasible
						extendedLabel.vertex = a.tail_vertex_id;
						extendedLabel.nextArc = a.id;
						updateNodesToProcess(extendedLabel);
					}
				}
			}
		}

		long totalTime = System.currentTimeMillis()-startTime;
		dataModel.heuristicPricingTime+=totalTime;
		if (dataModel.print_log) logger.debug("Time solving (heuristically) the pricing problem (s): " + getTimeInSeconds(totalTime)); 
	}


	/**
	 * Selects a set of labels to process (the one with the most remaining load)
	 */
	public ArrayList<Label> labelsToProcessNext(){

		ArrayList<Label> labelsToProcessNext = new ArrayList<Label>();
		PPVertex currentVertex = nodesToProcess.poll();
		byte vertex_type = currentVertex.vertex_type;
		while(true) {
			Label currentLabel = currentVertex.unprocessedLabels.poll();
			if(labelsToProcessNext.isEmpty()) labelsToProcessNext.add(currentLabel);
			else {
				boolean isDominated = false;
				for(Label L2: labelsToProcessNext) {
					
					if (vertex_type <= C1) isDominated = isDominatedRouting(currentLabel, L2, vertex_type);
					else isDominated = isDominatedCharging(currentLabel, L2);
					if(isDominated) break;
				}
				if(!isDominated) labelsToProcessNext.add(currentLabel);
			}
			if(currentVertex.unprocessedLabels.isEmpty() || (vertex_type<=C1 && currentVertex.unprocessedLabels.peek().remainingLoad<currentLabel.remainingLoad)) break;
		}

		if(!currentVertex.unprocessedLabels.isEmpty()) nodesToProcess.add(currentVertex);
		return labelsToProcessNext;
	}


	/**
	 * Given a new (non-dominated) label, updates the nodes to be processed
	 */
	public void updateNodesToProcess(Label extendedLabel) {
		PPVertex currentVertex = PPvertices[extendedLabel.vertex];
		if(currentVertex.vertex_type == Source) PPvertices[extendedLabel.vertex].unprocessedLabels.add(extendedLabel);
		else if(currentVertex.unprocessedLabels.isEmpty()) {currentVertex.unprocessedLabels.add(extendedLabel); nodesToProcess.add(currentVertex);}
		else currentVertex.unprocessedLabels.add(extendedLabel);
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
				if(head <= dataModel.C) {
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
	public Label extendLabel(Label currentLabel, Arc routing_arc, byte arc_type, double modifiedCost) {

		int source = routing_arc.tail;
		if (currentLabel.unreachable[source-1] || currentLabel.ng_path[source-1]) return null;

		// Update the remaining time and check feasibility
		int remainingTime = currentLabel.remainingTime-routing_arc.time;
		if (arc_type == AR1 && remainingTime < vertices[source].open_tw_nonfirst) return null;
		if(remainingTime>vertices[source].closing_tw) remainingTime = vertices[source].closing_tw;

		double reducedCost = currentLabel.reducedCost+modifiedCost;
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

		// Only negative reduced cost labels at the depot
		if (arc_type == AR0 && reducedCost >= pricingProblem.reducedCostThreshold-dataModel.precision) return null;
		
		int[] remainingEnergy = new int[dataModel.gamma + 1];
		remainingEnergy[0] = currentLabel.remainingEnergy[0]-routing_arc.energy; if (remainingEnergy[0] < 0) return null;
		for (int gam = 1; gam <= dataModel.gamma; gam++){
			if (currentLabel.remainingEnergy[gam-1] - routing_arc.energy_deviation < currentLabel.remainingEnergy[gam]){ remainingEnergy[gam] = currentLabel.remainingEnergy[gam-1] - routing_arc.energy - routing_arc.energy_deviation; }
			else { remainingEnergy[gam] = currentLabel.remainingEnergy[gam] - routing_arc.energy; }
			if (remainingEnergy[gam] < 0) return null;
		}
		
		// If the arc connects to a vertex i0 \in C0, then the time and energy resources must account for the depot-i arc
		if (arc_type == AR0){
			Arc depotArc = dataModel.graph.getEdge(0, source);
			remainingTime -= depotArc.time;
			remainingEnergy[0] -= depotArc.energy; if (remainingEnergy[0] < 0) return null;
			for (int gam = 1; gam <= dataModel.gamma; gam++){
				if (remainingEnergy[gam-1] - depotArc.energy_deviation < remainingEnergy[gam]){ remainingEnergy[gam] = remainingEnergy[gam-1] - depotArc.energy - depotArc.energy_deviation; }
				else { remainingEnergy[gam] -= depotArc.energy; }
				if (remainingEnergy[gam] < 0) return null;
			}
		}
		
		// Update charging time and check if it's feasible
		int chargingTime = dataModel.f_inverse[dataModel.E-remainingEnergy[dataModel.gamma]];
		if (chargingTime >= (int) (remainingTime/10)) return null;
		
		// After confirming that the label is feasible, update the remaining load
		int remainingLoad = currentLabel.remainingLoad-vertices[source].load;

		////////////////////////////////////////////
		/// Bounding Procedure
		////////////////////////////////////////////
		
		if (arc_type == AR0 && reducedCost + pricingProblem.charging_bounds.get(chargingTime).get((int)(remainingTime/10)) >= -dataModel.precision) return null;

		boolean[] unreachable = Arrays.copyOf(currentLabel.unreachable.clone(), currentLabel.unreachable.length);
		boolean[] ng_path = new boolean[dataModel.C];
		
		// Mark unreachable customers and ng-path cycling restrictions
		if(arc_type == AR1) {
			ng_path[source-1] = true;
			for (Arc c: dataModel.graph.incomingEdgesOf(source)) {
				if(c.tail==0 || unreachable[c.tail-1]) continue;
				//unreachable
				if (remainingLoad-vertices[c.tail].load<0 || remainingTime-c.min_time<vertices[c.tail].opening_tw || remainingEnergy[dataModel.gamma] - c.min_energy - dataModel.graph.getEdge(0, c.tail).min_energy < 0) {
					unreachable[c.tail-1] = true; }

				//ng-path
				if (currentLabel.ng_path[c.tail-1] && vertices[source].neighbors.contains(c.tail)) ng_path[c.tail-1] = true;
				else ng_path[c.tail-1] = false;
			}
		} else {
			ng_path = Arrays.copyOf(currentLabel.ng_path, currentLabel.ng_path.length);
			// We re-scale the remaining time, which represents the departure time of the route
			remainingTime = (int) (remainingTime/10);
		}

		Label extendedLabel = new Label(currentLabel.index, reducedCost, remainingLoad, remainingTime, remainingEnergy, chargingTime,unreachable, ng_path, eta, srcIndices);
		return extendedLabel;

	}

	/**
	 * Label extension procedure
	 */
	public Label extendLabelChargingTime(Label currentLabel, int t, byte arc_type, double modifiedCost) {

		// If t is a finishing charging time period and is not a value between b_r and d_r-1, it's not feasible
		if(arc_type == AC1 && (t < currentLabel.chargingTime || t >= currentLabel.remainingTime)) return null;
		// If the label is extended to the source node but has not charged enough, it's not feasible
		if(arc_type == AC3 && currentLabel.chargingTime>0 ) return null;

		double reducedCost = currentLabel.reducedCost+modifiedCost;
		reducedCost = Math.floor(reducedCost*10000)/10000;
		int chargingTime = currentLabel.chargingTime;
		if(arc_type < AC3) {
			chargingTime -= 1;
			if(chargingTime<0) return null; // If the label is extended through consecutive charging time periods and is charging more than necessary, deem it infeasible
		} else {
			if (reducedCost < this.bestReducedCost - dataModel.precision) this.bestReducedCost = reducedCost;
			if (reducedCost > -dataModel.precision) return null; // Only negative reduced costs labels will get to the source node
		}

		Label extendedLabel = new Label(currentLabel.index, reducedCost, currentLabel.remainingLoad, currentLabel.remainingTime, currentLabel.remainingEnergy, chargingTime , currentLabel.unreachable, currentLabel.ng_path, currentLabel.eta, currentLabel.srcIndices);
		return extendedLabel;
	}


	/**
	 * When the CG procedure terminates, the close function is invoked. 
	 */
	@Override
	public void close() {

		for (int i = 0; i < PPvertices.length; i++) {
			PPvertices[i].processedLabels = new ArrayList<Label>(dataModel.numArcs);
			PPvertices[i].unprocessedLabels =  new PriorityQueue<Label>(dataModel.numArcs, new Label.SortLabels()); }
		if (!this.pricingProblemInfeasible) for (int i = 0; i < vertices.length; i++) vertices[i].SRCIndices = new ArrayList<>(); 
		this.nodesToProcess = new PriorityQueue<PPVertex>(new SortVertices());
	}

	/**
	 * When the CG procedure terminates, the close function is invoked. 
	 */
	public void restart() {
		for (int i = 0; i < PPvertices.length; i++) {
			PPvertices[i].processedLabels = new ArrayList<Label>(dataModel.numArcs);
			PPvertices[i].unprocessedLabels =  new PriorityQueue<Label>(dataModel.numArcs, new Label.SortLabels());
		}
		this.nodesToProcess = new PriorityQueue<PPVertex>(new SortVertices());
	}

	/**
	 * This method produces zero or more columns. 
	 */
	@Override
	protected List<Route> generateNewColumns() {

		//Solve the problem and check the solution
		List<Route> newRoutes = new ArrayList<>(this.numCols);  			//list of routes
		List<Route> nonElementaryRoutes = new ArrayList<>(this.numCols);  //list of nonelementary routes
		
		this.runLabeling(); 											//runs the labeling algorithm

		if(PPvertices[0].unprocessedLabels.isEmpty()) {
			this.pricingProblemInfeasible=true; this.objective=Double.MAX_VALUE;
		} else {
			this.pricingProblemInfeasible=false;
			for (Label label: PPvertices[0].unprocessedLabels) {
				if (label.reducedCost<=-dataModel.precision) {		//generate new column if it has negative reduced cost
					
					int load = dataModel.Q - label.remainingLoad;
					int energy = dataModel.E-label.remainingEnergy[dataModel.gamma]; double reducedCost = label.reducedCost;
					int departureTime = label.remainingTime;
					
					// Retrieves the charging schedule
					int initialChargingTime = PPvertices[dataModel.PParcs[label.nextArc].head_vertex_id].node_number; 
					int chargingTime = 0;
					PPArc nextArc = dataModel.PParcs[label.nextArc];
					while(nextArc.arc_type>=AC2) {
						chargingTime++;
						
						label = PPvertices[nextArc.head_vertex_id].processedLabels.get(label.nextLabelIndex);
						nextArc = dataModel.PParcs[label.nextArc];
					}

					// Save the (last_t - 0j) Arc
					ArrayList<Integer> PParcs = new ArrayList<Integer>(dataModel.C);
					PParcs.add(nextArc.id);
					
					// Retrieve the route
					HashMap<Integer, Integer> route = new HashMap<Integer, Integer>(dataModel.C);
					ArrayList<Integer> arcs = new ArrayList<Integer>(dataModel.C);
					boolean isElementary = true;
					
					int j = PPvertices[nextArc.head_vertex_id].routing_vertex.node_id;
					route.put(j, 1);
					Arc routing_arc = dataModel.graph.getEdge(0,j); int cost = routing_arc.cost;
					arcs.add(routing_arc.id);
					
					label = PPvertices[nextArc.head_vertex_id].processedLabels.get(label.nextLabelIndex);
					while(label.vertex != depotID) {
						
						int i = j;
						nextArc = dataModel.PParcs[label.nextArc];
						j = PPvertices[nextArc.head_vertex_id].routing_vertex.node_id;
						
						if (route.containsKey(j)) {route.replace(j, route.get(j)+1); isElementary = false; } 
						else route.put(j, 1);
						routing_arc = dataModel.graph.getEdge(i,j); cost += routing_arc.cost;

						arcs.add(routing_arc.id); PParcs.add(nextArc.id);
						label = PPvertices[nextArc.head_vertex_id].processedLabels.get(label.nextLabelIndex);
						
					}

					// Retrieves the route sequence (of customers)
					int[] routeSequence = new int[arcs.size()-1];
					int counter = 0;
					for(Integer arcID: arcs) {
						if(counter>=routeSequence.length) break;
						routeSequence[counter] = dataModel.arcs[arcID].head;
						counter++;
					}

					Route column = new Route("heuristicLabeling", false, route, routeSequence, pricingProblem, cost, departureTime, energy, load, reducedCost, arcs, PParcs, initialChargingTime+chargingTime-1, chargingTime);
					if (isElementary) {newRoutes.add(column);}
					else {nonElementaryRoutes.add(column);}
				}
			}
			
		}

		if (dataModel.print_log) {
				logger.debug("Finished heuristic pricing: "+PPvertices[0].processedLabels.size()+" processed, "+PPvertices[0].unprocessedLabels.size()+" unprocessed.");
				logger.debug("Found " + newRoutes.size() + " columns");
		}
		
		close();
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

		pricingProblem.reducedCostThreshold = 0.0;
		pricingProblem.bestReducedCost = -Double.MAX_VALUE;
		// Update the objective function with the new dual values
		DirectedWeightedMultigraph<Integer, PPArc> PPgraph = dataModel.PPgraph;

		// Routing Arcs
		for (int i = 1; i <= dataModel.C; i++){
			int vertex_id = dataModel.C0_startID+i;
			for (PPArc arc: PPgraph.outgoingEdgesOf(vertex_id)){
				arc.modifiedCost = dataModel.graph.getEdge(0,i).cost + arc.routing_arc.cost - pricingProblem.dualCosts[i-1];
			}

			vertex_id = dataModel.C1_startID+i;
			for (PPArc arc: PPgraph.outgoingEdgesOf(vertex_id)){
				arc.modifiedCost = arc.routing_arc.cost - pricingProblem.dualCosts[i-1];
			}
		}

		// Charging Scheduling Arcs
		for (int t = 1; t <= dataModel.last_charging_period; t++){
			int vertex_id = dataModel.T_startID+t;
			for (PPArc arc: PPgraph.outgoingEdgesOf(vertex_id)){
				arc.modifiedCost = -pricingProblem.dualCosts[dataModel.C+t-1];
			}
		}
	
	}

	public boolean checkDominance(Label newLabel){

		byte arc_type = dataModel.PParcs[newLabel.nextArc].arc_type;
		if (arc_type <= AR1) return checkDominanceRouting(newLabel, arc_type);
		else return checkDominanceCharging(newLabel);
	}

	/**
	 * Verifies if a label is dominated. Returns true if it is, false otherwise.
	 * If the label is dominated it is discarded
	 * If the label is not dominated, the existing labels dominated by the label is discarded
	 * @param label to which check dominance
	 */
	public boolean checkDominanceRouting(Label newLabel, byte depot) {

		/* // DELETE BLOCK LATER
		int[] lookup_route = new int[]{0,12,9,3,20,10,1}; // DELETE LATER
		int[] nl_sequence = get_route_sequence(newLabel); // DELETE LATER
		boolean is_nl_subset = false; // DELETE LATER
		if (nl_sequence.length <= lookup_route.length){
			is_nl_subset = sequence_is_subset(nl_sequence, lookup_route);
		} */

		PPVertex currentVertex = PPvertices[newLabel.vertex];
		ArrayList<Label> labelsToDelete = new ArrayList<Label>();
		for(Label existingLabel: currentVertex.unprocessedLabels) {

			/* // DELETE BLOCK LATER
			boolean existing_is_discarded = false; 
			int[] el_sequence = get_route_sequence(existingLabel); 
			boolean is_el_subset = false;
			if (el_sequence.length <= lookup_route.length){ 
				is_el_subset = sequence_is_subset(el_sequence, lookup_route);
			} */

			if(isDominatedRouting(existingLabel, newLabel, depot)) {
				//existing_is_discarded = true; // DELETE LATER
				labelsToDelete.add(existingLabel);
			}
			
		}
		currentVertex.unprocessedLabels.removeAll(labelsToDelete);
		if(currentVertex.unprocessedLabels.isEmpty()) nodesToProcess.remove(currentVertex);

		//boolean new_is_discarded = false; // DELETE LATER
		for(Label existingLabel: currentVertex.processedLabels) {
			//int[] el_sequence = get_route_sequence(existingLabel); // DELETE LATER
			if(isDominatedRouting(newLabel, existingLabel, depot)) return true;
			
		}

		return false;
	}

	public boolean checkDominanceCharging(Label newLabel) {

		PPVertex currentVertex = PPvertices[newLabel.vertex];
		ArrayList<Label> labelsToDelete = new ArrayList<Label>();
		for(Label existingLabel: currentVertex.unprocessedLabels) { if(isDominatedCharging(existingLabel, newLabel)) labelsToDelete.add(existingLabel); }
		currentVertex.unprocessedLabels.removeAll(labelsToDelete);
		
		if(currentVertex.unprocessedLabels.isEmpty()) nodesToProcess.remove(currentVertex);

		for(Label existingLabel: currentVertex.processedLabels) { if(isDominatedCharging(newLabel, existingLabel)) return true; }

		return false;
	}


	/**
	 * Verifies if L1 is (strongly) dominated by L2
	 * @param L1, L2 labels
	 */
	public boolean isDominatedRouting(Label L1, Label L2, byte depot) {
		
		/* int[] nl_sequence = get_route_sequence(L1); // DELETE LATER
		int[] el_sequence = get_route_sequence(L2); // DELETE LATER */

		if (depot==AR1 && L2.remainingLoad<L1.remainingLoad) return false; 	//load
		if (L2.reducedCost-L1.reducedCost>dataModel.precision) return false; 	//reduced cost
		if (L2.remainingTime<L1.remainingTime) return false; 					//time
		
		for (int gam=0; gam<=dataModel.gamma; gam++){
			if (L2.remainingEnergy[gam]<L1.remainingEnergy[gam]) return false; //energy
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

		return true;
	}

	public boolean isDominatedCharging(Label L1, Label L2) {
		
		if (L2.chargingTime>L1.chargingTime) return false;
		if (L2.reducedCost-L1.reducedCost>dataModel.precision) return false;
		return true;

	}


	public int[] get_route_sequence(Label label) {

		Label new_label = label.clone();

		ArrayList<Integer> arcs = new ArrayList<Integer>(dataModel.C);
		while(new_label.vertex != depotID) {	
			
			PPArc nextArc = dataModel.PParcs[new_label.nextArc];
			int i = PPvertices[nextArc.tail_vertex_id].routing_vertex.node_id;
			int j = PPvertices[nextArc.head_vertex_id].routing_vertex.node_id;
			
			Arc routing_arc = dataModel.graph.getEdge(i,j);

			arcs.add(routing_arc.id);
			new_label = PPvertices[nextArc.head_vertex_id].processedLabels.get(new_label.nextLabelIndex);
			
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
	public class SortVertices implements Comparator<PPVertex> {

		@Override
		public int compare(PPVertex vertex1, PPVertex vertex2) {
			
			// If one of the vertices belongs to C1 and the other to C0, then give priority to the one in C1
			if(vertex2.vertex_type==C0 && vertex1.vertex_type==C1) return -1;
			if(vertex1.vertex_type==C0 && vertex2.vertex_type==C1) return 1;

			
			// If both vertices are charging vertices, give priority to the one with largest t
			if(vertex1.vertex_type==Tt && vertex2.vertex_type==Tt) { 
				if(vertex1.node_number>vertex2.node_number) return -1;
				else return 1;
			}

			// If one of the vertices is a routing vertex (C0 or C1) and the other is a charging vertex, give priority to the routing one
			if(vertex2.vertex_type==Tt) return -1;
			if(vertex1.vertex_type==Tt) return 1;
			
			// If both vertices are C0 or if both vertices are C1, choose according the current unprocessed labels
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
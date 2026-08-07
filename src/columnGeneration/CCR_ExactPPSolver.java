package columnGeneration;

import model.EVRPTW;
import model.EVRPTW.Arc;
import model.EVRPTW.PPArc;
import model.EVRPTW.Vertex;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Set;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;


/**
 * This class provides a heuristic solver for the ng-SPPRC pricing problem
 * It uses a relaxed dominance rule (but an exact dominance rule)
 */
public final class CCR_ExactPPSolver extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem> {

	public Vertex[] vertices = dataModel.vertices; 			//vertices of the instance
	public PriorityQueue<Vertex> nodesToProcess; 			//labels that need be processed

	public double bestReducedCost;
	private int Gamma;
	private int depotID;

	public boolean canTriggerRollback;
	public int nLabels;
	public int rollbackThreshold;

	/**
	 * Labeling algorithm to solve the ng-SPPRC
	 */
	public CCR_ExactPPSolver(EVRPTW dataModel, PricingProblem pricingProblem) {
		super(dataModel, pricingProblem);
		this.name="ExactLabelingSolver"; //Set a name for the solver
		this.nodesToProcess = new PriorityQueue<Vertex>(dataModel.V, new SortVertices());
		this.Gamma = dataModel.gamma;
		this.depotID = dataModel.T_startID;
	}

	/**
	 * Runs the labeling algorithm
	 */
	public void runLabeling() {

		dataModel.rollbackTrigger = false;
		dataModel.exactPricing = true;
		this.bestReducedCost = Double.MAX_VALUE;

		// Initialization
		int[] remain_energy = new int[Gamma + 1]; Arrays.fill( remain_energy, dataModel.E);
		CCRLabel initialLabel = new CCRLabel(0, -pricingProblem.dualCost, dataModel.Q, vertices[dataModel.C+1].closing_tw, remain_energy, 0,new boolean[dataModel.C], new boolean[dataModel.C], new boolean[pricingProblem.subsetRowCuts.size()], new HashSet<Integer>(pricingProblem.subsetRowCuts.size()));
		initialLabel.index = 0; initialLabel.vertex = dataModel.C+1; initialLabel.nextArc = 0;
		this.nodesToProcess.add(vertices[dataModel.C+1]);
		vertices[dataModel.C+1].unprocessedLabels.add(initialLabel);
		
		////////////////////////////////////////////
		/// Routing Labeling
		////////////////////////////////////////////
		
		this.nLabels = 0;
		long startTime = System.currentTimeMillis();
		while (!nodesToProcess.isEmpty() && System.currentTimeMillis()<timeLimit && (!canTriggerRollback || nLabels < rollbackThreshold)) {
			ArrayList<CCRLabel> labelsToProcessNext = routingLabelsToProcessNext();
			Set<Arc> incomingArcs = new HashSet<Arc>(dataModel.graph.incomingEdgesOf(labelsToProcessNext.get(0).vertex));
			
			for (CCRLabel currentLabel: labelsToProcessNext) {
				
				boolean isDominated = checkRoutingDominance(currentLabel);
				if(isDominated) continue;
				
				currentLabel.index = vertices[currentLabel.vertex].processedLabels.size();
				vertices[currentLabel.vertex].processedLabels.add(currentLabel);
				
				for(Arc a: incomingArcs) extendLabel(currentLabel, a);

			}
		}

		////////////////////////////////////////////
		/// Depot Labels
		////////////////////////////////////////////
		
		if (System.currentTimeMillis()>=timeLimit) vertices[0].unprocessedLabels.clear();
		
		while (!vertices[0].unprocessedLabels.isEmpty()){

			CCRLabel currentLabel = vertices[0].unprocessedLabels.poll();
			boolean isDominated = checkDepotDominance(currentLabel);
			if(isDominated) continue;
			
			currentLabel.index = vertices[currentLabel.vertex].processedLabels.size();
			vertices[currentLabel.vertex].processedLabels.add(currentLabel);
			vertices[0].processedLabels.add(currentLabel);
			
			for (Arc a: dataModel.graph.incomingEdgesOf(currentLabel.vertex)) extendLabelChargingTime(currentLabel, a);

		}

		/////////////////////////////////////
		/// Charging Scheduling Labeling
		////////////////////////////////////
		
		while (!nodesToProcess.isEmpty() && System.currentTimeMillis()<timeLimit) {
			ArrayList<CCRLabel> labelsToProcessNext = chargingLabelsToProcessNext();
			
			for (CCRLabel currentLabel: labelsToProcessNext) {
				
				boolean isDominated = checkChargingDominance(currentLabel);
				if(isDominated) continue;
				
				currentLabel.index = vertices[currentLabel.vertex].processedLabels.size();
				vertices[currentLabel.vertex].processedLabels.add(currentLabel);
				
				for (Arc a: dataModel.graph.incomingEdgesOf(currentLabel.vertex)) extendLabelChargingTime(currentLabel, a);
			}
		}

		long totalTime = System.currentTimeMillis()-startTime;
		dataModel.heuristicPricingTime+=totalTime;
		if (dataModel.print_log) logger.debug("Time solving (exactly) the pricing problem (s): " + getTimeInSeconds(totalTime)); 
	}

	/**
	 * Selects a set of labels to process (the one with the most remaining load)
	 */
	public ArrayList<CCRLabel> routingLabelsToProcessNext(){

		ArrayList<CCRLabel> labelsToProcessNext = new ArrayList<CCRLabel>();
		Vertex currentVertex = nodesToProcess.poll();
		
		while(true) {
			CCRLabel currentLabel = currentVertex.unprocessedLabels.poll();
			if(labelsToProcessNext.isEmpty()) labelsToProcessNext.add(currentLabel);
			else {
				boolean isDominated = false;
				for(CCRLabel L2: labelsToProcessNext) {
					
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

	public ArrayList<CCRLabel> chargingLabelsToProcessNext(){

		ArrayList<CCRLabel> labelsToProcessNext = new ArrayList<CCRLabel>();
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
	 * CCRLabel extension procedure
	 */
	public CCRLabel extendLabel(CCRLabel currentLabel, Arc arc) {

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
		boolean is_energy_feasible = update_worst_case_energy_resource(remainingEnergy, currentLabel.remainingEnergy, arc);
		if (!is_energy_feasible) return null;
		
		int chargingTime = dataModel.f_inverse[dataModel.E-remainingEnergy[dataModel.gamma]];

		//Quick check
		if(source>0 && remainingTime-dataModel.graph.getEdge(0, source).time<vertices[0].opening_tw) return null;
		//Check whether the extension is actually feasible
		if(remainingTime<vertices[source].opening_tw || chargingTime>= (int) (remainingTime/10)) return null;

		// Unreachable resources
		boolean[] unreachable = null;
		CCRLabel extendedLabel = null;

		//Mark unreachable customers and ng-path cycling restrictions
		if(source > 0) {

			// Mark unreachable customers and ng-path cycling restrictions
			unreachable = Arrays.copyOf(currentLabel.unreachable.clone(), currentLabel.unreachable.length);

			unreachable[source-1] = true;
			for (Arc c: dataModel.graph.incomingEdgesOf(source)) {
				if(c.tail==0 || unreachable[c.tail-1]) continue;
				//unreachable
				if (remainingLoad-vertices[c.tail].load<0 || remainingTime-c.min_time<vertices[c.tail].opening_tw || remainingEnergy[dataModel.gamma]-c.min_energy < 0 || 
					Math.min(remainingTime-c.min_time, vertices[c.tail].closing_tw)-dataModel.graph.getEdge(0, c.tail).min_time<vertices[0].opening_tw ||
					remainingEnergy[dataModel.gamma]-c.min_energy - dataModel.graph.getEdge(0, c.tail).min_energy<0) {
					unreachable[c.tail-1] = true;
				}
			}

			extendedLabel = new CCRLabel(currentLabel.index, reducedCost, remainingLoad, remainingTime, remainingEnergy, chargingTime,unreachable, currentLabel.ng_path, eta, srcIndices);
			extendedLabel.vertex = source;
			extendedLabel.nextArc = arc.id;
			vertices[extendedLabel.vertex].unprocessedLabels.add(extendedLabel);
			if (vertices[extendedLabel.vertex].unprocessedLabels.size() == 1) nodesToProcess.add(vertices[extendedLabel.vertex]);

		} else {
			
			// We re-scale the remaining time, which represents the departure time of the route
			remainingTime = (int) (remainingTime/10);

			////////////////////////////////////////////
			/// Bounding Procedure
			////////////////////////////////////////////

			extendedLabel = new CCRLabel(currentLabel.index, reducedCost, remainingLoad, remainingTime, remainingEnergy, chargingTime,unreachable, currentLabel.ng_path, eta, srcIndices);
			extendedLabel.vertex = source;
			extendedLabel.nextArc = arc.id;

			pricingProblem.frcCCRLabels.add(extendedLabel);
			double min_col_rc = reducedCost + pricingProblem.charging_bounds.get(chargingTime).get(remainingTime);
			if (min_col_rc < this.bestReducedCost - dataModel.precision) this.bestReducedCost = min_col_rc;
			if (min_col_rc >= -dataModel.precision) return null;

			vertices[0].unprocessedLabels.add(extendedLabel);
		
		}

		return extendedLabel;

	}

	/**
	 * CCRLabel extension procedure
	 */
	public CCRLabel extendLabelChargingTime(CCRLabel currentLabel, Arc arc) {

		int source = arc.tail;

		if(arc.head == 0 && (source - dataModel.V < currentLabel.chargingTime || source - dataModel.V >= currentLabel.remainingTime)) return null;

		double reducedCost = currentLabel.reducedCost + arc.modifiedCost;
		reducedCost = Math.floor(reducedCost*10000)/10000;
		if (reducedCost >= -dataModel.precision) return null; // Only negative reduced costs labels will get to the source node
		
		int chargingTime = currentLabel.chargingTime;
		CCRLabel extendedLabel = null;
		if(source != dataModel.V) {
			chargingTime -= 1;
			if (chargingTime < 0) return null;

			extendedLabel = new CCRLabel(currentLabel.index, reducedCost, currentLabel.remainingLoad, currentLabel.remainingTime, currentLabel.remainingEnergy, chargingTime , currentLabel.unreachable, currentLabel.ng_path, currentLabel.eta, currentLabel.srcIndices);
			vertices[source].unprocessedLabels.add(extendedLabel);
			if (vertices[source].unprocessedLabels.size() == 1) nodesToProcess.add(vertices[source]);

		} else {
			if (chargingTime > 0) return null;

			extendedLabel = new CCRLabel(currentLabel.index, reducedCost, currentLabel.remainingLoad, currentLabel.remainingTime, currentLabel.remainingEnergy, chargingTime , currentLabel.unreachable, currentLabel.ng_path, currentLabel.eta, currentLabel.srcIndices);
			vertices[source].unprocessedLabels.add(extendedLabel);
		
		}

		extendedLabel.vertex = source;
		extendedLabel.nextArc = arc.id;

		return extendedLabel;
	}

	private boolean update_worst_case_energy_resource(int[] remainingEnergy, int[] currentEnergy, Arc routing_arc){

		int gamma_change = Gamma + 1;
		for (int gam = 1; gam <= Gamma; gam++)
			if (currentEnergy[gam-1] - routing_arc.energy_deviation < currentEnergy[gam]) { gamma_change = gam; break; }
		if (gamma_change <= Gamma){
			remainingEnergy[Gamma] = currentEnergy[Gamma-1] - routing_arc.energy - routing_arc.energy_deviation;
			if (remainingEnergy[Gamma] < 0) return false; }
		for (int gam = Gamma-1; gam >= gamma_change; gam--)
			remainingEnergy[gam] = currentEnergy[gam-1] - routing_arc.energy - routing_arc.energy_deviation;
		for (int gam = 0; gam < gamma_change; gam ++){
			remainingEnergy[gam] = currentEnergy[gam] - routing_arc.energy;
			if (remainingEnergy[gam] < 0) return false; }

		return true;
	}


	/**
	 * When the CG procedure terminates, the close function is invoked. 
	 */
	@Override
	public void close() {
		
		// Save information of the routing subgraph for Fixing by Reduced Cost
		for (int i = 1; i <= dataModel.C; i++) pricingProblem.SRCIndices.add(new ArrayList<>(vertices[i].SRCIndices));

		for (int i = 0; i < vertices.length; i++) {
			vertices[i].processedLabels = new ArrayList<CCRLabel>(dataModel.numArcs);
			vertices[i].unprocessedLabels =  new PriorityQueue<CCRLabel>(dataModel.numArcs, new CCRLabel.SortLabels(dataModel.V)); }

		if (!this.pricingProblemInfeasible) for (int i = 0; i < vertices.length; i++) vertices[i].SRCIndices = new ArrayList<>(); 
		this.nodesToProcess = new PriorityQueue<Vertex>(new SortVertices());
	}

	/**
	 * When the CG procedure terminates, the close function is invoked. 
	 */
	public void restart() {
		for (int i = 0; i < vertices.length; i++) {
			vertices[i].processedLabels = new ArrayList<CCRLabel>(dataModel.numArcs);
			vertices[i].unprocessedLabels =  new PriorityQueue<CCRLabel>(dataModel.numArcs, new CCRLabel.SortLabels(dataModel.V));
		}
		this.nodesToProcess = new PriorityQueue<Vertex>(new SortVertices());
		pricingProblem.frcCCRLabels.clear();
	}

	/**
	 * This method produces zero or more columns. 
	 */
	@Override
	protected List<Route> generateNewColumns() {

		this.canTriggerRollback = dataModel.cut_iterations > 1;
		this.rollbackThreshold = dataModel.rollbackBaseLine*dataModel.rollbackFactor;

		//Solve the problem and check the solution
		boolean existsElementaryRoute=false;
		boolean maxNeighborhoodSize=false;
		List<Route> newRoutes = new ArrayList<>();  			//list of routes
		List<Route> nonElementaryRoutes = new ArrayList<>();  //list of nonelementary routes

		while (!existsElementaryRoute && !maxNeighborhoodSize){
			this.runLabeling();

			dataModel.rollbackExplosion = nLabels;
			if (canTriggerRollback && this.nLabels >= this.rollbackThreshold) { // If the rollback is triggered, return an empty list of columns
				dataModel.rollbackTrigger = true;
				dataModel.exactPricing = false;
				this.close();
				return new ArrayList<Route>(); }

			if(vertices[dataModel.V].unprocessedLabels.isEmpty()) {
				existsElementaryRoute = true; pricingProblemInfeasible=false; this.objective=this.bestReducedCost;
				
			} else {
				this.pricingProblemInfeasible=false;
				for (CCRLabel label: vertices[dataModel.V].unprocessedLabels) {
					int departureTime = label.remainingTime;
					int load = dataModel.Q - label.remainingLoad;
					if (label.reducedCost<=-dataModel.precision) {		//generate new column if it has negative reduced cost
						boolean isElementary = true;
						double reducedCost = label.reducedCost; int energy = dataModel.E-label.remainingEnergy[dataModel.gamma];
					
					HashMap<Integer, Integer> route = new HashMap<Integer, Integer>(dataModel.C);
					ArrayList<Integer> arcs = new ArrayList<Integer>(dataModel.C);
					int initialChargingTime = dataModel.arcs[label.nextArc].head-dataModel.V; int chargingTime = 0;
					int i = label.vertex;
					while(i != 0) {
						Arc currentArc = dataModel.arcs[label.nextArc];
						int j = currentArc.head;
						if (i != dataModel.V) chargingTime++;

						label = vertices[j].processedLabels.get(label.nextLabelIndex); i = j;
					}

					int last_t = initialChargingTime + chargingTime - 1;
					ArrayList<Integer> PParcs = new ArrayList<Integer>();

					// Retrieve first customer visited in the route, and the corresponding EC2FC arc
					int cost = 0; 
					Arc currentArc = dataModel.arcs[label.nextArc];
					cost += currentArc.cost; int j = currentArc.head;
					PPArc currentPPArc = dataModel.PPgraph.getEdge(dataModel.T_startID+last_t, dataModel.C0_startID+j);
					arcs.add(currentArc.id); PParcs.add(currentPPArc.id);
					label = vertices[j].processedLabels.get(label.nextLabelIndex); i = j;

					// Retrieve second customer visited in the route, and the corresponding AR0 arc
					currentArc = dataModel.arcs[label.nextArc];
					cost += currentArc.cost; j = currentArc.head;
					int j_PPix = dataModel.C1_startID+j; if (j == dataModel.C+1) j_PPix = depotID;
					currentPPArc = dataModel.PPgraph.getEdge(dataModel.C0_startID+i, j_PPix);
					route.put(i, 1); arcs.add(currentArc.id); PParcs.add(currentPPArc.id);
					label = vertices[j].processedLabels.get(label.nextLabelIndex); i = j;

					while(i != dataModel.C+1) {
						currentArc = dataModel.arcs[label.nextArc];
						cost += currentArc.cost; j = currentArc.head;
						j_PPix = dataModel.C1_startID+j; if (j == dataModel.C+1) j_PPix = depotID;
						currentPPArc = dataModel.PPgraph.getEdge(dataModel.C1_startID+i, j_PPix);
						if(route.containsKey(i)) {route.replace(i, route.get(i)+1); isElementary = false; } 
						else route.put(i, 1);
						arcs.add(currentArc.id); PParcs.add(currentPPArc.id);

						label = vertices[j].processedLabels.get(label.nextLabelIndex); i = j;
					}

					//Gets the route sequence (of customers)
					int[] routeSequence = new int[arcs.size()-1];
					int counter = 0;
					for(Integer arc: arcs) {
						if(counter>=routeSequence.length) break;
						routeSequence[counter] = dataModel.arcs[arc].head;
						counter++;
					}
					Route column = new Route("exactLabeling", false, route, routeSequence, pricingProblem, cost, departureTime, energy, load, reducedCost, arcs, PParcs, initialChargingTime+chargingTime-1, chargingTime);
					
					if (isElementary) {existsElementaryRoute = true; newRoutes.add(column);}
						else {nonElementaryRoutes.add(column);}
					}
				}
				if (!existsElementaryRoute) {
					maxNeighborhoodSize = !enlargeNeighborhoods(nonElementaryRoutes); 
					if(!maxNeighborhoodSize) {
						nonElementaryRoutes = new ArrayList<Route>();newRoutes=new ArrayList<>();
						restart();} //restart //run again
					else {newRoutes = nonElementaryRoutes; existsElementaryRoute = true;}
				}
			}
		}

		pricingProblem.bestReducedCost = this.bestReducedCost;

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
		//Already done by the heuristic labeling (must be invoked first)
	}

	public boolean checkRoutingDominance(CCRLabel newLabel) {
		
		Vertex currentVertex = vertices[newLabel.vertex];

		ArrayList<CCRLabel> labelsToDelete = new ArrayList<CCRLabel>();
		for(CCRLabel existingLabel: currentVertex.unprocessedLabels) {

			if(isDominatedRouting(existingLabel, newLabel)) {
				labelsToDelete.add(existingLabel);
			}
			
		}
		currentVertex.unprocessedLabels.removeAll(labelsToDelete);
		if(currentVertex.unprocessedLabels.isEmpty()) nodesToProcess.remove(currentVertex);

		for(CCRLabel existingLabel: currentVertex.processedLabels) {
			if(isDominatedRouting(newLabel, existingLabel)) return true;
		}

		return false;
	}

	public boolean checkDepotDominance(CCRLabel newLabel) {
		
		Vertex currentVertex = vertices[0];

		ArrayList<CCRLabel> labelsToDelete = new ArrayList<CCRLabel>();
		for(CCRLabel existingLabel: currentVertex.unprocessedLabels)  if(isDominatedDepot(existingLabel, newLabel))  labelsToDelete.add(existingLabel);
		currentVertex.unprocessedLabels.removeAll(labelsToDelete);
		if(currentVertex.unprocessedLabels.isEmpty()) nodesToProcess.remove(currentVertex);

		for(CCRLabel existingLabel: currentVertex.processedLabels) if(isDominatedDepot(newLabel, existingLabel)) return true;

		return false;
	}

	public boolean checkChargingDominance(CCRLabel newLabel) {
		
		Vertex currentVertex = vertices[newLabel.vertex];

		ArrayList<CCRLabel> labelsToDelete = new ArrayList<CCRLabel>();
		for(CCRLabel existingLabel: currentVertex.unprocessedLabels) if(isDominatedCharging(existingLabel, newLabel)) labelsToDelete.add(existingLabel);
		currentVertex.unprocessedLabels.removeAll(labelsToDelete);
		if(currentVertex.unprocessedLabels.isEmpty()) nodesToProcess.remove(currentVertex);

		for(CCRLabel existingLabel: currentVertex.processedLabels) if(isDominatedCharging(newLabel, existingLabel)) return true;

		return false;
	}

	/**
	 * Verifies if L1 is (strongly) dominated by L2
	 * @param L1, L2 labels
	 */
	public boolean isDominatedDepot(CCRLabel L1, CCRLabel L2) {
		
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
	public boolean isDominatedRouting(CCRLabel L1, CCRLabel L2) {
		
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

		// Ng-paths and unreachable resources
		Vertex currentVertex = vertices[L1.vertex];
		for(int i: currentVertex.neighbors) {
			
			//boolean check_binaries = (L2.ng_path[i-1] || L2.unreachable[i-1]) && !(L1.ng_path[i-1] || L1.unreachable[i-1]);
			boolean other_way = L2.ng_path[i-1] && (!L1.unreachable[i-1] && !L1.ng_path[i-1]); // Dani's way
			if (other_way) {
				return false;
			}
		}

		return true;
	}

	public boolean isDominatedCharging(CCRLabel L1, CCRLabel L2) {
		
		if (L2.chargingTime>L1.chargingTime) return false;
		if (L2.reducedCost-L1.reducedCost>dataModel.precision) return false;
		return true;

	}

	public int[] get_route_sequence(CCRLabel label) {

		CCRLabel new_label = label.clone();

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
		
	}

	/**
	 * When the Branch-and-Price algorithm backtracks, branching decisions are reversed.
	 * @param bd BranchingDecision
	 */
	@Override
	public void branchingDecisionReversed(BranchingDecision bd) {
		
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

			if(vertex2.node_id==0 && (vertex1.node_id>0 && vertex1.node_id<=dataModel.C)) return -1;
			if(vertex1.node_id==0 && (vertex2.node_id>0 && vertex2.node_id<=dataModel.C)) return 1;

			if(vertex1.node_id<dataModel.V && vertex2.node_id>=dataModel.V) return -1;
			if(vertex1.node_id>=dataModel.V && vertex2.node_id<dataModel.V) return 1;
			if(vertex1.node_id>=dataModel.V && vertex2.node_id>=dataModel.V) {
				if(vertex1.node_id>vertex2.node_id) return -1;
				else return 1;
			}

			CCRLabel L1 = vertex1.unprocessedLabels.peek();
			CCRLabel L2 = vertex2.unprocessedLabels.peek();
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
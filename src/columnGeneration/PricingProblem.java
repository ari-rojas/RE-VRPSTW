package columnGeneration;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Comparator;
import java.util.Map;
import java.util.HashMap;
import java.util.HashSet;

import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblem;
import model.EVRPTW;
import model.EVRPTW.Arc;
import model.EVRPTW.PPArc;
import model.EVRPTW.Vertex;
import model.EVRPTW.PPVertex;

/**
 * This class defines the pricing problem. 
 * We simply extend the pricing problem included in the framework (there is no need for modification, only one pricing problem)
 */
public final class PricingProblem extends AbstractPricingProblem<EVRPTW> {

	public ArrayList<SubsetRowInequality> subsetRowCuts; 				//subset row cuts considered
	public double bestReducedCost = -Double.MAX_VALUE; 					//best reduced cost found by the exact labeling
	public double reducedCostThreshold = 0; 							//minimum reduced cost when arriving at the depot source

	public Map<Integer, Map<Integer,Double>> charging_bounds;

	// Information for Fixing by Reduced Costs procedure
	public ArrayList<ArrayList<Label>> bwLabels = new ArrayList<>();
	public ArrayList<ArrayList<Integer>> SRCIndices = new ArrayList<>();
	public int[] infeasiblePPArcs;
	public PPVertex[] PPvertices;
	public Vertex[] vertices;

	// General information
	private int Gamma = dataModel.gamma;
	private int depotID = dataModel.T_startID;

	// Identifiers for the differnt types of vertices and arcs in the Pricing Problem Routing SubGraph
	public static final byte C0 = EVRPTW.C0; 	 	// Customer depot nodes
	public static final byte C1 = EVRPTW.C1; 		// Non-first customer nodes
	public static final byte AR0 = EVRPTW.AR0;		// Routing arcs between customer depot nodes and non-first customer nodes
	public static final byte AR1 = EVRPTW.AR1;		// Routing arcs between non-first customer nodes

	public PricingProblem(EVRPTW modelData, String name) {
		super(modelData, name);
	}

	public void compute_charging_bounds(){

		// 1. Get the minimum and maximum possible departure times
		int minT = (int) (dataModel.vertices[0].opening_tw/10);
		int maxT = dataModel.last_charging_period + 1;

		///////////////////////////////////////////////////////////////////////////
		/// Preprocess the dual information
		///////////////////////////////////////////////////////////////////////////

		// Precompute fixed sums of the charging capacity dual variables
		double[] S = new double[maxT]; S[0] = 0.0;
		for (int t = 1; t < maxT; t++) {
			double dual = this.dualCosts[dataModel.C + t - 1];
			S[t] = S[t - 1] + dual;
		}

		///////////////////////////////////////////////////////////////////////////////////
		/// Compute the bounds for every combination of chargingTime b and departureTime d
		///////////////////////////////////////////////////////////////////////////////////

		this.charging_bounds = new HashMap<>();
		for (int b = 1; b <= dataModel.f_inverse[dataModel.E]; b++){

			double min_rc = Double.MAX_VALUE;
			int first_departure = Math.max(b+1, minT);
			for (int last_t = b; last_t < first_departure-1; last_t ++){
				double rc = - (S[last_t] - S[last_t-b]);
				if (rc < min_rc - dataModel.precision) min_rc = rc;
			}

			Map<Integer, Double> boundsMap = new HashMap<>();
			for (int d = first_departure; d <= maxT; d++){
				
				double rc = - (S[d-1] - S[d-b-1]);
				if (rc < min_rc - dataModel.precision) min_rc = rc;
				
				boundsMap.put(d, min_rc);
			}
			
			this.charging_bounds.put(b, boundsMap);
			
		}

	}

	public Map<Integer, Double> fixByReducedCosts(long timeLimit, double UB, double LB){
		
		double FRC_gap = UB - LB;
		this.PPvertices = dataModel.PPvertices;
		this.vertices = dataModel.vertices;

		cleanBackwardLabels();
		Map<Integer, Double> arcsToRemove = new HashMap<Integer, Double>();
		
		long startTime = System.currentTimeMillis();

		// j \in C1 Vertices: Non-first customer visits nodes
		// Last j corresponds to the Returning Depot vertex
		for (int j = 1; j <= dataModel.C+1; j++){

			ArrayList<Label> backwardLabels = this.bwLabels.get(j);
			for (PPArc arc: dataModel.PPgraph.incomingEdgesOf(dataModel.C1_startID+j)){

				if (System.currentTimeMillis() > timeLimit) break;
				if (infeasiblePPArcs[arc.id] > 0) continue;

				double min_rc = Double.POSITIVE_INFINITY;
				if (arc.arc_type == AR0) min_rc = computeArcRCProxy_Depot(arc, backwardLabels, this.bwLabels.get(0));
				else {
					int i = PPvertices[arc.tail_vertex_id].node_number;
					min_rc = computeArcRCProxy_Routing(arc, backwardLabels, this.bwLabels.get(i));
				}

				if (!Double.isInfinite(min_rc) && min_rc - bestReducedCost > FRC_gap + dataModel.precision) {
					arcsToRemove.put(arc.id, min_rc); this.infeasiblePPArcs[arc.id] ++; }

			}
		}

		long totalTime = System.currentTimeMillis()-startTime;
		dataModel.exactPricingTime+=totalTime;
		if (dataModel.print_log) {
			logger.debug("Total time on Variable Fixing by Reduced Cost: " + getTimeInSeconds(totalTime));
		}

		return arcsToRemove;

	}

	public double computeArcRCProxy_Routing(PPArc arc, ArrayList<Label> labels_head, ArrayList<Label> labels_source) {

		Arc routing_arc = arc.routing_arc;
		double mod_cost = arc.modifiedCost;

		double bestCandidate = Double.POSITIVE_INFINITY;
		int max_q = labels_source.size();

		// For each backward label at j
		for (Label currentLabel : labels_head) {

			Label extendedLabel = null;
			extendedLabel = extendLabel(currentLabel, routing_arc, AR1, mod_cost);

			if (extendedLabel == null) continue; // Skip if the extension if not feasible

			// ---------- Compute reduced-cost proxy ----------
			
			double newReducedCost = currentLabel.reducedCost;
			BitSet ng_reachable = new BitSet();
			for (int i = 0; i < dataModel.C; i++) if (!extendedLabel.ng_path[i]) ng_reachable.set(i);
			
			boolean foundDominatingSet = false;
			for (int q = 0; q < max_q; q++) {

				Label l2 = labels_source.get(q);
				if (l2.reducedCost-extendedLabel.reducedCost>dataModel.precision) break; // Early stop when l2.reducedCost is worse than that of extendedLabel
				
				boolean dominated = isDominatedRoutingResources(extendedLabel, l2);
				if (!dominated) continue; // Skip source labels that don't dominate the extended label

				for (int c = ng_reachable.nextSetBit(0); c >= 0; c = ng_reachable.nextSetBit(c + 1)) {
					if (!l2.ng_path[c])  ng_reachable.clear(c); }

				if (ng_reachable.isEmpty()) {
					newReducedCost -= l2.reducedCost;
					foundDominatingSet = true; break; }
			}

			if (!foundDominatingSet) continue; // Skip head labels for which there is no dominating set

			// ---------- 4. Update best bound ----------
			if (newReducedCost < bestCandidate - dataModel.precision)
				bestCandidate = newReducedCost;
		}

		return bestCandidate;
	}

	public double computeArcRCProxy_Depot(PPArc arc, ArrayList<Label> labels_head, ArrayList<Label> labels_source) {

		Arc routing_arc = arc.routing_arc;
		double mod_cost = arc.modifiedCost;

		double bestCandidate = Double.POSITIVE_INFINITY;
		int max_q = labels_source.size();

		// For each backward label at j
		for (Label currentLabel : labels_head) {

			Label extendedLabel = null;
			extendedLabel = extendLabel(currentLabel, routing_arc, AR0, mod_cost);

			if (extendedLabel == null) continue; // Skip if the extension if not feasible

			// ---------- Compute reduced-cost proxy ----------
			
			double newReducedCost = currentLabel.reducedCost;
			
			boolean foundDominatingSet = false;
			for (int q = 0; q < max_q; q++) {

				Label l2 = labels_source.get(q);
				if (l2.reducedCost-extendedLabel.reducedCost>dataModel.precision) break; // Early stop when l2.reducedCost is worse than that of extendedLabel
				
				boolean dominated = isDominatedDepot(extendedLabel, l2);
				if (!dominated) continue; // Skip source labels that don't dominate the extended label
				else {
					newReducedCost -= l2.reducedCost;
					foundDominatingSet = true; break; }
			}

			if (!foundDominatingSet) continue; // Skip head labels for which there is no dominating set

			// ---------- 4. Update best bound ----------
			if (newReducedCost < bestCandidate - dataModel.precision)
				bestCandidate = newReducedCost;
		}

		return bestCandidate;
	}

	private void cleanBackwardLabels() {

		// (Super) Depot labels
		ArrayList<Label> labels = this.bwLabels.get(0);
		ArrayList<Label> labels_to_remove = new ArrayList<Label>();
		for (int ix = 0; ix < labels.size(); ix++){
			Label l1 = labels.get(ix); boolean dominated = false;
			for (int ix2 = ix+1; ix2 < labels.size(); ix2++) if (isDominatedDepot(l1, labels.get(ix2))) { dominated = true; break; }
			if (dominated) labels_to_remove.add(l1);
		} labels.removeAll(labels_to_remove);
		labels.sort(Comparator.comparing(l -> l.reducedCost));

		// C1 PP vertices labels
		for (int i = 1; i <= dataModel.C; i++){
			labels = this.bwLabels.get(i); labels_to_remove = new ArrayList<Label>();
			for (int ix = 0; ix < labels.size(); ix++){
				Label l1 = labels.get(ix); boolean dominated = false;
				for (int ix2 = ix+1; ix2 < labels.size(); ix2++) if (isDominatedRouting(l1, labels.get(ix2))) { dominated = true; break; }
				if (dominated) labels_to_remove.add(l1);
			} labels.removeAll(labels_to_remove);
			labels.sort(Comparator.comparing(l -> l.reducedCost));
		}

	}

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
				reducedCost-=this.dualCosts[dualIndex];
				srcIndices.remove(srcIndex);
			}
			else {eta[srcIndex]=true; srcIndices.add(srcIndex);}
		}
		reducedCost = Math.floor(reducedCost*10000)/10000;
		
		int[] remainingEnergy = new int[Gamma + 1];
		boolean is_energy_feasible = update_worst_case_energy_resource(remainingEnergy, currentLabel.remainingEnergy, routing_arc);
		if (!is_energy_feasible) return null;
		
		// If the arc connects to a vertex i0 \in C0, then the time and energy resources must account for the depot-i arc
		if (arc_type == AR0){
			Arc depotArc = dataModel.graph.getEdge(0, source);
			remainingTime -= depotArc.time;
			
			is_energy_feasible = update_worst_case_energy_resource(remainingEnergy, remainingEnergy, depotArc);
			if (!is_energy_feasible) return null;
		}
		
		// Update charging time and check if it's feasible
		int chargingTime = dataModel.f_inverse[dataModel.E-remainingEnergy[Gamma]];
		if (chargingTime >= (int) (remainingTime/10)) return null;
		
		// After confirming that the label is feasible, update the remaining load
		int remainingLoad = currentLabel.remainingLoad-vertices[source].load;

		// Unreachable resources
		boolean[] unreachable = Arrays.copyOf(currentLabel.unreachable.clone(), currentLabel.unreachable.length);
		boolean[] ng_path = new boolean[dataModel.C];
		
		// Mark unreachable customers and ng-path cycling restrictions
		if(arc_type == AR1) {
			ng_path[source-1] = true;
			for (Arc c: dataModel.graph.incomingEdgesOf(source)) {
				if(c.tail==0 || unreachable[c.tail-1]) continue;
				//unreachable
				if (remainingLoad-vertices[c.tail].load<0 || remainingTime-c.min_time<vertices[c.tail].opening_tw || remainingEnergy[Gamma] - c.min_energy - dataModel.graph.getEdge(0, c.tail).min_energy < 0) {
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

	public boolean isDominatedDepot(Label L1, Label L2) {

		if (L2.reducedCost-L1.reducedCost>dataModel.precision) return false; 	//reduced cost
		if (L2.remainingTime<L1.remainingTime) return false; 					//departure time
		if (L2.chargingTime>L1.chargingTime) return false;						//charging time

		return true;
	}

	public boolean isDominatedRouting(Label L1, Label L2) {

		if (L2.remainingLoad<L1.remainingLoad) return false; 					//load
		if (L2.reducedCost-L1.reducedCost>dataModel.precision) return false; 	//reduced cost
		if (L2.remainingTime<L1.remainingTime) return false; 					//time
		
		for (int gam=0; gam<=Gamma; gam++){
			if (L2.remainingEnergy[gam]<L1.remainingEnergy[gam]) return false; 	//energy
		}
		
		//reducedCost
		double reducedCostL2 = 0;
		for(int i: L2.srcIndices) {
			if(!L1.eta[i]) {
				SubsetRowInequality src = this.subsetRowCuts.get(i);
				if(!L2.unreachable[src.cutSet[0]-1] || !L2.unreachable[src.cutSet[1]-1] || !L2.unreachable[src.cutSet[2]-1]) {
					int dualIndex = dataModel.C+dataModel.last_charging_period+i;
					reducedCostL2+=this.dualCosts[dualIndex];
				}
			}
			if (L2.reducedCost-reducedCostL2-L1.reducedCost>dataModel.precision) return false;
		}

		if (L2.reducedCost-reducedCostL2-L1.reducedCost>dataModel.precision) return false;

		// Ng-paths and unreachable resources
		Vertex currentVertex = PPvertices[L1.vertex].routing_vertex;
		for(int i: vertices[currentVertex.node_id].neighbors) {
			
			//boolean check_binaries = (L2.ng_path[i-1] || L2.unreachable[i-1]) && !(L1.ng_path[i-1] || L1.unreachable[i-1]);
			boolean other_way = L2.ng_path[i-1] && (!L1.unreachable[i-1] && !L1.ng_path[i-1]); // Dani's way
			if (other_way)  return false;
		}

		return true;
	}

	public boolean isDominatedRoutingResources(Label L1, Label L2) {

		if (L2.remainingLoad<L1.remainingLoad) return false; 					//load
		if (L2.reducedCost-L1.reducedCost>dataModel.precision) return false; 	//reduced cost
		if (L2.remainingTime<L1.remainingTime) return false; 					//time
		
		for (int gam=0; gam<=Gamma; gam++){
			if (L2.remainingEnergy[gam]<L1.remainingEnergy[gam]) return false; 	//energy
		}
		
		//reducedCost
		double reducedCostL2 = 0;
		for(int i: L2.srcIndices) {
			if(!L1.eta[i]) {
				SubsetRowInequality src = this.subsetRowCuts.get(i);
				if(!L2.unreachable[src.cutSet[0]-1] || !L2.unreachable[src.cutSet[1]-1] || !L2.unreachable[src.cutSet[2]-1]) {
					int dualIndex = dataModel.C+dataModel.last_charging_period+i;
					reducedCostL2+=this.dualCosts[dualIndex];
				}
			}
			if (L2.reducedCost-reducedCostL2-L1.reducedCost>dataModel.precision) return false;
		}

		if (L2.reducedCost-reducedCostL2-L1.reducedCost>dataModel.precision) return false;

		// Unreachable resources
		int vertex = PPvertices[L1.vertex].routing_vertex.node_id;
		for(Arc arc: dataModel.graph.incomingEdgesOf(vertex)) {
			int i = arc.tail;
			if (!L1.unreachable[i-1] && L2.unreachable[i-1])  return false;
		}

		return true;
	}

	public double getTimeInSeconds(double time) {
		double realTime = time*0.001;
		realTime = Math.floor(realTime*100)/100; //two decimals
		return realTime;
	}

}
package columnGeneration;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.Map;
import java.util.HashMap;
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
	public int[] infeasibleArcs;
	public PPVertex[] PPvertices;
	public Vertex[] vertices;

	// General information
	private int Gamma = dataModel.gamma;

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

		Map<Integer, Double> arcsToRemove = new HashMap<Integer, Double>();
		ArrayList<ArrayList<PartialBackwardSequence>> bwSequences = getBackwardSequences(); // get *sorted* backward sub-paths
		
		long startTime = System.currentTimeMillis();
		for (int c = 1; c <= dataModel.C+1; c++){

			if (System.currentTimeMillis()>timeLimit) break;
			ArrayList<PartialBackwardSequence> bwSeqs_head = bwSequences.get(c);

			for (Arc arc: dataModel.graph.incomingEdgesOf(c)){

				if (System.currentTimeMillis()>timeLimit) break;
				if (infeasibleArcs[arc.id] == 0 && arc.tail <= dataModel.C+1 && arc.head <= dataModel.C+1){ // only routing arcs
					
					double min_rc = monodirectionalFRC(arc, bwSeqs_head, bwSequences.get(arc.tail));
					if (!Double.isInfinite(min_rc) && min_rc - bestReducedCost > FRC_gap + dataModel.precision) {
						arcsToRemove.put(arc.id, min_rc); this.infeasibleArcs[arc.id] ++; 
					}
					
				}
			}

		}

		bwSequences.clear();

		long totalTime = System.currentTimeMillis()-startTime;
		dataModel.exactPricingTime+=totalTime;
		if (dataModel.print_log) {
			logger.debug("Time merging forward and backward labels: " + getTimeInSeconds(totalTime));
		}

		return arcsToRemove;

	}

	public double monodirectionalFRC(Arc arc, ArrayList<PartialBackwardSequence> bwSeqs_head, ArrayList<PartialBackwardSequence> bwSeqs_source) {

		int source = arc.tail;
		double bestCandidate = Double.POSITIVE_INFINITY;
		int gamma = dataModel.gamma;

		int max_q = bwSeqs_source.size();

		// For each backward seq at node j
		for (PartialBackwardSequence seq : bwSeqs_head) {

			// ---------- 1. Feasibility of extending through arc (i,j) ----------

			// Elementarity
			if (source > 0 && (seq.unreachable.get(source-1) || seq.ng_path.get(source-1))) continue;

			// Load feasibility
			int remLoad = seq.remainingLoad - dataModel.vertices[source].load;
			if (remLoad < 0) continue;

			// Time feasibility
			int remTime = seq.remainingTime - arc.time;
			if (remTime > dataModel.vertices[source].closing_tw) remTime = dataModel.vertices[source].closing_tw;
			else if (remTime < dataModel.vertices[source].opening_tw) continue;

			// Energy feasibility (worst-case logic preserved)
			int remainingEnergy = seq.worstEnergy - arc.energy;
			if (gamma > 0 && seq.remainingEnergy_gminus1 - arc.energy_deviation < seq.worstEnergy) {
				remainingEnergy = seq.remainingEnergy_gminus1 - arc.energy - arc.energy_deviation;
			}
			if (remainingEnergy < 0) continue;

			// ---------- 2. Build Ω(θ) = reachable customers after extension ----------

			BitSet unreachable = new BitSet(); unreachable.or(seq.unreachable);
			BitSet ng_path = new BitSet();
			if (source > 0) ng_path.set(source-1);
			else ng_path.or(seq.ng_path);

			//Mark unreachable customers and ng-path cycling restrictions
			if(source>0) {
				
				for (Arc c: dataModel.graph.incomingEdgesOf(source)) {
					if(c.tail == 0 || unreachable.get(c.tail-1)) continue;
					//unreachable
					if (remLoad-vertices[c.tail].load<0 || remTime-c.time<vertices[c.tail].opening_tw || remainingEnergy-c.min_energy < 0 || 
						Math.min(remTime-c.time, vertices[c.tail].closing_tw)-dataModel.graph.getEdge(0, c.tail).time<vertices[0].opening_tw ||
						remainingEnergy-c.min_energy-dataModel.graph.getEdge(0, c.tail).min_energy<0) {
						unreachable.set(c.tail-1);
					}

					//ng-path
					if (seq.ng_path.get(c.tail-1) && vertices[source].neighbors.contains(c.tail)) ng_path.set(c.tail-1);
				}
			}

			BitSet reachables = new BitSet(dataModel.C);
			for (int i = 1; i <= dataModel.C; i++) {
				if (!unreachable.get(i-1) && !ng_path.get(i-1)) reachables.set(i-1);
			}

			// ---------- 3. Compute reduced-cost proxy via Ω-coverage ----------

			double newReducedCost = seq.reducedCost + arc.modifiedCost;
			boolean foundDominatingSet = false;

			for (int q = 1; q <= max_q; q++) {

				PartialBackwardSequence l2 = bwSeqs_source.get(q-1);

				for (int c = reachables.nextSetBit(0); c >= 0; c = reachables.nextSetBit(c + 1)) {
					if (!l2.unreachable.get(c) && !l2.ng_path.get(c))  reachables.clear(c);
				}

				if (reachables.isEmpty()) {
					newReducedCost -= l2.reducedCost;
					foundDominatingSet = true;
					break;
				}
			}

			if (!foundDominatingSet) continue;

			// ---------- 4. Update best bound ----------
			if (newReducedCost < bestCandidate - dataModel.precision)
				bestCandidate = newReducedCost;
		}

		return bestCandidate;
	}

	private ArrayList<ArrayList<PartialBackwardSequence>> getBackwardSequences() {

		ArrayList<ArrayList<PartialBackwardSequence>> bwSequences = new ArrayList<>();

		// (Super) Depot labels
		ArrayList<Label> labels_to_clean = this.bwLabels.get(0);
		ArrayList<PartialBackwardSequence> bwSeq = new ArrayList<>();
		for (int ix = 0; ix < labels_to_clean.size(); ix++){
			Label l1 = labels_to_clean.get(ix); boolean dominated = false;
			for (int ix2 = ix+1; ix2 < labels_to_clean.size(); ix2++) if (isDominatedDepot(l1, labels_to_clean.get(ix2))) { dominated = true; break; }
			if (dominated) continue;
			else bwSeq.add(new PartialBackwardSequence(l1.reducedCost, l1.remainingEnergy, l1.remainingLoad, l1.remainingTime, l1.unreachable, l1.ng_path));
		} bwSeq.sort(Comparator.comparing(l -> l.reducedCost));
		bwSequences.add(bwSeq);

		// C1 PP vertices labels
		for (int i = 1; i <= dataModel.C; i++){
			labels_to_clean = this.bwLabels.get(i); bwSeq = new ArrayList<>();
			for (int ix = 0; ix < labels_to_clean.size(); ix++){
				Label l1 = labels_to_clean.get(ix); boolean dominated = false;
				for (int ix2 = ix+1; ix2 < labels_to_clean.size(); ix2++) if (isDominatedRouting(l1, labels_to_clean.get(ix2))) { dominated = true; break; }
				if (dominated) continue;
				else bwSeq.add(new PartialBackwardSequence(l1.reducedCost, l1.remainingEnergy, l1.remainingLoad, l1.remainingTime, l1.unreachable, l1.ng_path));
			} bwSeq.sort(Comparator.comparing(l -> l.reducedCost));
			bwSequences.add(bwSeq);
		}

		// Returning depot (unique) label
		bwSeq = new ArrayList<>();
		Label l1 = this.bwLabels.get(dataModel.C+1).get(0);
		bwSeq.add(new PartialBackwardSequence(l1.reducedCost, l1.remainingEnergy, l1.remainingLoad, l1.remainingTime, l1.unreachable, l1.ng_path));
		bwSequences.add(bwSeq);

		this.bwLabels.clear();

		return bwSequences;
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



	private final class PartialBackwardSequence {

		public double reducedCost;
		public int worstEnergy;
		public int remainingEnergy_gminus1;
		public int remainingLoad;
		public int remainingTime;
		public BitSet unreachable;
		public BitSet ng_path;

		private PartialBackwardSequence(double rc, int[] remEnergy, int remLoad, int remTime, boolean[] unreach, boolean[] ngp){
			this.reducedCost = rc;
			this.remainingLoad = remLoad;
			this.remainingTime = remTime;
			
			this.worstEnergy = remEnergy[dataModel.gamma];
			if (dataModel.gamma > 0) this.remainingEnergy_gminus1 = remEnergy[dataModel.gamma-1];
			else this.remainingEnergy_gminus1 = 0;

			this.unreachable = new BitSet(); this.ng_path = new BitSet();
			for (int i = 0; i < dataModel.C; i++) {
				if (unreach[i]) this.unreachable.set(i);
				if (ngp[i]) this.ng_path.set(i);
			}

		}
	}
}
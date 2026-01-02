package columnGeneration;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.HashMap;
import java.util.Map;
import java.util.HashSet;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.TreeSet;

import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblem;

import branchAndPrice.ChargingTimeInequality;
import model.EVRPTW;
import model.EVRPTW.Arc;
import model.EVRPTW.Vertex;

/**
 * This class defines the pricing problem. 
 * We simply extend the pricing problem included in the framework (there is no need for modification, only one pricing problem)
 */
public final class PricingProblem extends AbstractPricingProblem<EVRPTW> {

	public ArrayList<SubsetRowInequality> subsetRowCuts; 				//subset row cuts considered
	public Set<ChargingTimeInequality> branchesOnChargingTimes;			//branching on charging times
	public double bestReducedCost = -Double.MAX_VALUE; 					//best reduced cost found by the exact labeling
	public double reducedCostThreshold = 0; 							//minimum reduced cost when arriving at the depot source

	// Information for Fixing by Reduced Costs procedure
	public ArrayList<ArrayList<PartialForwardSequence>> fwSequences = new ArrayList<>();
	public ArrayList<ArrayList<Label>> bwLabels = new ArrayList<>();
	public ArrayList<ArrayList<Integer>> SRCIndices = new ArrayList<>();
	public int[] infeasibleArcs;
	public Vertex[] vertices;

	// Charging pricing information
	public Map<Integer, Map<Integer, Double>> charging_bounds;

	public PricingProblem(EVRPTW modelData, String name) {
		super(modelData, name);
	}

	public Map<Integer, Double> fixByReducedCosts(long timeLimit){
		
		double FRC_gap = dataModel.UB_FRC - dataModel.LB_FRC;
		this.vertices = dataModel.vertices;

		Map<Integer, Double> arcsToRemove = new HashMap<Integer, Double>();
		FixByReducedCostSolver FRC = new FixByReducedCostSolver(dataModel, timeLimit);

		ArrayList<ArrayList<PartialBackwardSequence>> bwSequences = FRC.getBackwardSequences(); // get *sorted* backward sub-paths
		this.fwSequences = FRC.runForwardLabeling(); // get *sorted* forward sub-paths
		this.compute_charging_bounds(); // Compute the charging bounds for ALL charging and departure times
		
		long startTime = System.currentTimeMillis();
		for (int c = 1; c <= dataModel.C+1; c++){
			ArrayList<PartialBackwardSequence> backwardSequences = bwSequences.get(c);

			for (Arc arc: dataModel.graph.incomingEdgesOf(c)){
				if (infeasibleArcs[arc.id] == 0 && arc.tail <= dataModel.C+1 && arc.head <= dataModel.C+1){ // only routing arcs
					
					ArrayList<PartialBackwardSequence> filteredBw = new ArrayList<PartialBackwardSequence>();
					for (PartialBackwardSequence bwS: backwardSequences) { if ((arc.tail == 0) || (arc.tail > 0 && !bwS.unreachable.get(arc.tail-1) && !bwS.ng_path.get(arc.tail-1))) filteredBw.add(bwS); }
					ArrayList<PartialForwardSequence> forwardSequences = this.fwSequences.get(arc.tail);

					if (!filteredBw.isEmpty()){
						double max_rc = arc.modifiedCost + forwardSequences.get(forwardSequences.size()-1).reducedCost + filteredBw.get(filteredBw.size()-1).reducedCost;

						if (max_rc > FRC_gap + dataModel.precision){ // This is an optimistic bound of the worst reduced cost, missing the charging bound, not guaranteed to be the actual worst
							double min_rc = findMinimumRCPath(filteredBw, forwardSequences, arc);
							if (!Double.isInfinite(min_rc) && min_rc - bestReducedCost > FRC_gap + dataModel.precision) {
								arcsToRemove.put(arc.id, min_rc); this.infeasibleArcs[arc.id] ++; 
							}
							if (min_rc < bestReducedCost - dataModel.precision) logger.debug("!!! Arc {} has a merged label with a reduced cost of {}", new Object[]{arc.toString(), min_rc});
						}

					}
				}
			}

			backwardSequences.clear();
		}

		this.fwSequences.clear(); bwSequences.clear(); this.charging_bounds.clear();

		long totalTime = System.currentTimeMillis()-startTime;
		dataModel.exactPricingTime+=totalTime;
		if (dataModel.print_log) {
			logger.debug("Time merging forward and backward labels: " + FRC.getTimeInSeconds(totalTime));
		}

		return arcsToRemove;

	}

	private double findMinimumRCPath(ArrayList<PartialBackwardSequence> bwSequences, ArrayList<PartialForwardSequence> fwSequences, Arc arc) {

        int nFw = fwSequences.size();
        int nBw = bwSequences.size();

        PriorityQueue<MergeState> pq = new PriorityQueue<>( (s1, s2) -> Double.compare(s1.rc, s2.rc) );
        pq.add(new MergeState(0, 0, arc.modifiedCost + bwSequences.get(0).reducedCost + fwSequences.get(0).reducedCost));

        HashSet<Long> visited = new HashSet<>();
        visited.add(key(0, 0));

		double bestReducedCost = Double.POSITIVE_INFINITY;
        while (!pq.isEmpty()) {
            MergeState current = pq.poll();
            
			int fw = current.f;
            PartialForwardSequence fwSeq = fwSequences.get(fw);
            
			int bw = current.b;
            PartialBackwardSequence bwSeq = bwSequences.get(bw);

			if (!Double.isInfinite(bestReducedCost) && bestReducedCost <= current.rc - dataModel.precision){ return bestReducedCost; }
			else {
				MergedSequence mergedPath = mergeLabel(fwSeq, bwSeq, arc);
				if (mergedPath != null){
					double chBound = this.charging_bounds.get(mergedPath.chargingTime).get(mergedPath.departureTime);
					double complete_rc = mergedPath.reducedCost + chBound;
					
					if (complete_rc < bestReducedCost - dataModel.precision){
						bestReducedCost = complete_rc;
						if (chBound < dataModel.precision){ return bestReducedCost; } // if the charging bound of the new best label is 0, is optimal
					}
				}
			}

            // neighbor: (i+1, j)
            if (fw + 1 < nFw) {
                long k = key(fw + 1, bw);
                if (visited.add(k)) {
                    pq.add(new MergeState(fw + 1, bw, arc.modifiedCost + fwSequences.get(fw+1).reducedCost + bwSeq.reducedCost));
                }
            }

            // neighbor: (i, j+1)
            if (bw + 1 < nBw) {
                long k = key(fw, bw + 1);
                if (visited.add(k)) {
                    pq.add(new MergeState(fw, bw + 1, arc.modifiedCost + fwSeq.reducedCost + bwSequences.get(bw + 1).reducedCost));
                }
            }
        }

        return bestReducedCost;
    }

	private long key(int fw, int bw) {
		return (((long) fw) << 32) | (bw & 0xffffffffL);
	}

	private MergedSequence mergeLabel(PartialForwardSequence fwSequence, PartialBackwardSequence bwSeq, Arc arc){
		
		/////////////////////////////////
		/// MERGE FEASIBILITY ASSESSMENT
		/////////////////////////////////
		
		ArrayList<Integer> arcExtensions = fwSequence.arcsSequence;

		// Initialize current_ng_path as BitSet from bwSeq.ng_path
		BitSet currentNg = new BitSet(dataModel.C);
		currentNg.or(bwSeq.ng_path);

		// (ng-path) Elementarity assessment
		for (int arcID : arcExtensions) {

			Arc arcExt = dataModel.arcs[arcID]; int source = arcExt.tail;
			if (source > 0 && bwSeq.unreachable.get(source - 1) || currentNg.get(source - 1)) return null;

			// Build next ng-path state
			BitSet nextNg = new BitSet(dataModel.C);

			if (source > 0) {
				nextNg.set(source - 1); // Visiting source marks it as visited

				for (Arc c : dataModel.graph.incomingEdgesOf(source)) {
					int tail = c.tail;
					if (tail == 0 || bwSeq.unreachable.get(tail - 1)) continue;

					// ng-path rule
					if (currentNg.get(tail - 1) && vertices[source].neighbors.contains(tail)) nextNg.set(tail - 1);
				}

			} else nextNg.or(currentNg);

			// Move forward
			currentNg = nextNg;
		}
		
		if (fwSequence.remainingLoad + bwSeq.remainingLoad - dataModel.Q < 0) return null; // Load feasibility
		if (fwSequence.cumulativeTime + arc.time > bwSeq.remainingTime) return null; // Time feasibility
		
		// Worst-case energy feasibility
		int remainingEnergy = bwSeq.nominalEnergy - arc.energy - fwSequence.nominalEnergy; if (remainingEnergy < 0) return null; // Nominal energy consumption

		if (dataModel.gamma > 0){
			PriorityQueue<Integer> energy_deviations = new PriorityQueue<>(bwSeq.worstEnergyDevs);
			energy_deviations.add(arc.energy_deviation); energy_deviations.addAll(fwSequence.worstEnergyDevs);
			
			for (int g = 1; g <= dataModel.gamma; g++) { // Worst-case energy consumption
				int dev = energy_deviations.poll();
				remainingEnergy -= dev;
				if (remainingEnergy < 0) return null;
			}
		}

		int chargingTime = dataModel.f_inverse[dataModel.E-remainingEnergy];
		
		/////////////////////////////////
		/// RESOURCES UPDATE
		/////////////////////////////////
		
		double reducedCost = fwSequence.reducedCost + bwSeq.reducedCost + arc.modifiedCost;
		reducedCost = Math.floor(reducedCost*10000)/10000;
		
		int remainingTime = bwSeq.remainingTime - arc.time;
		if(remainingTime > vertices[arc.tail].closing_tw) remainingTime = vertices[arc.tail].closing_tw;
		
		for (int arcID: arcExtensions){
			Arc extArc = dataModel.arcs[arcID];
			int source = extArc.tail;
			
			remainingTime -= extArc.time;
			if(remainingTime > vertices[source].closing_tw) remainingTime = vertices[source].closing_tw;
		}
		
		int departure = (int)(remainingTime/10);
		if (chargingTime >= departure) return null; // Charging interval feasibility

		return new MergedSequence(reducedCost, chargingTime, departure);
	}

	private void compute_charging_bounds(){

		this.charging_bounds = new HashMap<>();

		// 1. Group by charging time b
		Map<Integer, TreeSet<Integer>> charging_times = new HashMap<>();
		int minT = (int) (dataModel.vertices[0].opening_tw/10);
		int maxT = dataModel.last_charging_period + 1;
		for (int b = 1; b <= dataModel.f_inverse[dataModel.E]; b++){
			TreeSet<Integer> departures = new TreeSet<>();
			for (int d = Math.max(b+1, minT); d <= maxT; d++) { departures.add(d); }
			charging_times.put(b, departures);
		}

		// Precompute fixed sums of the charging dual variables
		double[] S = new double[maxT]; S[0] = 0.0;
		for (int t = 1; t < maxT; t++) { S[t] = S[t - 1] + this.dualCosts[dataModel.C + t - 1]; }

		process_charging_times_map(charging_times, S);

	}

	private void process_charging_times_map(Map<Integer, TreeSet<Integer>> charging_times, double[] S){

		for (Map.Entry<Integer, TreeSet<Integer>> e : charging_times.entrySet()){

			int b = e.getKey();
        	TreeSet<Integer> departures = e.getValue();

			int initial_t = 1;
			double rc = - (S[b]-S[0]); double min_rc = rc;

			Map<Integer, Double> boundsMap = new LinkedHashMap<>();
			
			for (int d: departures){
				if (d <= b) continue; // skip if departure does not allow for sufficient charging

				for (int t=initial_t; t<=d-b-1; t++){
					rc = - (S[t+b] - S[t]);
					if (rc < min_rc - dataModel.precision) min_rc = rc;
				}
				boundsMap.put(d, min_rc);
				initial_t = d-b;
			}
			
			this.charging_bounds.put(b, boundsMap);
			
		}

	}

	private static final class MergeState {
        final int f, b;
        final double rc;
        MergeState(int f, int b, double rc) { this.f = f; this.b = b; this.rc = rc; }
    }

	private final class PartialForwardSequence {

		public double reducedCost;
		public ArrayList<Integer> arcsSequence;
		public int nominalEnergy;
		public ArrayList<Integer> worstEnergyDevs;
		public int remainingLoad;
		public int cumulativeTime;

		private PartialForwardSequence(double rc, ArrayList<Integer> aSeq, int nomEnergy, ArrayList<Integer> devs, int remLoad, int cumTime){
			this.reducedCost = rc;
			this.arcsSequence = aSeq;
			this.nominalEnergy = nomEnergy;
			this.worstEnergyDevs = devs;
			this.remainingLoad = remLoad;
			this.cumulativeTime = cumTime;
		}
	}

	private final class PartialBackwardSequence {

		public double reducedCost;
		public int nominalEnergy;
		public PriorityQueue<Integer> worstEnergyDevs;
		public int remainingLoad;
		public int remainingTime;
		public BitSet unreachable;
		public BitSet ng_path;

		private PartialBackwardSequence(double rc, int[] remEnergy, int remLoad, int remTime, boolean[] unreach, boolean[] ngp){
			this.reducedCost = rc;
			this.remainingLoad = remLoad;
			this.remainingTime = remTime;
			
			this.nominalEnergy = remEnergy[0];
			this.worstEnergyDevs = new PriorityQueue<>(Comparator.reverseOrder());
			if (dataModel.gamma > 0) {
				for (int g = 1; g <= dataModel.gamma; g++) this.worstEnergyDevs.add(remEnergy[g-1] - remEnergy[g]);
			}

			this.unreachable = new BitSet(); this.ng_path = new BitSet();
			for (int i = 0; i < dataModel.C; i++) {
				if (unreach[i]) this.unreachable.set(i);
				if (ngp[i]) this.ng_path.set(i);
			}

		}
	}

	private final class MergedSequence {

		public double reducedCost;
		public int chargingTime;
		public int departureTime;

		private MergedSequence(double rc, int b, int d){
			this.reducedCost = rc;
			this.chargingTime = b;
			this.departureTime = d;
		}
	}

	public final class FixByReducedCostSolver {

		public EVRPTW dataModel;
		public PricingProblem pricingProblem;
		public Vertex[] vertices; 			//vertices of the instance
		public PriorityQueue<Vertex> nodesToProcess; 			//labels that need be processed
		public long timeLimit;

		/**
		 * Labeling algorithm to solve the ng-SPPRC
		 */
		public FixByReducedCostSolver(EVRPTW dataModel, long timeLimit) {
			this.dataModel = dataModel;
			this.pricingProblem = PricingProblem.this;
			this.vertices = dataModel.vertices;
			this.nodesToProcess = new PriorityQueue<Vertex>(dataModel.V, new SortVertices());
			this.timeLimit = timeLimit;
		}

		public ArrayList<ArrayList<PartialBackwardSequence>> getBackwardSequences() {

			ArrayList<ArrayList<PartialBackwardSequence>> bwSequences = new ArrayList<>();
			for (int i = 0; i <= dataModel.C+1; i++){
				ArrayList<PartialBackwardSequence> bwSeq = new ArrayList<>();
				for (Label l: pricingProblem.bwLabels.get(i)){ bwSeq.add(new PartialBackwardSequence(l.reducedCost, l.remainingEnergy, l.remainingLoad, l.remainingTime, l.unreachable, l.ng_path)); }
				bwSeq.sort( Comparator.comparing(l -> l.reducedCost) );
				bwSequences.add(bwSeq);
			}

			pricingProblem.bwLabels.clear();

			return bwSequences;
		}

		public ArrayList<ArrayList<PartialForwardSequence>> runForwardLabeling() {

			//Initialization
			int[] remain_energy = new int[dataModel.gamma + 1]; Arrays.fill( remain_energy, dataModel.E);
			Label initialLabel = new Label(0, 0, 0, 0, dataModel.Q, vertices[0].opening_tw, remain_energy, 0, new boolean[dataModel.C], new boolean[dataModel.C], new boolean[pricingProblem.subsetRowCuts.size()], new HashSet<Integer>(pricingProblem.subsetRowCuts.size()));
			this.nodesToProcess.add(vertices[0]);
			initialLabel.index = 0;
			vertices[0].unprocessedLabels.add(initialLabel);

			//Labeling algorithm
			long startTime = System.currentTimeMillis();
			while (!nodesToProcess.isEmpty() && System.currentTimeMillis()<timeLimit) {
				ArrayList<Label> labelsToProcessNext = labelsToProcessNext();
				for(Label currentLabel: labelsToProcessNext) {
					boolean isDominated = checkForwardDominance(currentLabel);
					if(isDominated) continue;
					else {currentLabel.index = vertices[currentLabel.vertex].processedLabels.size(); vertices[currentLabel.vertex].processedLabels.add(currentLabel);}
					for(Arc a: dataModel.graph.outgoingEdgesOf(currentLabel.vertex)) {
						if(this.pricingProblem.infeasibleArcs[a.id] > 0) continue;
						Label extendedLabel;
						extendedLabel = extendForwardLabel(currentLabel, a);
						if (extendedLabel!=null) { //verifies if the extension is feasible
							updateNodesToProcess(extendedLabel);
						}
					}
				}
			}

			ArrayList<ArrayList<PartialForwardSequence>> fwSequences = new ArrayList<>();
			for (int i = 0; i <= dataModel.C; i++){
				ArrayList<PartialForwardSequence> allSequences = new ArrayList<>();
				for (Label label: vertices[i].processedLabels){ allSequences.add(get_forward_sequence(label)); }
				allSequences.sort( Comparator.comparing(l -> l.reducedCost) );
				fwSequences.add(allSequences);
			}

			for (int i = 0; i <= dataModel.C+1; i++) {
				vertices[i].processedLabels = new ArrayList<Label>(dataModel.numArcs);
				vertices[i].unprocessedLabels =  new PriorityQueue<Label>(dataModel.numArcs, new Label.SortLabels());
			}
			vertices[dataModel.C+1].processedLabels = new ArrayList<Label>(dataModel.numArcs);
			vertices[dataModel.C+1].unprocessedLabels =  new PriorityQueue<Label>(dataModel.numArcs, new Label.SortLabels());

			long totalTime = System.currentTimeMillis()-startTime;
			dataModel.exactPricingTime+=totalTime;
			if (dataModel.print_log) logger.debug("Time running forward routing labeling algorithm: " + getTimeInSeconds(totalTime));

			return fwSequences;
		}

		private PartialForwardSequence get_forward_sequence(Label fwL){

			ArrayList<Integer> aSeq = new ArrayList<>();
			boolean[] route = new boolean[dataModel.C+1];
			
			Label currentLabel = fwL.clone();
			int currentVertex = currentLabel.vertex;
			while(currentVertex!=0) {
				route[currentVertex] = true;

				Arc currentArc = dataModel.arcs[currentLabel.nextArc];
				int nextVertex = currentArc.tail;

				currentLabel = vertices[nextVertex].processedLabels.get(currentLabel.nextLabelIndex);
				if(currentArc.tail>=0 && currentArc.tail<=dataModel.C) aSeq.add(currentArc.id);
				currentVertex = nextVertex;
			}

			ArrayList<Integer> worst_energy_deviations = new ArrayList<>();
			for (int g=0; g<dataModel.gamma; g++){ worst_energy_deviations.add(fwL.remainingEnergy[g] - fwL.remainingEnergy[g+1]); }

			return new PartialForwardSequence(fwL.reducedCost, aSeq, dataModel.E - fwL.remainingEnergy[0], worst_energy_deviations, fwL.remainingLoad, fwL.remainingTime);
		}

		public Label extendForwardLabel(Label currentLabel, Arc arc) {

			int head = arc.head; int depot = dataModel.C+1;
			if (head < depot && (currentLabel.unreachable[head-1] || currentLabel.ng_path[head-1])) return null;

			double reducedCost = currentLabel.reducedCost+arc.modifiedCost;
			reducedCost = Math.floor(reducedCost*10000)/10000;

			int remainingLoad = currentLabel.remainingLoad-vertices[head].load;
			int remainingTime = currentLabel.remainingTime+arc.time; // For forward labeling time resource must be computed in ascending order
			if(remainingTime<vertices[head].opening_tw) remainingTime = vertices[head].opening_tw;

			int[] remainingEnergy = new int[dataModel.gamma + 1];
			remainingEnergy[0] = currentLabel.remainingEnergy[0]-arc.energy; if (remainingEnergy[0] < 0) return null;
			for (int gam = 1; gam <= dataModel.gamma; gam++){
				if (currentLabel.remainingEnergy[gam-1] - arc.energy_deviation < currentLabel.remainingEnergy[gam]){ remainingEnergy[gam] = currentLabel.remainingEnergy[gam-1] - arc.energy - arc.energy_deviation; }
				else { remainingEnergy[gam] = currentLabel.remainingEnergy[gam] - arc.energy; }
				if (remainingEnergy[gam] < 0) return null;
			}
			
			int chargingTime = dataModel.f_inverse[dataModel.E-remainingEnergy[dataModel.gamma]];

			//Quick check
			if (head < depot){
				if(remainingTime+dataModel.graph.getEdge(head, depot).time > vertices[depot].closing_tw) return null;
			}

			//Check whether the extension is actually feasible
			if(remainingTime>vertices[head].closing_tw || chargingTime > dataModel.f_inverse[dataModel.E]) return null;

			boolean[] unreachable = Arrays.copyOf(currentLabel.unreachable.clone(), currentLabel.unreachable.length);
			boolean[] ng_path = new boolean[dataModel.C];
			if (head<depot) ng_path[head-1] = true;
			else ng_path = Arrays.copyOf(currentLabel.ng_path, currentLabel.ng_path.length);

			//Mark unreachable customers and ng-path cycling restrictions
			for (Arc c: dataModel.graph.outgoingEdgesOf(head)) {
				if(c.head==depot || unreachable[c.head-1]) continue;
				//unreachable
				if (remainingLoad-vertices[c.head].load<0 || remainingTime+c.time>vertices[c.head].closing_tw || remainingEnergy[dataModel.gamma]-c.min_energy < 0 || 
					Math.max(remainingTime+c.time, vertices[c.head].opening_tw)+dataModel.graph.getEdge(c.head, depot).time>vertices[depot].closing_tw ||
					remainingEnergy[dataModel.gamma]-c.min_energy - dataModel.graph.getEdge(c.head, depot).min_energy<0) {
					unreachable[c.head-1] = true;
				}

				//ng-path
				if (currentLabel.ng_path[c.head-1] && vertices[head].neighbors.contains(c.head)) ng_path[c.head-1] = true;
				else ng_path[c.head-1] = false;
			}
			
			Label extendedLabel = new Label(head, arc.id, currentLabel.index, reducedCost, remainingLoad, remainingTime, remainingEnergy, chargingTime,unreachable, ng_path, currentLabel.eta.clone(), new HashSet<Integer>(currentLabel.srcIndices));
			return extendedLabel;

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
						isDominated = isForwardDominated(currentLabel, L2);
						if(isDominated) break;
					}
					if(!isDominated) labelsToProcessNext.add(currentLabel);
				}
				if(currentVertex.unprocessedLabels.isEmpty() || (currentVertex.id<=dataModel.C && currentVertex.unprocessedLabels.peek().remainingLoad<currentLabel.remainingLoad)) break;
			}

			if(!currentVertex.unprocessedLabels.isEmpty()) nodesToProcess.add(currentVertex);
			return labelsToProcessNext;
		}

		/**
		 * Given a new (non-dominated) label, updates the nodes to be processed
		 */
		public void updateNodesToProcess(Label extendedLabel) {
			Vertex currentVertex = vertices[extendedLabel.vertex];
			if(currentVertex.id == dataModel.V) vertices[extendedLabel.vertex].unprocessedLabels.add(extendedLabel);
			else if(currentVertex.unprocessedLabels.isEmpty()) {currentVertex.unprocessedLabels.add(extendedLabel); nodesToProcess.add(currentVertex);}
			else currentVertex.unprocessedLabels.add(extendedLabel);
		}

		public boolean checkForwardDominance(Label newLabel) {

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

				if(isForwardDominated(existingLabel, newLabel)) {
					//existing_is_discarded = true; // DELETE LATER
					labelsToDelete.add(existingLabel);
				}
				
			}
			currentVertex.unprocessedLabels.removeAll(labelsToDelete);
			if(currentVertex.unprocessedLabels.isEmpty()) nodesToProcess.remove(currentVertex);

			//boolean new_is_discarded = false; // DELETE LATER
			for(Label existingLabel: currentVertex.processedLabels) {
				//int[] el_sequence = get_route_sequence(existingLabel); // DELETE LATER
				if(isForwardDominated(newLabel, existingLabel)) return true;
				
			}

			return false;
		}

		public boolean isForwardDominated(Label L1, Label L2) {
			
			/* int[] nl_sequence = get_route_sequence(L1); // DELETE LATER
			int[] el_sequence = get_route_sequence(L2); // DELETE LATER */

			if (L1.vertex<dataModel.C+1 && L2.remainingLoad<L1.remainingLoad) return false; 	//load
			if (L2.reducedCost-L1.reducedCost>dataModel.precision) return false; 	//reduced cost
			if (L1.remainingTime<L2.remainingTime) return false; 					//time
			
			for (int gam = 0; gam <= dataModel.gamma; gam ++){ if (L2.remainingEnergy[gam]<L1.remainingEnergy[gam]) return false; } // energy

			// Ng-paths and unreachable resources
			Vertex currentVertex = vertices[L1.vertex];
			if (currentVertex.id > 0) {
				for(int i=1; i<=dataModel.C; i++) {
					
					//boolean check_binaries = (L2.ng_path[i-1] || L2.unreachable[i-1]) && !(L1.ng_path[i-1] || L1.unreachable[i-1]);
					boolean other_way = L2.ng_path[i-1] && (!L1.unreachable[i-1] && !L1.ng_path[i-1]); // Dani's way
					if (other_way) {
						return false;
					}
				}
			}

			return true;
			
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

}
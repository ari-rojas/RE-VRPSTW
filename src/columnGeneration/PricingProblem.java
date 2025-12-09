package columnGeneration;

import java.util.Collections;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.HashMap;
import java.util.Map;
import java.util.HashSet;
import java.util.List;
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
	public ArrayList<ArrayList<PartialSequence>> fwSequences = new ArrayList<>();
	public ArrayList<ArrayList<Label>> bwLabels = new ArrayList<>();
	public ArrayList<ArrayList<Integer>> SRCIndices = new ArrayList<>();
	public int[] infeasibleArcs;

	// Charging pricing information
	public double[] S;
	public boolean[] negative_charging_duals;
	public Map<Integer, Map<Integer, Double>> charging_bounds;
	public Map<Integer, Map<Integer, Double>> charging_reducedCosts;
	public int maxT;

	public PricingProblem(EVRPTW modelData, String name) {
		super(modelData, name);
		
	}

	public Map<Integer, Double> fixByReducedCosts(long timeLimit){

		FixByReducedCostSolver FRC = new FixByReducedCostSolver(dataModel, timeLimit);
		this.fwSequences = FRC.runForwardLabeling();

		Map<Integer, List<Integer>> mergedMap = new HashMap<>();
		ArrayList<Label> mergedLabels = new ArrayList<>();

		long startTime = System.currentTimeMillis(); int cont = 0; int arcCont = 0;
		for (Arc arc: dataModel.graph.edgeSet()){
			
			if (arc.minCostAlternative && infeasibleArcs[arc.id] == 0 && arc.tail <= dataModel.C+1 && arc.head <= dataModel.C+1){ // only routing arcs
				arcCont ++;

				ArrayList<PartialSequence> forwardSequences = this.fwSequences.get(arc.tail);
				ArrayList<Label> backwardLabels = this.bwLabels.get(arc.head);

				for (PartialSequence fwL: forwardSequences){ // Attempt to merge all the forward and backward labels that use this arc
					for (Label bwL: backwardLabels){

						Label newLabel = mergeLabel(fwL, bwL, arc);
						if (newLabel != null){
							newLabel.index = cont;
							mergedLabels.add(newLabel);
							mergedMap.computeIfAbsent(arc.id, k -> new ArrayList<Integer>()).add(cont);
							cont ++;
						}

					}
				}
			}
			
		}

		////////////////////////////////////////////
		/// DEBUG
		///////////////////////////////////////////
		
		logger.debug("Total of arcs evaluated: " + arcCont);
		logger.debug("Arcs with feasible merged labels:");
		for (Map.Entry<Integer, List<Integer>> entry: mergedMap.entrySet()){
			logger.debug("Arc "+ dataModel.arcs[entry.getKey()].toString() + " , with "+ entry.getValue().size() + " labels");
		}

		int maxMerged = 0; int arcMaxMerged = -1;
		for (Map.Entry<Integer, List<Integer>> entry: mergedMap.entrySet()) { if (entry.getValue().size() > maxMerged) { maxMerged = entry.getValue().size(); arcMaxMerged = entry.getKey(); } }
		logger.debug("Arc with the most merged labels: {} with {} labels", new Object[]{dataModel.arcs[arcMaxMerged].toString(), maxMerged});
		for (int labelIx: mergedMap.get(arcMaxMerged)) { logger.debug(mergedLabels.get(labelIx).toString()); }

		////////////////////////////////////////////
		/// END DEBUG
		///////////////////////////////////////////


		// Just for the charging times bound computation
		for (Label label: bwLabels.get(0)){ mergedLabels.add(label.clone()); }

		long totalTime = System.currentTimeMillis()-startTime;
		dataModel.exactPricingTime+=totalTime;
		if (dataModel.print_log) {
			logger.debug("Time merging forward and backward labels: " + FRC.getTimeInSeconds(totalTime));
			logger.debug("Found "+mergedLabels.size()+" merged labels");
		}

		startTime = System.currentTimeMillis();
		this.compute_charging_bounds(mergedLabels);
		totalTime = System.currentTimeMillis()-startTime;
		if (dataModel.print_log) {
			logger.debug("Time computing charging bounds: " + FRC.getTimeInSeconds(totalTime));
		}

		/// Computing the best reduced cost of the pricing problem
		bestReducedCost = Double.MAX_VALUE;
		for (Label label: bwLabels.get(0)) {
			double rc = label.reducedCost + this.charging_bounds.get(label.chargingTime).get((int)(label.remainingTime/10));
			if (rc < bestReducedCost - dataModel.precision) bestReducedCost = rc; 
		}

		//////////////////////////////////////////////////////////////////////////////////////////
		/* logger.debug("Printing the non-dominated backward labels at the depot (0):");
		for (Label label: this.bwLabels.get(0)){
			int departure = (int)(label.remainingTime/10);
			logger.debug("Label: {}, Bound: {}",new Object[]{label.toString(), charging_bounds.get(label.chargingTime).get(departure)});
		} */
		///////////////////////////////////////////////////////////////////////////////////////////

		Map<Integer, Double> arcsToRemove = new HashMap<Integer, Double>();
		for (Map.Entry<Integer, List<Integer>> entry: mergedMap.entrySet()){

			int arcID = entry.getKey(); List<Integer> labelIDs = entry.getValue();

			double min_rc = Double.MAX_VALUE;
			for (int labID: labelIDs){
				Label label = mergedLabels.get(labID);
				double rc = label.reducedCost + this.charging_bounds.get(label.chargingTime).get((int)(label.remainingTime/10));
				if (rc < min_rc - dataModel.precision) min_rc = rc;
			}

			if (min_rc - bestReducedCost > dataModel.UB_FRC - dataModel.LB_FRC + dataModel.precision) arcsToRemove.put(arcID, min_rc);

		}

		return arcsToRemove;

	}

	private Label mergeLabel(PartialSequence fwSequence, Label bwL, Arc arc){
		
		/////////////////////////////////
		/// MERGE FEASIBILITY ASSESSMENT
		/////////////////////////////////
		
		boolean[] fwRoute = fwSequence.route;
		ArrayList<Integer> arcExtensions = fwSequence.arcsSequence;
		
		// Elementarity assessment
		if ( arc.tail > 0 && (bwL.unreachable[arc.tail - 1] || bwL.ng_path[arc.tail - 1])) return null;
		for (int i=1; i<=dataModel.C; i++){ 
			if (fwRoute[i] && (bwL.unreachable[i-1] || bwL.ng_path[i-1])) return null;
		}
		
		if (fwSequence.remainingLoad + bwL.remainingLoad - dataModel.Q < 0) return null; // Load feasibility
		if (fwSequence.cumulativeTime + arc.time > bwL.remainingTime) return null; // Time feasibility
		
		// Worst-case energy feasibility
		int[] remainingEnergy = new int[dataModel.gamma + 1];
		remainingEnergy[0] = bwL.remainingEnergy[0] - arc.energy - fwSequence.nominalEnergy; if (remainingEnergy[0] < 0) return null; // Nominal energy consumption

		ArrayList<Integer> energy_deviations = new ArrayList<>();
		energy_deviations.add(arc.energy_deviation); energy_deviations.addAll(fwSequence.worstEnergyDevs);
		for (int g=0; g<dataModel.gamma; g++){ energy_deviations.add(bwL.remainingEnergy[g] - bwL.remainingEnergy[g+1]); }
		
		energy_deviations.sort(Comparator.reverseOrder()); int cont = 0;
		for (Integer e_dev: energy_deviations){ cont ++; remainingEnergy[cont] = remainingEnergy[cont-1]-e_dev; if (cont == dataModel.gamma) break; }
		if (remainingEnergy[dataModel.gamma] < 0) return null; // Worst-case energy consumption

		int chargingTime = dataModel.f_inverse[dataModel.E-remainingEnergy[dataModel.gamma]];
		
		/////////////////////////////////
		/// RESOURCES UPDATE
		/////////////////////////////////
		
		Label updatedLabel = new Label(0, arc.id, -1, 0, 0, 0, new int[dataModel.gamma + 1], 0, new boolean[dataModel.C], new boolean[dataModel.C], new boolean[1], new HashSet<>() );

		for (int g = 0; g<= dataModel.gamma; g++) updatedLabel.remainingEnergy[g] = remainingEnergy[g];
		updatedLabel.chargingTime = chargingTime;
		updatedLabel.remainingLoad = fwSequence.remainingLoad + bwL.remainingLoad - dataModel.Q;
		
		double reducedCost = bwL.reducedCost + arc.modifiedCost;
		int remainingTime = bwL.remainingTime - arc.time;
		if(remainingTime > dataModel.vertices[arc.tail].closing_tw) remainingTime = dataModel.vertices[arc.tail].closing_tw;
		
		for (int arcID: arcExtensions){
			Arc extArc = dataModel.arcs[arcID];
			int source = extArc.tail;
			
			reducedCost += extArc.modifiedCost;
			remainingTime -= extArc.time;
			if(remainingTime> dataModel.vertices[source].closing_tw) remainingTime = dataModel.vertices[source].closing_tw;
		}
		reducedCost = Math.floor(reducedCost*10000)/10000;
		
		updatedLabel.reducedCost = reducedCost;
		
		if (chargingTime >= (int)(remainingTime/10)) return null; // Charging interval feasibility
		updatedLabel.remainingTime = remainingTime;

		updatedLabel.vertex = 0;

		return updatedLabel;
	}

	public void compute_charging_bounds(ArrayList<Label> labels){

		this.charging_reducedCosts = new HashMap<>();
		this.charging_bounds = new HashMap<>();

		// 1. Group labels by charging time b
		Map<Integer, TreeSet<Integer>> charging_times = new HashMap<>();
		for (Label label : labels) charging_times.computeIfAbsent(label.chargingTime, k -> new TreeSet<Integer>()).add((int)(label.remainingTime / 10));

		// Precompute fixed sums of the charging dual variables
		this.maxT = 0; for (TreeSet<Integer> set : charging_times.values()) this.maxT = Math.max(this.maxT, set.last()-1);
		this.S = new double[dataModel.last_charging_period + 1]; this.S[0] = 0.0;
		this.negative_charging_duals = new boolean[dataModel.last_charging_period + 1];
		for (int t = 1; t <= this.maxT; t++) {
			double dual = this.dualCosts[dataModel.C + t - 1];
			this.negative_charging_duals[t] = dual < -dataModel.precision;
			this.S[t] = this.S[t - 1] + dual;
		}

		/* logger.debug("Identifying time periods with negative dual");
		for (int t = 1; t <= this.maxT; t++){
			if (this.negative_charging_duals[t]) logger.debug("Period "+t+": Yes");
		} */

		if (!charging_times.isEmpty()) process_charging_times_map(charging_times);

	}

	public void process_charging_times_map(Map<Integer, TreeSet<Integer>> charging_times){

		for (Map.Entry<Integer, TreeSet<Integer>> e : charging_times.entrySet()){

			int b = e.getKey();
        	TreeSet<Integer> departures = e.getValue();

			int initial_t = 1;
			double rc = - (this.S[b]-this.S[0]); double min_rc = rc;

			Map<Integer, Double> reducedCostsMap = new LinkedHashMap<>(); reducedCostsMap.put(b, rc);
			Map<Integer, Double> boundsMap = new LinkedHashMap<>();
			
			for (int d: departures){
				if (d <= b) continue; // skip if departure does not allow for sufficient charging

				for (int t=initial_t; t<=d-b-1; t++){
					rc = - (this.S[t+b] - this.S[t]);
					reducedCostsMap.put(t+b, rc);
					if (rc < min_rc - dataModel.precision) min_rc = rc;
				}
				boundsMap.put(d, min_rc);
				initial_t = d-b;
			}
				
			this.charging_reducedCosts.put(b, reducedCostsMap);
			this.charging_bounds.put(b, boundsMap);
			
		}

	}

	private final class PartialSequence {

		public ArrayList<Integer> arcsSequence;
		public boolean[] route;
		public int nominalEnergy;
		public ArrayList<Integer> worstEnergyDevs;
		public int remainingLoad;
		public int cumulativeTime;

		private PartialSequence(ArrayList<Integer> aSeq, boolean[] route, int nomEnergy, ArrayList<Integer> devs, int remLoad, int cumTime){
			this.arcsSequence = aSeq;
			this.route = route;
			this.nominalEnergy = nomEnergy;
			this.worstEnergyDevs = devs;
			this.remainingLoad = remLoad;
			this.cumulativeTime = cumTime;
		}
	}

	public final class FixByReducedCostSolver {

		public EVRPTW dataModel;
		public PricingProblem pricingProblem;
		public Vertex[] vertices; 			//vertices of the instance
		public PriorityQueue<Vertex> nodesToProcess; 			//labels that need be processed
		public long timeLimit;

		
		public ArrayList<ArrayList<Label>> fwLabels = new ArrayList<>();

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

		public ArrayList<ArrayList<PartialSequence>> runForwardLabeling() {

			this.fwLabels = new ArrayList<>();

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
						if(a.head<=dataModel.C+1 && !a.minCostAlternative) continue;
						if(this.pricingProblem.infeasibleArcs[a.id] > 0) continue;
						Label extendedLabel;
						extendedLabel = extendForwardLabel(currentLabel, a);
						if (extendedLabel!=null) { //verifies if the extension is feasible
							updateNodesToProcess(extendedLabel);
						}
					}
				}
			}

			for (int i = 0; i <= dataModel.C+1; i++) { this.fwLabels.add(vertices[i].processedLabels); }

			ArrayList<ArrayList<PartialSequence>> fwSequences = new ArrayList<>();
			for (int i = 0; i <= dataModel.C; i++){
				ArrayList<PartialSequence> allSequences = new ArrayList<>();
				for (Label label: vertices[i].processedLabels){ allSequences.add(get_forward_sequence(label)); }
				fwSequences.add(allSequences);
				//logger.debug("Forward Labels at vertex {}: {}. Forward Sequences: {}", new Object[]{i, this.fwLabels.get(i).size(), fwSequences.get(i).size()});
			}

			long totalTime = System.currentTimeMillis()-startTime;
			dataModel.exactPricingTime+=totalTime;
			if (dataModel.print_log) logger.debug("Time running forward routing labeling algorithm: " + getTimeInSeconds(totalTime)); 

			return fwSequences;
		}

		private PartialSequence get_forward_sequence(Label fwL){

			ArrayList<Integer> aSeq = new ArrayList<>();
			boolean[] route = new boolean[dataModel.C+1];
			
			Label currentLabel = fwL.clone();
			int currentVertex = currentLabel.vertex;
			while(currentVertex!=0) {
				route[currentVertex] = true;

				Arc currentArc = dataModel.arcs[currentLabel.nextArc];
				int nextVertex = currentArc.tail;

				currentLabel = this.fwLabels.get(nextVertex).get(currentLabel.nextLabelIndex);
				if(currentArc.tail>=0 && currentArc.tail<=dataModel.C) aSeq.add(currentArc.id);
				currentVertex = nextVertex;
			}

			ArrayList<Integer> worst_energy_deviations = new ArrayList<>();
			for (int g=0; g<dataModel.gamma; g++){ worst_energy_deviations.add(fwL.remainingEnergy[g] - fwL.remainingEnergy[g+1]); }

			return new PartialSequence(aSeq, route, dataModel.E - fwL.remainingEnergy[0], worst_energy_deviations, fwL.remainingLoad, fwL.remainingTime);
		}

		public Label extendBackwardLabel(Label currentLabel, Arc arc) {

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
			if(source>0 && remainingTime-dataModel.graph.getEdge(0, source).minimumTime<vertices[0].opening_tw) return null;
			if (source>0 && remainingEnergy[dataModel.gamma] - dataModel.graph.getEdge(0, source).minimumEnergy < 0) return null;

			//Check whether the extension is actually feasible
			if(remainingTime<vertices[source].opening_tw || chargingTime>= (int) (remainingTime/10)) return null;

			boolean[] unreachable = Arrays.copyOf(currentLabel.unreachable.clone(), currentLabel.unreachable.length);
			boolean[] ng_path = new boolean[dataModel.C];
			if(source>0) ng_path[source-1] = true;
			else ng_path = Arrays.copyOf(currentLabel.ng_path, currentLabel.ng_path.length);

			//Mark unreachable customers and ng-path cycling restrictions
			if(source>0) {
				
				int lastTail = -1;
				for (Arc c: dataModel.graph.incomingEdgesOf(source)) {
					if(c.tail==lastTail || c.tail==0 || unreachable[c.tail-1]) continue;
					//unreachable
					if (remainingLoad-vertices[c.tail].load<0 || remainingTime-c.minimumTime<vertices[c.tail].opening_tw || 
							remainingEnergy[dataModel.gamma]-c.minimumEnergy<0 || 
							Math.min(remainingTime-c.minimumTime, vertices[c.tail].closing_tw)-dataModel.graph.getEdge(0, c.tail).minimumTime<vertices[0].opening_tw
							|| remainingEnergy[dataModel.gamma]-c.minimumEnergy - dataModel.graph.getEdge(0, c.tail).minimumEnergy<0) {
						unreachable[c.tail-1] = true;
					}
					//ng-path
					if (currentLabel.ng_path[c.tail-1] && vertices[source].neighbors.contains(c.tail)) ng_path[c.tail-1] = true;
					else ng_path[c.tail-1] = false;
					lastTail = c.tail;
				}
			}
			Label extendedLabel = new Label(source, arc.id, currentLabel.index, reducedCost, remainingLoad, remainingTime, remainingEnergy, chargingTime,unreachable, ng_path, eta, srcIndices);
			return extendedLabel;

		}

		public Label extendForwardLabel(Label currentLabel, Arc arc) {

			int head = arc.head; int depot = dataModel.C+1;
			if (head < depot && (currentLabel.unreachable[head-1] || currentLabel.ng_path[head-1])) return null;

			double reducedCost = currentLabel.reducedCost+arc.modifiedCost;
			boolean[] eta = currentLabel.eta.clone();
			HashSet<Integer> srcIndices = new HashSet<Integer>(currentLabel.srcIndices);
			for(int srcIndex: this.pricingProblem.SRCIndices.get(head)) {
				if(currentLabel.eta[srcIndex]) {
					eta[srcIndex] = false;
					int dualIndex = dataModel.C+dataModel.last_charging_period+srcIndex;
					reducedCost-=pricingProblem.dualCosts[dualIndex];
					srcIndices.remove(srcIndex);
				}
				else {eta[srcIndex]=true; srcIndices.add(srcIndex);}
			}
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
				if(remainingTime+dataModel.graph.getEdge(head, depot).minimumTime > vertices[depot].closing_tw) return null;
				if (remainingEnergy[dataModel.gamma] - dataModel.graph.getEdge(head, depot).minimumEnergy < 0) return null;
			}

			//Check whether the extension is actually feasible
			if(remainingTime>vertices[head].closing_tw || chargingTime > dataModel.f_inverse[dataModel.E]) return null;

			boolean[] unreachable = Arrays.copyOf(currentLabel.unreachable.clone(), currentLabel.unreachable.length);
			boolean[] ng_path = new boolean[dataModel.C];
			if (head < depot) ng_path[head-1] = true;

			//Mark unreachable customers and ng-path cycling restrictions
			int lastHead = -1;
			for (Arc c: dataModel.graph.outgoingEdgesOf(head)) {
				if(c.head==lastHead || c.head==depot || unreachable[c.head-1]) continue;
				//unreachable
				if (remainingLoad-vertices[c.head].load<0 || remainingTime+c.minimumTime>vertices[c.head].closing_tw || 
						remainingEnergy[dataModel.gamma]-c.minimumEnergy<0 || 
						Math.max(remainingTime+c.minimumTime, vertices[c.head].opening_tw)+dataModel.graph.getEdge(c.head, depot).minimumTime>vertices[depot].closing_tw
						|| remainingEnergy[dataModel.gamma]-c.minimumEnergy - dataModel.graph.getEdge(0, c.head).minimumEnergy<0) {
					unreachable[c.head-1] = true;
				}
				//ng-path
				if (currentLabel.ng_path[c.head-1] && vertices[head].neighbors.contains(c.head)) ng_path[c.head-1] = true;
				else ng_path[c.head-1] = false;
				lastHead = c.head;
			}
			
			Label extendedLabel = new Label(head, arc.id, currentLabel.index, reducedCost, remainingLoad, remainingTime, remainingEnergy, chargingTime,unreachable, ng_path, eta, srcIndices);
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

			if (L1.vertex>0 && L2.remainingLoad<L1.remainingLoad) return false; 	//load
			if (L2.reducedCost-L1.reducedCost>dataModel.precision) return false; 	//reduced cost
			if (L1.remainingTime<L2.remainingTime) return false; 					//time
			
			for (int gam = 0; gam <= dataModel.gamma; gam ++){ if (L2.remainingEnergy[gam]<L1.remainingEnergy[gam]) return false; } // energy
			
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
			if (L2.reducedCost-L1.reducedCost>dataModel.precision) return false; 	//reduced cost

			
			// Ng-paths and unreachable resources
			Vertex currentVertex = vertices[L1.vertex];
			if (currentVertex.id > 0) {
				for(int i: currentVertex.neighbors) {
					
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
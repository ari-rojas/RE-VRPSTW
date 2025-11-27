package columnGeneration;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Set;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblem;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;

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

	public ArrayList<ArrayList<Label>> fwLabels = new ArrayList<>();
	public ArrayList<ArrayList<Label>> bwLabels = new ArrayList<>();
	public ArrayList<ArrayList<Integer>> SRCIndices = new ArrayList<>();
	public int[] infeasibleArcs;

	public PricingProblem(EVRPTW modelData, String name) {
		super(modelData, name);
		
	}

	public void fixByReducedCosts(long timeLimit){

		FixByReducedCostSolver FRC = new FixByReducedCostSolver(dataModel, timeLimit);
		this.fwLabels = FRC.runForwardLabeling();

		Map<Integer, List<Integer>> mergedMap = new HashMap<>();
		ArrayList<Label> mergedLabels = new ArrayList<>();

		int cont = 0;
		for (Arc arc: dataModel.arcs){
			
			if (arc.head <= dataModel.C+1){ // only routing arcs

				ArrayList<Label> forwardLabels = this.fwLabels.get(arc.tail);
				ArrayList<Label> backwardLabels = this.bwLabels.get(arc.head);

				for (Label fwL: forwardLabels){
					for (Label bwL: backwardLabels){

						Label newLabel = mergeLabel(FRC, fwL, bwL, arc);
						if (newLabel != null){
							newLabel.index = cont;
							mergedLabels.add(newLabel);
							if (mergedMap.containsKey(arc.id)) mergedMap.get(arc.id).add(cont);
							else mergedMap.put(arc.id, new ArrayList<>(cont));
							cont ++;
						}

					}
				}
			}
			
		}

		

	}

	private Label mergeLabel(FixByReducedCostSolver FRC, Label fwL, Label bwL, Arc arc){

		PartialSequence fwSequence = get_forward_sequence(fwL, arc);
		boolean[] fwRoute = fwSequence.route;
		ArrayList<Integer> arcExtensions = fwSequence.arcsSequence;

		for (int i=1; i<=dataModel.C; i++){ // Elementarity assessment
			bwL.unreachable[i-1] = bwL.unreachable[i-1] || bwL.ng_path[i-1];
			if (fwRoute[i] && bwL.unreachable[i-1]) return null;
		}
		
		Arc currentArc = arc;
		if (fwL.remainingLoad + bwL.remainingLoad - dataModel.Q < 0) return null; // Load feasibility
		if (fwL.remainingTime + dataModel.arcs[currentArc.id].time > bwL.remainingTime) return null; // Time feasibility
		
		Label updatedLabel = bwL.clone();
		for (int arcID: arcExtensions){ // Worst-case energy feasibility
			updatedLabel = FRC.extendBackwardLabel(updatedLabel, dataModel.arcs[arcID]);
			if (updatedLabel==null) return null;
		}

		int chargingTime = dataModel.f_inverse[dataModel.E-updatedLabel.remainingEnergy[dataModel.gamma]];
		if (chargingTime > dataModel.f_inverse[dataModel.E] || chargingTime >= (int) updatedLabel.remainingTime/10) return null; // Charging interval feasibility

		// Reduced cost computing (SRCs missing from the label extensions)
		double reducedCost = updatedLabel.reducedCost;
		for (int arcID: arcExtensions){
			
			int source = dataModel.arcs[arcID].tail;
			HashSet<Integer> srcIndices = new HashSet<Integer>(updatedLabel.srcIndices);
			boolean[] eta = updatedLabel.eta.clone();

			for(int srcIndex: dataModel.vertices[source].SRCIndices) {
				if (updatedLabel.eta[srcIndex]) {
					eta[srcIndex] = false;
					int dualIndex = dataModel.C+dataModel.last_charging_period+srcIndex;
					reducedCost -= this.dualCosts[dualIndex];
					srcIndices.remove(srcIndex);
				} else {eta[srcIndex]=true; srcIndices.add(srcIndex);}
			}

			updatedLabel.eta = eta; updatedLabel.srcIndices = srcIndices;
		}
		reducedCost = Math.floor(reducedCost*10000)/10000; updatedLabel.reducedCost = reducedCost;

		return updatedLabel;
	}

	private PartialSequence get_forward_sequence(Label fwL, Arc initialArc){

		ArrayList<Integer> aSeq = new ArrayList<>(initialArc.id);
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

		return new PartialSequence(aSeq, route);
	}

	private final class PartialSequence {

		public ArrayList<Integer> arcsSequence;
		public boolean[] route;

		private PartialSequence(ArrayList<Integer> aSeq, boolean[] route){
			this.arcsSequence = aSeq;
			this.route = route;
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

		public ArrayList<ArrayList<Label>> runForwardLabeling() {

			ArrayList<ArrayList<Label>> fwLabels = new ArrayList<>();

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

			for (int i = 0; i <= dataModel.C+1; i++){ fwLabels.add(vertices[i].processedLabels); }

			long totalTime = System.currentTimeMillis()-startTime;
			dataModel.exactPricingTime+=totalTime;
			if (dataModel.print_log) logger.debug("Time running forward routing labeling algorithm: " + getTimeInSeconds(totalTime)); 

			return fwLabels;
		}

		public Label extendBackwardLabel(Label currentLabel, Arc arc) {

			int source = arc.tail;

			double reducedCost = currentLabel.reducedCost+arc.modifiedCost;
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

			// Stronger check
			if (source>0 && remainingEnergy[dataModel.gamma] - dataModel.graph.getEdge(0, source).minimumEnergy < 0) return null;

			// Mark current customer as unreachable (elementarity)
			boolean[] unreachable = Arrays.copyOf(currentLabel.unreachable.clone(), currentLabel.unreachable.length);
			if(source>0) unreachable[source-1] = true;
			
			Label extendedLabel = new Label(source, arc.id, currentLabel.index, reducedCost, remainingLoad, remainingTime, remainingEnergy, currentLabel.chargingTime+0,unreachable, currentLabel.ng_path.clone(), currentLabel.eta.clone(), new HashSet<Integer>(currentLabel.srcIndices));
			return extendedLabel;

		}

		public Label extendForwardLabel(Label currentLabel, Arc arc) {

			int head = arc.head; int depot = dataModel.C+1;
			if (head < depot && (currentLabel.unreachable[head-1] || currentLabel.ng_path[head-1])) return null;

			double reducedCost = currentLabel.reducedCost+arc.modifiedCost;
			boolean[] eta = currentLabel.eta.clone();
			HashSet<Integer> srcIndices = new HashSet<Integer>(currentLabel.srcIndices);
			for(int srcIndex: vertices[head].SRCIndices) {
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
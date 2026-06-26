package columnGeneration;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Comparator;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.HashMap;
import java.util.HashSet;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblem;

import branchAndPrice.FixArc;
import branchAndPrice.RemoveArc;
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
	public double bestReducedCost; 										//best reduced cost found by the exact labeling
	public double reducedCostThreshold = 0; 							//minimum reduced cost when arriving at the depot source

	public Map<Integer, Map<Integer,Double>> charging_bounds;

	// Information for Fixing by Reduced Costs procedure
	public ArrayList<ArrayList<Label>> bwLabels;
	public ArrayList<ArrayList<PartialBackwardSequence>> bwSequences;
	public ArrayList<ArrayList<PartialForwardSequence>> fwDepotSequences;
	public ArrayList<ArrayList<PartialForwardSequence>> fwC1Sequences;
	public double[] bwBounds;

	public BitSet nonFixablePPArcs;
	public ArrayList<Label> frcRouteLabels;

	public ArrayList<ArrayList<Integer>> SRCIndices = new ArrayList<>();
	public int[] infeasiblePPArcs;
	public PPVertex[] PPvertices = dataModel.PPvertices;
	public Vertex[] vertices = dataModel.vertices;
	public PriorityQueue<PPVertex> nodesToProcess;

	// General information
	private int Gamma = dataModel.gamma;
	private int depotID = dataModel.T_startID;

	// Identifiers for the differnt types of vertices and arcs in the Pricing Problem Routing SubGraph
	public static final byte C0 = EVRPTW.C0; 	 	// Customer depot nodes
	public static final byte C1 = EVRPTW.C1; 		// Non-first customer nodes
	public static final byte AR0 = EVRPTW.AR0;		// Routing arcs between customer depot nodes and non-first customer nodes
	public static final byte AR1 = EVRPTW.AR1;		// Routing arcs between non-first customer nodes

	public double FRC_gap;
	private int problematic_arc = 1271;

	public PricingProblem(EVRPTW modelData, String name) {
		super(modelData, name);
		this.infeasiblePPArcs = new int[dataModel.numArcs];
		this.nonFixablePPArcs = new BitSet();
		this.frcRouteLabels = new ArrayList<Label>();
		this.nodesToProcess = new PriorityQueue<PPVertex>(dataModel.PPvertices.length-dataModel.C, new SortVertices());
	}

	public Map<Integer, Double> fixByReducedCosts(long timeLimit, double UB, double LB){
		
		this.FRC_gap = UB-LB;
		this.bwSequences = new ArrayList<ArrayList<PartialBackwardSequence>>();
		this.fwDepotSequences = new ArrayList<ArrayList<PartialForwardSequence>>();
		this.fwC1Sequences = new ArrayList<ArrayList<PartialForwardSequence>>();

		Map<Integer, Double> arcsToRemove = new HashMap<Integer, Double>();

		long startTime = System.currentTimeMillis();
		for (Label label: this.frcRouteLabels){

			Label nextLabel = label;
			while(nextLabel.vertex != depotID) {
				this.nonFixablePPArcs.set(nextLabel.nextArc);
				PPArc nextArc = dataModel.PParcs[nextLabel.nextArc];
				int j = PPvertices[nextArc.head_vertex_id].node_number;
				nextLabel = this.bwLabels.get(j).get(nextLabel.nextLabelIndex);
			}
		}

		long totalTime = System.currentTimeMillis()-startTime;
		if (dataModel.print_log) {
			logger.debug("Identified " + this.nonFixablePPArcs.cardinality()+"/"+(dataModel.lenAR0+dataModel.lenAR1) + " non-fixable PP Arcs");
			logger.debug("Time identifying first non-fixable PP Arcs: " + getTimeInSeconds(totalTime));
		}

		this.cleanBackwardLabels();
		this.runForwardLabeling(timeLimit);
		
		startTime = System.currentTimeMillis();
		for (int j = 1; j <= dataModel.C+1; j++){
			ArrayList<PartialBackwardSequence> backwardSequences = bwSequences.get(j);

			for (PPArc arc: dataModel.PPgraph.incomingEdgesOf(PPvertices[dataModel.C1_startID+j].id)){

				if (System.currentTimeMillis()>timeLimit) break;
				if (infeasiblePPArcs[arc.id] > 0 || nonFixablePPArcs.get(arc.id)) continue;
				
				ArrayList<PartialForwardSequence> forwardSequences;
				int i = PPvertices[arc.tail_vertex_id].node_number;
				if (arc.arc_type == AR0) forwardSequences = this.fwDepotSequences.get(i);
				else forwardSequences = this.fwC1Sequences.get(i);

				//if (arc.id == this.problematic_arc) for (PartialForwardSequence fww: forwardSequences) logger.debug(fww.toString());
				
				//if (arc.id == this.problematic_arc) for (PartialForwardSequence ffw: forwardSequences) logger.debug(ffw.reducedCost+" "+ffw.routingArcsSequence.toString());
				
				double min_rc = findMinimumRCPath_acc(backwardSequences, forwardSequences, arc.routing_arc, arc.modifiedCost, arc.id);

				if (min_rc - bestReducedCost - dataModel.precision > FRC_gap) arcsToRemove.put(arc.id, min_rc);
				if (min_rc < bestReducedCost - 1) logger.debug("!!! Arc {} has a merged label with a reduced cost of {}", new Object[]{arc.toString(), min_rc});

			}

			backwardSequences.clear();
		}

		this.bwSequences.clear(); this.fwC1Sequences.clear(); this.fwDepotSequences.clear();

		totalTime = System.currentTimeMillis()-startTime;
		dataModel.exactPricingTime+=totalTime;
		if (dataModel.print_log) {
			logger.debug("Time merging forward and backward labels: " + getTimeInSeconds(totalTime));
		}

		return arcsToRemove;

	}

	private double findMinimumRCPath_acc(ArrayList<PartialBackwardSequence> bwSequences, ArrayList<PartialForwardSequence> fwSequences, Arc arc, double modifiedCost, int ppid) {

		int nFw = fwSequences.size(); int nBw = bwSequences.size();
		PriorityQueue<MergeState> pq = new PriorityQueue<>( (s1, s2) -> Double.compare(s1.rc, s2.rc) );
		
		for (int ixFw = 0; ixFw < nFw; ixFw++){
			PartialForwardSequence fwSeq = fwSequences.get(ixFw);

			//fwSeq.routingArcsSequence.equals(new ArrayList<>(List.of(551, 611, 583, 627, 675, 43)));

			for (int ixBw = 0; ixBw < nBw; ixBw++){
				PartialBackwardSequence bwSeq = bwSequences.get(ixBw);
				
				double route_rc = fwSeq.reducedCost + modifiedCost + bwSeq.reducedCost;
				route_rc = Math.floor(route_rc*10000)/10000;
				if (route_rc - 1 - bestReducedCost > this.FRC_gap) continue;
				
				if (bwSeq.ng.intersects(fwSeq.ng)) continue;										// ng-Elementarity
				if (bwSeq.remainingTime - arc.time < fwSeq.cumulativeTime) continue; 				// Time feasibility
				if (bwSeq.worstRemainEnergy - arc.energy <  fwSeq.nominalEnergy) continue; 			// Worst-case Energy of the backwards - rest of nominal energy
				if (bwSeq.remainingLoad < fwSeq.cumulativeLoad) continue; 							// Load feasibility
        		
				route_rc += getMergeSRCs_RC(fwSeq.eta, bwSeq.eta);
				route_rc = Math.floor(route_rc*10000)/10000;
				if (route_rc - 1 - bestReducedCost > this.FRC_gap) continue;
				pq.add(new MergeState(ixFw, ixBw, Math.floor(route_rc*10000)/10000));
			}
		}

		double min_merged_rc = this.FRC_gap+1+2*dataModel.precision;
        while (!pq.isEmpty()) {
			
            MergeState current = pq.poll();
            
			int fw = current.f; PartialForwardSequence fwSeq = fwSequences.get(fw);
			int bw = current.b; PartialBackwardSequence bwSeq = bwSequences.get(bw);

			if (min_merged_rc <= current.rc + dataModel.precision){ return min_merged_rc; }
			else {
				MergedSequence mergedPath = mergeLabel_acc(fwSeq, bwSeq, arc, current.rc, ppid);
				if (mergedPath != null){ // If found a feasible merged label
					double chBound = this.charging_bounds.get(mergedPath.chargingTime).get(mergedPath.departureTime);
					double complete_rc = mergedPath.reducedCost + chBound;
					
					if (complete_rc < min_merged_rc - dataModel.precision){
						min_merged_rc = complete_rc;
						// If found a feasible column with lower RC than the gap, the arc won't be fixed
						// If the charging bound is 0, the column's reduced cost is optimal for the FRC expression
						if (min_merged_rc - bestReducedCost <= this.FRC_gap + dataModel.precision || chBound <= dataModel.precision)  return min_merged_rc; 
					}
				}
			}
        }

        return 1e5;
    }

	private double getMergeSRCs_RC(boolean[] etaFw, boolean[] etaBw){

		double additional_rc = 0;
		for (int ix = 0; ix < etaFw.length; ix++){
			if (etaFw[ix] && etaBw[ix]) additional_rc -= this.dualCosts[dataModel.C+dataModel.last_charging_period+ix];
		}

		return additional_rc;
	}

	private MergedSequence mergeLabel_acc(PartialForwardSequence fwSequence, PartialBackwardSequence bwSeq, Arc routing_arc, double reducedCost, int ppid){
		
		///////////////////////////////////
		/// MERGE FEASIBILITY ASSESSMENT
		///////////////////////////////////
		
		// Worst-case energy feasibility
		int remainingEnergy = bwSeq.nominalEnergy - routing_arc.energy - fwSequence.nominalEnergy; // Nominal energy consumption
		
		ArrayList<Integer> fwDevs = fwSequence.worstEnergyDevs;
		ArrayList<Integer> bwDevs = new ArrayList<>(bwSeq.worstEnergyDevs);
		for (int g = 0; g < Gamma; g++) if (routing_arc.energy_deviation >= bwDevs.get(g)) { bwDevs.add(g, routing_arc.energy_deviation); break; }

		int ixFw = 0; int ixBw = 0;
		for (int g = 1; g <= Gamma; g++){
			if (fwDevs.get(ixFw) >= bwDevs.get(ixBw)) { remainingEnergy -= fwDevs.get(ixFw); ixFw ++; }
			else { remainingEnergy -= bwDevs.get(ixBw); ixBw ++; }
			if (remainingEnergy < 0) return null;
		}
		
		int chargingTime = dataModel.f_inverse[dataModel.E-remainingEnergy];
		
		/////////////////////////////////
		/// LATEST DEPARTURE TIME
		/////////////////////////////////
		
		int source = routing_arc.tail;
		int remainingTime = bwSeq.remainingTime - routing_arc.time;
		if (remainingTime > vertices[source].closing_tw) remainingTime = vertices[source].closing_tw;
		
		ArrayList<Integer> fwRoutingArcs = fwSequence.routingArcsSequence;
		for (int arcID: fwRoutingArcs){ // Loop over the arc extensions leading up to a C0 vertex
			Arc routeArc = dataModel.arcs[arcID];
			source = routeArc.tail;
			
			remainingTime -= routeArc.time;
			if (remainingTime > vertices[source].closing_tw) remainingTime = vertices[source].closing_tw;
		}
		
		int departure = (int)(remainingTime/10);
		if (chargingTime >= departure) return null; 	// Charging interval feasibility

		return new MergedSequence(reducedCost, chargingTime, departure);
	}

	private void cleanBackwardLabels() {

		//////////////////////////////////
		/// C1 vertices labels
		//////////////////////////////////
		
		this.bwSequences.add(null);

		for (int i = 1; i <= dataModel.C+1; i++){
			ArrayList<Label> labels = this.bwLabels.get(i); ArrayList<Label> labels_to_remove = new ArrayList<Label>();
			for (int ix = 0; ix < labels.size(); ix++){
				Label l1 = labels.get(ix); boolean dominated = false;
				for (int ix2 = ix+1; ix2 < labels.size(); ix2++) if (isDominatedBackwardRouting(l1, labels.get(ix2))) { dominated = true; break; }
				if (dominated) labels_to_remove.add(l1);
			} labels.removeAll(labels_to_remove);
			labels.sort(Comparator.comparingDouble(l -> l.reducedCost));

			ArrayList<PartialBackwardSequence> allSequences = new ArrayList<PartialBackwardSequence>();
			for (Label label: labels) allSequences.add(new PartialBackwardSequence(label.reducedCost, label.remainingEnergy, label.remainingLoad, label.remainingTime, label.ng_path, label.eta));
			this.bwSequences.add(allSequences);
		}

		this.bwLabels.clear();

		this.bwBounds = new double[2*dataModel.C+2];
		for (int i = 1; i <= dataModel.C; i++){ // Bound for depot-customer nodes
			double min_label_rc = Double.MAX_VALUE;
			for (PPArc arc: dataModel.PPgraph.outgoingEdgesOf(dataModel.C0_startID+i)){
				if (infeasiblePPArcs[arc.id] > 0) continue;
				int j = PPvertices[arc.head_vertex_id].node_number;
				if (j == 0) j = dataModel.C+1;
				double rc = this.bwSequences.get(j).get(0).reducedCost + arc.modifiedCost;
				if (rc < min_label_rc - dataModel.precision) min_label_rc = rc;
			}
			this.bwBounds[i] = min_label_rc;
		}

		for (int i = 1; i <= dataModel.C; i++){ // Bound for non-first customer nodes
			double min_label_rc = Double.MAX_VALUE;
			for (PPArc arc: dataModel.PPgraph.outgoingEdgesOf(dataModel.C1_startID+i)){
				if (infeasiblePPArcs[arc.id] > 0) continue;
				int j = PPvertices[arc.head_vertex_id].node_number;
				if (j == 0) j = dataModel.C+1;
				double rc = this.bwSequences.get(j).get(0).reducedCost + arc.modifiedCost;
				if (rc < min_label_rc - dataModel.precision) min_label_rc = rc;
			}
			this.bwBounds[dataModel.C+i] = min_label_rc;
		}

		this.bwBounds[2*dataModel.C+1] = 0; // Bound for returning depot node

	}

	private ArrayList<PartialBackwardSequence> get_specific_bwSequence(){

		ArrayList<PartialBackwardSequence> lista = new ArrayList<>();
		for (PartialBackwardSequence bwS: this.bwSequences.get(1)){
			if (Math.abs(bwS.reducedCost-1671)<1e-2) lista.add(bwS);
		}

		return lista;
	}

	private ArrayList<Integer> get_backward_arcs_sequence(Label bwL){

		ArrayList<Integer> aSeq = new ArrayList<>();

		Label currentLabel = bwL.clone();
		while (true) {
			PPArc nextArc = dataModel.PParcs[currentLabel.nextArc];
			aSeq.add(nextArc.routing_arc.id);
			if (nextArc.head_vertex_id == depotID) break;
			currentLabel = this.bwLabels.get(PPvertices[nextArc.head_vertex_id].node_number).get(currentLabel.nextLabelIndex);
		}

		return aSeq;
	}

	//////////////////////////////////////////////////
	/// FORWARD LABELING
	//////////////////////////////////////////////////
	
	public void runForwardLabeling(long timeLimit) {

		// Initialization
		int[] remain_energy = new int[dataModel.gamma + 1]; Arrays.fill(remain_energy, dataModel.E);
		ForwardLabel l0 = new ForwardLabel(0, 0, 0, 0, vertices[0].opening_tw, vertices[0].closing_tw, 0, remain_energy, 0, new boolean[dataModel.C], new boolean[dataModel.C], new boolean[this.subsetRowCuts.size()], new HashSet<Integer>(this.subsetRowCuts.size()) );
		for (int i = 1; i <= dataModel.C; i++){
			PPVertex depotVertex = PPvertices[dataModel.C0_startID+i];
			Arc routingArc = dataModel.graph.getEdge(0,i);
			
			ForwardLabel initialLabel = extendForwardLabel(l0, routingArc, 0, depotVertex.id);
			if (initialLabel != null){
				initialLabel.vertex = depotVertex.id;
				depotVertex.unprocessedForwardLabels.add(initialLabel);
				this.nodesToProcess.add(depotVertex);
			}
		}

		////////////////////////////////////////////
		/// Routing Labeling
		////////////////////////////////////////////
		
		long startTime = System.currentTimeMillis();
		while (!nodesToProcess.isEmpty() && System.currentTimeMillis()<timeLimit) {
			ArrayList<ForwardLabel> labelsToProcessNext = routingLabelsToProcessNext();
			Set<PPArc> outgoingArcs = new HashSet<PPArc>(dataModel.PPgraph.outgoingEdgesOf(labelsToProcessNext.get(0).vertex));
			outgoingArcs.removeIf(arc -> infeasiblePPArcs[arc.id] > 0 || arc.head_vertex_id == depotID); // No need to extend to the returning depot

			for (ForwardLabel currentLabel: labelsToProcessNext) {
				
				boolean isDominated = checkDominance(currentLabel);
				if (isDominated) continue;
				
				currentLabel.index = PPvertices[currentLabel.vertex].processedForwardLabels.size();
				PPvertices[currentLabel.vertex].processedForwardLabels.add(currentLabel);
				
				for (PPArc a: outgoingArcs) {
					ForwardLabel extendedLabel = extendForwardLabel(currentLabel, a.routing_arc, a.modifiedCost, a.head_vertex_id);
					if (extendedLabel!=null) { // Verifies if the extension is feasible
						extendedLabel.vertex = a.head_vertex_id;
						extendedLabel.previousArc = a.id;
						
						PPvertices[extendedLabel.vertex].unprocessedForwardLabels.add(extendedLabel);
						if (PPvertices[a.head_vertex_id].unprocessedForwardLabels.size() == 1) nodesToProcess.add(PPvertices[a.head_vertex_id]);
						
					}
				}
			}
		}

		///////////////////////////////
		/// Labels cleanse
		///////////////////////////////

		if (System.currentTimeMillis()<timeLimit){
			
			this.fwDepotSequences.add(null);
			for (int i = 1; i <= dataModel.C; i++) {
				ArrayList<PartialForwardSequence> allSequences = new ArrayList<PartialForwardSequence>();
				if (PPvertices[dataModel.C0_startID+i].processedForwardLabels.size() > 0){
					ForwardLabel label = PPvertices[dataModel.C0_startID+i].processedForwardLabels.get(0);
					allSequences.add(get_forward_sequence(label));
				}
				this.fwDepotSequences.add(allSequences);
			}
			
			this.fwC1Sequences.add(null);
			for (int i = 1; i <= dataModel.C; i++){
				ArrayList<ForwardLabel> labels = new ArrayList<ForwardLabel>(PPvertices[dataModel.C1_startID+i].processedForwardLabels);
				ArrayList<ForwardLabel> labels_to_remove = new ArrayList<>();
				for (int ix = 0; ix < labels.size(); ix++){
					ForwardLabel l1 = labels.get(ix);
					for (int ix2 = ix+1; ix2 < labels.size(); ix2++) if (isDominatedRouting(l1, labels.get(ix2))) {
						labels_to_remove.add(l1); break; }
				} labels.removeAll(labels_to_remove);

				ArrayList<PartialForwardSequence> allSequences = new ArrayList<PartialForwardSequence>();
				for (ForwardLabel label: labels){
					PartialForwardSequence debugFW = get_forward_sequence(label);
					allSequences.add(debugFW);
				}
				//allSequences.sort( Comparator.comparing(l -> l.reducedCost) );
				this.fwC1Sequences.add(allSequences);
			}

			for (int i = 0; i < PPvertices.length; i++) {
				PPvertices[i].processedForwardLabels = new ArrayList<ForwardLabel>(dataModel.numArcs);
				PPvertices[i].unprocessedForwardLabels =  new PriorityQueue<ForwardLabel>(dataModel.numArcs, new ForwardLabel.SortForwardLabels()); }

			long totalTime = System.currentTimeMillis()-startTime;
			dataModel.exactPricingTime+=totalTime;
			if (dataModel.print_log) logger.debug("Time running forward routing labeling algorithm: " + getTimeInSeconds(totalTime));
				
		} else {
			if (dataModel.print_log) logger.debug("Caught timeout while running the forward labeling algorithm");
		}

	}

	private PartialForwardSequence get_forward_sequence(ForwardLabel fwL){

		ArrayList<Integer> aSeq = new ArrayList<>();
		
		ForwardLabel currentLabel = fwL.clone();
		while(PPvertices[currentLabel.vertex].vertex_type > C0) {
			PPArc previousArc = dataModel.PParcs[currentLabel.previousArc];
			aSeq.add(previousArc.routing_arc.id);
			currentLabel = PPvertices[previousArc.tail_vertex_id].processedForwardLabels.get(currentLabel.previousLabelIndex);
		}
		aSeq.add(dataModel.graph.getEdge(0,PPvertices[currentLabel.vertex].node_number).id);

		ArrayList<Integer> worst_energy_deviations = new ArrayList<>();
		for (int g = 0; g < Gamma; g++){ worst_energy_deviations.add(fwL.remainingEnergy[g] - fwL.remainingEnergy[g+1]); }

		BitSet ng = new BitSet();
		for (int i = 0; i < dataModel.C; i++) if (fwL.ng_path[i]) ng.set(i);

		return new PartialForwardSequence(fwL.reducedCost, aSeq, dataModel.E - fwL.remainingEnergy[0], worst_energy_deviations, fwL.cumulativeLoad, fwL.cumulativeTime, ng, fwL.eta);
	}

	public ArrayList<ForwardLabel> routingLabelsToProcessNext(){

		ArrayList<ForwardLabel> labelsToProcessNext = new ArrayList<ForwardLabel>();
		PPVertex currentVertex = nodesToProcess.poll();
		
		while(true) {
			ForwardLabel currentLabel = currentVertex.unprocessedForwardLabels.poll();
			if(labelsToProcessNext.isEmpty()) labelsToProcessNext.add(currentLabel);
			else {
				boolean isDominated = false;
				for(ForwardLabel L2: labelsToProcessNext) {
					
					isDominated = isDominatedRouting(currentLabel, L2);
					if(isDominated) break;
				}
				if (!isDominated) labelsToProcessNext.add(currentLabel);
			}
			if(currentVertex.unprocessedForwardLabels.isEmpty() || (currentVertex.unprocessedForwardLabels.peek().cumulativeLoad>currentLabel.cumulativeLoad)) break;
		}

		if(!currentVertex.unprocessedForwardLabels.isEmpty()) nodesToProcess.add(currentVertex);
		return labelsToProcessNext;
	}

	public boolean checkDominance(ForwardLabel newLabel) {
		
		PPVertex currentVertex = PPvertices[newLabel.vertex];

		ArrayList<ForwardLabel> labelsToDelete = new ArrayList<ForwardLabel>();
		for(ForwardLabel existingLabel: currentVertex.unprocessedForwardLabels) {
			boolean isDominated = isDominatedRouting(existingLabel, newLabel);
			if (isDominated) {
				labelsToDelete.add(existingLabel);
			}
		}
		currentVertex.unprocessedForwardLabels.removeAll(labelsToDelete);
		if(currentVertex.unprocessedForwardLabels.isEmpty()) nodesToProcess.remove(currentVertex);

		for(ForwardLabel existingLabel: currentVertex.processedForwardLabels) if(isDominatedRouting(newLabel, existingLabel)) return true;

		return false;
	}

	public ForwardLabel extendForwardLabel(ForwardLabel currentLabel, Arc routing_arc, double modifiedCost, int pp_head_ix) {

		int head = routing_arc.head;
		if (currentLabel.unreachable[head-1] || currentLabel.ng_path[head-1]) return null;

		// Update the remaining time and check feasibility
		int cumulativeTime = currentLabel.cumulativeTime+routing_arc.time;
		if (cumulativeTime<vertices[head].opening_tw) cumulativeTime = vertices[head].opening_tw;
		if (cumulativeTime>vertices[head].closing_tw) return null;
		
		int travelTimes = currentLabel.travelTimes+routing_arc.time;
		int latestDeparture = currentLabel.latestDeparture;
		if ((int)((vertices[head].closing_tw - travelTimes)/10) < latestDeparture) latestDeparture = (int)((vertices[head].closing_tw - travelTimes)/10);
		
		double reducedCost = currentLabel.reducedCost+modifiedCost;
		boolean[] eta = currentLabel.eta.clone();
		HashSet<Integer> srcIndices = new HashSet<Integer>(currentLabel.srcIndices);
		for(int srcIndex: this.SRCIndices.get(head)) {
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

		if (remainingEnergy[Gamma]-dataModel.graph.getEdge(head, dataModel.C+1).min_energy < 0) return null;

		// Charging Time
		int minEnergyRoute_remEn = remainingEnergy[0] - dataModel.graph.getEdge(head, dataModel.C+1).min_energy; // Nominal Energy
		int[] bwDevs = vertices[head].minEnergy_DepotPath_Devs;
		int ixFw = 1; int ixBw = 0;
		for (int g = 1; g <= Gamma; g++){
			if (remainingEnergy[ixFw-1]-remainingEnergy[ixFw] >= bwDevs[ixBw]) { minEnergyRoute_remEn -= (remainingEnergy[ixFw-1]-remainingEnergy[ixFw]); ixFw ++; }
			else { minEnergyRoute_remEn -= bwDevs[ixBw]; ixBw ++; }
		} if (minEnergyRoute_remEn < 0) return null;
		
		int chargingTime = dataModel.f_inverse[dataModel.E-minEnergyRoute_remEn];
		if (chargingTime >= latestDeparture) return null;
		double chargingBound = Math.floor(this.charging_bounds.get(chargingTime).get(latestDeparture)*10000)/10000;
		
		/////////////////////////////////////////
		/// Bounding Procedure
		/////////////////////////////////////////
		
		if (reducedCost + this.bwBounds[pp_head_ix] + chargingBound > this.FRC_gap + dataModel.precision) return null;

		// After confirming that the label is feasible, update the remaining load
		int cumulativeLoad = currentLabel.cumulativeLoad+vertices[head].load;
		
		// Unreachable resources
		boolean[] unreachable = Arrays.copyOf(currentLabel.unreachable.clone(), currentLabel.unreachable.length);
		boolean[] ng_path = new boolean[dataModel.C];
		
		// Mark unreachable customers and ng-path cycling restrictions
		ng_path[head-1] = true;
		for (Arc c: dataModel.graph.outgoingEdgesOf(head)) {
			if (c.head == dataModel.C+1 || unreachable[c.head-1]) continue;
			//unreachable
			if (cumulativeLoad+vertices[c.head].load>dataModel.Q || cumulativeTime+c.min_time>vertices[c.head].closing_tw || remainingEnergy[Gamma] - c.min_energy - dataModel.graph.getEdge(c.head, dataModel.C+1).min_energy < 0) {
				unreachable[c.head-1] = true; }
				
			//ng-path
			if (currentLabel.ng_path[c.head-1] && vertices[head].neighbors.contains(c.head)) ng_path[c.head-1] = true;
			else ng_path[c.head-1] = false;
		}
			
		ForwardLabel extendedLabel = new ForwardLabel(currentLabel.index, reducedCost, chargingBound, cumulativeLoad, cumulativeTime, latestDeparture, travelTimes, remainingEnergy, chargingTime, unreachable, ng_path, eta, srcIndices);

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

	public boolean isDominatedRouting(ForwardLabel L1, ForwardLabel L2) {

		if (L2.cumulativeLoad>L1.cumulativeLoad) return false; 													//load
		if (L2.reducedCost-L1.reducedCost>dataModel.precision) return false; 	//reduced cost
		if (L2.cumulativeTime>L1.cumulativeTime) return false;
		if (L2.latestDeparture<L1.latestDeparture) return false;
		if (L2.travelTimes>L1.travelTimes) return false; 													//time
		
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

	public boolean isDominatedBackwardRouting(Label L1, Label L2) {

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

	///////////////////////////////////////////////////
	/// FRC CLASSES
	///////////////////////////////////////////////////
	
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

	private final class PartialForwardSequence {

		public final double reducedCost;
		public final ArrayList<Integer> routingArcsSequence;
		public final int nominalEnergy;
		public final ArrayList<Integer> worstEnergyDevs;
		public final int cumulativeLoad;
		public final int cumulativeTime;
		public final BitSet ng;
		public final boolean[] eta;

		private PartialForwardSequence(double rc, ArrayList<Integer> aSeq, int nomEnergy, ArrayList<Integer> devs, int cumLoad, int cumTime, BitSet ng, boolean[] eta){
			this.reducedCost = rc;
			this.routingArcsSequence = aSeq;
			this.nominalEnergy = nomEnergy;
			this.worstEnergyDevs = devs;
			this.cumulativeLoad = cumLoad; //iykyk
			this.cumulativeTime = cumTime;
			this.ng = ng;
			this.eta = eta;
		}

		@Override
		public String toString(){
			return "Forward Label. Reduced Cost: "+reducedCost+" Route Arcs: "+routingArcsSequence.toString()+" Nominal Energy: "+nominalEnergy+ "Worst Energy Devs: "+worstEnergyDevs.toString()+" Cumulative Load: "+cumulativeLoad+" Cumulative Time: "+cumulativeTime+" Ng-Elementarity: "+ng.toString()+" Eta SRCs: "+Arrays.toString(eta);
		}
	}

	private final class PartialBackwardSequence {

		public final double reducedCost;
		public final int nominalEnergy;
		public final int worstRemainEnergy;
		public final ArrayList<Integer> worstEnergyDevs;
		public final int remainingLoad;
		public final int remainingTime;
		public final BitSet ng;
		public final boolean[] eta;

		private PartialBackwardSequence(double rc, int[] remEnergy, int remLoad, int remTime, boolean[] ngp, boolean[] eta){
			this.reducedCost = rc;
			this.remainingLoad = remLoad;
			this.remainingTime = remTime;
			
			this.nominalEnergy = remEnergy[0];
			this.worstRemainEnergy = remEnergy[Gamma];
			this.worstEnergyDevs = new ArrayList<Integer>();
			for (int g = 0; g < Gamma; g++){ worstEnergyDevs.add(remEnergy[g] - remEnergy[g+1]); }

			this.ng = new BitSet();
			for (int i = 0; i < dataModel.C; i++)  if (ngp[i]) this.ng.set(i);

			this.eta = eta;
		}

		@Override
		public String toString(){
			return "Backward Label. Reduced Cost: "+reducedCost+" Remaining Load: "+remainingLoad+" Remaining Time: "+remainingTime+" Remaining Nominal Energy: "+nominalEnergy+" Remaining Worst Energy: "+worstRemainEnergy+" Worst Energy Devs: "+worstEnergyDevs.toString()+" Ng-Elementarity: "+ng.toString()+" Eta SRCs: "+Arrays.toString(eta);
		}
	}

	private static final class MergeState {
        final int f, b;
        final double rc;
        MergeState(int f, int b, double rc) { this.f = f; this.b = b; this.rc = rc; }
    }

	///////////////////////////////////////////////////
	/// DON'T TOUCH
	///////////////////////////////////////////////////

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

	public double getTimeInSeconds(double time) {
		double realTime = time*0.001;
		realTime = Math.floor(realTime*100)/100; //two decimals
		return realTime;
	}

	/**
	 * Listen to branching decisions. The pricing problem is changed by the branching decisions.
	 * @param bd BranchingDecision
	 */
	@Override
	public void branchingDecisionPerformed(BranchingDecision bd) {
		if(bd instanceof FixArc) { 			//Fixing one arc
			FixArc fixArcDecision = (FixArc) bd;
			for(int infeasibleArc: fixArcDecision.infeasiblePPArcs) this.infeasiblePPArcs[infeasibleArc] ++;
		}else if(bd instanceof RemoveArc) {//Removing one arc
			RemoveArc removeArcDecision= (RemoveArc) bd;
			this.infeasiblePPArcs[removeArcDecision.arcID] ++;
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
			for(int infeasibleArc: fixArcDecision.infeasiblePPArcs) this.infeasiblePPArcs[infeasibleArc] --;
		}else if(bd instanceof RemoveArc) {//Removing one arc
			RemoveArc removeArcDecision= (RemoveArc) bd;
			this.infeasiblePPArcs[removeArcDecision.arcID] --;
		}
	}

	public class SortVertices implements Comparator<PPVertex> {

		@Override
		public int compare(PPVertex vertex1, PPVertex vertex2) {
			
			ForwardLabel L1 = vertex1.unprocessedForwardLabels.peek();
			ForwardLabel L2 = vertex2.unprocessedForwardLabels.peek();
			
			if (vertex1.vertex_type == C0 && vertex2.vertex_type == C0){
				if (L1.reducedCost<L2.reducedCost) return -1;
				else return 1;}
			if (vertex1.vertex_type == C0) return -1;
			if (vertex2.vertex_type == C0) return 1;

			// Choose according the current unprocessed labels
			if(L1.cumulativeLoad<L2.cumulativeLoad) return -1;
			if(L1.cumulativeLoad>L2.cumulativeLoad) return 1;
			if(L1.remainingEnergy[Gamma]>L2.remainingEnergy[Gamma]) return -1;
			if(L1.remainingEnergy[Gamma]<L2.remainingEnergy[Gamma]) return 1;
			if(L1.cumulativeTime<L2.cumulativeTime) return -1;
			if(L1.cumulativeTime>L2.cumulativeTime) return 1;
			if(L1.reducedCost+L1.chargingBound<L2.reducedCost+L2.chargingBound) return -1;
			if(L1.reducedCost+L1.chargingBound>L2.reducedCost+L2.chargingBound) return 1;
			
			return 0;
		}
	}

}
package columnGeneration;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.TreeSet;

import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblem;
import branchAndPrice.ChargingTimeInequality;
import model.EVRPTW;
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

	//Charging pricing information
	private BitSet negative_charging_duals;
	private Map<Integer, Map<Integer, Double>> charging_bounds;
	public Map<Integer, Map<Integer, Double>> charging_reducedCosts;

	private HashMap<Integer, Integer> nonDominatedT;
	public Map<Integer, BitSet> last_charging_periods;

	public PricingProblem(EVRPTW modelData, String name) {
		super(modelData, name);
	}

	public void compute_charging_bounds(){

		this.charging_reducedCosts = new HashMap<>();
		this.charging_bounds = new HashMap<>();
		this.negative_charging_duals = new BitSet();

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
		for (int t = 1; t < maxT; t++) {
			double dual = this.dualCosts[dataModel.C + t - 1];
			this.negative_charging_duals.set(t, dual < -dataModel.precision);
			S[t] = S[t - 1] + dual;
		}

		for (Map.Entry<Integer, TreeSet<Integer>> e : charging_times.entrySet()){

			int b = e.getKey();
        	TreeSet<Integer> departures = e.getValue();

			int initial_t = 1;
			double rc = - (S[b]-S[0]); double min_rc = rc;

			Map<Integer, Double> reducedCostsMap = new LinkedHashMap<>(); reducedCostsMap.put(b, rc);
			Map<Integer, Double> boundsMap = new LinkedHashMap<>();
			
			for (int d: departures){
				if (d <= b) continue; // skip if departure does not allow for sufficient charging

				for (int t=initial_t; t<=d-b-1; t++){
					rc = - (S[t+b] - S[t]);
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

	public ArrayList<Label> charging_pricing_filtering(ArrayList<Label> labels){

		this.nonDominatedT = new HashMap<>();

		//////////////////////////////////////////////////////////
		/// 1. Bounding Procedure
		//////////////////////////////////////////////////////////

		// To avoid constantly recomputing the departure times of the labels, we save them in the vertex field
		for (Label label: labels) label.vertex = (int)(label.remainingTime/10);

		// Only labels that will generate at least one column with negative reduced cost are accounted for
		ArrayList<Label> filtered_labels = new ArrayList<>();
		for (Label l: labels){ if (l.reducedCost + this.charging_bounds.get(l.chargingTime).get(l.vertex) < -dataModel.precision) filtered_labels.add(l); }
		
		//////////////////////////////////////////////////////////
		/// 2. Labels dominance
		//////////////////////////////////////////////////////////

		// 1. Group labels by charging time b
		Map<Integer, List<Label>> labelsByB = new HashMap<>();
		for (Label label : filtered_labels) labelsByB.computeIfAbsent(label.chargingTime, k -> new ArrayList<>()).add(label);

		// 2. Dominance between labels of same chargingTime b
		// For each chargingTime b, get the last-charging-time-periods that have a non-dominated column, and their corresponding label
		Map<Integer, BitSet> columnsIndicator = new HashMap<>();
		Map<Long, Label> columnsMap = new HashMap<>(); BitSet t_set = new BitSet();
		for (Map.Entry<Integer, List<Label>> entry : labelsByB.entrySet()) { filter_labels_same_chargingTime(entry, t_set, columnsMap, columnsIndicator); }

		// 3. Dominance between columns of same lastChargingTimePeriod t and different chargingTime b
		/* for (int t = t_set.nextSetBit(0); t >= 1; t = t_set.nextSetBit(t + 1)){
			BitSet colsIndicator = columnsIndicator.get(t);
			for (int b = colsIndicator.nextSetBit(1); b >= 1; b = colsIndicator.nextSetBit(b + 1)){
				
				Label currentLabel = columnsMap.get(pack(b,t));
				double current_rc = currentLabel.reducedCost;

				for (int b2 = colsIndicator.nextSetBit(b + 1); b2 >= 1; b2 = colsIndicator.nextSetBit(b2+1)){
					Label otherLabel = columnsMap.get(pack(b2,t));
					if (current_rc < otherLabel.reducedCost - dataModel.precision){
						colsIndicator.clear(b2);
					}
				}
			}
		} */

		// 4. Dominance between labels of different lastChargingTimePeriod t and chargingTime b
		for (int t = t_set.previousSetBit(t_set.length()-1); t >= 1; t = t_set.previousSetBit(t - 1)){
			// We get the FULL BitSet, to evaluate the dominance of each column, even if it has already been deemed dominated
			BitSet colsIndicator = columnsIndicator.get(t);

			for (int b = colsIndicator.previousSetBit(colsIndicator.length() - 1); b >= 1; b = colsIndicator.previousSetBit(b - 1)){

				Label currentLabel = columnsMap.get(pack(b,t));
				double current_rc = currentLabel.reducedCost;

				// At t, currentLabel can be dominated by columns with b2 < b
				boolean dominated = false;
				for (int b2 = colsIndicator.previousSetBit(b - 1); b2 >= 1; b2 = colsIndicator.previousSetBit(b2 - 1)){
					Label otherLabel = columnsMap.get(pack(b2, t));
					if (otherLabel.reducedCost < current_rc - dataModel.precision){ // currentLabel is dominated
						colsIndicator.clear(b);
						dominated = true; break;
					}
				}

				if (dominated) continue;

				int previous_t = t;
				for (int t2 = t_set.previousSetBit(t - 1); t2 >= t-b+1 ; t2 = t_set.previousSetBit(t2 - 1)){

					// Update the reduced costs (as the labels are "extended")
					for (int tt = previous_t; tt > t2; tt--){
						double dual = this.dualCosts[dataModel.C + tt - 1];
						current_rc -= dual;
					}

					BitSet colsIndicator_t2 = columnsIndicator.get(t2);

					// For the columns above the diagonal, currentLabel has a b strictly greater than their corresponding b2, thus
					// currentLabel can be dominated by them.
					for (int b2 = colsIndicator_t2.previousSetBit(b - (t-t2)); b2 >= 1; b2 = colsIndicator_t2.previousSetBit(b2 - 1)){
						Label otherLabel = columnsMap.get(pack(b2, t2));
						if (otherLabel.reducedCost < current_rc - dataModel.precision){ // currentLabel is dominated
							colsIndicator.clear(b);
							dominated = true; break;
						}
					}
					
					// For the columns below the diagonal, currentLabel has a b strictly lower than their corresponding b2, thus
					// currentLabel can dominate them. For the column in the diagonal, currentLabel and otherLabel have the same b.
					
					for (int b2 = colsIndicator_t2.previousSetBit(b - 1); b2 >= b - (t-t2) ; b2 = colsIndicator_t2.previousSetBit(b2 - 1)){
						Label otherLabel = columnsMap.get(pack(b2,t2));
						if (current_rc < otherLabel.reducedCost - dataModel.precision){ // otherLabel is dominated
							colsIndicator_t2.clear(b2);
						}
					}
					
					if (t2 >= this.nonDominatedT.get(currentLabel.index) && !this.negative_charging_duals.get(t)){
						colsIndicator_t2.clear(b);
					}

					for (int b2 = colsIndicator_t2.nextSetBit(b + 1); b2 >= 1; b2 = colsIndicator_t2.nextSetBit(b2 + 1)){
						Label otherLabel = columnsMap.get(pack(b2,t2));
						if (current_rc < otherLabel.reducedCost - dataModel.precision){ // otherLabel is dominated
							colsIndicator_t2.clear(b2);
						}
					}

					// If the currentLabel is dominated, it won't be extended to previous time periods.
					if (dominated) break;
					previous_t = t2;
				}

			}
		
		}

		// Retrieve the resulting non-dominated columns
		this.last_charging_periods = new HashMap<>();
		filtered_labels.clear();
		for (Map.Entry<Integer, BitSet> entry: columnsIndicator.entrySet()){

			int t = entry.getKey();
			BitSet colsIndicator = entry.getValue();

			for (int b = colsIndicator.nextSetBit(1); b >= 0; b = colsIndicator.nextSetBit(b + 1)){

				Label label = columnsMap.get(pack(b,t));
				if (this.last_charging_periods.containsKey(label.index)) this.last_charging_periods.get(label.index).set(t);
				else {
					BitSet newBit = new BitSet(); newBit.set(t);
					this.last_charging_periods.put(label.index, newBit);
					filtered_labels.add(label);
				}
			}
		}

		return filtered_labels;
	}

	private void filter_labels_same_chargingTime(Map.Entry<Integer, List<Label>> entry, BitSet t_set, Map<Long, Label> colsMap, Map<Integer, BitSet> colsInd){

		int b = entry.getKey();
		List<Label> labels_group = entry.getValue();
		List<Label> labels_to_remove = new ArrayList<>();

		//////////////////////////////////////////////////////////////////////////////
		///  1. Dominance between labels of same chargingTime and same departureTime
		//////////////////////////////////////////////////////////////////////////////

		// For every unique departure time d keep only the label with minimum reduced cost
		Map<Integer, Label> bestPerDeparture = new HashMap<>();
		for (Label l : labels_group) { // detects dominated labels of same departure time
			int d = l.vertex; double rc = l.reducedCost;

			if (!bestPerDeparture.containsKey(d)) bestPerDeparture.put(d, l);
			else {
				Label best = bestPerDeparture.get(d);
				if (rc < best.reducedCost - dataModel.precision) {
					labels_to_remove.add(best); // new best, old fully dominated
					bestPerDeparture.put(d, l);
				} else { labels_to_remove.add(l); } // the new one is fully dominates
			}
		}

		//////////////////////////////////////////////////////////////////////////////
		///  2. Dominance between labels of same chargingTime but diff departureTime
		//////////////////////////////////////////////////////////////////////////////

		

		// 2.a. Sort in descending order of departure time
		List<Label> sorted = new ArrayList<>(bestPerDeparture.values());
		sorted.sort(Comparator.comparingInt((Label l) -> l.vertex).reversed());

		// 2.b. Sweep to detect dominance in time periods
		Label bestLabel = sorted.get(0);
		double bestRC = bestLabel.reducedCost;
		
		int ix = 1;
		while (ix < sorted.size()) {
			Label l = sorted.get(ix); double rc = l.reducedCost;
			if (rc < bestRC - dataModel.precision) { // the current best label is partially dominated on t = 1 ... d
				int d = l.vertex;
				this.nonDominatedT.put(bestLabel.index, d);
				bestLabel = l; bestRC = rc; // update sweep front
			}  else  { labels_to_remove.add(l); }

			ix ++;
		}
		this.nonDominatedT.put(bestLabel.index, b);

		// Remove all the fully dominated labels
		labels_group.removeAll(labels_to_remove);

		//////////////////////////////////////////////////////////////////////////////
		///  3. Dominance between columns of the same label
		//////////////////////////////////////////////////////////////////////////////

		

		// For each non-fully dominated label, find all the last charging time periods for which they have a non-dominated column
		for (Label label: labels_group){

			int t = label.vertex-1; // Starting at departureTime - 1
			while (t >= this.nonDominatedT.get(label.index)){

				colsMap.put(pack(b, t), label);
				BitSet colsIndicator = colsInd.get(t);
				if (colsIndicator == null){ t_set.set(t); BitSet cI = new BitSet(); cI.set(b); colsInd.put(t,cI); }
				else colsIndicator.set(b);

				/* int next_t = t-b; // find the next non-dominated time period
				for (int tt = t; tt >= t-b+1; tt--){
					if (this.negative_charging_duals.get(tt)){
						next_t = tt - 1;
						break;
					}
				}

				t = next_t; */
				t--;
			}
		}

	}

	/**
	 * Verifies if L1 is (strongly) dominated by L2
	 * @param L1, L2 labels
	 */
	public boolean isDominated(Label L1, Label L2) {

		if(L1.vertex>dataModel.C) { //charging time vertices
			if (L2.chargingTime>L1.chargingTime) return false;
			if (L2.reducedCost-L1.reducedCost>dataModel.precision) return false;
			return true;

		}else { //customer vertices

			if (L1.vertex>0 && L2.remainingLoad<L1.remainingLoad) return false; 	//load
			if (L2.reducedCost-L1.reducedCost>dataModel.precision) return false; 	//reduced cost
			if (L2.remainingTime<L1.remainingTime) return false; 					//time
			
			for (int gam = 0; gam <= dataModel.gamma; gam ++){
				if (L2.remainingEnergy[gam]<L1.remainingEnergy[gam]) return false;				 //energy
			}

			// Ng-paths and unreachable resources
			Vertex currentVertex = dataModel.vertices[L1.vertex];
			if (currentVertex.id > 0) {
				for(int i: dataModel.vertices[currentVertex.id].neighbors) {
					
					//boolean check_binaries = (L2.ng_path[i-1] || L2.unreachable[i-1]) && !(L1.ng_path[i-1] || L1.unreachable[i-1]);
					boolean other_way = L2.ng_path[i-1] && (!L1.unreachable[i-1] && !L1.ng_path[i-1]); // Dani's way
					if (other_way) {
						return false;
					}
				}
			}

			return true;
		}
	}

	static long pack(int b, int t) {
    	return (((long) b) << 32) | (t & 0xffffffffL);
	}

}
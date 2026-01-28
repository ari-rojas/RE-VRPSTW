package columnGeneration;

import java.util.ArrayList;
import java.util.BitSet;
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

/**
 * This class defines the pricing problem. 
 * We simply extend the pricing problem included in the framework (there is no need for modification, only one pricing problem)
 */
public final class PricingProblem extends AbstractPricingProblem<EVRPTW> {

	public ArrayList<SubsetRowInequality> subsetRowCuts; 				//subset row cuts considered
	public Set<ChargingTimeInequality> branchesOnChargingTimes;			//branching on charging times
	public double bestReducedCost = -Double.MAX_VALUE; 					//best reduced cost found by the exact labeling
	public double reducedCostThreshold = 0; 							//minimum reduced cost when arriving at the depot source

	public int maxCols = 800;
	public boolean isExact;

	//Charging pricing information
	private BitSet negative_charging_duals;
	public Map<Integer, Map<Integer, Double>> charging_reducedCosts;

	public Map<Integer, Map<Integer,Double>> charging_bounds;
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

		// Include the charging branching information for the bounds
		int i=0;
		for(ChargingTimeInequality branching: this.branchesOnChargingTimes) {
			this.dualCosts[dataModel.C + branching.timestep - 1] += this.dualCosts[dataModel.C+dataModel.last_charging_period+this.subsetRowCuts.size()+i]; i++;
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
		this.isExact = true;

		//////////////////////////////////////////////////////////
		/// 1. Bounding Procedure and Preprocessing
		//////////////////////////////////////////////////////////

		// To avoid constantly recomputing the departure times of the labels, we save them in the vertex field
		int ix = 0;
		for (Label label: labels) { label.index = ix; label.vertex = (int)(label.remainingTime/10); ix ++; }

		// Only labels that will generate at least one column with negative reduced cost are accounted for
		ArrayList<Label> filtered_labels = new ArrayList<>();
		for (Label l: labels){ if (l.reducedCost < -dataModel.precision) filtered_labels.add(l); }

		//////////////////////////////////////////////////////////
		/// 2. Labels dominance
		//////////////////////////////////////////////////////////

		// 1. Group labels by charging time b
		Map<Integer, List<Label>> labelsByB = new HashMap<>();
		for (Label label : filtered_labels) labelsByB.computeIfAbsent(label.chargingTime, k -> new ArrayList<>()).add(label);

		// 2. Dominance between labels of same chargingTime b
		// For each chargingTime b, get the last-charging-time-periods that have a non-dominated column, and their corresponding label
		BitSet b_set = new BitSet();
		Map<Integer, BitSet> columnsIndicator = new HashMap<>();
		for (Map.Entry<Integer, List<Label>> entry : labelsByB.entrySet()) {
			int b = entry.getKey(); b_set.set(b);
			columnsIndicator.put(b, new BitSet());
			filter_labels_same_chargingTime(entry);
		}

		Map<Long, Label> columnsMap = new HashMap<>();
		this.last_charging_periods = new HashMap<>();
		filtered_labels.sort( Comparator.comparingDouble( l -> l.reducedCost) );
		
		// 3. Dominance between labels of different chargingTime b and different last_t t
		for (Label currentLabel: filtered_labels){
			
			int b = currentLabel.chargingTime;
			int index = currentLabel.index;

			int t = (int)(currentLabel.remainingTime/10) - 1; // starting at departureTime - 1
			while (t >= this.nonDominatedT.get(index)){

				boolean dominated = false;

				// For chargingTimes b2 > b, the label is extended forward
				for (int b2 = b_set.nextSetBit(b + 1); b2 >= 0; b2 = b_set.nextSetBit(b2 + 1)){

					double current_rc = currentLabel.reducedCost;
					BitSet colsIndicator = columnsIndicator.get(b2); int previous_t = t;
					
					for (int t2 = colsIndicator.previousSetBit(t + b2 - 1); t2 >= t + (b2-b); t2 = colsIndicator.previousSetBit(t2 - 1)){
						
						// Update the reduced cost of currentLabel (as the labels would be "extended")
						for (int tt = previous_t+1; tt <= t2; tt++){ current_rc += this.dualCosts[dataModel.C + tt - 1]; }

						// Check whether currentLabel is dominated or not
						Label otherLabel = columnsMap.get(pack(b2,t2));
						if (otherLabel.reducedCost < current_rc - dataModel.precision){  dominated = true; break; }

						previous_t = t2;
					}

					if (dominated) break;
				}

				if (dominated){ t = t - 1; continue; }
				
				// For chargingTimes b2 < b, the procedure is separated between time periods t2 < t and t2 >= t
				for (int b2 = b_set.previousSetBit(b - 1); b2 >= 1; b2 = b_set.previousSetBit(b2 - 1)){
					
					double current_rc = currentLabel.reducedCost;
					BitSet colsIndicator = columnsIndicator.get(b2); int previous_t = t;
					
					// For t2 >= t, the label is extended forward
					for (int t2 = colsIndicator.previousSetBit(t + b2 - 1); t2 >= t; t2 = colsIndicator.previousSetBit(t2 - 1)){
						// Update the reduced cost of currentLabel (as the labels would be "extended")
						for (int tt = previous_t+1; tt <= t2; tt++){ current_rc += this.dualCosts[dataModel.C + tt - 1]; }

						// Check whether currentLabel is dominated or not
						Label otherLabel = columnsMap.get(pack(b2,t2));
						if (otherLabel.reducedCost < current_rc - dataModel.precision){  dominated = true; break; }

						previous_t = t2;
					}

					if (dominated) break;

					current_rc = currentLabel.reducedCost;

					// For t2 < t, the label is extended backwards
					for (int t2 = colsIndicator.previousSetBit(t - 1); t2 >= t + (b2-b); t2 = colsIndicator.previousSetBit(t2 - 1)){
						// Update the reduced cost of currentLabel (as the labels would be "extended")
						for (int tt = previous_t; tt > t2; tt--){ current_rc -= this.dualCosts[dataModel.C + tt - 1]; }

						// Check whether currentLabel is dominated or not
						Label otherLabel = columnsMap.get(pack(b2,t2));
						if (otherLabel.reducedCost < current_rc - dataModel.precision){  dominated = true; break; }

						previous_t = t2;
					}

					if (dominated) break;
				}

				if (dominated){ t = t - 1; continue; }

				// If the column with last_t t of currentLabel is NOT dominated, map the column
				
				if (currentLabel.reducedCost + this.charging_reducedCosts.get(b).get(t) < - dataModel.precision){
					columnsMap.put(pack(b,t), currentLabel);
					columnsIndicator.get(b).set(t);
					if (this.last_charging_periods.containsKey(index)) this.last_charging_periods.get(index).set(t);
					else {
						BitSet newBit = new BitSet(); newBit.set(t);
						this.last_charging_periods.put(index, newBit);
					}
				}
				
				// and find the next non-dominated time period
				int next_t = t-b; double current_rc = currentLabel.reducedCost;
				for (int tt = t; tt >= t-b+2; tt--){
					current_rc -= this.dualCosts[dataModel.C + tt - 1];
					if (currentLabel.reducedCost < current_rc - dataModel.precision){
						next_t = tt - 1;
						break;
					}
				}

				t = next_t;
			}

			//if (columnsMap.size() > this.maxCols) { this.isExact = false; break;}
		}
		
		ArrayList<Label> to_remove = new ArrayList<>();
		for (Label l: filtered_labels) if (!last_charging_periods.containsKey(l.index)) to_remove.add(l);
		filtered_labels.removeAll(to_remove);

		return filtered_labels;
	}

	private void filter_labels_same_chargingTime(Map.Entry<Integer, List<Label>> entry){

		int b = entry.getKey();
		List<Label> labels_group = entry.getValue();

		//////////////////////////////////////////////////////////////////////////////
		///  1. Dominance between labels of same chargingTime but diff departureTime
		//////////////////////////////////////////////////////////////////////////////

		// 1.a. Sort in descending order of departure time
		PriorityQueue<Label> sorted = new PriorityQueue<>((l1, l2) -> Integer.compare(l2.vertex, l1.vertex)); // sort the labels by descending departure time
		sorted.addAll(labels_group);

		// 1.b. Sweep to detect dominance in time periods
		Label bestLabel = sorted.poll();
		while (!sorted.isEmpty()) {

			// the current bestLabel is partially dominated on t = 1 ... d by l
			Label l = sorted.poll();
			int d = l.vertex;
			this.nonDominatedT.put(bestLabel.index, d);

			bestLabel = l; // update sweep front
		}
		this.nonDominatedT.put(bestLabel.index, b);

	}

	/**
	 * Verifies if L1 is (strongly) dominated by L2
	 * @param L1, L2 labels
	 */
	public boolean isDominated(Label L1, Label L2) {

		//customer vertices
		if (L2.remainingLoad<L1.remainingLoad) return false; 	//load
		if (L2.reducedCost-L1.reducedCost>dataModel.precision) return false; 	//reduced cost
		if (L2.remainingTime<L1.remainingTime) return false; 					//time
		
		for (int gam = 0; gam <= dataModel.gamma; gam ++){
			if (L2.remainingEnergy[gam]<L1.remainingEnergy[gam]) return false;				 //energy
		}

		return true;
		
	}

	static long pack(int b, int t) {
    	return (((long) b) << 32) | (t & 0xffffffffL);
	}

}
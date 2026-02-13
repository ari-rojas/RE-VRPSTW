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

	public int maxCols = 700;
	public boolean isExact;

	//Charging pricing information
	private BitSet negative_charging_duals;
	public Map<Integer, Map<Integer, Double>> charging_reducedCosts;

	private HashMap<Integer, Integer> nonDominatedT;
	public Map<Integer, BitSet> last_charging_periods;

	public double[] last_charging_branch_duals;
	public double[] initial_charging_branch_duals;

	public PricingProblem(EVRPTW modelData, String name) {
		super(modelData, name);
	}

	public void compute_charging_bounds(){

		this.negative_charging_duals = new BitSet();

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
			this.negative_charging_duals.set(t, dual < -dataModel.precision);
			S[t] = S[t - 1] + dual;
		}

		// Include the charging branching dual information
		this.last_charging_branch_duals = new double[maxT];
		this.initial_charging_branch_duals = new double[maxT];
		int i=0;
		for(ChargingTimeInequality branching: this.branchesOnChargingTimes) {
			double dual = this.dualCosts[dataModel.C+dataModel.last_charging_period+this.subsetRowCuts.size()+i];
			if (branching.startCharging) this.initial_charging_branch_duals[branching.timestep] = dual;
			else this.last_charging_branch_duals[branching.timestep] = dual;
			i++;
		}

		///////////////////////////////////////////////////////////////////////////////////
		/// Compute the bounds for every combination of chargingTime b and departureTime d
		///////////////////////////////////////////////////////////////////////////////////

		this.charging_reducedCosts = new HashMap<>();
		for (int b = 1; b <= dataModel.f_inverse[dataModel.E]; b++){

			Map<Integer, Double> reducedCostsMap = new HashMap<>();
			int first_departure = Math.max(b+1, minT);
			for (int last_t = b; last_t < first_departure-1; last_t ++){
				double rc = - (S[last_t] - S[last_t-b]) - this.last_charging_branch_duals[last_t] - this.initial_charging_branch_duals[last_t-b+1];
				reducedCostsMap.put(last_t, rc);
			}
			
			for (int d = first_departure; d <= maxT; d++){
				
				double rc = - (S[d-1] - S[d-b-1]) - this.last_charging_branch_duals[d-1] - this.initial_charging_branch_duals[d-b];
				reducedCostsMap.put(d-1, rc);
				
			}
			
			this.charging_reducedCosts.put(b, reducedCostsMap);
			
		}

	}

	public ArrayList<Label> charging_pricing_filtering(ArrayList<Label> labels){

		this.nonDominatedT = new HashMap<>();
		this.last_charging_periods = new HashMap<>();
		this.isExact = true;

		//////////////////////////////////////////////////////////
		/// 1. Bounding Procedure and Preprocessing
		//////////////////////////////////////////////////////////

		// To avoid constantly recomputing the departure times of the labels, we save them in the vertex field
		ArrayList<Label> filtered_labels = new ArrayList<>(); int ix = 0;
		for (Label label: labels) {
			label.index = ix; label.vertex = (int)(label.remainingTime/10);
			filtered_labels.add(label); ix ++;
		}

		//////////////////////////////////////////////////////////
		/// 2. Labels dominance
		//////////////////////////////////////////////////////////

		// 1. Group labels by charging time b
		Map<Integer, List<Label>> labelsByB = new HashMap<>();
		for (Label label : filtered_labels) labelsByB.computeIfAbsent(label.chargingTime, k -> new ArrayList<>()).add(label);

		// 2. Dominance between labels of same chargingTime b
		// For each chargingTime b, get the last-charging-time-periods that have a non-dominated column, and their corresponding label
		PriorityQueue<RouteColumn> columnsQueue = new PriorityQueue<>(
			Comparator
			.comparingDouble((RouteColumn r) -> r.reducedCost)						// 1) lowest reducedCost first
			.thenComparingInt(r -> r.b)                         					// 2) lowest b first
			.thenComparing((r1, r2) -> Integer.compare(r2.last_t, r1.last_t)) 		// 3) highest last_t first
		);
		for (Map.Entry<Integer, List<Label>> entry : labelsByB.entrySet()) exhaustive_filter_labels_same_chargingTime(entry, columnsQueue);
		
		// The "processed" columns are mapped
		Map<Long, Label> columnsMap = new HashMap<>();
		int maxB = dataModel.f_inverse[dataModel.E]; int maxT = dataModel.last_charging_period;
		BitSet t_set = new BitSet();
		BitSet[] columnsIndicator = new BitSet[maxT+1];
		for (int t = 1; t <= maxT; t++) columnsIndicator[t] = new BitSet();

		while (!columnsQueue.isEmpty()){

			RouteColumn column = columnsQueue.poll();
			int b = column.b; int t = column.last_t;
			Label currentLabel = column.routeLabel;

			boolean dominated = false;

			// See if there are any columns of the same route that might dominate it
			for (int t2 = t_set.nextSetBit(t+1); t2 > 0 && t2 <= t+b-1; t2 = t_set.nextSetBit(t2 + 1)){
				if (columnsIndicator[t2].get(b)){
					Label otherLabel = columnsMap.get(pack(b,t2));
					if (otherLabel.index == currentLabel.index && !this.negative_charging_duals.get(t+1)) dominated = true;
					break; // Only needs to check dominance with respect to the first found column
				}
			}

			if (dominated) continue;

			// For last_t t2 >= t
			double current_rc = currentLabel.reducedCost; int previous_t = t;
			for (int t2 = t_set.nextSetBit(t); t2 > 0 && t2 <= t+maxB-1; t2 = t_set.nextSetBit(t2 + 1)){

				BitSet colsIndicator = columnsIndicator[t2];

				// Update the reduced cost of currentLabel (as the labels would be "extended")
				for (int tt = previous_t+1; tt <= t2; tt++){ current_rc += this.dualCosts[dataModel.C + tt - 1]; }

				// For b2 < b
				for (int b2 = colsIndicator.previousSetBit(b-1); b2 >= t2-t+1; b2 = colsIndicator.previousSetBit(b2 - 1)){
					Label otherLabel = columnsMap.get(pack(b2,t2));
					if (otherLabel.reducedCost < current_rc - dataModel.precision) {dominated = true; break;}
				}

				if (dominated) break;

				// For b2 > b
				int b_low = Math.max(b+1, t2-t+1);
				for (int b2 = colsIndicator.previousSetBit(b+t2-t); b2 >= b_low; b2 = colsIndicator.previousSetBit(b2 - 1)){
					Label otherLabel = columnsMap.get(pack(b2,t2));
					if (otherLabel.reducedCost < current_rc - dataModel.precision) {dominated = true; break;}
				}

				if (dominated) break;
				previous_t = t2;
			}

			if (dominated) continue;

			// For last_t t2 < t
			current_rc = currentLabel.reducedCost; previous_t = t;
			int t_low = Math.max(1, t-b+1);
			for (int t2 = t_set.previousSetBit(t-1);  t2 >= t_low; t2 = t_set.previousSetBit(t2 - 1)){

				BitSet colsIndicator = columnsIndicator[t2];

				// Update the reduced cost of currentLabel (as the labels would be "extended")
				for (int tt = previous_t; tt > t2; tt--){ current_rc -= this.dualCosts[dataModel.C + tt - 1]; }

				for (int b2 = colsIndicator.previousSetBit(b+t2-t); b2 >= 1; b2 = colsIndicator.previousSetBit(b2-1)){
					Label otherLabel = columnsMap.get(pack(b2,t2));
					if (otherLabel.reducedCost < current_rc - dataModel.precision) {dominated = true; break;}
				}

				if (dominated) break;
				previous_t = t2;

			}

			if (dominated) continue;

			// If the column is NOT dominated, map the column
			int index = currentLabel.index;
			columnsMap.put(pack(b,t), currentLabel);
			t_set.set(t); columnsIndicator[t].set(b);
			if (this.last_charging_periods.containsKey(index)) this.last_charging_periods.get(index).set(t);
			else {
				BitSet newBit = new BitSet(); newBit.set(t);
				this.last_charging_periods.put(index, newBit);
			}

			//if (columnsMap.size() > this.maxCols) break;
			
		}

		ArrayList<Label> to_remove = new ArrayList<>();
		for (Label l: filtered_labels) if (!last_charging_periods.containsKey(l.index)) to_remove.add(l);
		filtered_labels.removeAll(to_remove);

		return filtered_labels;
	}

	private void filter_labels_same_chargingTime(Map.Entry<Integer, List<Label>> entry, PriorityQueue<RouteColumn> colsQueue){

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
			
			Label l = sorted.poll();
			int d = l.vertex;

			// The current bestLabel is partially dominated from t = 1 to t = d(l)
			// therefore has a column for every last_charging_time_period t \in {d(l), ..., d-1}
			for (int t = d; t < bestLabel.vertex; t++){
				double col_rc = bestLabel.reducedCost + this.charging_reducedCosts.get(b).get(t);
				if (col_rc < - dataModel.precision) colsQueue.add(new RouteColumn(b,t, bestLabel.reducedCost, bestLabel));
			}

			bestLabel = l; // update sweep front
		}
		for (int t = b; t < bestLabel.vertex; t++){
			double col_rc = bestLabel.reducedCost + this.charging_reducedCosts.get(b).get(t);
			if (col_rc < - dataModel.precision) colsQueue.add(new RouteColumn(b,t, bestLabel.reducedCost, bestLabel));
		}

	}

	private void exhaustive_filter_labels_same_chargingTime(Map.Entry<Integer, List<Label>> entry, PriorityQueue<RouteColumn> colsQueue){

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
			
			Label l = sorted.poll();
			int d = l.vertex;

			// The current bestLabel is partially dominated from t = 1 to t = d(l)
			// therefore has a column for every last_charging_time_period t \in {d(l), ..., d-1}
			for (int t = d; t < bestLabel.vertex; t++){
				double col_rc = bestLabel.reducedCost + this.charging_reducedCosts.get(b).get(t);
				if (col_rc < - dataModel.precision) colsQueue.add(new RouteColumn(b,t, col_rc, bestLabel));
			}

			bestLabel = l; // update sweep front
		}
		for (int t = b; t < bestLabel.vertex; t++){
			double col_rc = bestLabel.reducedCost + this.charging_reducedCosts.get(b).get(t);
			if (col_rc < - dataModel.precision) colsQueue.add(new RouteColumn(b,t, col_rc, bestLabel));
		}

	}

	private class RouteColumn{

		int b;
		int last_t;
		double reducedCost;
		Label routeLabel;

		private RouteColumn(int b, int last_t, double rc, Label routeLabel){
			this.b = b;
			this.last_t = last_t;
			this.reducedCost = rc;
			this.routeLabel = routeLabel;
		}
	}

	static long pack(int b, int t) {
    	return (((long) b) << 32) | (t & 0xffffffffL);
	}

}
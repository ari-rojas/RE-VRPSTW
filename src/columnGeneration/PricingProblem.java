package columnGeneration;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
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

	//Charging pricing information
	private BitSet negative_charging_duals;
	private Map<Integer, Map<Integer, Double>> charging_bounds;
	public Map<Integer, Map<Integer, Double>> charging_reducedCosts;

	public Map<Integer, Integer> nonDominatedT;
	public Map<Integer, ArrayList<Integer>> last_charging_periods;

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
		this.last_charging_periods = new HashMap<>();

		ArrayList<Label> filtered_labels = new ArrayList<>();

		// To avoid constantly recomputing the departure times of the labels, we save them in the vertex field
		for (Label label: labels) label.vertex = (int)(label.remainingTime/10);

		// 1. Group labels by charging time b
		Map<Integer, List<Label>> labelsByB = new HashMap<>();
		for (Label label : labels) labelsByB.computeIfAbsent(label.chargingTime, k -> new ArrayList<>()).add(label);

		// 2. Dominance filtering by charging time b
		for (Map.Entry<Integer, List<Label>> entry : labelsByB.entrySet()) {
			filter_labels_same_chargingTime(entry, last_charging_periods); // returns the departures with non-fully-dominated labels
			filtered_labels.addAll(entry.getValue());
		}

		// 3. Remove the labels that won't have any negative reduced cost column
		List<Label> labels_to_remove = new ArrayList<>();
		for (Label l: filtered_labels) {
			if (l.reducedCost + this.charging_bounds.get(l.chargingTime).get(l.vertex) > -dataModel.precision) labels_to_remove.add(l);
		}
		labels.removeAll(labels_to_remove);

		return filtered_labels;
	}

	private void filter_labels_same_chargingTime(Map.Entry<Integer, List<Label>> entry, Map<Integer, ArrayList<Integer>> nonDom_last_charg_periods){

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
			ArrayList<Integer> nonDom_last_ts = new ArrayList<>();

			int t = label.vertex-1; // Starting at departureTime - 1
			while (t >= this.nonDominatedT.get(label.index)){

				nonDom_last_ts.add(t);

				int next_t = t-b; // find the next non-dominated time period
				for (int tt = t; tt >= t-b+1; tt--){
					if (this.negative_charging_duals.get(tt)){
						next_t = tt - 1;
						break;
					}
				}

				t = next_t;
			}

			nonDom_last_charg_periods.put(label.index, nonDom_last_ts);
		}
	}
}
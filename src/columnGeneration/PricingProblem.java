package columnGeneration;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
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
	public double[] S;
	public boolean[] negative_charging_duals;
	public Map<Integer, Map<Integer, Double>> charging_reducedCosts;
	public int maxT;

	public Map<Integer, Boolean> fullyDominated;
	public Map<Integer, Integer> nonDominatedT;

	private Map<Integer, TreeSet<Integer>> nonDominatedD;
	private Map<Integer, List<Double>> nonDominatedRC;
	private Label fakeLabel;

	public PricingProblem(EVRPTW modelData, String name) {
		super(modelData, name);

		this.fakeLabel = new Label(dataModel.V, 0, 0, 0, 0, 0, new int[0], 0, new boolean[0], new boolean[0], new boolean[0], new HashSet<>());
		this.fakeLabel.index = -1;
	}

	public void compute_charging_bounds(ArrayList<Label> labels){

		this.fullyDominated = new HashMap<>();
		this.nonDominatedT = new HashMap<>();
		this.nonDominatedRC = new HashMap<>();
		this.nonDominatedD = new HashMap<>();

		this.charging_reducedCosts = new HashMap<>();

		// 1. Group labels by charging time b
		Map<Integer, List<Label>> labelsByB = new HashMap<>();
		for (Label label : labels) labelsByB.computeIfAbsent(label.chargingTime, k -> new ArrayList<>()).add(label);

		// 2. Dominance filtering per charging time b
		for (Map.Entry<Integer, List<Label>> entry : labelsByB.entrySet()) {
			TreeSet<Integer> departures = filter_labels_same_chargingTime(entry); // returns the departures with non-fully-dominated labels
			if (!departures.isEmpty()) this.nonDominatedD.put(entry.getKey(), departures);
		}

		// Precompute fixed sums of the charging dual variables
		this.maxT = 0; for (TreeSet<Integer> set : this.nonDominatedD.values()) this.maxT = Math.max(this.maxT, set.last()-1);
		this.S = new double[dataModel.last_charging_period + 1]; this.S[0] = 0.0;
		this.negative_charging_duals = new boolean[dataModel.last_charging_period + 1];
		for (int t = 1; t <= this.maxT; t++) {
			double dual = this.dualCosts[dataModel.C + t - 1];
			this.negative_charging_duals[t] = dual < -dataModel.precision;
			this.S[t] = this.S[t - 1] + dual;
		}

		if (!this.nonDominatedD.isEmpty()) process_charging_times_map(this.nonDominatedD);

	}

	public void update_charging_bounds(ArrayList<Label> labels){

		this.fullyDominated = new HashMap<>();
		this.nonDominatedT = new HashMap<>();
		int prevT = this.maxT;

		// 1. Group labels by charging time b
		Map<Integer, List<Label>> labelsByB = new HashMap<>();
		for (Label label : labels) labelsByB.computeIfAbsent(label.chargingTime, k -> new ArrayList<>()).add(label);

		// 2. Dominance filtering per charging time b
		for (Map.Entry<Integer, List<Label>> entry : labelsByB.entrySet()){
			TreeSet<Integer> departures = filter_labels_same_chargingTime(entry);
			if (!departures.isEmpty()) this.nonDominatedD.put(entry.getKey(), departures);
		}

		// Update precomputed fixed sums of the charging dual variables
		for (TreeSet<Integer> set : this.nonDominatedD.values()) this.maxT = Math.max(this.maxT, set.last()-1);
		for (int t = prevT+1; t <= maxT; t++) {
			double dual = this.dualCosts[dataModel.C + t - 1];
			this.negative_charging_duals[t] = dual < -dataModel.precision;
			this.S[t] = this.S[t - 1] + dual;
		}

		if (!this.nonDominatedD.isEmpty()) process_charging_times_map(this.nonDominatedD);
		
	}

	public void process_charging_times_map(Map<Integer, TreeSet<Integer>> charging_times){

		for (Map.Entry<Integer, TreeSet<Integer>> e : charging_times.entrySet()){

			int b = e.getKey();
        	TreeSet<Integer> departures = e.getValue();

			int initial_t = 1;
			double rc = - (this.S[b]-this.S[0]);

			Map<Integer, Double> reducedCostsMap = new LinkedHashMap<>(); reducedCostsMap.put(b, rc);
			
			for (int d: departures){
				if (d <= b) continue; // skip if departure does not allow for sufficient charging

				for (int t=initial_t; t<=d-b-1; t++){
					rc = - (this.S[t+b] - this.S[t]);
					reducedCostsMap.put(t+b, rc);
				}
				initial_t = d-b;
			}
				
			this.charging_reducedCosts.put(b, reducedCostsMap);
			
		}

	}

	private TreeSet<Integer> filter_labels_same_chargingTime(Map.Entry<Integer, List<Label>> entry){

		int b = entry.getKey();
		List<Label> group = entry.getValue();

		// (a) group by departure time d and keep only the labels with minimum reduced cost
		Map<Integer, Label> bestPerDeparture = new HashMap<>();
		if (this.nonDominatedD.containsKey(b)) {
			List<Integer> departures = new ArrayList<>(this.nonDominatedD.get(b)); Collections.reverse(departures); // get both lists in descending order
			List<Double> nonDomRC = this.nonDominatedRC.get(b);

			for (int i = 0; i<departures.size(); i++) {
				int d = departures.get(i);
				double rc = nonDomRC.get(i);

				Label newFake = fakeLabel.clone(); newFake.reducedCost = rc; newFake.remainingTime = d*10; newFake.index = -1;
				bestPerDeparture.put(d, newFake);
			}
		}

		for (Label l : group) { // detects dominated labels of same departure time
			int d = (int) (l.remainingTime / 10);
			double rc = l.reducedCost;

			if (!bestPerDeparture.containsKey(d)) bestPerDeparture.put(d, l);
			else {
				Label best = bestPerDeparture.get(d);
				if (rc < best.reducedCost - dataModel.precision) {
					// new best, old fully dominated
					fullyDominated.put(best.index, true);
					bestPerDeparture.put(d, l);
				} else { fullyDominated.put(l.index, true); } // the new one is fully dominatws
			}
		}

		// (b) sort in descending order of departure time
		List<Label> sorted = new ArrayList<>(bestPerDeparture.values());
		sorted.sort(Comparator.comparingInt((Label l) -> (int) (l.remainingTime / 10)).reversed());

		// (c) sweep to detect dominance in time periods
		TreeSet<Integer> departuresForB = new TreeSet<>();

		Label bestLabel = sorted.get(0);
		double bestRC = bestLabel.reducedCost; int bestD = (int) (bestLabel.remainingTime / 10);
		departuresForB.add(bestD); fullyDominated.put(bestLabel.index, false);
		List<Double> nonDomRC = new ArrayList<>(Arrays.asList(bestRC));
		
		int ix = 1; // Loops over the Departure Times in descending order
		while (ix < sorted.size()) {

			Label l = sorted.get(ix);
			double rc = l.reducedCost;

			if (rc < bestRC - dataModel.precision) { // the current best label is partially dominated on t = 1 ... d
				
				int d = (int) (l.remainingTime/10);
				this.nonDominatedT.put(bestLabel.index, d);
				
				// update sweep front
				bestLabel = l; bestRC = rc; bestD = d;
				departuresForB.add(d); fullyDominated.put(l.index, false);
				nonDomRC.add(rc);
			}  else  { fullyDominated.put(l.index, true); }

			ix ++;
		}
		this.nonDominatedT.put(bestLabel.index, b);
		
		this.nonDominatedRC.put(b,nonDomRC);

		return departuresForB;
	}
}
package columnGeneration;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
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
	public Map<Integer, Map<Integer, Double>> charging_bounds;
	public Map<Integer, Map<Integer, Double>> charging_reducedCosts;
	public int maxT;

	public PricingProblem(EVRPTW modelData, String name) {
		super(modelData, name);
	}

	public void compute_charging_bounds(ArrayList<Label> labels){

		Map<Integer, TreeSet<Integer>> charging_times = new HashMap<>();

		for (Label label: labels){
			int departureTime = (int) (label.remainingTime/10);
			charging_times.computeIfAbsent(label.chargingTime, k -> new TreeSet<>()).add(departureTime);
		}

		// Precompute fixed sums of the charging dual variables
		this.maxT = 0; for (TreeSet<Integer> set : charging_times.values()) this.maxT = Math.max(this.maxT, set.last());
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

		this.charging_reducedCosts = new HashMap<>();
		this.charging_bounds = new HashMap<>();
		if (!charging_times.isEmpty()) process_charging_times_map(charging_times);

	}

	public void update_charging_bounds(ArrayList<Label> labels){

		int prevT = this.maxT;
		Map<Integer, TreeSet<Integer>> charging_times = new HashMap<>(); // new charging times that have not been tracked before
		Map<Integer, TreeSet<Integer>> existing_charging_times = new HashMap<>();

		for (Label label: labels){
			int d = (int) (label.remainingTime/10);
			int b = label.chargingTime;
			if (this.charging_bounds.containsKey(b)) existing_charging_times.computeIfAbsent(b, k -> new TreeSet<>()).add(d);
			else charging_times.computeIfAbsent(b, k -> new TreeSet<>()).add(d);
		}

		// Update precomputed fixed sums of the charging dual variables
		for (TreeSet<Integer> set : existing_charging_times.values()) this.maxT = Math.max(this.maxT, set.last());
		for (TreeSet<Integer> set : charging_times.values()) this.maxT = Math.max(this.maxT, set.last());
		for (int t = prevT+1; t <= maxT; t++) {
			double dual = this.dualCosts[dataModel.C + t - 1];
			this.negative_charging_duals[t] = dual < -dataModel.precision;
			this.S[t] = this.S[t - 1] + dual;
		}

		if (!charging_times.isEmpty()) process_charging_times_map(charging_times);
		if (!existing_charging_times.isEmpty()) update_charging_times_map(existing_charging_times);
		
	}

	public void process_charging_times_map(Map<Integer, TreeSet<Integer>> charging_times){

		for (Map.Entry<Integer, TreeSet<Integer>> e : charging_times.entrySet()){

			int b = e.getKey();
        	TreeSet<Integer> departures = e.getValue();

			int initial_t = 1;
			double rc = - (this.S[b]-this.S[0]); double min_rc = rc;

			Map<Integer, Double> reducedCostsMap = new HashMap<>(); reducedCostsMap.put(b, rc);
			Map<Integer, Double> boundsMap = new HashMap<>();
			
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

	public void update_charging_times_map(Map<Integer, TreeSet<Integer>> charging_times){

		for (Map.Entry<Integer, TreeSet<Integer>> e : charging_times.entrySet()){
			
			int b = e.getKey();
			Map<Integer, Double> reducedCostsMap = this.charging_reducedCosts.get(b);
			Map<Integer, Double> boundsMap = this.charging_bounds.get(b);
			
        	TreeSet<Integer> departures = e.getValue(); int lastD = Collections.max(boundsMap.keySet());
			departures.addAll(boundsMap.keySet());
			
			int initial_t = 1;
			double rc = reducedCostsMap.get(b); double min_rc = rc;
			for (int d: departures){
				if (d <= b) continue; // skip if departure does not allow for sufficient charging

				if (!boundsMap.containsKey(d)){ // if this is a new departure time
					if (d <= lastD){ // if the departure is less than the previous maximum departure of the same charging time, then its reducedCost is already in the object
						for (int t=initial_t; t<=d-b-1; t++){
							rc = reducedCostsMap.get(t+b);
							if (rc < min_rc - dataModel.precision) min_rc = rc;
						}

					} else { // if not, the value needs to be added
						for (int t=initial_t; t<=d-b-1; t++){
							rc = - (this.S[t+b] - this.S[t]);
							reducedCostsMap.put(t+b, rc);
							if (rc < min_rc - dataModel.precision) min_rc = rc;
						}

					}
					boundsMap.put(d, min_rc);
				} else { min_rc = boundsMap.get(d); }

				initial_t = d-b;
			}
			
		}
	}
}
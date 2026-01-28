package columnGeneration;

import java.util.ArrayList;
import java.util.Set;
import java.util.Map;
import java.util.HashMap;
import java.util.LinkedHashMap;
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

	public Map<Integer, Map<Integer,Double>> charging_bounds;

	public PricingProblem(EVRPTW modelData, String name) {
		super(modelData, name);
	}

	public void compute_charging_bounds(){

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
		for (int t = 1; t < maxT; t++) {
			double dual = this.dualCosts[dataModel.C + t - 1];
			S[t] = S[t - 1] + dual;
		}

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
}
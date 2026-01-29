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
	public double[] last_charging_branch_duals;

	public PricingProblem(EVRPTW modelData, String name) {
		super(modelData, name);
	}

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

		// Include the charging branching dual information
		this.last_charging_branch_duals = new double[maxT]; int i=0;
		for(ChargingTimeInequality branching: this.branchesOnChargingTimes) {
			this.last_charging_branch_duals[branching.timestep] = this.dualCosts[dataModel.C+dataModel.last_charging_period+this.subsetRowCuts.size()+i]; i++;
		}

		///////////////////////////////////////////////////////////////////////////////////
		/// Compute the bounds for every combination of chargingTime b and departureTime d
		///////////////////////////////////////////////////////////////////////////////////

		this.charging_bounds = new HashMap<>();
		for (int b = 1; b <= dataModel.f_inverse[dataModel.E]; b++){

			double min_rc = Double.MAX_VALUE;
			Map<Integer, Double> boundsMap = new LinkedHashMap<>();
			for (int d = Math.max(b+1, minT); d <= maxT; d++){
				
				double rc = - (S[d-1] - S[d-b-1]) - this.last_charging_branch_duals[d-1];
				if (rc < min_rc - dataModel.precision) min_rc = rc;
				
				boundsMap.put(d, min_rc);
			}
			
			this.charging_bounds.put(b, boundsMap);
			
		}

	}
}
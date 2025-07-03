package columnGeneration;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import org.jorlib.frameworks.columnGeneration.master.MasterData;
import org.jorlib.frameworks.columnGeneration.util.OrderedBiMap;

import branchAndPrice.NumberVehiclesInequalities;
import ilog.concert.IloNumVar;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import model.EVRPTW;

/**
 * Class that stores the information for the Master Problem (MP).
 */
public final class VRPMasterData extends MasterData<EVRPTW, Route, PricingProblem, IloNumVar>{

	public final IloCplex cplex;												//CPLEX instance
	public final PricingProblem pricingProblem;									//list of pricing problems
	public HashMap<SubsetRowInequality, IloRange> subsetRowInequalities;		//mapping of the Subset row inequalities to constraints in the CPLEX model
	public Map<NumberVehiclesInequalities, IloRange> branchingNumberOfVehicles;	//mapping of branching decisions on the number of vehicles
	public Map<Route, Double> routeValueMap;									//routes used (only non-zero routes are considered) 

	public VRPMasterData(IloCplex cplex, PricingProblem pricingProblem, Map<PricingProblem, OrderedBiMap<Route, IloNumVar>> varMap) {
		super(varMap);
		this.cplex = cplex;
		this.pricingProblem = pricingProblem;
		this.subsetRowInequalities = new LinkedHashMap<>();
		this.routeValueMap = new HashMap<>();
		this.branchingNumberOfVehicles = new HashMap<NumberVehiclesInequalities, IloRange>();
	}

	public Map<NumberVehiclesInequalities, IloRange> getBranchingNumberOfVehicles(){
		return new HashMap<NumberVehiclesInequalities, IloRange>(this.branchingNumberOfVehicles);
	}

	public void setBranchingNumberOfVehicles(Map<NumberVehiclesInequalities, IloRange> map){
		this.branchingNumberOfVehicles = map;
	}

	@Override
	public void addColumn(Route column, IloNumVar variable) {
		if (!( ((OrderedBiMap)this.varMap.get(column.associatedPricingProblem)).containsKey(column) )) {
			((OrderedBiMap)this.varMap.get(column.associatedPricingProblem)).put(column, variable);
		}
   	}

}
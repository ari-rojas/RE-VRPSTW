package columnGeneration;

import ilog.concert.IloColumn;
import ilog.concert.IloConstraint;
import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import model.EVRPTW;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.io.TimeLimitExceededException;
import org.jorlib.frameworks.columnGeneration.master.AbstractMaster;
import org.jorlib.frameworks.columnGeneration.master.OptimizationSense;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractInequality;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.CutHandler;
import org.jorlib.frameworks.columnGeneration.util.OrderedBiMap;
import branchAndPrice.BranchEndChargingTimeDown;
import branchAndPrice.BranchEndChargingTimeUp;
import branchAndPrice.BranchInitialChargingTimeDown;
import branchAndPrice.BranchInitialChargingTimeUp;
import branchAndPrice.BranchVehiclesDown;
import branchAndPrice.BranchVehiclesUp;
import branchAndPrice.ChargingTimeInequality;
import branchAndPrice.FixArc;
import branchAndPrice.NumberVehiclesInequalities;
import branchAndPrice.RemoveArc;
import java.io.PrintWriter;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Implementation of the Master Problem (MP). 
 * It is handled by CPLEX.
 */
public final class Master extends AbstractMaster<EVRPTW, Route, PricingProblem, VRPMasterData> {

	private IloObjective obj; 						//objective
	private IloRange[] visitCustomerConstraints; 	//partitioning constraints
	private IloRange[] chargersCapacityConstraints; //capacity constraints
	private IloRange roundedCapacityInequality; 	//(weak) rounded capacity inequality
	private int minimumNumberOfVehicles; 			//for the weak rounded capacity inequality
	private List<Route> solutionKeeper; 			//stores the solution found

	// Lexicographic constraints
	private IloRange costLexicoInequality;
	private List<IloRange> uniqueCustomerRouteInequalities;

	public Master(EVRPTW modelData, PricingProblem pricingProblem, CutHandler<EVRPTW, VRPMasterData> cutHandler) {
		super(modelData, pricingProblem, cutHandler, OptimizationSense.MINIMIZE);
	}

	/** Builds the CPLEX problem (a linear program). */
	@Override
	protected VRPMasterData buildModel() {
		IloCplex cplex=null;
		try {
			cplex =new IloCplex(); 			//create CPLEX instance
			cplex.setOut(null); 			//disable CPLEX output
			//			System.out.println(cplex.getVersion());
			cplex.setParam(IloCplex.Param.RootAlgorithm, IloCplex.Algorithm.Primal); //Primal Simplex
			cplex.setParam(IloCplex.Param.Simplex.Tolerances.Feasibility, 1e-9);
			cplex.setParam(IloCplex.Param.RandomSeed, 30);

			obj= cplex.addMinimize();		//objective
			//Partitioning constraints
			int totalLoad = 0;
			visitCustomerConstraints=new IloRange[dataModel.C];
			for(int i=0; i< dataModel.C; i++) {
				visitCustomerConstraints[i] = cplex.addEq(cplex.linearNumExpr(), 1, "visitCustomer_"+(i+1));
				totalLoad+=dataModel.vertices[i].load;
			}

			//Chargers capacity constraints
			chargersCapacityConstraints = new IloRange[dataModel.last_charging_period];
			for (int t = 0; t < dataModel.last_charging_period; t++)
				chargersCapacityConstraints[t] = cplex.addLe(cplex.linearNumExpr(), dataModel.B, "capacity_"+(t+1));

			//Rounded capacity constraint
			this.minimumNumberOfVehicles =  (int) Math.ceil((double) totalLoad/dataModel.Q);
			roundedCapacityInequality = cplex.addGe(cplex.linearNumExpr(), minimumNumberOfVehicles, "capacity inequality");

			//Lexicographic constraints
			costLexicoInequality = null;
			uniqueCustomerRouteInequalities = new ArrayList<>();

		} catch (IloException e) {
			e.printStackTrace();
		}

		//Container for the variables
		Map<PricingProblem,OrderedBiMap<Route, IloNumVar>> varMap=new LinkedHashMap<>();
		varMap.put(pricingProblems.get(0),new OrderedBiMap<>());

		//New data object which will hold data from the Master Problem (including the optimization engine).
		return new VRPMasterData(cplex, pricingProblems.get(0), varMap);
	}

	/** Solves the MP problem (through CPLEX) and returns whether it was solved to optimality. */
	@Override
	protected boolean solveMasterProblem(long timeLimit) throws TimeLimitExceededException {
		try {
			solutionKeeper = new ArrayList<Route>();
			double timeRemaining=Math.max(1,(timeLimit-System.currentTimeMillis())/1000.0);
			masterData.cplex.setParam(IloCplex.DoubleParam.TiLim, timeRemaining); 				//set time limit in seconds
			//solve the model
			if(!masterData.cplex.solve() || masterData.cplex.getStatus()!=IloCplex.Status.Optimal){
				if(masterData.cplex.getCplexStatus()==IloCplex.CplexStatus.AbortTimeLim) 		//Aborted due to time limit
					throw new TimeLimitExceededException();
				else {
					masterData.cplex.exportModel("./results/log/"+dataModel.algorithm+"/"+dataModel.experiment+"/model.lp");

					List<IloRange> constraints = new ArrayList<>();
					constraints.addAll(Arrays.asList(visitCustomerConstraints));
					constraints.addAll(Arrays.asList(chargersCapacityConstraints));
					constraints.addAll(masterData.subsetRowInequalities.values());
					constraints.addAll(masterData.branchingNumberOfVehicles.values());
					constraints.addAll(masterData.branchingChargingTimes.values());
					constraints.add(roundedCapacityInequality);
					if (costLexicoInequality != null){
						constraints.add(costLexicoInequality);
					}
					IloConstraint[] constraintArray = constraints.toArray(new IloConstraint[0]);
					double[] prefs = new double[constraints.size()];
					Arrays.fill(prefs, 1.0);  // Set all preferences to 1.0

					// Now refine the conflict
					masterData.cplex.refineConflict(constraintArray, prefs);
					IloCplex.ConflictStatus[] status = masterData.cplex.getConflict(constraintArray);
					try (PrintWriter writer = new PrintWriter(new FileWriter("./results/log/"+dataModel.algorithm+"/"+dataModel.experiment+"/Infeasibility_report.txt"))) {
						writer.println("Infeasibility Conflict Report:");
						writer.println("--------------------------------");

						for (int i = 0; i < constraintArray.length; i++) { writer.println("Constraint: " + constraintArray[i].getName() + " | Status: " + status[i]); }

						writer.println("--------------------------------");
						writer.println("End of report.");
					} catch (IOException e) { System.err.println("Failed to write conflict report: " + e.getMessage()); }

					throw new RuntimeException(dataModel.algorithm + " - " + dataModel.experiment + " - " + dataModel.instanceName + ". Master problem solve failed! Status: "+ masterData.cplex.getStatus());
				}

			}else{
				masterData.objectiveValue= masterData.cplex.getObjValue();
				//Print solution
				List<Route> solution=getSolution();
				if (dataModel.print_log) {
					logger.debug("Objective: "+ masterData.objectiveValue);
					logger.debug("Number of columns: " + masterData.getNrColumns() + " Number of SRC separated: " + masterData.subsetRowInequalities.size());
					logger.debug("Number of vehicle branches: " + masterData.branchingNumberOfVehicles.size() + " Number of charging time branches: " + masterData.branchingChargingTimes.size());
					logger.debug("Columns (only non-zero columns are returned):");
					for(Route route: solution)
						logger.debug(route.toString());
					
					/* logger.debug("Printing dual variables");
					logger.debug("Capacity constraint: " + String.valueOf(masterData.cplex.getDual(roundedCapacityInequality)));
					logger.debug("Customer constraints: " + masterData.cplex.getDuals(visitCustomerConstraints).toString());
					for(int i=0; i< dataModel.C; i++) {
						logger.debug("Customer " + String.valueOf(i+1) + ": " +masterData.cplex.getDual(visitCustomerConstraints[i]));
					}
					logger.debug("Charging periods constraints:" + masterData.cplex.getDuals(chargersCapacityConstraints).toString());
					for (int t = 0; t < dataModel.last_charging_period; t++) {
						logger.debug("Period " + String.valueOf(t+1) + ": " + masterData.cplex.getDual(chargersCapacityConstraints[t]));
					}
					logger.debug("SRCs:");
					for(SubsetRowInequality subsetRowInequality: masterData.subsetRowInequalities.keySet()) {
						double dual = masterData.cplex.getDual(masterData.subsetRowInequalities.get(subsetRowInequality));
						logger.debug(subsetRowInequality.toString() + ": " + dual);
					} */
					
				}
			}
		} catch (IloException e) {
			e.printStackTrace();
		}
		return true;
	}

	/**
	 * We store the dual information in the pricing problem object.
	 * This method is invoked after the MP has been solved. 
	 */
	@Override
	public void initializePricingProblem(PricingProblem pricingProblem){
		try {

			pricingProblem.branchesOnChargingTimes = masterData.branchingChargingTimes.keySet();
			double[] dualsPartition= masterData.cplex.getDuals(visitCustomerConstraints);
			double[] dualsCapacity = masterData.cplex.getDuals(chargersCapacityConstraints);
			double[] dualsSRC = new double[masterData.subsetRowInequalities.size()];

			ArrayList<SubsetRowInequality> SRCToConsider = new ArrayList<SubsetRowInequality>();
			int s = 0;
			for(SubsetRowInequality subsetRowInequality: masterData.subsetRowInequalities.keySet()) {
				double dual = masterData.cplex.getDual(masterData.subsetRowInequalities.get(subsetRowInequality));
				if(dual<0) {
					SRCToConsider.add(subsetRowInequality);
					dualsSRC[s] = dual;
					for(int i: subsetRowInequality.cutSet) dataModel.vertices[i].SRCIndices.add(s);
					s++;
				}
			}

			double [] duals = new double[dualsPartition.length + dualsCapacity.length+ s + pricingProblem.branchesOnChargingTimes.size()];  //resultant array of size first array and second array  
			for (int i = 0; i < dualsPartition.length; i++) 
				duals[i] = dualsPartition[i];

			for (int i = 0; i < dualsCapacity.length; i++) 
				duals[dualsPartition.length+i] = dualsCapacity[i];

			for (int i = 0; i < s; i++)
				duals[dualsPartition.length+dualsCapacity.length+i] = dualsSRC[i];

			pricingProblem.subsetRowCuts = SRCToConsider;

			int i = 0;
			for(IloRange branching: masterData.branchingChargingTimes.values()) {
				duals[dualsPartition.length+dualsCapacity.length+s+i] = masterData.cplex.getDual(branching);
				i++;
			}

			double dualConstant = 0; //constant dual values (not depending on the arc)
			dualConstant+=masterData.cplex.getDual(roundedCapacityInequality);

			// branching on vehicles duals
			for(NumberVehiclesInequalities branching: masterData.branchingNumberOfVehicles.keySet())
				dualConstant+=masterData.cplex.getDual(masterData.branchingNumberOfVehicles.get(branching));

			pricingProblem.initPricingProblem(duals, dualConstant);

		} catch (IloException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Function that adds a new column to the CPLEX problem.
	 * This method is invoked when a Pricing Problem generated a new column.
	 */
	@Override
	public void addColumn(Route column) {
		try {

			// register column with objective
			IloColumn iloColumn= masterData.cplex.column(obj,column.cost);

			// register column with partitioning constraint
			for(int i: column.route.keySet())
				iloColumn=iloColumn.and(masterData.cplex.column(visitCustomerConstraints[i-1], column.route.get(i)));

			// register column with chargers capacity constraints
			for (int t = column.initialChargingTime; t <= (column.initialChargingTime+ column.chargingTime-1); t++)
				iloColumn=iloColumn.and(masterData.cplex.column(chargersCapacityConstraints[t-1], 1));

			// register (artificial) column with rounded capacity inequality and branching decisions (vehicles)
			if(column.isArtificialColumn) {
				iloColumn=iloColumn.and(masterData.cplex.column(roundedCapacityInequality, this.minimumNumberOfVehicles));
				for (NumberVehiclesInequalities branch: masterData.branchingNumberOfVehicles.keySet()) {
					IloRange branchConstraint = masterData.branchingNumberOfVehicles.get(branch);
					if(!branch.lessThanOrEqual) iloColumn = iloColumn.and(masterData.cplex.column(branchConstraint,branch.coefficient));
				}
			}

			if(!column.isArtificialColumn) {

				// register column with rounded capacity inequality
				iloColumn=iloColumn.and(masterData.cplex.column(roundedCapacityInequality, 1));


				// register the column with Subset Row Inequalities Constraints
				for(SubsetRowInequality subsetRowInequality: masterData.subsetRowInequalities.keySet()) {
					// check the number of visits to the customers in the triplet
					int coeff = getCoefficient(column, subsetRowInequality);
					if(coeff>0){
						IloRange subsetRowInequalityConstraint=masterData.subsetRowInequalities.get(subsetRowInequality);
						iloColumn = iloColumn.and(masterData.cplex.column(subsetRowInequalityConstraint, coeff));
					}
				}

				// register the column with the branching decision (number of vehicles)
				for (NumberVehiclesInequalities branch: masterData.branchingNumberOfVehicles.keySet()) {
					IloRange branchConstraint = masterData.branchingNumberOfVehicles.get(branch);
					iloColumn = iloColumn.and(masterData.cplex.column(branchConstraint,1));
				}

				// register the column with branching decision (charging time)
				for (ChargingTimeInequality branch: masterData.branchingChargingTimes.keySet()) {
					IloRange branchConstraint = masterData.branchingChargingTimes.get(branch);
					if(branch.startCharging && column.initialChargingTime==branch.timestep) {
						iloColumn = iloColumn.and(masterData.cplex.column(branchConstraint, 1));
					}else if(!branch.startCharging && (column.initialChargingTime+column.chargingTime-1)==branch.timestep) {
						iloColumn = iloColumn.and(masterData.cplex.column(branchConstraint, 1));
					}
				}
			}

			// create the variable and store it
			IloNumVar var= masterData.cplex.numVar(iloColumn, 0, Double.MAX_VALUE, "x_"+masterData.getNrColumns());
			masterData.cplex.add(var);
			masterData.addColumn(column, var);
		} catch (IloException e) {
			e.printStackTrace();
		}
	}


	/** If a violated inequality has been found add it to the MP. */
	private void addCut(SubsetRowInequality subsetRowInequality){

		if(!(subsetRowInequality instanceof SubsetRowInequality))
			throw new RuntimeException("Error only Subset-Row Cuts are considered");

		if(masterData.subsetRowInequalities.containsKey(subsetRowInequality))
			throw new RuntimeException("Error, duplicate subset-row cut is being generated! This cut should already exist in the master problem: "+subsetRowInequality);
		// create the inequality in CPLEX
		try {
			IloLinearNumExpr expr=masterData.cplex.linearNumExpr();
			// register the columns with this constraint.
			for(Route route: masterData.getColumnsForPricingProblemAsList(masterData.pricingProblem)){
				if(route.isArtificialColumn) continue;
				// check the number of visits to the customers in the triplet
				int coeff = getCoefficient(route, subsetRowInequality);
				if(coeff>0){
					IloNumVar var=masterData.getVar(masterData.pricingProblem,route);
					expr.addTerm(coeff, var);
				}
			}
			IloRange subsetRowConstraint = masterData.cplex.addLe(expr, 1, "subsetRow_"+Arrays.toString(subsetRowInequality.cutSet));
			masterData.subsetRowInequalities.put(subsetRowInequality, subsetRowConstraint);

		} catch (IloException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Computes the coefficient of a route in SRC.
	 * @param route for which the coefficient is calculated.
	 * @param subsetRowInequality considered.
	 */
	public int getCoefficient(Route route, SubsetRowInequality subsetRowInequality) {
		int visits = 0;
		for(int i: subsetRowInequality.cutSet)
			visits+=route.route.getOrDefault(i, 0);
		return (int) Math.floor(0.5*visits);
	}

	/** Returns the solution, i.e columns with non-zero values in the CPLEX problem. */
	@Override
	public List<Route> getSolution() {
		List<Route> solution=new ArrayList<>();
		if(!solutionKeeper.isEmpty()) return solutionKeeper;
		try {
			Route[] routes=masterData.getVarMap().getKeysAsArray(new Route[masterData.getNrColumns()]);
			IloNumVar[] vars=masterData.getVarMap().getValuesAsArray(new IloNumVar[masterData.getNrColumns()]);
			double[] values= masterData.cplex.getValues(vars);

			// iterate over each column and add it to the solution if it has a non-zero value
			for(int i=0; i<routes.length; i++){
				Route clone_route = routes[i].clone();
				clone_route.value=values[i];
				if(values[i]>=config.PRECISION){
					solution.add(clone_route);
				}
			}
		} catch (IloException e) { //the MP has been modified and needs to be solved again
			// TODO Auto-generated catch block
			try {
				masterData.cplex.solve();
				Route[] routes=masterData.getVarMap().getKeysAsArray(new Route[masterData.getNrColumns()]);
				IloNumVar[] vars=masterData.getVarMap().getValuesAsArray(new IloNumVar[masterData.getNrColumns()]);
				double[] values= masterData.cplex.getValues(vars);

				// iterate over each column and add it to the solution if it has a non-zero value
				for(int i=0; i<routes.length; i++){
					Route clone_route = routes[i].clone();
					clone_route.value=values[i];
					if(values[i]>=config.PRECISION){
						solution.add(clone_route);
					}
				}
			} catch (IloException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		} 
		solutionKeeper = solution;
		return solution;
	}

	/** Close the CPLEX problem. */
	@Override
	public void close() {
		masterData.cplex.end();
	}

	/** Print the solution. */
	@Override
	public void printSolution() {
		System.out.println("Master solution:");
		for(Route r : this.getSolution())
			System.out.println(r);
	}

	/**
	 * Checks whether there are any violated SR inequalities, thereby invoking the cut handler.
	 * @return true if violated inequalities have been found (and added to the master problem).
	 */
	@Override
	public boolean hasNewCuts(){
		masterData.routeValueMap = new HashMap<Route, Double>();
		boolean isInteger = true;
		for(Route route : this.getSolution()) {
			if(route.isArtificialColumn) return false;
			masterData.routeValueMap.put(route, route.value);
			if(route.value>0+dataModel.precision && route.value<1-dataModel.precision) isInteger=false;
		}
		if(isInteger) return false;

		return super.hasNewCuts();
	}

	/**
	 * Listen to branching decisions
	 */
	@Override
	public void branchingDecisionPerformed(BranchingDecision bd) {
		// for simplicity, we simply destroy the master problem and rebuild it. Of course, something more sophisticated may be used which retains the master problem.
		Set<NumberVehiclesInequalities> vehiclesInequalities = masterData.branchingNumberOfVehicles.keySet(); 	//keep branching decisions
		Set<ChargingTimeInequality> chargingInequalities = masterData.branchingChargingTimes.keySet(); 			//keep branching decisions

		this.close(); 																							//close the old CPLEX model
		masterData=this.buildModel(); 																			//create a new model without any columns
		cutHandler.setMasterData(masterData); 																	//inform the cutHandler about the new master model
		for(NumberVehiclesInequalities inequality: vehiclesInequalities) addBranchingOnVehichlesInequality(inequality);
		for(ChargingTimeInequality inequality: chargingInequalities) addChargingTimeInequality(inequality);


		if (bd instanceof BranchVehiclesDown) {
			BranchVehiclesDown branching = (BranchVehiclesDown) bd;
			addBranchingOnVehichlesInequality(branching.inequality);
			for(AbstractInequality src: branching.poolOfCuts) addCut((SubsetRowInequality) src);
		}
		else if (bd instanceof BranchVehiclesUp) {
			BranchVehiclesUp branching = (BranchVehiclesUp) bd;
			addBranchingOnVehichlesInequality(branching.inequality);
			for(AbstractInequality src: branching.poolOfCuts) addCut((SubsetRowInequality) src);
		}
		else if(bd instanceof FixArc) {
			FixArc fixArcDecision = (FixArc) bd;
			for(AbstractInequality src: fixArcDecision.poolOfCuts) addCut((SubsetRowInequality) src);
		}
		else if(bd instanceof RemoveArc) {
			RemoveArc removeArcDecision= (RemoveArc) bd;
			for(AbstractInequality src: removeArcDecision.poolOfCuts) addCut((SubsetRowInequality) src);
		}
		else if (bd instanceof BranchInitialChargingTimeDown) {
			BranchInitialChargingTimeDown branching = (BranchInitialChargingTimeDown) bd;
			addChargingTimeInequality(branching.inequality);
			for(AbstractInequality src: branching.poolOfCuts) addCut((SubsetRowInequality) src);
		}
		else if (bd instanceof BranchInitialChargingTimeUp) {
			BranchInitialChargingTimeUp branching = (BranchInitialChargingTimeUp) bd;
			addChargingTimeInequality(branching.inequality);
			for(AbstractInequality src: branching.poolOfCuts) addCut((SubsetRowInequality) src);
		}
		else if (bd instanceof BranchEndChargingTimeDown) {
			BranchEndChargingTimeDown branching = (BranchEndChargingTimeDown) bd;
			addChargingTimeInequality(branching.inequality);
			for(AbstractInequality src: branching.poolOfCuts) addCut((SubsetRowInequality) src);
		}
		else if (bd instanceof BranchEndChargingTimeUp) {
			BranchEndChargingTimeUp branching = (BranchEndChargingTimeUp) bd;
			addChargingTimeInequality(branching.inequality);
			for(AbstractInequality src: branching.poolOfCuts) addCut((SubsetRowInequality) src);
		}
	}

	/**
	 * Undo branching decisions during backtracking in the Branch-and-Price tree
	 */
	@Override
	public void branchingDecisionReversed(BranchingDecision bd) {
		if (bd instanceof BranchVehiclesDown) {
			BranchVehiclesDown branching = (BranchVehiclesDown) bd;
			masterData.branchingNumberOfVehicles.remove(branching.inequality);
		}
		else if (bd instanceof BranchVehiclesUp) {
			BranchVehiclesUp branching = (BranchVehiclesUp) bd;
			masterData.branchingNumberOfVehicles.remove(branching.inequality);
		}
		else if (bd instanceof BranchInitialChargingTimeDown) { 
			BranchInitialChargingTimeDown branching = (BranchInitialChargingTimeDown) bd;
			masterData.branchingChargingTimes.remove(branching.inequality);
		}
		else if (bd instanceof BranchInitialChargingTimeUp) {
			BranchInitialChargingTimeUp branching = (BranchInitialChargingTimeUp) bd;
			masterData.branchingChargingTimes.remove(branching.inequality);
		}
		else if (bd instanceof BranchEndChargingTimeDown) {
			BranchEndChargingTimeDown branching = (BranchEndChargingTimeDown) bd;
			masterData.branchingChargingTimes.remove(branching.inequality);
		}
		else if (bd instanceof BranchEndChargingTimeUp) {
			BranchEndChargingTimeUp branching = (BranchEndChargingTimeUp) bd;
			masterData.branchingChargingTimes.remove(branching.inequality);
		}
	}

	/**
	 * Creates a branching decision constraint (when branching on the number of vehicles used)
	 */
	public void addBranchingOnVehichlesInequality(NumberVehiclesInequalities inequality) {
		try {
			IloRange branchingConstraint;
			if (inequality.lessThanOrEqual) branchingConstraint = masterData.cplex.addLe(masterData.cplex.linearNumExpr(), inequality.coefficient, "branching_"+inequality.toString());
			else branchingConstraint = masterData.cplex.addGe(masterData.cplex.linearNumExpr(), inequality.coefficient, "branching_"+inequality.toString());
			masterData.branchingNumberOfVehicles.put(inequality, branchingConstraint);
		}
		catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} //new constraint
	}

	/**
	 * Creates a branching decision constraint (when branching on the number of vehicles used)
	 */
	public void addChargingTimeInequality(ChargingTimeInequality inequality) {
		try {
			IloRange branchingConstraint;
			if (inequality.lessThanOrEqual) branchingConstraint = masterData.cplex.addLe(masterData.cplex.linearNumExpr(), inequality.coefficient, "branching_"+inequality.toString());
			else branchingConstraint = masterData.cplex.addGe(masterData.cplex.linearNumExpr(), inequality.coefficient, "branching_"+inequality.toString());
			masterData.branchingChargingTimes.put(inequality, branchingConstraint);
		}
		catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} //new constraint
	}

	/**
	 * Obtains a lower bound on the MP problem. 
	 * You need an upper bound U on the number of vehicles used in an optimal solution because LB can be computed as: 
	 * z_RMP + U * z_PP, where z_RMP is the objective value of the current RMP solution
	 */
	@Override
	public double getBoundComponent() {
		if(pricingProblems.get(0).bestReducedCost==-Double.MAX_VALUE) { //do nothing
			return 0;
		}else {
			int maxK = dataModel.C;
			for (NumberVehiclesInequalities branching: masterData.branchingNumberOfVehicles.keySet())
				if(branching.lessThanOrEqual && maxK>branching.coefficient) maxK=branching.coefficient;

			double boundOnObjective = this.getObjective()+this.pricingProblems.get(0).bestReducedCost*maxK;
			if (dataModel.print_log) {
				if(boundOnObjective>0) logger.debug("Computed LB: " + boundOnObjective + " z_RMP="+ this.getObjective() + " U=" +maxK +" z_PP=" +this.pricingProblems.get(0).bestReducedCost);
			}
			pricingProblems.get(0).bestReducedCost = -Double.MAX_VALUE;
			return boundOnObjective;
		}
	}

	public VRPMasterData getMasterData(){
		return masterData;
	}

	public Master copy(){

		return new Master(this.dataModel, this.pricingProblems.get(0), this.cutHandler);
	}

	public void addColumnsDepletion(Map<int[], List<Route>> cols){
		
		for (int[] route: cols.keySet()) {
			for (Route column: cols.get(route)){
				try {

					// register column with objective
					IloColumn iloColumn= masterData.cplex.column(obj, column.departureTime-(column.initialChargingTime+column.chargingTime));
		
					// register column with partitioning constraint
					for(int i: column.route.keySet())
						iloColumn=iloColumn.and(masterData.cplex.column(visitCustomerConstraints[i-1], column.route.get(i)));
		
					// register column with chargers capacity constraints
					for (int t = column.initialChargingTime; t <= (column.initialChargingTime+ column.chargingTime-1); t++)
						iloColumn=iloColumn.and(masterData.cplex.column(chargersCapacityConstraints[t-1], 1));
		
					// register (artificial) column with rounded capacity inequality and branching decisions (vehicles)
					if(column.isArtificialColumn) {
						iloColumn=iloColumn.and(masterData.cplex.column(roundedCapacityInequality, this.minimumNumberOfVehicles));
						for (NumberVehiclesInequalities branch: masterData.branchingNumberOfVehicles.keySet()) {
							IloRange branchConstraint = masterData.branchingNumberOfVehicles.get(branch);
							if(!branch.lessThanOrEqual) iloColumn = iloColumn.and(masterData.cplex.column(branchConstraint,branch.coefficient));
						}
					}
		
					if(!column.isArtificialColumn) {
		
						// register column with rounded capacity inequality
						iloColumn=iloColumn.and(masterData.cplex.column(roundedCapacityInequality, 1));
		
						// register the column with Subset Row Inequalities Constraints
						for(SubsetRowInequality subsetRowInequality: masterData.subsetRowInequalities.keySet()) {
							// check the number of visits to the customers in the triplet
							int coeff = getCoefficient(column, subsetRowInequality);
							if(coeff>0){
								IloRange subsetRowInequalityConstraint=masterData.subsetRowInequalities.get(subsetRowInequality);
								iloColumn = iloColumn.and(masterData.cplex.column(subsetRowInequalityConstraint, coeff));
							}
						}
		
						// register the column with the branching decision (number of vehicles)
						for (NumberVehiclesInequalities branch: masterData.branchingNumberOfVehicles.keySet()) {
							IloRange branchConstraint = masterData.branchingNumberOfVehicles.get(branch);
							iloColumn = iloColumn.and(masterData.cplex.column(branchConstraint,1));
						}
		
						// register the column with branching decision (charging time)
						for (ChargingTimeInequality branch: masterData.branchingChargingTimes.keySet()) {
							IloRange branchConstraint = masterData.branchingChargingTimes.get(branch);
							if(branch.startCharging && column.initialChargingTime==branch.timestep) {
								iloColumn = iloColumn.and(masterData.cplex.column(branchConstraint, 1));
							}else if(!branch.startCharging && (column.initialChargingTime+column.chargingTime-1)==branch.timestep) {
								iloColumn = iloColumn.and(masterData.cplex.column(branchConstraint, 1));
							}
						}
					}
		
					// create the variable and store it
					IloNumVar var= masterData.cplex.numVar(iloColumn, 0, Double.MAX_VALUE, "x_"+masterData.getNrColumns());
					masterData.cplex.add(var);
					masterData.addColumn(column, var);
				} catch (IloException e) {
					e.printStackTrace();
				}

			}
		}
	}

	public double minimizeBatteryDepletion(List<Route> solution, long timeLimit, List<Route> cols, List<AbstractInequality> src_list, Map<NumberVehiclesInequalities, IloRange> vehic_branches_map, Map<ChargingTimeInequality, IloRange> time_branches_map, double minCost){

		Set<NumberVehiclesInequalities> vehiclesInequalities = vehic_branches_map.keySet(); 	//keep branching decisions
		Set<ChargingTimeInequality> chargingInequalities = time_branches_map.keySet();
		
		cutHandler.setMasterData(masterData);
		
		// Add all constraints added throughout the BPC root path
		for(NumberVehiclesInequalities inequality: vehiclesInequalities) addBranchingOnVehichlesInequality(inequality);
		for(ChargingTimeInequality inequality: chargingInequalities) addChargingTimeInequality(inequality);
		for(AbstractInequality src: src_list) addCut((SubsetRowInequality) src);

		// Add necessary columns
		List<int[]> unique_customer_routes = retrieve_unique_customer_routes(solution);
		Map<int[], List<Route>> columns_to_add = filter_columns_by_customer_routes(unique_customer_routes, cols);
		this.addColumnsDepletion(columns_to_add);

		IloLinearNumExpr expr;
		double new_cost = 0;

		try {

			expr = masterData.cplex.linearNumExpr();
			for(Route route: cols){
				IloNumVar var=masterData.getVar(masterData.pricingProblem, route);
				expr.addTerm(route.cost, var);
			}
			//logger.debug("MP obj inside minimizeBatteryDepletion: "+minCost+" - "+Math.round(minCost));
			costLexicoInequality = masterData.cplex.addLe(expr, Math.round(minCost), "minCost");
			//logger.debug("Cost constraint before solving: "+"<="+ cost_constraint.getUB());

			//Impose use of unique customer routes
			int ix = 0;
			for (int[] route: unique_customer_routes){
				IloLinearNumExpr lhs = masterData.cplex.linearNumExpr();
				for(Route column: columns_to_add.get(route)){
					IloNumVar var=masterData.getVar(masterData.pricingProblem, column);
					lhs.addTerm(1, var);
				}
				uniqueCustomerRouteInequalities.add(masterData.cplex.addGe(lhs, 1, "uniqueRoute"+ix));
				ix += 1;
			}
			
			masterData.cplex.setParam(IloCplex.Param.Simplex.Tolerances.Feasibility, 1e-6);
			this.masterData.optimal = this.solveMasterProblem(timeLimit);
			new_cost = masterData.cplex.getValue(expr);
			
			masterData.cplex.exportModel("./results/log/"+dataModel.algorithm+"/"+dataModel.experiment+"/model.lp");
			/* double lhs = masterData.cplex.getValue(cost_constraint.getExpr());
			masterData.cplex.writeSolution("./results/log/"+dataModel.algorithm+"/"+dataModel.experiment+"/solution"+lhs+".lp");
			logger.debug("Master optimal: "+((boolean)(masterData.cplex.getStatus()==IloCplex.Status.Optimal)));
			logger.debug("Cost constraint after solving: "+lhs+"<="+cost_constraint.getUB()); */
	

		} catch (TimeLimitExceededException e) {
			System.out.println("Time limit exceeded: " + e.getMessage());
		} catch (IloException e) {
			System.out.println("CPLEX encountered an error: " + e.getMessage());
		}

		return Math.round(new_cost);
		
	}

	private List<int[]> retrieve_unique_customer_routes(List<Route> solution){

		List<int[]> unique_routes = new ArrayList<>();
		Set<List<Integer>> seen = new HashSet<>();

		for (Route column: solution) {
			int[] arr = column.routeSequence;
			List<Integer> asList = Arrays.stream(arr).boxed().collect(Collectors.toList());
			if (seen.add(asList)) {
				unique_routes.add(arr);
			}
		}

		return unique_routes;
	}

	private Map<int[], List<Route>> filter_columns_by_customer_routes(List<int[]> routes, List<Route> cols){

		Map<int[], List<Route>> columns_to_add = new LinkedHashMap<>();
		for (Route column: cols){
			for (int[] route: routes){
				if (Arrays.equals(column.routeSequence, route)) {
					
					List<Route> cols_subset = columns_to_add.get(route);
					if (cols_subset == null) columns_to_add.put(route, new ArrayList<>(Collections.singletonList(column)));
					else cols_subset.add(column);
				
				}
			}
		}

		return columns_to_add;
	}

}
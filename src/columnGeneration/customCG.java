package columnGeneration;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.colgenMain.ColGen;
import org.jorlib.frameworks.columnGeneration.io.TimeLimitExceededException;
import org.jorlib.frameworks.columnGeneration.master.AbstractMaster;
import org.jorlib.frameworks.columnGeneration.master.MasterData;
import org.jorlib.frameworks.columnGeneration.master.OptimizationSense;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractInequality;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;
import org.jorlib.frameworks.columnGeneration.pricing.DefaultPricingProblemSolverFactory;
import org.jorlib.frameworks.columnGeneration.pricing.PricingProblemBundle;
import org.jorlib.frameworks.columnGeneration.pricing.PricingProblemManager;

import branchAndPrice.RemoveArc;
import branchAndPrice.ExtendBAPNotifier;
import branchAndPrice.NumberVehiclesInequalities;
import ilog.concert.IloColumn;
import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import model.EVRPTW;

/**
 * This class is a custom implementation of the ColGen class (jORLib)
 * It is implemented to compute a lower bound on the master problem
 */
public class customCG extends ColGen<EVRPTW, Route, PricingProblem> {

	private final ExtendBAPNotifier extendedNotifier;
	public ArrayList<Route> incumbentSolution = new ArrayList<Route>(); 	//stores the incumbent solution found throughout the CG
	public int incumbentSolutionObjective = (int) Double.MAX_VALUE; 		// stores the incumbent solution objective found throughout the CG
	
	public int BBnodeID;
	public boolean needsChargingBranchingPricing;
	public int contExact;

	public OptimalSolutionMemory solutionMemory;
	private boolean masterSolutionIsInteger;
	private boolean hasExceededPricingSoftThreshold = false;
	private double gapReduction = 0;

	public List<BranchingDecision> branchingFRC;

	private final Map<Class<? extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem>>, PricingProblemBundle<EVRPTW, Route, PricingProblem>> pricingProblemBundles;
	private final customPricingProblemManager cPricingProblemManager;

	private static final Map<Class<? extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem>>, Boolean> solverCapabilities = new HashMap<>();
	static {
		solverCapabilities.put(HeuristicLabelingThirdPricingProblemSolver.class, false);
		solverCapabilities.put(HeuristicLabelingPricingProblemSolver.class, false);
		solverCapabilities.put(HeuristicLabelingSecondPricingProblemSolver.class, false);
		solverCapabilities.put(HeuristicMinCostLabelingPricingProblemSolver.class, false);

		solverCapabilities.put(CBHeuristicSecondPricingProblemSolver.class, true);
		solverCapabilities.put(CBHeuristicMinCostPricingProblemSolver.class, true);
	}

	public customCG(EVRPTW dataModel, AbstractMaster<EVRPTW, Route, PricingProblem, ? extends MasterData> master,
			PricingProblem pricingProblem,
			List<Class<? extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem>>> solvers,
			List<Route> initSolution, int cutoffValue, double boundOnMasterObjective, int nodeID, boolean needsCB, ExtendBAPNotifier notifier) {
		super(dataModel, master, pricingProblem, solvers, initSolution, cutoffValue, boundOnMasterObjective);
		this.BBnodeID = nodeID;
		this.needsChargingBranchingPricing = needsCB;
		this.extendedNotifier = notifier;
		this.branchingFRC = new ArrayList<BranchingDecision>();
		
		this.pricingProblemManager.close();
		this.pricingProblemBundles = new HashMap<Class<? extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem>>, PricingProblemBundle<EVRPTW, Route, PricingProblem>>();

		for(Class<? extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem>> solverClass : solvers) {
			DefaultPricingProblemSolverFactory<EVRPTW, Route, PricingProblem> factory = new DefaultPricingProblemSolverFactory<EVRPTW, Route, PricingProblem>(solverClass, dataModel);
			PricingProblemBundle<EVRPTW, Route, PricingProblem> bundle = new PricingProblemBundle<EVRPTW, Route, PricingProblem>(solverClass, pricingProblems, factory);
			pricingProblemBundles.put(solverClass, bundle);
		}

		this.cPricingProblemManager = new customPricingProblemManager(pricingProblems, pricingProblemBundles);
	}

	public customCG(EVRPTW dataModel, AbstractMaster<EVRPTW, Route, PricingProblem, ? extends MasterData> master,
			List<PricingProblem> pricingProblems,
			List<Class<? extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem>>> solvers,
			PricingProblemManager<EVRPTW, Route, PricingProblem> pricingProblemManager, List<Route> initSolution,
			int cutoffValue, double boundOnMasterObjective, int nodeID, boolean needsCB, ExtendBAPNotifier notifier) {
		super(dataModel, master, pricingProblems, solvers, pricingProblemManager, initSolution, cutoffValue, boundOnMasterObjective);
		this.BBnodeID = nodeID;
		this.needsChargingBranchingPricing = needsCB;
		this.extendedNotifier = notifier;
		this.branchingFRC = new ArrayList<BranchingDecision>();

		this.pricingProblemManager.close();
		this.pricingProblemBundles = new HashMap<Class<? extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem>>, PricingProblemBundle<EVRPTW, Route, PricingProblem>>();

		for(Class<? extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem>> solverClass : solvers) {
			DefaultPricingProblemSolverFactory<EVRPTW, Route, PricingProblem> factory = new DefaultPricingProblemSolverFactory<EVRPTW, Route, PricingProblem>(solverClass, dataModel);
			PricingProblemBundle<EVRPTW, Route, PricingProblem> bundle = new PricingProblemBundle<EVRPTW, Route, PricingProblem>(solverClass, pricingProblems, factory);
			pricingProblemBundles.put(solverClass, bundle);
		}

		this.cPricingProblemManager = new customPricingProblemManager(pricingProblems, pricingProblemBundles);
	}

	public customCG(EVRPTW arg0, AbstractMaster<EVRPTW, Route, PricingProblem, ? extends MasterData> arg1,
			List<PricingProblem> arg2,
			List<Class<? extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem>>> arg3, List<Route> arg4,
			int arg5, double arg6, int arg7, boolean arg8, ExtendBAPNotifier arg9) {
		super(arg0, arg1, arg2, arg3, arg4, arg5, arg6);
		this.BBnodeID = arg7;
		this.needsChargingBranchingPricing = arg8;
		this.extendedNotifier = arg9;
		this.branchingFRC = new ArrayList<BranchingDecision>();
		
		this.pricingProblemManager.close();
		this.pricingProblemBundles =  new HashMap<Class<? extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem>>, PricingProblemBundle<EVRPTW, Route, PricingProblem>>();

		for(Class<? extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem>> solverClass : solvers) {
			DefaultPricingProblemSolverFactory<EVRPTW, Route, PricingProblem> factory = new DefaultPricingProblemSolverFactory<EVRPTW, Route, PricingProblem>(solverClass, dataModel);
			PricingProblemBundle<EVRPTW, Route, PricingProblem> bundle = new PricingProblemBundle<EVRPTW, Route, PricingProblem>(solverClass, pricingProblems, factory);
			pricingProblemBundles.put(solverClass, bundle);
		}

		this.cPricingProblemManager = new customPricingProblemManager(pricingProblems, pricingProblemBundles);
	}

	@Override
	protected double calculateBoundOnMasterObjective(
			Class<? extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem>> solver) {
		// TODO Auto-generated method stub
		return master.getBoundComponent(); //obtains a lower bound on the master problem objective
	}


	/**
	 * Solve the Column Generation problem. First the master problem is solved. Next the pricing problems(s) is (are) solved. To solve the pricing problems, the pricing
	 * solvers are invoked one by one in a hierarchical fashion. First the first solver is invoked to solve the pricing problems. Any new columns generated are immediately returned.
	 * If it fails to find columns, the next solver is invoked and so on. If the pricing problem discovers new columns, they are added to the master problem and the method continues
	 * with the next column generation iteration.<br>
	 * If no new columns are found, the method checks for violated inequalities. If there are violated inequalities, they are added to the master problem and the method continues with the
	 * next column generation iteration.<br>
	 * The solve procedure terminates under any of the following conditions:
	 * <ol>
	 * <li>the solver could not identify new columns</li>
	 * <li>Time limit exceeded</li>
	 * <li>The bound on the best attainable solution to the master problem is worse than the cutoff value. Assuming that the master is a minimization problem, the Colgen procedure is terminated if {@code ceil(boundOnMasterObjective) >= cutoffValue}</li>
	 * <li>The solution to the master problem is provable optimal, i.e the bound on the best attainable solution to the master problem equals the solution of the master problem.</li>
	 * </ol>
	 * @param timeLimit Future point in time (ms) by which the procedure should be finished. Should be defined as: {@code System.currentTimeMilis()+<desired runtime>}
	 * @throws TimeLimitExceededException Exception is thrown when time limit is exceeded
	 */
	@Override
	public void solve(long timeLimit) throws TimeLimitExceededException{
		//Set time limit pricing problems
		cPricingProblemManager.setTimeLimit(timeLimit);
		colGenSolveTime=System.currentTimeMillis();
		this.incumbentSolutionObjective = this.cutoffValue;

		boolean foundNewColumns=false; 				//identify whether the pricing problem generated new columns
		boolean hasNewCuts; 						//identify whether the master problem violates any valid inequalities
		notifier.fireStartCGEvent();

		int cg_iterations = 0;

		if (this.BBnodeID == 0) { this.contExact = 0; dataModel.rollbackBaseLine = 0; }
		dataModel.rollbackTrigger = false;
		dataModel.cut_iterations = 1;

		do{
			cg_iterations ++;

			nrOfColGenIterations++;
			hasNewCuts=false;

			//Solve the master
			this.invokeMaster(timeLimit);
			if (objectiveMasterProblem<boundOnMasterObjective-dataModel.precision) {
				throw new UnsupportedOperationException("Problem with LB");
			}

			//We can stop when the optimality gap is closed. We still need to check for violated inequalities though.
			if (Math.abs(objectiveMasterProblem - boundOnMasterObjective)<config.PRECISION){
				//Check whether there are inequalities. Otherwise potentially an infeasible integer solution (e.g. TSP solution with subtours) might be returned.
				if (dataModel.CUTSENABLED){
					long time=System.currentTimeMillis();
					hasNewCuts=master.hasNewCuts();
					masterSolveTime+=(System.currentTimeMillis()-time); //Generating inequalities is considered part of the master problem
					if (hasNewCuts) continue;
					else break;
				} else break;
			}

			//Solve the pricing problem and possibly update the bound on the master problem objective
			List<Route> newColumns=this.invokePricingProblems(timeLimit); //List containing new columns generated by the pricing problem
			foundNewColumns=!newColumns.isEmpty();

			//Check whether the boundOnMasterObjective exceeds the cutoff value
			if (dataModel.rollbackTrigger){
				this.perform_rollback(solutionMemory); break; }
			else if (boundOnMasterExceedsCutoffValue()) break;
			else if (System.currentTimeMillis() >= timeLimit){ 			//check whether we are still within the timeLimit
				notifier.fireTimeLimitExceededEvent();
				throw new TimeLimitExceededException(); }
			else if (dataModel.CUTSENABLED && !foundNewColumns){ 		//check for inequalities. This can only be done if the master problem hasn't changed (no columns can be added).

				// Check if the gap reduction was enough.
				// In case the reduction was bad, break the Column and Cut Generation to branch directly
				if (dataModel.cut_iterations > 1 && this.hasExceededPricingSoftThreshold && this.gapReduction < dataModel.gapReductionRequirement){
					extendedNotifier.fireGapReductionEvent(this.gapReduction);
					break; }

				dataModel.cut_iterations ++;
				long time = System.currentTimeMillis();
				hasNewCuts = master.hasNewCuts();
				masterSolveTime += (System.currentTimeMillis()-time);	//generating inequalities is considered part of the master problem
				
				this.update_solution_memory();
				dataModel.cleanSRCs();	// MODIFICATION
			}

		} while (foundNewColumns || hasNewCuts);
		
		colGenSolveTime = System.currentTimeMillis() - colGenSolveTime;
		notifier.fireFinishCGEvent();

	}

	/**
	 * Invokes the solve method of the Master Problem, fires corresponding events and queries the results.
	 * @param timeLimit Future point in time by which the Master Problem must be finished
	 * @throws TimeLimitExceededException TimeLimitExceededException
	 */
	protected void invokeMaster(long timeLimit) throws TimeLimitExceededException {
		notifier.fireStartMasterEvent();
		long time=System.currentTimeMillis();
		master.solve(timeLimit);
		objectiveMasterProblem =master.getObjective();
		masterSolveTime+=(System.currentTimeMillis()-time);

		//Check if we have found an integer solution
		masterSolutionIsInteger = true;
		for(Route route: master.getSolution())
			if(route.value>0+config.PRECISION && route.value<1-config.PRECISION) {masterSolutionIsInteger = false; break;}

		//Update incumbent solution
		if(masterSolutionIsInteger && this.cutoffValue>objectiveMasterProblem) {
			this.incumbentSolution = new ArrayList<>();
			this.cutoffValue = (int) (master.getObjective()+0.5);
			this.incumbentSolutionObjective = this.cutoffValue;
			for(Route route: master.getSolution()) {
				Route newRoute = route.clone();
				newRoute.value = route.value;
				this.incumbentSolution.add(newRoute);
			}
		}
		notifier.fireFinishMasterEvent();
	}

	/**
	 * Invokes the solve methods of the algorithms which solve the Pricing Problem. In addition, after solving the Pricing Problems
	 * and before any new columns are added to the Master Problem, this method invokes the {@link #calculateBoundOnMasterObjective(Class solver) calculateBoundOnMasterObjective} method.
	 * @param timeLimit Future point in time by which the Pricing Problem must be finished
	 * @return list of new columns which have to be added to the Master Problem, or an empty list if no columns could be identified
	 * @throws TimeLimitExceededException TimeLimitExceededException
	 */
	@Override
	protected List<Route> invokePricingProblems(long timeLimit) throws TimeLimitExceededException {
		
		//Solve the pricing problem
		List<Route> newColumns=new ArrayList<Route>();
		long time=System.currentTimeMillis();

		//Update data in pricing problems
		for(PricingProblem pricingProblem : pricingProblems){
			master.initializePricingProblem(pricingProblem);
		}

		//Solve pricing problems in the order of the pricing algorithms
		notifier.fireStartPricingEvent();
		cPricingProblemManager.setTimeLimit(timeLimit);
		((PricingProblem) pricingProblems.get(0)).compute_charging_bounds();
		boolean exact = false;
		for(Class<? extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem>> solver : solvers){
			
			if (needsChargingBranchingPricing == solverCapabilities.get(solver)) {
				newColumns = cPricingProblemManager.solvePricingProblems(solver);
				if (dataModel.rollbackTrigger) break;
			}

			//Stop when we found new columns
			if(!newColumns.isEmpty()){
				break;
			}
			exact = true;
		}
		
		notifier.fireFinishPricingEvent(newColumns);

		pricingSolveTime+=(System.currentTimeMillis()-time);
		nrGeneratedColumns+=newColumns.size();
		
		// Update of Lower Bound
		if(exact) 
			this.hasExceededPricingSoftThreshold = this.hasExceededPricingSoftThreshold || (dataModel.rollbackExplosion >= dataModel.pricingSoftFactor*dataModel.rollbackBaseLine);
			
			if (!newColumns.isEmpty()) this.boundOnMasterObjective = (optimizationSenseMaster == OptimizationSense.MINIMIZE ? Math.max(boundOnMasterObjective,this.calculateBoundOnMasterObjective(solvers.get(1))) : Math.min(boundOnMasterObjective,this.calculateBoundOnMasterObjective(solvers.get(1))));
			else { // The RMP bound is optimal
				
				// Look for an integer solution if
				// i) the current MP solution is NOT integer, and
				// ii) the current gap is greater than 5%
				if (!masterSolutionIsInteger && (1-this.boundOnMasterObjective/this.cutoffValue) > 0.05){
					
					Master mMaster = (Master) master;
					VRPMasterData masterData = mMaster.getMasterData();
					
					double IPtime = System.currentTimeMillis();
					extendedNotifier.fireIPSolutionEvent();
					try { this.solveIP(master.getColumns(pricingProblems.get(0)), masterData.subsetRowInequalities.keySet(), masterData.branchingNumberOfVehicles.keySet()); } 
					catch (IloException e) { e.printStackTrace(); logger.debug(e.getMessage()); }
					extendedNotifier.fireFinishIPSolutionEvent(System.currentTimeMillis() - IPtime);
				
				}

				if (!masterSolutionIsInteger &&  (1-this.boundOnMasterObjective/this.cutoffValue) <= 0.05){
					perform_fixing_by_reduced_cost(timeLimit);}

				// If the IP found a better integer solution, the gap reduction is computed using the newly updated Upper Bound
				this.gapReduction = (master.getObjective()-this.boundOnMasterObjective)/(this.cutoffValue-this.boundOnMasterObjective);
				this.boundOnMasterObjective = master.getObjective(); // Update the Bound before adding cuts
			}
			

			if (this.BBnodeID == 0 && dataModel.cut_iterations == 1){
				this.contExact ++;
				dataModel.rollbackBaseLine = (dataModel.rollbackBaseLine*(contExact-1)+dataModel.rollbackExplosion)/contExact;
			}

		// Add columns to the master problem
		if(!newColumns.isEmpty()){
			for(Route column : newColumns){
				master.addColumn(column);
			}
		}
		return newColumns;
	}

	private void perform_rollback(OptimalSolutionMemory memory){

		this.extendedNotifier.fireRollbackEvent(dataModel.rollbackBaseLine, dataModel.rollbackExplosion);
		this.objectiveMasterProblem = memory.previousMPObjective;
		this.boundOnMasterObjective = memory.previousMPBound;

		Master mMaster = (Master) master;
		mMaster.rollbackReconstruction(memory.previousColumns, memory.previousCuts, memory.previousMPSolution, memory.previousMPObjective);

		VRPMasterData mData = mMaster.getMasterData();
		this.extendedNotifier.fireFinishRollbackEvent(mMaster.getSolution(), mData.objectiveValue, mData.getNrColumns(), mData.subsetRowInequalities.size(), mData.branchingNumberOfVehicles.size());
	
	}

	protected void perform_fixing_by_reduced_cost(long timeLimit){

		/////////////////////////// PERFORM VARIABLE FIXING BY REDUCED COST //////////////////////////////

		extendedNotifier.fireFixingByReducedCostEvent(this.cutoffValue, this.boundOnMasterObjective);
		PricingProblem pricingProblem = (PricingProblem)pricingProblems.get(0);
		Map<Integer, Double> arcsToRemove = pricingProblem.fixByReducedCosts(timeLimit, this.cutoffValue, this.boundOnMasterObjective);
		
		// Deleting columns containing the eliminated arcs
		Master mMaster = (Master) master;
		Set<Route> columns = mMaster.getColumns(pricingProblem); List<Route> filtered_columns = new ArrayList<>();
		for (Route col: columns) if (!col.PParcs.stream().anyMatch(arcsToRemove.keySet()::contains)) filtered_columns.add(col);
		
		// Creating fake branching decisions
		List<BranchingDecision> removals = new ArrayList<BranchingDecision>();
		List<AbstractInequality> cuts = new ArrayList<>(mMaster.getMasterData().subsetRowInequalities.keySet());
		for (int arcID: arcsToRemove.keySet()){ removals.add(new RemoveArc(pricingProblem, arcID, dataModel.PParcs[arcID].arc_type, dataModel, cuts, 0));}
		this.branchingFRC.addAll(removals);

		// Updating the Master and Pricing Problem using the fake branches
		for (BranchingDecision bd: removals) { master.branchingDecisionPerformed(bd); pricingProblem.branchingDecisionPerformed(bd); }
		master.addColumns(filtered_columns);

		// Updating the Pricing Problem Solvers using the fake branches
		for (Class<? extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem>> solver: solvers){
			PricingProblemBundle<EVRPTW, Route, PricingProblem> bundle =  this.pricingProblemBundles.get(solver);
			AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem> solverInstance = bundle.solverInstances.get(0);
			for (BranchingDecision bd: removals) solverInstance.branchingDecisionPerformed(bd);
		}

		extendedNotifier.fireFinishFixingByReducedCostEvent(arcsToRemove, pricingProblem.bestReducedCost);
		
	}

	public void solveIP(Set<Route> columns, Set<SubsetRowInequality> subsetRowInequalities, Set<NumberVehiclesInequalities> vehiclesInequalities) throws IloException {

		Map<Route, IloIntVar> solution = new HashMap<Route, IloIntVar>();
		IloCplex cplex = new IloCplex();
		cplex.setOut(null); 			//disable CPLEX output
		cplex.setParam(IloCplex.Param.RandomSeed, 30);
		cplex.setParam(IloCplex.Param.Threads, 1);
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1e-4);

		// Define the objective
		IloObjective obj= cplex.addMinimize();

		// Routing set partitioning constraints
		IloRange[] visitCustomerConstraints=new IloRange[dataModel.C];
		for(int i=0; i< dataModel.C; i++)
			visitCustomerConstraints[i] = cplex.addEq(cplex.linearNumExpr(), 1, "visitCustomer_"+(i+1));

		// Charging capacity constraints
		IloRange[]  chargersCapacityConstraints = new IloRange[dataModel.last_charging_period];
		for (int t = 0; t < dataModel.last_charging_period; t++)
			chargersCapacityConstraints[t] = cplex.addLe(cplex.linearIntExpr(), dataModel.B, "capacity_"+(t+1));

		// Subset Row Cuts
		IloRange[] SRCs = new IloRange[subsetRowInequalities.size()]; int ix = 0;
		for (SubsetRowInequality subsetRowInequality: subsetRowInequalities){
			SRCs[ix] = cplex.addLe(cplex.linearNumExpr(), 1, "src_"+Arrays.toString(subsetRowInequality.cutSet));
			ix ++; }
		
		// Number of vehicles branches
		IloRange[] NumVehiclesBranches = new IloRange[vehiclesInequalities.size()]; ix = 0;
		for (NumberVehiclesInequalities branch: vehiclesInequalities){
			if (branch.lessThanOrEqual) NumVehiclesBranches[ix] = cplex.addLe(cplex.linearNumExpr(), branch.coefficient, "branching_"+branch.toString());
			else NumVehiclesBranches[ix] = cplex.addGe(cplex.linearNumExpr(), branch.coefficient, "branching_"+branch.toString());
			ix ++; }

		for (Route route: columns) {

			if (route.isArtificialColumn) continue;
			Route column = route.clone();
			IloColumn iloColumn = cplex.column(obj,column.cost);

			for(int i: route.route.keySet())
				iloColumn = iloColumn.and(cplex.column(visitCustomerConstraints[i-1], column.route.get(i)));

			for (int t = column.lastChargingTime; t >= (column.lastChargingTime-column.chargingTime+1); t--)
				iloColumn = iloColumn.and(cplex.column(chargersCapacityConstraints[t-1], 1));

			ix = 0;
			for (SubsetRowInequality subsetRowInequality: subsetRowInequalities) {
				iloColumn = iloColumn.and(cplex.column(SRCs[ix], getSRCCoefficient(column, subsetRowInequality)));
				ix ++; }

			ix = 0;
			for (NumberVehiclesInequalities branch: vehiclesInequalities) {
				iloColumn = iloColumn.and(cplex.column(NumVehiclesBranches[ix],1));
				ix ++; }

			//Create the variable and store it
			IloIntVar var= cplex.intVar(iloColumn, 0, Integer.MAX_VALUE);
			cplex.add(var);
			solution.put(column, var);
		}

		//Set time limit
		cplex.setParam(IloCplex.Param.TimeLimit, 60.0); //set time limit in seconds (in this case 10 seconds)
		if(cplex.solve() && cplex.getStatus()==IloCplex.Status.Optimal){
			int objVal = (int) (cplex.getObjValue()+0.05);

			if (objVal < this.cutoffValue){ // Found a better solution, update current incumbent
				this.cutoffValue = objVal;
				this.incumbentSolutionObjective = objVal;
				
				ArrayList<Route> optimalSolution = new ArrayList<Route>();
				if (dataModel.print_log) logger.debug("Found integer solution. Objective: "+this.incumbentSolutionObjective);
				for (Route route: solution.keySet()) {
					double value = cplex.getValue(solution.get(route));
					if(value > 0.5){
						Route newRoute = route.clone();
						newRoute.value = 1;
						optimalSolution.add(newRoute);

						if (dataModel.print_log) logger.debug(newRoute.toString());
					}
				}

				this.incumbentSolution = optimalSolution;
			}
		} else {
			if (dataModel.print_log) logger.debug("Did not find an integer solution");
		}
		cplex.close();
		cplex.end();
	}

	public int getSRCCoefficient(Route route, SubsetRowInequality subsetRowInequality) {
		int visits = 0;
		for(int i: subsetRowInequality.cutSet) visits+=route.route.getOrDefault(i, 0);
		return (int) Math.floor(0.5*visits);
	}

	public void update_solution_memory(){

		// Saves the current CG state in case it is needed in the future for a rollback
		List<Route> memoryColumns = new ArrayList<>();
		for (Route column: master.getColumns(pricingProblems.get(0))){
			Route newCol = column.clone(); newCol.BBnode = column.BBnode;
			memoryColumns.add(newCol);
		}

		List<SubsetRowInequality> memorySRCs = new ArrayList<>();
		for (SubsetRowInequality src: ((Master)master).getMasterData().subsetRowInequalities.keySet()) memorySRCs.add(src);

		List<Route> memoryIncumbent = new ArrayList<>();
		for (Route column: this.incumbentSolution){
			Route newCol = column.clone(); newCol.value = column.value;
			memoryIncumbent.add(newCol);
		}
		
		List<Route> memorySolution = master.getSolution();

		this.solutionMemory = new OptimalSolutionMemory(memoryColumns, memorySRCs, this.boundOnMasterObjective, memorySolution, this.objectiveMasterProblem, memoryIncumbent, this.incumbentSolutionObjective);

	}

	@Override
	public void close() {
    	this.master.close();
		this.cPricingProblemManager.close();
   	}

	public class OptimalSolutionMemory{

		public List<Route> previousColumns;
		public List<SubsetRowInequality> previousCuts;
		public double previousMPBound;
		public List<Route> previousMPSolution;
		public double previousMPObjective;
		public List<Route> previousIncumbentSolution;
		public int previousIncumbentObjective;

		public OptimalSolutionMemory(List<Route> previousColumns, List<SubsetRowInequality> previousCuts, double previousNodeBound, List<Route> previousMPSolution,
			double previousMPObjective, List<Route> previousIncumbentSolution, int previousIncumbentObjective){

			this.previousColumns = previousColumns;
			this.previousCuts = previousCuts;
			this.previousMPBound = previousNodeBound;
			this.previousMPSolution = previousMPSolution;
			this.previousMPObjective = previousMPObjective;
			this.previousIncumbentSolution = previousIncumbentSolution;
			this.previousIncumbentObjective = previousIncumbentObjective;

		}
	}
}
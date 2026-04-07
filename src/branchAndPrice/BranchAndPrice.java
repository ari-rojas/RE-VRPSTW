package branchAndPrice;

import java.util.Collection;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.Iterator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchAndPrice;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.CGListener;
import org.jorlib.frameworks.columnGeneration.io.TimeLimitExceededException;
import org.jorlib.frameworks.columnGeneration.master.OptimizationSense;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;
import org.jorlib.frameworks.columnGeneration.util.MathProgrammingUtil;
import org.jorlib.frameworks.columnGeneration.util.OrderedBiMap;

import columnGeneration.Master;
import columnGeneration.PricingProblem;
import columnGeneration.Route;
import columnGeneration.VRPMasterData;
import columnGeneration.customCG;
import ilog.concert.IloColumn;
import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import model.EVRPTW;


/**
 * Branch-and-Price class
 */
public final class BranchAndPrice extends AbstractBranchAndPrice<EVRPTW,Route,PricingProblem> {

	PricingProblem pricingProblem; 					//pricing problem
	public static final double PRECISION=0.001; 	//precision considered for the fractional solutions (nodes)
	private final ExtendBAPNotifier extendedNotifier;

	private List<Integer> chargingNodes = new ArrayList<Integer>();
	private List<Integer> arcFlowNodes = new ArrayList<Integer>();
	private long timeChargingBranching = 0;

	public BranchAndPrice(EVRPTW modelData, Master master, PricingProblem pricingProblem,
			List<Class<? extends AbstractPricingProblemSolver<EVRPTW,Route,PricingProblem>>> solvers,
			List<? extends AbstractBranchCreator<EVRPTW,Route,PricingProblem>> branchCreators,
			int objectiveInitialSolution,
			List<Route> initialSolution){
		
		super(modelData, master, pricingProblem, solvers, branchCreators, 0, objectiveInitialSolution);
		this.warmStart(objectiveInitialSolution, initialSolution);
		this.pricingProblem = pricingProblem;

		this.extendedNotifier = new ExtendBAPNotifier(this);

		this.setNodeOrdering(new Comparator<BAPNode>() {

			@Override
			public int compare(BAPNode node1, BAPNode node2) {
				if(node1.getBound()<=node2.getBound()) return -1;
				else return 1;
			}
			
		}); //Best Node First (BNF)
		//		this.setNodeOrdering(new BFSbapNodeComparator()); //Breadth-First Search (BFS)
		//		this.setNodeOrdering(new DFSbapNodeComparator()); //Depth-First Search (DFS)
	}

	/**
	 * Generates an artificial solution. Columns in the artificial solution are of high cost such that they never end up in the final solution
	 * if a feasible solution exists, since any feasible solution is assumed to be cheaper than the artificial solution. The artificial solution is used
	 * to guarantee that the master problem has a feasible solution.
	 * @return artificial solution
	 */
	@Override
	protected List<Route> generateInitialFeasibleSolution(BAPNode<EVRPTW,Route> node) {	
		//Dummy (artificial) routes to identify infeasibility
		HashMap<Integer, Integer> route=new HashMap<Integer, Integer>(dataModel.C);
		int[] routeSequence = new int[dataModel.C];
		for(int i=0; i< dataModel.C; i++) {route.put(i+1, 1); routeSequence[i] = i+1;}
		return Collections.singletonList(new Route("initSolution", true, route, routeSequence, pricingProblem, (int) Math.pow(10, 20), 0, 0, 0, 0.0, new ArrayList<Integer>(), 0, 0)); //dummy 
	}

	/**
	 * Checks whether the given node is integer
	 * @param node Node in the Branch-and-Price tree
	 * @return true if the solution is an integer solution
	 */
	@Override
	protected boolean isIntegerNode(BAPNode<EVRPTW, Route> node) {

		if(node.nodeID == 0) { //stores the information for the root node
			dataModel.columnsRootNode=master.getColumns(this.pricingProblem).size();
			dataModel.cutsRootNode=node.getInequalities().size();
		}

		return isIntegerSolution(node.getSolution());
	}

	protected boolean isIntegerSolution(List<Route> solution){

		for(Route route: solution)
			if(route.value>0+PRECISION && route.value<1-PRECISION) {return false;}

		return true;
	}

	protected void updateNodeGeneratedColumns(BAPNode<EVRPTW, Route> bapNode){

		//Inherit the routes generated
		List<Route> routesToAdd = new ArrayList<Route>();
		for(Route column: master.getColumns(this.pricingProblem)) {
			if(column.BBnode==-1) {
				column.BBnode=bapNode.nodeID;
				routesToAdd.add(column);
			}
		}
		bapNode.addInitialColumns(routesToAdd);
		//Inherit the cuts generated (not necessary)

		//Solve MIP at root node (optional)
		if(bapNode.nodeID == 0) {
			try {solveIPAtRootNode(bapNode);} 
			catch (IloException e) {e.printStackTrace();}
		}

	}

	/**
	 * To have a stronger upper bound, we solve the MIP at the root (with the generated columns)
	 */
	public void solveIPAtRootNode(BAPNode<EVRPTW, Route> node) throws IloException {

		Map<Route, IloIntVar> solution = new HashMap<Route, IloIntVar>();
		IloCplex cplex =new IloCplex(); 									//create CPLEX instance
		cplex.setOut(null);													//disable CPLEX output
		cplex.setParam(IloCplex.IntParam.Threads, config.MAXTHREADS); 		//set number of threads that may be used by the cplex

		//Define the objective
		IloObjective obj= cplex.addMinimize();
		//Define partitioning constraints
		IloRange[] visitCustomerConstraints=new IloRange[dataModel.C];
		for(int i=0; i< dataModel.C; i++)
			visitCustomerConstraints[i] = cplex.addEq(cplex.linearNumExpr(), 1, "visitCustomer_"+(i+1));

		//define constrains (capacitated station)
		IloRange[]  chargersCapacityConstraints = new IloRange[dataModel.last_charging_period];
		for (int t = 0; t < dataModel.last_charging_period; t++)
			chargersCapacityConstraints[t] = cplex.addLe(cplex.linearIntExpr(), dataModel.B, "capacity_"+(t+1));

		for(Route route: node.getInitialColumns()) {

			Route column = route.clone();
			//Register column with objective
			IloColumn iloColumn= cplex.column(obj,column.cost);

			//Register column with partitioning constraint
			for(int i: route.route.keySet())
				iloColumn=iloColumn.and(cplex.column(visitCustomerConstraints[i-1], column.route.get(i)));

			//Register column with chargers capacity constraints
			for (int t = column.initialChargingTime; t <= (column.initialChargingTime+ column.chargingTime-1); t++)
				iloColumn=iloColumn.and(cplex.column(chargersCapacityConstraints[t-1], 1));


			//Create the variable and store it
			IloIntVar var= cplex.intVar(iloColumn, 0, 1);
			cplex.add(var);
			solution.put(column, var);
		}

		//Set time limit
		cplex.setParam(IloCplex.DoubleParam.TiLim, 10.0); //set time limit in seconds (in this case 10 seconds)
		if(cplex.solve() && cplex.getStatus()==IloCplex.Status.Optimal && cplex.getCplexTime()<10){
			objectiveIncumbentSolution = (int) (cplex.getObjValue()+0.05);
			upperBoundOnObjective = objectiveIncumbentSolution;
			//retrieve solution
			List<Route> optimalSolution = new ArrayList<Route>();
			for (Route route: solution.keySet()) {
				double value = cplex.getValue(solution.get(route));
				if(value>=config.PRECISION){
					Route newRoute = route.clone();
					newRoute.value = value;
					optimalSolution.add(newRoute);
				}
			}
			incumbentSolution = optimalSolution;
		}
		cplex.close();
		cplex.end();
	}

	/**
	 * Solve a given Branch-and-Price node
	 * @param bapNode node in Branch-and-Price tree
	 * @param timeLimit future point in time by which the method must be finished
	 * @throws TimeLimitExceededException TimeLimitExceededException
	 */
	@Override
	protected void solveBAPNode(BAPNode<EVRPTW,Route> bapNode, long timeLimit) throws TimeLimitExceededException {
		// I want to override it with my method
	}

	protected CGResult solveNode(BAPNode<EVRPTW,Route> bapNode, long timeLimit) throws TimeLimitExceededException {
		
		customCG cg=null;
		try {
			dataModel.cleanSRCs(); // MODIFICATION
			cg = new customCG(dataModel, master, pricingProblems, solvers, pricingProblemManager, bapNode.getInitialColumns(), objectiveIncumbentSolution, bapNode.getBound()); //Solve the node
			for(CGListener listener : columnGenerationEventListeners) cg.addCGEventListener(listener);
			cg.solve(timeLimit);
		} finally {
			//Update statistics
			if(cg != null) {
				timeSolvingMaster += cg.getMasterSolveTime();
				timeSolvingPricing += cg.getPricingSolveTime();
				totalNrIterations += cg.getNumberOfIterations();
				totalGeneratedColumns += cg.getNrGeneratedColumns();
				//				if(cg.incumbentSolutionObjective<=this.objectiveIncumbentSolution) {this.objectiveIncumbentSolution = cg.incumbentSolutionObjective; this.incumbentSolution=cg.incumbentSolution;}
				notifier.fireFinishCGEvent(bapNode, cg.getBound(), cg.getObjective(), cg.getNumberOfIterations(), cg.getMasterSolveTime(), cg.getPricingSolveTime(), cg.getNrGeneratedColumns());
			}
		}

		ArrayList<Route> solution = new ArrayList<Route>(cg.getSolution().size()); //if not, it overwrites the value
		for(Route route: cg.getSolution()) {Route newRoute = route.clone(); newRoute.value = route.value; solution.add(newRoute);}
		bapNode.storeSolution(cg.getObjective(), cg.getBound(), solution, cg.getCuts());

		return new CGResult(cg.incumbentSolution, cg.incumbentSolutionObjective);
	}

	protected void processIntegerNode(BAPNode<EVRPTW, Route> bapNode){

		int integerObjective = MathProgrammingUtil.doubleToInt(bapNode.getObjective());
		this.notifier.fireNodeIsIntegerEvent(bapNode, bapNode.getBound(), integerObjective);
		this.objectiveIncumbentSolution = integerObjective;
		this.upperBoundOnObjective = (double)integerObjective;
		this.incumbentSolution = bapNode.getSolution();
	}

	protected boolean findIntegerSolution(BAPNode<EVRPTW, Route> bapNode){
		
		long time=System.currentTimeMillis();
		this.extendedNotifier.fireLexicographicMasterEvent(bapNode);
		boolean integer_solution_exists = false;

		// Retrieve the unique routes from the fractional solution
		List<Route> solution = bapNode.getSolution();
		int[] charging_times = new int[solution.size()]; int[] departure_times = new int[solution.size()]; int n = 0;
		LinkedHashMap<ArrayList<Integer>, Route> unique_routes = new LinkedHashMap<ArrayList<Integer>, Route>();
		for (Route column: solution){
			if (unique_routes.containsKey(column.arcs)) continue;
			else {
				unique_routes.put(column.arcs, column);
				charging_times[n] = column.chargingTime;
				departure_times[n] = column.departureTime;
				n ++;
			}
		}
		logger.debug("There are "+n+" unique routes in the fractional solution.");

		// Solve the charging scheduling problem
		int maxT = dataModel.last_charging_period;
		try { integer_solution_exists = this.solveChargingScheduling(n, maxT, charging_times, departure_times, dataModel.B, unique_routes);}
		catch (IloException e) {e.printStackTrace();}

		// Retrieve and store the solution, if it exists
		if (integer_solution_exists){
			bapNode.storeSolution(bapNode.getObjective(), bapNode.getBound(), new ArrayList<>(unique_routes.values()), this.master.getCuts());
		}

		this.extendedNotifier.fireFinishLexicographicMasterEvent(bapNode, true);
		this.timeChargingBranching += System.currentTimeMillis() - time;

		return integer_solution_exists;
		
	}

	protected boolean solveChargingScheduling(int n, int maxT, int[] charging_times, int[] departure_times, double B, LinkedHashMap<ArrayList<Integer>,Route> unique_routes) throws IloException {

		boolean integer_solution = false;
		try {
			IloCplex cplex = new IloCplex();
			cplex.setOut(null); 			//disable CPLEX output
			cplex.setParam(IloCplex.Param.Threads, 1);

			// x[t] = x_(t,t+1), t = 0..maxT
			IloNumVar[] x = new IloNumVar[maxT + 1];
			for (int t = 0; t <= maxT; t++) {
				x[t] = cplex.numVar(0.0, Double.MAX_VALUE, "x_" + t + "_" + (t + 1));
			}

			// y[r][t] = y_t^r, r = 0..n-1, b_r <= t < d_r
			IloIntVar[][] y = new IloIntVar[n][maxT + 1];
			for (int r = 0; r < n; r++) {
				for (int t = charging_times[r]; t < departure_times[r]; t++) {
					y[r][t] = cplex.boolVar("y_" + t + "_" + r);
				}
			}

			// -----------------------------
			// Boundary constraints
			// x_(0,1) = 0
			// x_(T,T+1) = 0, with T = maxT
			// -----------------------------
			cplex.addEq(x[0], 0.0, "x_0_1_zero");
			cplex.addEq(x[maxT], 0.0, "x_" + maxT + "_" + (maxT + 1) + "_zero");

			// -----------------------------
			// 1) Charging capacity
			// x_(t,t+1) + sum_{r : b_r <= t < d_r} y_t^r <= B     for all t = 1..maxT
			// -----------------------------
			for (int t = 1; t <= maxT; t++) {
				IloLinearNumExpr lhs = cplex.linearNumExpr();

				lhs.addTerm(1.0, x[t]);
				for (int r = 0; r < n; r++) {
					if (charging_times[r] <= t && t < departure_times[r]) lhs.addTerm(1.0, y[r][t]);
				}

				cplex.addLe(lhs, B, "capacity_" + t);
			}

			// -----------------------------
			// 2) Flow conservation
			// x_(t,t+1) + sum_{r : b_r <= t < d_r} y_t^r
			//   = x_(t-1,t) + sum_{r : t+b_r-1 < d_r} y_{t+b_r-1}^r
			// for all t = 1..maxT
			// -----------------------------
			for (int t = 1; t <= maxT; t++) {
				IloLinearNumExpr lhs = cplex.linearNumExpr();
				IloLinearNumExpr rhs = cplex.linearNumExpr();

				// Left-hand side: x_(t,t+1) + sum_{r : b_r <= t < d_r} y_t^r
				lhs.addTerm(1.0, x[t]);
				for (int r = 0; r < n; r++) {
					if (charging_times[r] <= t && t < departure_times[r])  lhs.addTerm(1.0, y[r][t]);
				}

				// Right-hand side: x_(t-1,t) + sum_{r : t+b_r-1 < d_r} y_{t+b_r-1}^r
				rhs.addTerm(1.0, x[t - 1]);
				for (int r = 0; r < n; r++) {
					int tau = t + charging_times[r] - 1; // tau = t + b_r - 1
					if (tau < departure_times[r]) rhs.addTerm(1.0, y[r][tau]);
				}

				cplex.addEq(lhs, rhs, "flow_" + t);
			}

			// -----------------------------
			// 3) Each EV charges once
			// sum_{t : b_r <= t < d_r} y_t^r = 1     for all r = 0..n-1
			// -----------------------------
			for (int r = 0; r < n; r++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();

				for (int t = charging_times[r]; t < departure_times[r]; t++) expr.addTerm(1.0, y[r][t]);
				cplex.addEq(expr, 1.0, "charge_once_" + r);
			}

			// Solve the model, and if a feasible solution exists, retrieve it
			if (cplex.solve()) {
                integer_solution = true;

				int r = 0;
				for (Map.Entry<ArrayList<Integer>, Route> entry: unique_routes.entrySet()){

					Route new_column = entry.getValue().clone();
					new_column.value = 1;
					for (int t=departure_times[r]-1; t>=charging_times[r]; t--){
						if (cplex.getValue(y[r][t]) > 0.5) {new_column.initialChargingTime = t-charging_times[r]+1; break;}
					}

					unique_routes.put(entry.getKey(), new_column);
					r ++;
				}
				
			}
            cplex.end();

        } catch (IloException e) {
            e.printStackTrace();
        }

        return integer_solution;
    }

	/**
	 * Run the BAP algorithm
	 * @param timeLimit time limit for the algorithm
	 */
	@Override
	public void runBranchAndPrice(long timeLimit) {
		this.notifier.fireStartBAPEvent();
		this.runtime = System.currentTimeMillis();
		BAPNode<EVRPTW, Route> rootNode = (BAPNode<EVRPTW, Route>)this.queue.peek();
		if (rootNode.getInitialColumns().isEmpty()) {
		   rootNode.addInitialColumns(this.generateInitialFeasibleSolution(rootNode));
		}

		CGResult cgIncumbent;
  
		while(!this.queue.isEmpty()) {
			BAPNode<EVRPTW, Route> bapNode = (BAPNode<EVRPTW, Route>)this.queue.poll();
			this.notifier.fireNextNodeEvent(bapNode);
			if (this.nodeCanBePruned(bapNode)) { // If can be pruned by bound BEFORE solving it
				this.notifier.firePruneNodeEvent(bapNode, bapNode.getBound());
				++this.nodesProcessed;
			} else {
				this.graphManipulator.next(bapNode);
				if (bapNode.nodeID != 0) { bapNode.addInitialColumns(this.generateInitialFeasibleSolution(bapNode)); }
				
				long time = 0;
				try { // Try solving the node
					if (this.chargingNodes.contains(bapNode.nodeID)) { time = System.currentTimeMillis(); }//logger.debug("TIME BRANCHING - Starting to process node "+bapNode.nodeID);} // TIME BRANCHING
					if (this.arcFlowNodes.contains(bapNode.nodeID)) { dataModel.CUTSENABLED = true; } else { dataModel.CUTSENABLED = false; }
					cgIncumbent = this.solveNode(bapNode, timeLimit);
					if (this.chargingNodes.contains(bapNode.nodeID)) { timeChargingBranching += (System.currentTimeMillis()-time); }//logger.debug("TIME BRANCHING - Finished processing node "+bapNode.nodeID);} // TIME BRANCHING
				} catch (TimeLimitExceededException var8) { // Catch runtime exceeded exception
					this.queue.add(bapNode);
					this.notifier.fireTimeOutEvent(bapNode);
					break;
				} catch (UnsupportedOperationException e) { // MODIFICATION - PROBLEMS WITH THE LB
					this.extendedNotifier.fireCGProblemsLBEvent(bapNode);
					++nodesProcessed;
					break;
				} catch (RuntimeException e) { // MODIFICATION - IF THE MASTER PROBLEM IS INFEASIBLE DUE TO BRANCHING (THE PRICING IS NOT EVEN INVOKED) THE NODE IS PRUNED
					this.extendedNotifier.fireCGMasterIsInfeasibleEvent(bapNode);
					logger.debug(e.getMessage());
					++nodesProcessed; continue;
				}
	
				if (this.nodeCanBePruned(bapNode)) { // If can be pruned by bound AFTER solving the node
					this.notifier.firePruneNodeEvent(bapNode, bapNode.getBound());
					++this.nodesProcessed;
				} else if (this.isInfeasibleNode(bapNode)) { // If can be pruned by infeasibility DUE TO artificial columns in the solution
					this.notifier.fireNodeIsInfeasibleEvent(bapNode);
					++this.nodesProcessed;
				} else { // If it is either integer solution and hence pruned by optimality OR fractional and should branch
					if (this.isIntegerNode(bapNode)) { // If is integer, update incumbent
						this.processIntegerNode(bapNode);
					} else {

						// Update of the global Primal Bound in case the local Primal Bound of the node is better
						if (this.isIntegerSolution(cgIncumbent.cgIncumbentSolution) && (int) cgIncumbent.cgIncumbentObjective < objectiveIncumbentSolution) {
							int integerObjective = MathProgrammingUtil.doubleToInt(cgIncumbent.cgIncumbentObjective);
							this.objectiveIncumbentSolution = integerObjective;
							this.upperBoundOnObjective = cgIncumbent.cgIncumbentObjective;
							this.incumbentSolution = cgIncumbent.cgIncumbentSolution;
						}

						this.updateNodeGeneratedColumns(bapNode);
						List<BAPNode<EVRPTW, Route>> newBranches = new ArrayList();
						
						// Initialize Branch Creator
						BranchingRules bc = (BranchingRules)this.branchCreators.iterator().next();
						
						if (this.chargingNodes.contains(bapNode.nodeID)) { time = System.currentTimeMillis(); }
						// Look for Number of Vehicles or Customers Arc Flow branching
						boolean foundBranches = false;
						foundBranches = bc.canPerformFirstBranching(bapNode.getSolution());
						if (this.chargingNodes.contains(bapNode.nodeID)) { timeChargingBranching += (System.currentTimeMillis()-time); }
						if (foundBranches){
							if (this.chargingNodes.contains(bapNode.nodeID)) { time = System.currentTimeMillis(); }
							this.notifier.fireNodeIsFractionalEvent(bapNode, bapNode.getBound(), bapNode.getObjective());
							newBranches.addAll(bc.getFirstBranches(bapNode));

							if (bc.branchOnCustomerArcs || this.arcFlowNodes.contains(bapNode.nodeID)){
								this.arcFlowNodes.add(newBranches.get(0).nodeID);
								this.arcFlowNodes.add(newBranches.get(1).nodeID);
							}

							if (this.chargingNodes.contains(bapNode.nodeID)) { 
								timeChargingBranching += (System.currentTimeMillis()-time);
								this.chargingNodes.add(newBranches.get(0).nodeID);
								this.chargingNodes.add(newBranches.get(1).nodeID);
							}
						} else {
							
							time = System.currentTimeMillis();

							foundBranches = this.findIntegerSolution(bapNode);
							//foundBranches = false;

							if (foundBranches){
								
								this.processIntegerNode(bapNode);

							} else {
							
								foundBranches = bc.canPerformBranching(bapNode.getSolution());
								if (foundBranches){
									this.notifier.fireNodeIsFractionalEvent(bapNode, bapNode.getBound(), bapNode.getObjective());
									newBranches.addAll(bc.getBranches(bapNode));
								}

								if (this.arcFlowNodes.contains(bapNode.nodeID)){
									this.arcFlowNodes.add(newBranches.get(0).nodeID);
									this.arcFlowNodes.add(newBranches.get(1).nodeID);
								}

								timeChargingBranching += (System.currentTimeMillis()-time);
								this.chargingNodes.add(newBranches.get(0).nodeID);
								this.chargingNodes.add(newBranches.get(1).nodeID);

							}

						}
	
						if (!foundBranches) {
							throw new RuntimeException("BAP encountered fractional solution, but none of the BranchCreators produced any new branches?");
						}
						
						if (!newBranches.isEmpty()){
							this.queue.addAll(newBranches);
							this.notifier.fireBranchEvent(bapNode, Collections.unmodifiableList(newBranches));
						}
					}
	
				++this.nodesProcessed;
			  	}
		   	}
		}
  
		if (this.queue.isEmpty()) { // If all the BAP tree was explored, the incumbent solution is optimal
			this.isOptimal = true;
			if (this.optimizationSenseMaster == OptimizationSense.MINIMIZE) {
				this.lowerBoundOnObjective = (double)this.objectiveIncumbentSolution;
			} else {
				this.upperBoundOnObjective = (double)this.objectiveIncumbentSolution;
			}
		} else { // Else, cannot declare optimality
			this.isOptimal = false;
			Iterator var9;
			BAPNode bapNode;
			if (this.optimizationSenseMaster == OptimizationSense.MINIMIZE) {
				this.lowerBoundOnObjective = ((BAPNode)this.queue.peek()).getBound();
	
				for(var9 = this.queue.iterator(); var9.hasNext(); this.lowerBoundOnObjective = Math.min(this.lowerBoundOnObjective, bapNode.getBound())) {
					bapNode = (BAPNode)var9.next();
				}
			} else {
				this.upperBoundOnObjective = ((BAPNode)this.queue.peek()).getBound();
	
				for(var9 = this.queue.iterator(); var9.hasNext(); this.upperBoundOnObjective = Math.max(this.upperBoundOnObjective, bapNode.getBound())) {
					bapNode = (BAPNode)var9.next();
				}
			}
		}
  
		this.notifier.fireStopBAPEvent();
		this.runtime = System.currentTimeMillis() - this.runtime;

		double realTime = this.timeChargingBranching*0.001;
		realTime = Math.floor(realTime*100)/100;
		logger.debug("TIME BRANCHING - Total time is: "+realTime);

	}

	/**
	 * Test whether the given node can be pruned based on this bounds
	 * @param node node
	 * @return true if the node can be pruned
	 */
	@Override
	protected boolean nodeCanBePruned(BAPNode<EVRPTW,Route> node){
		//		System.out.println(Math.ceil(node.getBound()-config.PRECISION) + " >= " + this.objectiveIncumbentSolution);
		return Math.ceil(node.getBound()) >= (this.objectiveIncumbentSolution-config.PRECISION);
	}

	public void addExtendCGEventListener(ExtendBAPListener listener) {
		this.extendedNotifier.addExtendBAPListener(listener);
	}

	public void removeExtendCGEventListener(ExtendBAPListener listener) {
		this.extendedNotifier.removeExtendBAPListener(listener);
	}

	public class CGResult {

		public List<Route> cgIncumbentSolution;
		public double cgIncumbentObjective;

		public CGResult(List<Route> incumbentSolution, double primalBound){
			this.cgIncumbentSolution = incumbentSolution;
			this.cgIncumbentObjective = primalBound;
		}
	}
}
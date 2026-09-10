package branchAndPrice;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.Iterator;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchAndPrice;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.BAPListener;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.BranchEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.CGListener;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.StartEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.FinishEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.FinishProcessingNodeEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.NodeIsFractionalEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.NodeIsInfeasibleEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.NodeIsIntegerEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.ProcessingNextNodeEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.PruneNodeEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.TimeLimitExceededEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.bapNodeComparators.DFSbapNodeComparator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecisionListener;
import org.jorlib.frameworks.columnGeneration.io.TimeLimitExceededException;
import org.jorlib.frameworks.columnGeneration.master.OptimizationSense;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractInequality;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;
import org.jorlib.frameworks.columnGeneration.pricing.DefaultPricingProblemSolverFactory;
import org.jorlib.frameworks.columnGeneration.pricing.PricingProblemBundle;
import org.jorlib.frameworks.columnGeneration.util.MathProgrammingUtil;
import org.jorlib.frameworks.columnGeneration.util.Configuration;

import columnGeneration.Master;
import columnGeneration.PricingProblem;
import columnGeneration.Route;
import columnGeneration.SubsetRowInequality;
import columnGeneration.customCG;
import columnGeneration.customPricingProblemManager;
import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import model.EVRPTW;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Branch-and-Price class
 */
public final class BranchAndPrice {

	///////////////////////////////////////////
	/// Fields from AbstractBranchAndPrice
	///////////////////////////////////////////

	protected final Logger logger;
	protected final BranchAndPrice.BAPNotifier notifier;
	protected final Set<CGListener> columnGenerationEventListeners;
	protected final Configuration config;
	protected final EVRPTW dataModel;
	protected Master master;
	protected List<PricingProblem> pricingProblems;
	protected List<Class<? extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem>>> solvers;
	protected final OptimizationSense optimizationSenseMaster;
	protected int objectiveIncumbentSolution;
	protected List<Route> incumbentSolution;
	protected boolean isOptimal;
	protected Queue<BAPNode<EVRPTW, Route>> queue;
	protected int nodeCounter;
	protected BAPNode<EVRPTW, Route> rootNode;
	protected double upperBoundOnObjective;
	protected double lowerBoundOnObjective;
	protected int nodesProcessed;
	protected long timeSolvingMaster;
	protected long timeSolvingPricing;
	protected long runtime;
	protected int totalGeneratedColumns;
	protected int totalNrIterations;

	///////////////////////////////////////////
	/// CUSTOM FIELDS
	///////////////////////////////////////////

	private PricingProblem pricingProblem; 										//pricing problem
	public static final double PRECISION=0.001; 						//precision considered for the fractional solutions (nodes)
	private final ExtendBAPNotifier extendedNotifier;

	// Tracking of different features of the BPC tree nodes
	private List<Integer> chargingNodes = new ArrayList<Integer>();
	private List<Integer> arcFlowNodes = new ArrayList<Integer>();
	private long timeChargingBranching = 0;
	private Map<Integer,Boolean> comesFromRollback = new HashMap<Integer,Boolean>();

	// Tracking of the BPC tree nodes Root Path and Branching Decisions
	public final Map<Integer, ArrayList<Integer>> rootPaths = new HashMap<>();
	public final Map<Integer, ArrayList<BranchingDecision<EVRPTW, Route>>> branchingDecisions = new HashMap<>();

	private final Map<Class<? extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem>>, PricingProblemBundle<EVRPTW, Route, PricingProblem>> pricingProblemBundles;
	private final customPricingProblemManager<EVRPTW, Route, PricingProblem> cPricingProblemManager;
	private final customGraphManipulator cGraphManipulator;
	protected final BranchingRules branchCreator;

	public BranchAndPrice(EVRPTW modelData, Master master, List<PricingProblem> pricingProblems,
			List<Class<? extends AbstractPricingProblemSolver<EVRPTW,Route,PricingProblem>>> solvers,
			BranchingRules branchCreator, int objectiveInitialSolution, List<Route> initialSolution){
		
		this.logger = LoggerFactory.getLogger(AbstractBranchAndPrice.class);
		this.config = Configuration.getConfiguration();
		this.incumbentSolution = new ArrayList<Route>();
		this.isOptimal = false;
		this.nodeCounter = 0;
		this.upperBoundOnObjective = Double.MAX_VALUE;
		this.lowerBoundOnObjective = -Double.MAX_VALUE;
		this.nodesProcessed = 0;
		this.timeSolvingMaster = 0L;
		this.timeSolvingPricing = 0L;
		this.runtime = 0L;
		this.totalGeneratedColumns = 0;
		this.totalNrIterations = 0;
		this.dataModel = modelData;
		this.master = master;
		this.optimizationSenseMaster = master.getOptimizationSense();
		this.branchCreator = branchCreator;
		this.pricingProblems = pricingProblems;
		this.solvers = solvers;
		this.queue = new PriorityQueue<BAPNode<EVRPTW, Route>>(new DFSbapNodeComparator());
		this.objectiveIncumbentSolution = this.optimizationSenseMaster == OptimizationSense.MINIMIZE ? Integer.MAX_VALUE : -2147483647;
		this.lowerBoundOnObjective = 0;
		this.upperBoundOnObjective = objectiveInitialSolution;
		List<Integer> rootPath = new ArrayList<Integer>();
		int nodeID = this.nodeCounter++;
		rootPath.add(nodeID);
		if (this.optimizationSenseMaster == OptimizationSense.MINIMIZE) {
			this.rootNode = new BAPNode<EVRPTW, Route>(nodeID, rootPath, new ArrayList<Route>(), new ArrayList<AbstractInequality>(), lowerBoundOnObjective, Collections.emptyList());
		} else {
			this.rootNode = new BAPNode<EVRPTW, Route>(nodeID, rootPath, new ArrayList<Route>(), new ArrayList<AbstractInequality>(), upperBoundOnObjective, Collections.emptyList());
		}

		this.queue.add(this.rootNode);

		this.warmStart(objectiveInitialSolution, initialSolution);
		
		this.pricingProblem = pricingProblems.get(0);
		this.notifier = new BAPNotifier();
		this.extendedNotifier = new ExtendBAPNotifier(this);
      	this.columnGenerationEventListeners = new LinkedHashSet<CGListener>();
		
		// INITIALIZE TRACKERS WITH THE ROOT NODE INFORMATION

		ArrayList<Integer> rootP = new ArrayList<>(); rootP.add(0);
		this.rootPaths.put(0, rootP);
		this.branchingDecisions.put(0, new ArrayList<BranchingDecision<EVRPTW, Route>>());
		
		// CREATE CUSTOM GRAPH MANIPULATOR AND PRICING PROBLEM MANAGER

		this.cGraphManipulator = new customGraphManipulator(rootNode, rootPaths, branchingDecisions);
		this.pricingProblemBundles = new HashMap<Class<? extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem>>, PricingProblemBundle<EVRPTW, Route, PricingProblem>>();
		for(Class<? extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem>> solverClass : solvers) {
			DefaultPricingProblemSolverFactory<EVRPTW, Route, PricingProblem> factory = new DefaultPricingProblemSolverFactory<EVRPTW, Route, PricingProblem>(solverClass, dataModel);
			PricingProblemBundle<EVRPTW, Route, PricingProblem> bundle = new PricingProblemBundle<EVRPTW, Route, PricingProblem>(solverClass, pricingProblems, factory);
			pricingProblemBundles.put(solverClass, bundle);
		}

		this.cPricingProblemManager = new customPricingProblemManager<EVRPTW, Route, PricingProblem>(pricingProblems, pricingProblemBundles);
		
		// ADD THE BRANCHING DECISION LISTENERS
		
		this.addCBranchingDecisionListener(master);
		for(PricingProblem pProblem : pricingProblems) this.addCBranchingDecisionListener(pProblem);
		for(PricingProblemBundle<EVRPTW, Route, PricingProblem> bunddle : pricingProblemBundles.values()) {
			for(AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem> solverInstance : bunddle.solverInstances)  this.addCBranchingDecisionListener(solverInstance); }
		
		this.branchCreator.registerBAP(this);

		// NODE PRIORITY RULE
		this.setNodeOrdering(new Comparator<BAPNode<EVRPTW, Route>>() {

			@Override
			public int compare(BAPNode node1, BAPNode node2) {
				if(node1.getBound()<=node2.getBound()) return -1;
				else return 1;
			}
			
		});
	}

	public void warmStart(int objectiveInitialSolution, List<Route> initialSolution) {
		this.rootNode = (BAPNode<EVRPTW, Route>)this.queue.peek();
		if (this.rootNode.nodeID != 0) {
			throw new RuntimeException("This method can only be invoked at the start of the Branch-and-Price procedure, before runBranchAndPrice is invoked");
		} else {
			this.rootNode.addInitialColumns(initialSolution);
			this.objectiveIncumbentSolution = objectiveInitialSolution;
			this.incumbentSolution = new ArrayList<Route>(initialSolution);
			if (this.optimizationSenseMaster == OptimizationSense.MINIMIZE)  this.upperBoundOnObjective = (double)objectiveInitialSolution;
			else this.lowerBoundOnObjective = (double)objectiveInitialSolution;

		}
	}

	/**
	 * Generates an artificial solution. Columns in the artificial solution are of high cost such that they never end up in the final solution
	 * if a feasible solution exists, since any feasible solution is assumed to be cheaper than the artificial solution. The artificial solution is used
	 * to guarantee that the master problem has a feasible solution.
	 * @return artificial solution
	 */
	protected List<Route> generateInitialFeasibleSolution(BAPNode<EVRPTW,Route> node) {	
		//Dummy (artificial) routes to identify infeasibility
		HashMap<Integer, Integer> route=new HashMap<Integer, Integer>(dataModel.C);
		int[] routeSequence = new int[dataModel.C];
		for(int i=0; i< dataModel.C; i++) {route.put(i+1, 1); routeSequence[i] = i+1;}
		return Collections.singletonList(new Route("initSolution", true, route, routeSequence, pricingProblem, (int) Math.pow(10, 20), 0, 0, 0, 0.0, new ArrayList<Integer>(), new ArrayList<Integer>(), 0, 0)); //dummy 
	}

	protected BAPNode<EVRPTW, Route> updateNodeGeneratedColumns(BAPNode<EVRPTW, Route> bapNode, List<BranchingDecision<EVRPTW, Route>> removals){

		// Inherit the routes generated
		List<Route> columnsToAdd = new ArrayList<Route>();
		for(Route column: master.getColumns(this.pricingProblem)) {
			if(column.BBnode==-1)  column.BBnode=bapNode.nodeID;
			columnsToAdd.add(column);
		}

		List<Route> solution = new ArrayList<>();
		for (Route col: bapNode.getSolution()){
			Route newCol = col.clone(); newCol.value = col.value; newCol.BBnode = col.BBnode;
			solution.add(newCol);
		}

		List<AbstractInequality> inequalities = new ArrayList<>(bapNode.getInequalities());

		ArrayList<BranchingDecision<EVRPTW, Route>> brDecisions = branchingDecisions.get(bapNode.nodeID);
		brDecisions.addAll(removals);

		double bound = bapNode.getBound();
		
		// Destroy and re-create the rootNode
		bapNode = new BAPNode<EVRPTW, Route>(bapNode.nodeID, rootPaths.get(bapNode.nodeID), columnsToAdd, bapNode.getInitialInequalities(), bound, new ArrayList<BranchingDecision>(brDecisions));
		bapNode.storeSolution(bound, bound, solution, inequalities);

		if (!this.arcFlowNodes.contains(bapNode.nodeID) && !removals.isEmpty()) this.arcFlowNodes.add(bapNode.nodeID);

		return bapNode;

	}

	protected CGResult solveNode(BAPNode<EVRPTW,Route> bapNode, long timeLimit) throws TimeLimitExceededException {
		
		customCG cg=null;
		try {
			dataModel.cleanSRCs(); // MODIFICATION
			cg = new customCG(dataModel, master, pricingProblems, solvers, null, bapNode.getInitialColumns(), objectiveIncumbentSolution, bapNode.getBound(), bapNode.nodeID, this.chargingNodes.contains(bapNode.nodeID), this.extendedNotifier, this.pricingProblemBundles, this.cPricingProblemManager); //Solve the node
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

		dataModel.infeasiblePPArcs = pricingProblem.infeasiblePPArcs;

		return new CGResult(cg.incumbentSolution, cg.incumbentSolutionObjective, cg.branchingFRC);
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
		this.extendedNotifier.fireMIPMasterEvent(bapNode);
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

		this.extendedNotifier.fireFinishMIPMasterEvent(bapNode, true);
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
						if (cplex.getValue(y[r][t]) > 0.5) {new_column.lastChargingTime = t; break;}
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

	private void process_branching(BAPNode<EVRPTW, Route> bapNode, List<BAPNode<EVRPTW, Route>> newBranches, long time){

		// Initialize Branch Creator
		BranchingRules bc = this.branchCreator;
		
		if (this.chargingNodes.contains(bapNode.nodeID)) { time = System.currentTimeMillis(); }
		// Look for Number of Vehicles or Customers Arc Flow branching
		boolean foundBranches = false;
		foundBranches = bc.canPerformRoutingBranching(bapNode.getSolution());
		if (this.chargingNodes.contains(bapNode.nodeID)) { timeChargingBranching += (System.currentTimeMillis()-time); }
		if (foundBranches){
			if (this.chargingNodes.contains(bapNode.nodeID)) { time = System.currentTimeMillis(); }
			this.notifier.fireNodeIsFractionalEvent(bapNode, bapNode.getBound(), bapNode.getObjective());
			newBranches.addAll(bc.getRoutingBranches(bapNode));

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

			//foundBranches = this.findIntegerSolution(bapNode);
			foundBranches = false;

			if (foundBranches){
				
				this.processIntegerNode(bapNode);

			} else {
			
				foundBranches = bc.canPerformChargingBranching(bapNode.getSolution());
				if (foundBranches){
					this.notifier.fireNodeIsFractionalEvent(bapNode, bapNode.getBound(), bapNode.getObjective());
					newBranches.addAll(bc.getChargingBranches(bapNode));
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

		if (!newBranches.isEmpty()) {
			this.comesFromRollback.put(newBranches.get(0).nodeID, false);
			this.comesFromRollback.put(newBranches.get(1).nodeID,false);
		}

		if (!foundBranches) {
			throw new RuntimeException("BAP encountered fractional solution, but none of the BranchCreators produced any new branches?");
		}

	}

	/**
	 * Run the BAP algorithm
	 * @param timeLimit time limit for the algorithm
	 */
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
				this.cGraphManipulator.next(bapNode);
				if (bapNode.nodeID != 0) { bapNode.addInitialColumns(this.generateInitialFeasibleSolution(bapNode)); }
				
				long time = 0;
				try { // Try solving the node
					if (this.chargingNodes.contains(bapNode.nodeID)) { time = System.currentTimeMillis(); }//logger.debug("TIME BRANCHING - Starting to process node "+bapNode.nodeID);} // TIME BRANCHING
					if (this.arcFlowNodes.contains(bapNode.nodeID)) { dataModel.CUTSENABLED = true; } else { dataModel.CUTSENABLED = true; }
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
					bapNode = this.updateNodeGeneratedColumns(bapNode, cgIncumbent.branchingFRC);
					++this.nodesProcessed;
				} else if (this.isInfeasibleNode(bapNode)) { // If can be pruned by infeasibility DUE TO artificial columns in the solution
					this.notifier.fireNodeIsInfeasibleEvent(bapNode);
					++this.nodesProcessed;
				} else { // If it is either integer solution and hence pruned by optimality OR fractional and should branch
					if (this.isIntegerNode(bapNode)) { // If is integer, update incumbent
						this.processIntegerNode(bapNode);
						bapNode = this.updateNodeGeneratedColumns(bapNode, cgIncumbent.branchingFRC);
					} else {

						// Update of the global Primal Bound in case the local Primal Bound of the node is better
						if (this.isIntegerSolution(cgIncumbent.cgIncumbentSolution) && (int) cgIncumbent.cgIncumbentObjective < objectiveIncumbentSolution) {
							int integerObjective = MathProgrammingUtil.doubleToInt(cgIncumbent.cgIncumbentObjective);
							this.objectiveIncumbentSolution = integerObjective;
							this.upperBoundOnObjective = cgIncumbent.cgIncumbentObjective;
							this.incumbentSolution = cgIncumbent.cgIncumbentSolution;
						}

						bapNode = this.updateNodeGeneratedColumns(bapNode, cgIncumbent.branchingFRC);

						List<BAPNode<EVRPTW, Route>> newBranches = new ArrayList<>();
						this.process_branching(bapNode, newBranches, time);

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

	/////////////////////////////////////////////
	/// NODE SOLUTION CLASSIFIERS
	////////////////////////////////////////////

	/**
	 * Test whether the given node can be pruned based on this bounds
	 * @param node node
	 * @return true if the node can be pruned
	 */
	protected boolean nodeCanBePruned(BAPNode<EVRPTW,Route> node){
		//		System.out.println(Math.ceil(node.getBound()-config.PRECISION) + " >= " + this.objectiveIncumbentSolution);
		return Math.ceil(Math.floor(node.getBound()*10000)/10000) >= (this.objectiveIncumbentSolution-config.PRECISION);
	}

	protected boolean isInfeasibleNode(BAPNode<EVRPTW, Route> node) {
		for(Route column : node.getSolution()) {
			if (column.isArtificialColumn)  return true; }

		return false;
	}

	/**
	 * Checks whether the given node is integer
	 * @param node Node in the Branch-and-Price tree
	 * @return true if the solution is an integer solution
	 */
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

	public void setNodeOrdering(Comparator<BAPNode<EVRPTW, Route>> comparator) {
		Queue<BAPNode<EVRPTW, Route>> newQueue = new PriorityQueue<BAPNode<EVRPTW, Route>>(comparator);
		newQueue.addAll(this.queue);
		this.queue = newQueue;
	}

	public void close() {
		this.master.close();
		this.cPricingProblemManager.close();
   	}

	////////////////////////////////////////
	/// ADD AND REMOVE EVENT LISTENERS
	///////////////////////////////////////

	public void addCBranchingDecisionListener(BranchingDecisionListener listener) {
		this.cGraphManipulator.addBranchingDecisionListener(listener);
	}

	public void removeCBranchingDecisionListener(BranchingDecisionListener listener) {
		this.cGraphManipulator.removeBranchingDecisionListener(listener);
	}

	public void addBranchAndPriceEventListener(BAPListener listener) {
		this.notifier.addListener(listener);
	}

	public void removeBranchAndPriceEventListener(BAPListener listener) {
		this.notifier.removeListener(listener);
	}

	public void addExtendCGEventListener(ExtendBAPListener listener) {
		this.extendedNotifier.addExtendBAPListener(listener);
	}

	public void removeExtendCGEventListener(ExtendBAPListener listener) {
		this.extendedNotifier.removeExtendBAPListener(listener);
	}

	public void addColumnGenerationEventListener(CGListener listener) {
		this.columnGenerationEventListeners.add(listener);
	}

	public void removeColumnGenerationEventListener(CGListener listener) {
		this.columnGenerationEventListeners.add(listener);
	}

	public Integer getNewChildNodeID() {
		return this.nodeCounter++;
	}

	/**
	 * Computes the coefficient of a route in SRC.
	 * @param route for which the coefficient is calculated.
	 * @param subsetRowInequality considered.
	 */
	public int getSRCCoefficient(Route route, SubsetRowInequality subsetRowInequality) {
		int visits = 0;
		for(int i: subsetRowInequality.cutSet) visits+=route.route.getOrDefault(i, 0);
		return (int) Math.floor(0.5*visits);
	}

	public class CGResult {

		public List<Route> cgIncumbentSolution;
		public double cgIncumbentObjective;

		public List<BranchingDecision<EVRPTW, Route>> branchingFRC;

		public CGResult(List<Route> incumbentSolution, double primalBound, List<BranchingDecision<EVRPTW, Route>> branchingFRC){
			this.cgIncumbentSolution = incumbentSolution;
			this.cgIncumbentObjective = primalBound;
			this.branchingFRC = branchingFRC;
		}
	}

	////////////////////////////////////////////////
	/// GETTER METHODS
	////////////////////////////////////////////////
	
	protected int getUniqueNodeID() {
		return this.nodeCounter++;
	}

	public int getObjective() {
		return this.objectiveIncumbentSolution;
	}

	public double getBoundRootNode() {
		return this.rootNode.getBound();
	}

	public boolean hasSolution() {
		return !this.incumbentSolution.isEmpty();
	}

	public boolean isOptimal() {
		return this.isOptimal;
	}

	public double getBound() {
		return this.optimizationSenseMaster == OptimizationSense.MINIMIZE ? this.lowerBoundOnObjective : this.upperBoundOnObjective;
	}

	public int getNumberOfProcessedNodes() {
		return this.nodesProcessed;
	}

	public long getSolveTime() {
		return this.runtime;
	}

	public long getMasterSolveTime() {
		return this.timeSolvingMaster;
	}

	public long getPricingSolveTime() {
		return this.timeSolvingPricing;
	}

	public int getTotalGeneratedColumns() {
		return this.totalGeneratedColumns;
	}

	public int getTotalNrIterations() {
		return this.totalNrIterations;
	}

	public List<Route> getSolution() {
		return this.incumbentSolution;
	}

	/////////////////////////////////////
	/// BAP NOTIFIER
	////////////////////////////////////

	protected class BAPNotifier {
		private Set<BAPListener> listeners = new LinkedHashSet<BAPListener>();

		public BAPNotifier() {
		}

		public void addListener(BAPListener listener) {
			this.listeners.add(listener);
		}

		public void removeListener(BAPListener listener) {
			this.listeners.remove(listener);
		}

		public void fireStartBAPEvent() {
			StartEvent startEvent = null;

			for(BAPListener listener : this.listeners) {
				if (startEvent == null) {
				startEvent = new StartEvent(BranchAndPrice.this, BranchAndPrice.this.dataModel.getName(), BranchAndPrice.this.objectiveIncumbentSolution);
				}

				listener.startBAP(startEvent);
			}

		}

		public void fireStopBAPEvent() {
			FinishEvent finishEvent = null;

			for(BAPListener listener : this.listeners) {
				if (finishEvent == null) {
				finishEvent = new FinishEvent(BranchAndPrice.this);
				}

				listener.finishBAP(finishEvent);
			}

		}

		public void fireNodeIsFractionalEvent(BAPNode<EVRPTW, Route> node, double nodeBound, double nodeValue) {
			NodeIsFractionalEvent nodeIsFractionalEvent = null;

			for(BAPListener listener : this.listeners) {
				if (nodeIsFractionalEvent == null) {
				nodeIsFractionalEvent = new NodeIsFractionalEvent(BranchAndPrice.this, node, nodeBound, nodeValue);
				}

				listener.nodeIsFractional(nodeIsFractionalEvent);
			}

		}

		public void fireNodeIsIntegerEvent(BAPNode<EVRPTW, Route> node, double nodeBound, int nodeValue) {
			NodeIsIntegerEvent nodeIsIntegerEvent = null;

			for(BAPListener listener : this.listeners) {
				if (nodeIsIntegerEvent == null) {
				nodeIsIntegerEvent = new NodeIsIntegerEvent(BranchAndPrice.this, node, nodeBound, nodeValue);
				}

				listener.nodeIsInteger(nodeIsIntegerEvent);
			}

		}

		public void fireNodeIsInfeasibleEvent(BAPNode<EVRPTW, Route> node) {
			NodeIsInfeasibleEvent nodeIsInfeasibleEvent = null;

			for(BAPListener listener : this.listeners) {
				if (nodeIsInfeasibleEvent == null) {
				nodeIsInfeasibleEvent = new NodeIsInfeasibleEvent(BranchAndPrice.this, node);
				}

				listener.nodeIsInfeasible(nodeIsInfeasibleEvent);
			}

		}

		public void firePruneNodeEvent(BAPNode<EVRPTW, Route> node, double nodeBound) {
			PruneNodeEvent pruneNodeEvent = null;

			for(BAPListener listener : this.listeners) {
				if (pruneNodeEvent == null) {
				pruneNodeEvent = new PruneNodeEvent(BranchAndPrice.this, node, nodeBound, BranchAndPrice.this.objectiveIncumbentSolution);
				}

				listener.pruneNode(pruneNodeEvent);
			}

		}

		public void fireNextNodeEvent(BAPNode<EVRPTW, Route> node) {
			ProcessingNextNodeEvent processingNextNodeEvent = null;

			for(BAPListener listener : this.listeners) {
				if (processingNextNodeEvent == null) {
				processingNextNodeEvent = new ProcessingNextNodeEvent(BranchAndPrice.this, node, BranchAndPrice.this.queue.size(), BranchAndPrice.this.objectiveIncumbentSolution);
				}

				listener.processNextNode(processingNextNodeEvent);
			}

		}

		public void fireFinishCGEvent(BAPNode<EVRPTW, Route> node, double nodeBound, double nodeValue, int numberOfCGIterations, long masterSolveTime, long pricingSolveTime, int nrGeneratedColumns) {
			FinishProcessingNodeEvent finishProcessingNodeEvent = null;

			for(BAPListener listener : this.listeners) {
				if (finishProcessingNodeEvent == null) {
				finishProcessingNodeEvent = new FinishProcessingNodeEvent(BranchAndPrice.this, node, nodeBound, nodeValue, numberOfCGIterations, masterSolveTime, pricingSolveTime, nrGeneratedColumns);
				}

				listener.finishedColumnGenerationForNode(finishProcessingNodeEvent);
			}

		}

		public void fireBranchEvent(BAPNode<EVRPTW, Route> parentNode, List<BAPNode<EVRPTW, Route>> childNodes) {
			BranchEvent branchEvent = null;

			for(BAPListener listener : this.listeners) {
				if (branchEvent == null) {
				branchEvent = new BranchEvent(BranchAndPrice.this, childNodes.size(), parentNode, new ArrayList<BAPNode>(childNodes));
				}

				listener.branchCreated(branchEvent);
			}

		}

		public void fireTimeOutEvent(BAPNode<EVRPTW, Route> node) {
			TimeLimitExceededEvent timeLimitExceededEvent = null;

			for(BAPListener listener : this.listeners) {
				if (timeLimitExceededEvent == null) {
				timeLimitExceededEvent = new TimeLimitExceededEvent(BranchAndPrice.this, node);
				}

				listener.timeLimitExceeded(timeLimitExceededEvent);
			}

		}
	}
}
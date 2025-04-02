package branchAndPrice;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Iterator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchAndPrice;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.CGListener;
import org.jorlib.frameworks.columnGeneration.io.TimeLimitExceededException;
import org.jorlib.frameworks.columnGeneration.master.OptimizationSense;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;
import org.jorlib.frameworks.columnGeneration.util.MathProgrammingUtil;

import columnGeneration.Master;
import columnGeneration.PricingProblem;
import columnGeneration.Route;
import columnGeneration.customCG;
import ilog.concert.IloColumn;
import ilog.concert.IloException;
import ilog.concert.IloIntVar;
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

		boolean isInteger = true;
		List<Route> solution = node.getSolution();
		for(Route route: solution)
			if(route.value>0+PRECISION && route.value<1-PRECISION) {isInteger = false; break;}

		if(isInteger) return true;
		else {
			//Inherit the routes generated
			List<Route> routesToAdd = new ArrayList<Route>();
			for(Route column: master.getColumns(this.pricingProblem)) {
				if(column.BBnode==-1) {
					column.BBnode=node.nodeID;
					routesToAdd.add(column);
				}
			}
			node.addInitialColumns(routesToAdd);
			//Inherit the cuts generated (not necessary)

			//Solve MIP at root node (optional)
			if(node.nodeID == 0) {
				try {solveIPAtRootNode(node);} 
				catch (IloException e) {e.printStackTrace();}
			}
			return false;
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
		customCG cg=null;
		try {
			dataModel.cleanSRCs();
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
  
		while(!this.queue.isEmpty()) {
			BAPNode<EVRPTW, Route> bapNode = (BAPNode<EVRPTW, Route>)this.queue.poll();
			this.notifier.fireNextNodeEvent(bapNode);
			if (this.nodeCanBePruned(bapNode)) { // If can be pruned by bound BEFORE solving it
				this.notifier.firePruneNodeEvent(bapNode, bapNode.getBound());
				++this.nodesProcessed;
			} else {
				this.graphManipulator.next(bapNode);
				if (bapNode.nodeID != 0) { bapNode.addInitialColumns(this.generateInitialFeasibleSolution(bapNode)); }
	
				try { // Try solving the node
					this.solveBAPNode(bapNode, timeLimit);
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
					++nodesProcessed; continue;
				}
	
				if (this.nodeCanBePruned(bapNode)) { // If can be pruned by bound AFTER solving the node
					this.notifier.firePruneNodeEvent(bapNode, bapNode.getBound());
					++this.nodesProcessed;
				} else if (this.isInfeasibleNode(bapNode)) { // If can be pruned by infeasibility DUE TO artificial columns in the solution
					this.notifier.fireNodeIsInfeasibleEvent(bapNode);
					++this.nodesProcessed;
				} else { // If it is either integer solution and hence pruned by optimality OR fractional and should branch
					if (this.isIntegerNode(bapNode)) {
						int integerObjective = MathProgrammingUtil.doubleToInt(bapNode.getObjective());
						this.notifier.fireNodeIsIntegerEvent(bapNode, bapNode.getBound(), integerObjective);
						if (this.optimizationSenseMaster == OptimizationSense.MINIMIZE && (double)integerObjective < this.upperBoundOnObjective) {
						this.objectiveIncumbentSolution = integerObjective;
						this.upperBoundOnObjective = (double)integerObjective;
						this.incumbentSolution = bapNode.getSolution();
						} else if (this.optimizationSenseMaster == OptimizationSense.MAXIMIZE && (double)integerObjective > this.lowerBoundOnObjective) {
						this.objectiveIncumbentSolution = integerObjective;
						this.lowerBoundOnObjective = (double)integerObjective;
						this.incumbentSolution = bapNode.getSolution();
						}
					} else {

						this.notifier.fireNodeIsFractionalEvent(bapNode, bapNode.getBound(), bapNode.getObjective());
						List<BAPNode<EVRPTW, Route>> newBranches = new ArrayList();
						
						// Initialize branches boolean and Branch Creator
						BranchingRules bc = (BranchingRules)this.branchCreators.iterator().next();

						// Look for Number of Vehicles or Customers Arc Flow branching
						boolean firstBranches = false;
						firstBranches = bc.canPerformFirstBranching(bapNode.getSolution());
						if (firstBranches){ newBranches.addAll(bc.getFirstBranches(bapNode)); }

						if (!firstBranches) {

							this.extendedNotifier.fireLexicographicMasterEvent(bapNode);

							long time=System.currentTimeMillis();

							Master new_Master = ((Master)this.master).copy();
							new_Master.minimizeBatteryDepletion(timeLimit, new ArrayList<Route>(this.master.getColumns(this.pricingProblem)), bapNode.getInequalities(), this.master.getObjective());

							Double obj = master.getObjective();
							this.timeSolvingMaster += (System.currentTimeMillis()-time);

							this.extendedNotifier.fireFinishLexicographicMasterEvent(bapNode, obj);
							
							// Look for EndCT or InitialCT branching
							boolean otherBranches = false;
							otherBranches = bc.canPerformBranching(bapNode.getSolution());
							if (otherBranches){ newBranches.addAll(bc.getBranches(bapNode)); }

						}
	
						if (newBranches.isEmpty()) {
							throw new RuntimeException("BAP encountered fractional solution, but none of the BranchCreators produced any new branches?");
						}
	
						this.queue.addAll(newBranches);
						this.notifier.fireBranchEvent(bapNode, Collections.unmodifiableList(newBranches));
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
}
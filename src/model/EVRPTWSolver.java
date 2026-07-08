package model;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.StartEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.BAPListener;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.BranchEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.CGListener;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.CHListener;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.FinishEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.FinishGeneratingCutsEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.FinishMasterEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.FinishPricingEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.FinishProcessingNodeEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.NodeIsFractionalEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.NodeIsInfeasibleEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.NodeIsIntegerEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.ProcessingNextNodeEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.PruneNodeEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.StartGeneratingCutsEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.StartMasterEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.StartPricingEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.TimeLimitExceededEvent;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.colgenMain.AbstractColumn;
import org.jorlib.frameworks.columnGeneration.io.SimpleDebugger;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractCutGenerator;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.CutHandler;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;
import branchAndPrice.BranchAndPrice;
import branchAndPrice.BranchingRules;
import branchAndPrice.CGFinishFixingByReducedCostEvent;
import branchAndPrice.CGFixingByReducedCostEvent;
import branchAndPrice.CGMasterIsInfeasibleEvent;
import branchAndPrice.CGProblemsLBEvent;
import branchAndPrice.ExtendBAPListener;
import branchAndPrice.MIPMasterEvent;
import branchAndPrice.FinishMIPMasterEvent;
import branchAndPrice.FixArc;
import branchAndPrice.GapReductionEvent;
import branchAndPrice.RemoveArc;
import branchAndPrice.RollbackEvent;
import branchAndPrice.FinishRollbackEvent;
import branchAndPrice.IPSolutionEvent;
import branchAndPrice.FinishIPSolutionEvent;
import columnGeneration.HeuristicMinCostLabelingPricingProblemSolver;
import columnGeneration.HeuristicLabelingSecondPricingProblemSolver;
import columnGeneration.CBHeuristicSecondPricingProblemSolver;
import columnGeneration.CBHeuristicMinCostPricingProblemSolver;
import columnGeneration.Master;
import columnGeneration.PricingProblem;
import columnGeneration.Route;
import columnGeneration.SubsetRowInequality;
import columnGeneration.SubsetRowInequalityGenerator;
import columnGeneration.VRPMasterData;
import model.EVRPTW.Arc;
import model.EVRPTW.PPArc;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Solver class for the mE-VRSPTW (BPC algorithm).
 */
public final class EVRPTWSolver {

	private final EVRPTW dataModel;  		//information about the instance
	public Double upperBound; 				//upper bound on column generation solution (stronger is better).

	public BranchAndPrice bap;
	public CutHandler<EVRPTW, VRPMasterData> cutHandler;
	public boolean isOptimal;

	public EVRPTWSolver(EVRPTW dataModel, ArrayList<Route> initialColumns) throws FileNotFoundException{

		this.dataModel = dataModel;
		this.isOptimal = false;

		//Create a cutHandler, then create a SRC AbstractInequality Generator and add it to the handler
		this.cutHandler = new CutHandler<>();
		SubsetRowInequalityGenerator cutGen = new SubsetRowInequalityGenerator(dataModel);
		this.cutHandler.addCutGenerator(cutGen);

		//Create the pricing problem
		PricingProblem pricingProblem = new PricingProblem(dataModel, "EVRSPTWPricing");

		//Create the master problem
		Master master = new Master(dataModel, pricingProblem, this.cutHandler);

		//Define which solvers to use (one or more)
		List<Class<? extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem>>> solvers = new ArrayList<>(); // The solvers list of classes is restricted to subclasses of AbstractPricingProblemSolver with the specified parameters
	
		//solvers.add(HeuristicLabelingThirdPricingProblemSolver.class);
		//solvers.add(HeuristicLabelingPricingProblemSolver.class);
		solvers.add(HeuristicLabelingSecondPricingProblemSolver.class);
		solvers.add(HeuristicMinCostLabelingPricingProblemSolver.class);

		solvers.add(CBHeuristicSecondPricingProblemSolver.class);
		solvers.add(CBHeuristicMinCostPricingProblemSolver.class);
		
		//Create a set of initial columns and use it as an upper bound
		List<Route> initSolution = this.getInitialSolution(pricingProblem, initialColumns);

		//Define Branch creators
		BranchingRules branchCreator = new BranchingRules(dataModel, pricingProblem);

		//Create a Branch-and-Price instance
		this.bap = new BranchAndPrice(dataModel, master, Collections.singletonList(pricingProblem), solvers, branchCreator, upperBound.intValue(), initSolution);

		//OPTIONAL: Attach a debugger
		PersonalizedDebbuger debugger = new PersonalizedDebbuger(this.bap, this.cutHandler, false);
		this.bap.addExtendCGEventListener(debugger);
		this.bap.addBranchAndPriceEventListener(debugger);
		this.bap.addColumnGenerationEventListener(debugger);
		this.cutHandler.addCHEventListener(debugger);
	}

	public void solve(long timeLimit){

		//Solve the problem problem through Branch-and-Price
		this.bap.runBranchAndPrice(System.currentTimeMillis()+timeLimit);
	}

	public ArrayList<Route> close(){

		//Clean up:
		this.bap.close(); 		//close master and pricing problems
		this.cutHandler.close(); //close the cut handler. The close() call is propagated to all registered AbstractCutGenerator classes
		this.upperBound = getScaledObjective(this.bap.getObjective());
		this.isOptimal = this.bap.isOptimal();
		
		return new ArrayList<>(this.bap.getSolution());
	}

	/** Computes the charging schedule statistics for a given solution. */
	public double[] getChargingInformation(List<Route> solution) {

		//Min and max charging start and end time
		int minChargingTime = dataModel.last_charging_period;
		int maxChargingTime = 1;
		for(Route route:solution) {
			if(route.lastChargingTime-route.chargingTime+1<minChargingTime) minChargingTime = route.lastChargingTime-route.chargingTime+1;
			if(route.lastChargingTime>maxChargingTime) maxChargingTime = route.lastChargingTime;
		}


		//number of vehicles charging each timestep
		int[] vehiclesCharging = new int[maxChargingTime+1];
		for (int t = minChargingTime; t <= maxChargingTime; t++) {
			for(Route route: solution) {
				int charge = (route.lastChargingTime-route.chargingTime+1<=t && t<=route.lastChargingTime) ? 1 : 0;
				vehiclesCharging[t]+= charge;
			}
		}

		//number of timesteps chargers are in use and avg. number of vehicles charging
		double numberOfTimesInUse = 0;
		double totalVehiclesCharging = 0;
		double timestepsAtFullCapacity = 0;
		for (int t = minChargingTime; t <= maxChargingTime; t++) {
			if(vehiclesCharging[t]>0) {
				numberOfTimesInUse++;
				totalVehiclesCharging+=vehiclesCharging[t];
				if(vehiclesCharging[t]==dataModel.B) timestepsAtFullCapacity++;
			}
		}

		double averageTimeInUse = Math.floor(numberOfTimesInUse/(maxChargingTime-minChargingTime+1)*10000)/10000;
		double averageVehiclesCharging = totalVehiclesCharging/(maxChargingTime-minChargingTime+1);
		averageVehiclesCharging = Math.floor(averageVehiclesCharging/dataModel.B*10000)/10000;
		double averageTimeAtFullCapacity = timestepsAtFullCapacity/(maxChargingTime-minChargingTime+1);
		averageTimeAtFullCapacity = Math.floor(averageTimeAtFullCapacity*10000)/10000;
		return new double[] {averageTimeInUse, averageVehiclesCharging, averageTimeAtFullCapacity};
	}

	/**
	 * Create an initial solution for the mE-VRSPTW.
	 * Simple initial solution: visit each customer with a single vehicle and assign a charging schedule. 
	 * @return initial set of routes
	 */
	private List<Route> getInitialSolution(PricingProblem pricingProblem, ArrayList<Route> initialColumns){

		List<Route> initSolution = new ArrayList<>();
		List<ArrayList<Integer>> all_PParcs = new ArrayList<>();
		this.upperBound = 0.0;
		if (initialColumns.isEmpty()) this.upperBound = Math.pow(10, 20);
		for (Route init_col: initialColumns){
			Route new_route = new Route("initSolution", false, (HashMap<Integer, Integer>) init_col.route.clone(), (int[]) init_col.routeSequence.clone(), pricingProblem, init_col.cost, init_col.departureTime, init_col.energy, init_col.load, 0.0, (ArrayList<Integer>) init_col.arcs.clone(), (ArrayList<Integer>) init_col.PParcs.clone(), init_col.lastChargingTime, init_col.chargingTime);
			new_route.BBnode=0;
			initSolution.add(new_route);
			this.upperBound += new_route.cost;
			all_PParcs.add(init_col.PParcs);
		}

		// Dummy (artificial) column to identify infeasibility and initialize the CG
		HashMap<Integer, Integer> route=new HashMap<Integer, Integer>(dataModel.C);
		int[] routeSequence = new int[dataModel.C];
		for(int i=0; i< dataModel.C; i++) {route.put(i+1, 1); routeSequence[i] = i+1;}
		
		initSolution.add(new Route("initSolution", true, route, routeSequence, pricingProblem, (int) Math.pow(10, 20), 0, 0, 0, 0.0, new ArrayList<Integer>(), new ArrayList<>(), 0, 0)); //dummy 

		// Feasible direct shipping routes
		for(int i = 1; i <= dataModel.C; i++){ //a route for each customer
			route = new HashMap<Integer, Integer>();
			route.put(i, 1);
			
			Arc arc_0i = dataModel.graph.getEdge(0, i);
			PPArc arc = dataModel.PPgraph.getEdge(dataModel.C0_startID+i, dataModel.T_startID);
			Arc arc_i0 = arc.routing_arc;
				
			int cost = arc_0i.cost+arc_i0.cost;
			int energy = arc_0i.energy+arc_i0.energy;
			if (dataModel.gamma > 1) {
				energy += arc_0i.energy_deviation + arc_i0.energy_deviation;
			} else if (dataModel.gamma == 1) {
				if (arc_0i.energy_deviation > arc_i0.energy_deviation) energy += arc_0i.energy_deviation;
				else energy += arc_i0.energy_deviation;
			}

			if(energy>dataModel.E) continue;
			routeSequence = new int[] {i};
			ArrayList<Integer> arcs = new ArrayList<Integer>(dataModel.C);
			arcs.add(arc_0i.id);arcs.add(arc_i0.id);
			
			int latestDeparture = dataModel.vertices[i].closing_tw-arc_0i.time;
			latestDeparture = (int) (latestDeparture/10);
			int last_t = latestDeparture-1;
			
			ArrayList<Integer> PParcs = new ArrayList<Integer>(dataModel.C);
			PPArc lastT_0i_arc = dataModel.PPgraph.getEdge(dataModel.PPvertices[dataModel.T_startID+last_t].id, arc.tail_vertex_id);
			PParcs.add(lastT_0i_arc.id); PParcs.add(arc.id);
			if (all_PParcs.contains(PParcs)) continue; // no repeated columns

			//Add the route
			Route column = new Route("initSolution", false, route, routeSequence, pricingProblem, cost, latestDeparture, energy, dataModel.vertices[i].load, 0.0, arcs, PParcs, last_t, dataModel.f_inverse[energy]);
			column.BBnode=0;
			initSolution.add(column);
			
		}

		return initSolution;
	}

	/** Main class. Here the program starts.
	 *  The args[0] must be the instance name. 
	 *  The instance file needs to be in ./data
	 * */
	public static void main(String[] args) throws IOException{

		int gamma = Integer.parseInt(args[1]);
		String energy_deviation = args[3];
		if (!energy_deviation.equals("")) energy_deviation = "-"+energy_deviation;

		EVRPTW evrptw = new EVRPTW(args[0], gamma, 0, true, "FRC-1GR15-RB30"+energy_deviation, args[2], args[3]);
		//EVRPTW evrptw = new EVRPTW("RC105-50", 5, 0, true, "NF", "Debug", "");

		EVRPTWSolver Solver = new EVRPTWSolver(evrptw, new ArrayList<>());

		Solver.solve(10800000L); evrptw.fileOut.close();
		Solver.close();

	}

	public static void deleteStaticObject(Class<?> clazz, String fieldName) {
        try {
            Field field = clazz.getDeclaredField(fieldName);
            field.setAccessible(true); // Allow access to private or protected fields
            field.set(null, null);     // Set the static field to null
        } catch (NoSuchFieldException | IllegalAccessException e) {
            e.printStackTrace();
            throw new RuntimeException("Unable to delete the static object", e);
        }
    }

	/** Returns the real (double) objective (divided by 10). */
	public double getScaledObjective(double objective) {
		double realCost = objective;
		realCost = realCost*0.1+0.05;
		return Math.floor(realCost*10)/10;
	}

	/** Returns the time in seconds (and considering two decimals). */
	public double getTimeInSeconds(double time) {
		double realTime = time*0.001;
		realTime = Math.floor(realTime*100)/100; //two decimals
		return realTime;
	}


	/** Defines a personalized debbuger (for the log files). */
	public class PersonalizedDebbuger implements ExtendBAPListener, BAPListener, CGListener, CHListener{
		
		protected final Logger logger;
		protected String instanceName;
   		protected int bestIntegerSolution;

		public PersonalizedDebbuger(BranchAndPrice bap, CutHandler cutHandler, boolean captureColumnGenerationEventsBAP) {
			this.logger = LoggerFactory.getLogger(SimpleDebugger.class);
		}

		////////////////////////////////
		/// Extended BAP Listener
		////////////////////////////////

		public void IPRootNode(IPSolutionEvent IPRootEvent){
			if (dataModel.print_log) {
				logger.debug("=============== SOLVING IP ===============");
			}
		}

		public void finishIPRootNode(FinishIPSolutionEvent IPRootEvent){
			if (dataModel.print_log) {
				logger.debug("Time solving the IP: "+getTimeInSeconds(IPRootEvent.time));
			}
		}

		public void Rollback(RollbackEvent rollbackEvent){
			if (dataModel.print_log){
				logger.debug("=============== SRCs ROLLBACK ===============");
				logger.debug("Performing a rollback step due to a label explosion in the Pricing Problem.");
				logger.debug("Baseline: " + rollbackEvent.BL + ". Generated: " + rollbackEvent.explosion);
				logger.debug("Returning to the previous optimal RMP state before the last generation of cuts.");
			}
		}

		public void finishRollback(FinishRollbackEvent rollbackEvent){
			if (dataModel.print_log){
				logger.debug("Objective: "+ rollbackEvent.objective);
				logger.debug("Number of columns: " + rollbackEvent.ncols + " Number of SRC separated: " + rollbackEvent.nSRCs);
				logger.debug("Number of vehicle branches: " + rollbackEvent.nVehicleBranches);
				logger.debug("Columns (only non-zero columns are returned):");
				for(Route route: rollbackEvent.solution)
					logger.debug(route.toString());
			}
		}

		public void fixingByReducedCost(CGFixingByReducedCostEvent frcEvent){
			if (dataModel.print_log) {
				logger.debug("================ FIXING BY REDUCED COST ================");
				logger.debug("UB = {} LB = {} Gap = {} ({})",new Object[]{frcEvent.UB, frcEvent.LB, frcEvent.UB-frcEvent.LB, Math.round((1-frcEvent.LB/frcEvent.UB)*1e4)/1e4});
			}
		}

		public void finishFixingByReducedCost(CGFinishFixingByReducedCostEvent frcEvent){
			if (dataModel.print_log) {
				logger.debug("Removing {} arcs.",frcEvent.arcs.size());
				for(Map.Entry<Integer, Double> entry: frcEvent.arcs.entrySet()){
					int arcID = entry.getKey(); PPArc arc = dataModel.PParcs[arcID];
					double rc = entry.getValue();
					logger.debug("{}: RC - RC* = {} - {} = {}",new Object[]{arc.toString(), Math.round(rc*1e4)/1e4, Math.round(frcEvent.best_rc*1e4)/1e4, Math.round((rc-frcEvent.best_rc)*1e4)/1e4});
				}
			}
		}

		public void gapReduction(GapReductionEvent gapEvent){
			if (dataModel.print_log){
				logger.debug("Stopping to add more cuts in the current node.");
				logger.debug("Gap Reduction: "+(Math.floor(gapEvent.gapReduction*100)/100)+" < "+dataModel.gapReductionRequirement);
			}
		}

		public void CGProblemsLB(CGProblemsLBEvent cgProblemsLBEvent){
			if (dataModel.print_log) {
				logger.debug("Check problems with LB!");
			}
		}

		public void startMIPMaster(MIPMasterEvent mipEvent){
			if (dataModel.print_log) {

				logger.debug("================ MASTER - FINDING INTEGER SOLUTION ================");

			}
		}

		public void finishMIPMaster(FinishMIPMasterEvent MIPEvent){
			if (dataModel.print_log) {

				if (MIPEvent.found_integer_solution){
					logger.debug("Found integer solution:");
					List<Route> solution = MIPEvent.node.getSolution();
					for (Route column : solution){
						logger.debug(column.toString());
					}
				} else {
					logger.debug("Did not find integer solution, charging branching proceeds.");
				}

			}
		}

		public void CGMasterIsInfeasible(CGMasterIsInfeasibleEvent cgMasterIsInfeasibleEvent){
			if (dataModel.print_log) {
				logger.debug("CPLEX found the RMP to be infeasible.");
			}
		}

		////////////////////////////////
		/// BAP Listener
		////////////////////////////////

		public void startBAP(StartEvent startEvent) {
			this.instanceName = startEvent.instanceName;
			this.bestIntegerSolution = startEvent.objectiveIncumbentSolution;
			if (dataModel.print_log) this.logger.debug("BAP solving {} - Initial solution: {}", this.instanceName, startEvent.objectiveIncumbentSolution);
		}

		public void finishedColumnGenerationForNode(FinishProcessingNodeEvent finishProcessingNodeEvent) {}

		public void branchCreated(BranchEvent branchEvent) {
			if (dataModel.print_log) {
				logger.debug("================ BRANCHING ================");
				logger.debug("Branching - {} new nodes: ",branchEvent.nrBranches);
				for(BAPNode childNode : branchEvent.childNodes){

					BranchingDecision bd = childNode.getBranchingDecision();
					logger.debug("ChildNode {} - {}",childNode.nodeID, bd.toString());
					
					if (bd instanceof FixArc){
						FixArc fixArcBranch = (FixArc) bd;
						for (int arcID: fixArcBranch.infeasiblePPArcs) logger.debug(dataModel.PParcs[arcID].toString());
					} else if (bd instanceof RemoveArc) {
						RemoveArc removeArcBranch = (RemoveArc) bd;
						logger.debug(dataModel.PParcs[removeArcBranch.arcID].toString());

					}
					
				}
			}
		}

		public void processNextNode(ProcessingNextNodeEvent processingNextNodeEvent) {
			if (dataModel.print_log) {
				logger.debug("================ PROCESSING NODE {} ================", processingNextNodeEvent.node.nodeID);
				logger.debug("Nodes remaining in queue: {} - Node bound: {} - Incumbent solution: {}",new Object[]{processingNextNodeEvent.nodesInQueue, processingNextNodeEvent.node.getBound(), processingNextNodeEvent.objectiveIncumbentSolution});
			}
		}

		public void pruneNode(PruneNodeEvent pruneNodeEvent) {
			if (dataModel.print_log) {
				logger.debug("Pruning node {}. Bound: {}, best integer solution: {}", new Object[]{pruneNodeEvent.node.nodeID, pruneNodeEvent.nodeBound, pruneNodeEvent.bestIntegerSolution});
			}
		}

		public void nodeIsInfeasible(NodeIsInfeasibleEvent nodeIsInfeasibleEvent) {
			if (dataModel.print_log) logger.debug("Node {} is infeasible.", nodeIsInfeasibleEvent.node.nodeID);
		}

		public void nodeIsInteger(NodeIsIntegerEvent nodeIsIntegerEvent) {
			this.bestIntegerSolution = Math.min(this.bestIntegerSolution, nodeIsIntegerEvent.nodeValue);
			if (dataModel.print_log) logger.debug("Node {} is integer. Objective: {} (best integer solution: {})", new Object[]{nodeIsIntegerEvent.node.nodeID, nodeIsIntegerEvent.nodeValue, this.bestIntegerSolution});
		}

		public void nodeIsFractional(NodeIsFractionalEvent nodeIsFractionalEvent) {
			if (dataModel.print_log) logger.debug("Node {} is fractional. Objective: {}, bound: {}", new Object[]{nodeIsFractionalEvent.node.nodeID, nodeIsFractionalEvent.nodeValue, nodeIsFractionalEvent.nodeBound});
		}

		public void timeLimitExceeded(TimeLimitExceededEvent timeLimitExceededEvent) {
			if (timeLimitExceededEvent.node != null) {
				this.logger.debug("Caught timeout exception while processing node {}", timeLimitExceededEvent.node.nodeID);
			}
		}

		public void finishBAP(FinishEvent finishEvent) {

			if (dataModel.print_log) {
				logger.debug("================ SOLUTION BPC - " + instanceName +" ================");
				logger.debug("BAP terminated with objective: "+getScaledObjective(bap.getObjective()));
				logger.debug("Total Number of iterations: "+bap.getTotalNrIterations());
				logger.debug("Total Number of processed nodes: "+bap.getNumberOfProcessedNodes());
				logger.debug("Total Time spent on master problems (s): "+getTimeInSeconds(bap.getMasterSolveTime())+" Total time spent on pricing problems (s): "+getTimeInSeconds(bap.getPricingSolveTime()));
				logger.debug("Total running time (s): "+ getTimeInSeconds(System.currentTimeMillis()-bap.getSolveTime()));
				if(bap.hasSolution()) {
					logger.debug("Solution is optimal: "+bap.isOptimal());
					logger.debug("Columns (only non-zero columns are returned):");
					List<Route> solution = bap.getSolution();
					for (Route column : solution){
						logger.debug(column.toString());
					}
				}
			} else {

				logger.debug("B = "+dataModel.B+", Objective = "+getScaledObjective(bap.getObjective())+", K = "+bap.getSolution().size()+", Optimal: "+bap.isOptimal());
			}

		}
		
		////////////////////////////
		/// CG Listener
		////////////////////////////

		public void startCG(StartEvent startEvent) {
			if (dataModel.print_log) {
				this.instanceName = startEvent.instanceName;
				this.bestIntegerSolution = startEvent.objectiveIncumbentSolution;
				this.logger.debug("CG solving {} - Initial upper bound: {}", this.instanceName, startEvent.objectiveIncumbentSolution);
			}

		}

		public void finishCG(FinishEvent finishEvent) {
			if (dataModel.print_log) {
				this.logger.debug("Finished Column Generation for instance {}", this.instanceName);
			}

		}

		public void startMaster(StartMasterEvent startMasterEvent) {
      		if (dataModel.print_log) logger.debug("=============== MASTER {} ===============", startMasterEvent.columnGenerationIteration);
   		}

		public void finishMaster(FinishMasterEvent finishMasterEvent) {
			if (dataModel.print_log) {
				logger.debug("Finished master -> RMP objective: {}, LB: {}, UB: {}", new Object[]{finishMasterEvent.objective, finishMasterEvent.boundOnMasterObjective, finishMasterEvent.cutoffValue});
				logger.debug("Total running time (s): " + getTimeInSeconds(System.currentTimeMillis()-bap.getSolveTime()));
			}
		}

		public void startPricing(StartPricingEvent startPricing) {
			if (dataModel.print_log) logger.debug("=============== PRICING {} ===============", startPricing.columnGenerationIteration);
   		}

		public void finishPricing(FinishPricingEvent finishPricingEvent) {
			if (dataModel.print_log){
				String solver = (finishPricingEvent.columns.size()==0) ? "exactLabeling": finishPricingEvent.columns.get(0).creator;
				logger.debug("Finished pricing ({}, {} columns generated) -> CG objective: {}, CG bound: {}, CG cutoff: {}", new Object[]{solver, finishPricingEvent.columns.size(), finishPricingEvent.objective, finishPricingEvent.boundOnMasterObjective, finishPricingEvent.cutoffValue});
				for(AbstractColumn<?, ?> column : finishPricingEvent.columns){
					logger.debug(column.toString());
				}
			}
		}

		//////////////////////////////
		/// 
		/////////////////////////////
		
		public void startGeneratingCuts(StartGeneratingCutsEvent startGenerateCutsEvent) {
      		if (dataModel.print_log) this.logger.debug("=============== GENERATING CUTS ===============");
   		}

		public void finishGeneratingCuts(FinishGeneratingCutsEvent finishGenerateCutsEvent) {
			if (dataModel.print_log) {
				Map<AbstractCutGenerator, Integer> cutSummary = new LinkedHashMap();
				if (finishGenerateCutsEvent.separatedInequalities.isEmpty()) {
					this.logger.debug("No inequalities found!");
				} else {
					this.logger.debug("Cuts have been generated! Summary:");
					Iterator var3 = finishGenerateCutsEvent.separatedInequalities.iterator();

					while(var3.hasNext()) {
						SubsetRowInequality inequality = (SubsetRowInequality)var3.next();
						if (cutSummary.containsKey(inequality.maintainingGenerator)) {
						cutSummary.put(inequality.maintainingGenerator, (Integer)cutSummary.get(inequality.maintainingGenerator) + 1);
						} else {
						cutSummary.put(inequality.maintainingGenerator, 1);
						}
					}

					String summary = (String)cutSummary.keySet().stream().map((i) -> {
						return "-" + i.toString() + ": " + cutSummary.get(i);
					}).collect(Collectors.joining(", "));
					this.logger.debug(summary);
				}
			}
		}

	}
}
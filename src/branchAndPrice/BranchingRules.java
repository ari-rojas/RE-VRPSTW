package branchAndPrice;

import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;
import columnGeneration.PricingProblem;
import columnGeneration.Route;
import model.EVRPTW;
import model.EVRPTW.PPArc;

/**
 * Class which creates new branches in the Branch-and-Price tree. 
 * The class checks whether a fractional number of vehicles is used
 * The class checks whether there is a fractional arc in the solution 
 * The edge with a fractional value closest to 0.5 is selected for branching.
 * Two important methods are:
 * 	1. canPerformBranching that determines whether the particular branch creator can create the child nodes (there is a fractional arc to branch on)
 *  2. getBranches creates the actual branches
 */

public final class BranchingRules extends AbstractBranchCreator<EVRPTW, Route, PricingProblem>{

	private double vehiclesForBranching=0; 				//number of vehicles used in a solution
	public boolean branchingOnVehicles; 				//true if the branching is on the number of vehicles
	public boolean branchOnCustomerArcs; 				//true if the branching is performed on an arc between customers (or the depot)
	private int arcForBranching=-1; 					//arc to branch on
	private byte arcType;
	private double bestArcValue = 0; 					//current flow value of the arc to branch on
	private EVRPTW dataModel; 							//model data

	private double PRECISION = 0.001;
	private final int depotID;

	public BranchingRules(EVRPTW dataModel, PricingProblem pricingProblem){
		super(dataModel, pricingProblem);
		this.dataModel = dataModel;
		this.depotID = dataModel.T_startID;
	}

	/**
	 * Determine the next branching decision
	 * It can be the number of vehicles used
	 * Or if that is an integer number then a fractional arc.
	 * @param solution Fractional column generation solution
	 * @return true if a fractional number of vehicles is used or a fractional arc exists
	 */
	public boolean canPerformFirstBranching(List<Route> solution) {

		//Reset values
		this.vehiclesForBranching = 0;
		this.branchingOnVehicles = false;
		this.branchOnCustomerArcs = false;
		this.arcForBranching = -1;
		this.arcType = -1;
		this.bestArcValue = 0;

		// Aggregate route values
		for(Route route : solution){vehiclesForBranching+=route.value;}
		if(isFractional(vehiclesForBranching)) {branchingOnVehicles = true; return true;}

		// Determine whether there's a fractional routing arc for branching
		// The array separates the arcs by arc_type: 0 belonging to AR0, 1 belonging to AR1
		Map<Integer, Double>[] arcValues= (Map<Integer, Double>[]) new Map[2];
		for (int i = 0; i < EVRPTW.AC1; i++) arcValues[i] = new HashMap<>();

		//Aggregate edge values
		for(Route route : solution){
			if (route.value < 1-this.PRECISION) {
				ArrayList<Integer> PParcs = route.PParcs;
				for(int ix = 1; ix < PParcs.size(); ix++){ // Skip the first (ix = 0), as it is the AC1 arc of the column
					int arcID = PParcs.get(ix);
					PPArc arc = dataModel.PParcs[arcID];
					Double arcValue = arcValues[arc.arc_type].get(arcID);

					if (arcValue == null) arcValues[arc.arc_type].put(arcID, route.value);
					else arcValues[arc.arc_type].put(arcID, route.value+arcValue);
				}
			}
		}

		// Select the arc whose flow is closest to 0.5
		// Prioritize arcs belonging to AR1
		this.arcType = EVRPTW.AR1;
		for(int arcID : arcValues[1].keySet()){
			double value = arcValues[1].get(arcID);
			if(Math.abs(0.5-value) <= Math.abs(0.5- bestArcValue)){
				arcForBranching = arcID; bestArcValue = value;
				if (bestArcValue == 0.5 && dataModel.PParcs[arcID].head_vertex_id != depotID) { branchOnCustomerArcs = true; return true; }
			}
		}
		if (isFractional(bestArcValue)) { branchOnCustomerArcs = true; return true; }

		this.arcType = EVRPTW.AR0;
		// If no arcs in AR1 are fractional, check for fractional arcs belonging to AR0
		for(int arcID : arcValues[0].keySet()){
			double value = arcValues[0].get(arcID);
			if(Math.abs(0.5-value) <= Math.abs(0.5- bestArcValue)){
				arcForBranching = arcID; bestArcValue = value;
				if (bestArcValue == 0.5 && dataModel.PParcs[arcID].head_vertex_id != depotID) { branchOnCustomerArcs = true; return true; }
			}
		}
		if (isFractional(bestArcValue)) { branchOnCustomerArcs = true; return true; }

		return false;
	}

	@Override
	public boolean canPerformBranching(List<Route> solution) {

		this.arcForBranching = -1;
		this.arcType = 2;
		this.bestArcValue = 0;

		// Determine whether there's a fractional routing arc for branching
		// The array separates the arcs by arc_type: 0 belonging to AR0, 1 belonging to AR1
		Map<Integer, Double> arcValues = new HashMap<Integer, Double>();

		//Aggregate edge values
		for(Route route : solution){
			if (route.value < 1-this.PRECISION) {
				int arcID = route.PParcs.get(0);
				Double arcValue = arcValues.get(arcID);

				if (arcValue == null) arcValues.put(arcID, route.value);
				else arcValues.put(arcID, route.value+arcValue);
			}
		}

		// Select the arc whose flow is closest to 0.5
		this.arcType = EVRPTW.AC1;
		for(int arcID : arcValues.keySet()){
			double value = arcValues.get(arcID);
			if(Math.abs(0.5-value) <= Math.abs(0.5- bestArcValue)){
				arcForBranching = arcID; bestArcValue = value;
				if(bestArcValue == 0.5 && dataModel.PParcs[arcID].head_vertex_id != dataModel.T_startID) { branchOnCustomerArcs = true; return true; }
			}
		}
		if (isFractional(bestArcValue)) { branchOnCustomerArcs = true; return true; }

		return false;
	}

	public List<BAPNode<EVRPTW,Route>> getFirstBranches(BAPNode<EVRPTW,Route> parentNode) {
		
		BAPNode<EVRPTW,Route> node2; 		//one child node
		BAPNode<EVRPTW,Route> node1; 		//other child node

		if(branchingOnVehicles) {
			//Branch 1: number of vehicles down
			BranchVehiclesDown branchingDecision1 = new BranchVehiclesDown(this.pricingProblems.get(0), (int) Math.floor(vehiclesForBranching), parentNode.getInequalities());
			node1=this.createBranch(parentNode, branchingDecision1, parentNode.getInitialColumns(), parentNode.getInequalities());
			//Branch 2: number of vehicles up
			BranchVehiclesUp branchingDecision2 = new BranchVehiclesUp(this.pricingProblems.get(0), (int) Math.ceil(vehiclesForBranching), parentNode.getInequalities());
			node2=this.createBranch(parentNode, branchingDecision2, parentNode.getInitialColumns(), parentNode.getInequalities());
		} else {
			//Branch 1: remove the edge:
			RemoveArc branchingDecision1 = new RemoveArc(this.pricingProblems.get(0), arcForBranching, arcType, dataModel, parentNode.getInequalities(), bestArcValue);
			node2=this.createBranch(parentNode, branchingDecision1, parentNode.getInitialColumns(), parentNode.getInequalities());
			//Branch 2: fix the edge:
			FixArc branchingDecision2 = new FixArc(this.pricingProblems.get(0), arcForBranching, arcType, dataModel, parentNode.getInequalities(), bestArcValue);
			node1=this.createBranch(parentNode, branchingDecision2, parentNode.getInitialColumns(), parentNode.getInequalities());
		}
		
		return Arrays.asList(node1,node2);
	}

	/**
	 * Create the branches:
	 * branch 1: edge {@code edgeForBranching} must be used by {@code PricingProblem},</li>
	 * 	branch 2: edge {@code edgeForBranching} may NOT used by {@code PricingProblem},</li>
	 * @param parentNode Fractional node on which we branch
	 * @return List of child nodes
	 */
	@Override
	public List<BAPNode<EVRPTW,Route>> getBranches(BAPNode<EVRPTW,Route> parentNode) {
		
		BAPNode<EVRPTW,Route> node2; 		//one child node
		BAPNode<EVRPTW,Route> node1; 		//other child node
		
		//Branch 1: remove the edge:
		RemoveArc branchingDecision1 = new RemoveArc(this.pricingProblems.get(0), arcForBranching, arcType, dataModel, parentNode.getInequalities(), bestArcValue);
		node2=this.createBranch(parentNode, branchingDecision1, parentNode.getInitialColumns(), parentNode.getInequalities());
		//Branch 2: fix the edge:
		FixArc branchingDecision2 = new FixArc(this.pricingProblems.get(0), arcForBranching, arcType, dataModel, parentNode.getInequalities(), bestArcValue);
		node1=this.createBranch(parentNode, branchingDecision2, parentNode.getInitialColumns(), parentNode.getInequalities());
		
		return Arrays.asList(node1,node2);
	}

	private boolean isFractional(double value) {
		return Math.abs(value - (double)Math.round(value)) > this.PRECISION;
	}

}
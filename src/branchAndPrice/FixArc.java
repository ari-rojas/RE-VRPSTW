package branchAndPrice;

import java.util.ArrayList;
import java.util.List;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractInequality;
import columnGeneration.PricingProblem;
import columnGeneration.Route;
import model.EVRPTW;
import model.EVRPTW.PPVertex;
import model.EVRPTW.Arc;
import model.EVRPTW.PPArc;


/**
 * Ensure that an arc is used (Branching on the arc-flow variables >=)
 */
public final class FixArc implements BranchingDecision<EVRPTW,Route> {


	public final PricingProblem pricingProblem;				//pricing problem
	public final int arcID;									//arc on which we branch
	public final byte arc_type;
	public double flowValue;								//flow value of the arc on which we are branching
	public List<AbstractInequality> poolOfCuts;				//separated SRCs
	public EVRPTW dataModel;								//data model
	public ArrayList<Integer> infeasiblePPArcs;				//infeasible arcs by the branching decision

	public FixArc(PricingProblem pricingProblem, int arc, byte arc_type, EVRPTW dataModel, List<AbstractInequality> list, double flowValue){
		this.pricingProblem = pricingProblem;
		this.arcID = arc;
		this.arc_type = arc_type;
		this.dataModel = dataModel;
		this.poolOfCuts = list;
		this.infeasiblePPArcs = new ArrayList<Integer>();
		this.flowValue = flowValue;

		// Retrieve the arc (i,j)
		PPArc pparc = dataModel.PParcs[arc];

		// Remove all other incoming arcs of j
		for (PPArc other_arc: dataModel.PPgraph.incomingEdgesOf(pparc.head_vertex_id)) if (other_arc.id != arc) this.infeasiblePPArcs.add(other_arc.id);

		// Remove the node i^(1-\kappa_(i,j))
		// \kappa_(i,j): 1 if (i,j) \in AR0 \cup AC1, 0 otherwise
		int customer_depot = dataModel.PPvertices[pparc.tail_vertex_id].node_number;
		if (arc_type == EVRPTW.AC1) customer_depot = dataModel.PPvertices[pparc.head_vertex_id].node_number;

		int startID = dataModel.C1_startID;
		if (arc_type == EVRPTW.AR1) startID = dataModel.C0_startID;
		PPVertex vx_to_remove = dataModel.PPvertices[startID+customer_depot];
		remove_node(vx_to_remove);

		// If it's a routing arc and it does not end in the depot, remove all other outgoing arcs of i, and remove the j0 node
		if (arc_type <= EVRPTW.AR1 && pparc.head_vertex_id != dataModel.T_startID){
			for (PPArc other_arc: dataModel.PPgraph.outgoingEdgesOf(pparc.tail_vertex_id)) if (other_arc.id != arc) this.infeasiblePPArcs.add(other_arc.id);

			customer_depot = dataModel.PPvertices[pparc.head_vertex_id].node_number;
			vx_to_remove = dataModel.PPvertices[dataModel.C0_startID+customer_depot];
			remove_node(vx_to_remove);
		}
		
	}

	private void remove_node(PPVertex vx_to_remove){

		for (PPArc other_arc: dataModel.PPgraph.incomingEdgesOf(vx_to_remove.id)) this.infeasiblePPArcs.add(other_arc.id);
		for (PPArc other_arc: dataModel.PPgraph.outgoingEdgesOf(vx_to_remove.id)) this.infeasiblePPArcs.add(other_arc.id);
	}

	/**
	 * Determine whether the given inequality remains feasible for the child node
	 * @param inequality inequality
	 * @return true
	 */
	@Override
	public boolean inEqualityIsCompatibleWithBranchingDecision(AbstractInequality inequality) {
		return true;
	}

	/**
	 * Determine whether the given column remains feasible for the child node
	 * @param column column
	 * @return true if the column is compliant with the branching decision
	 */
	@Override
	public boolean columnIsCompatibleWithBranchingDecision(Route column) {
		
		if(column.associatedPricingProblem != this.pricingProblem) return false;
		if(column.isArtificialColumn) return true;

		if (this.arc_type <= EVRPTW.AR1){ // If the branching arc is a routing arc
			for (int edge: this.infeasiblePPArcs) if (column.PParcs.contains(edge)) return false; }
		else if (column.PParcs.get(0) == this.arcID) return false;

		return true;
	}

	@Override
	public String toString(){
		return "Fix: "+ dataModel.PParcs[arcID].toString() + " Current flow-value: " + this.flowValue;
	}
}
package branchAndPrice;

import java.util.List;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractInequality;
import columnGeneration.PricingProblem;
import columnGeneration.Route;
import model.EVRPTW;


/**
 * Prevent an arc from being used (Branching on the arc-flow variables <=)
 */
public final class RemoveArc implements BranchingDecision<EVRPTW, Route> {

	public final PricingProblem pricingProblem;				//pricing problem
	public final int arcID;									//arc on which we branch
	public final byte arc_type;									//arc on which we branch
	public double flowValue;								//flow value of the arc on which we are branching
	public List<AbstractInequality> poolOfCuts;				//separated SRCs
	public EVRPTW dataModel;								//data model

	public RemoveArc(PricingProblem pricingProblem, int arc, byte arc_type, EVRPTW dataModel, List<AbstractInequality> list, double flowValue){
		this.pricingProblem = pricingProblem;
		this.arcID = arc;
		this.arc_type = arc_type;
		this.dataModel = dataModel;
		this.poolOfCuts = list;
		this.flowValue = flowValue;
	}

	/**
	 * Determine whether the given inequality remains feasible for the child node
	 * @param inequality inequality
	 * @return true
	 */
	@Override
	public boolean inEqualityIsCompatibleWithBranchingDecision(AbstractInequality inequality) {
		return true;  //We only have SRC. They remain valid, independent of whether we remove an arc.
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
		
		if(column.PParcs.contains(arcID)) return false;
		
		return true;
	}

	@Override
	public String toString(){
		return "Remove: "+ dataModel.PParcs[arcID].toString() + " Current flow-value: " + this.flowValue;
	}
}
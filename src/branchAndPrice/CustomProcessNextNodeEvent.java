package branchAndPrice;

import java.util.EventObject;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;

public class CustomProcessNextNodeEvent extends EventObject{
    public final BAPNode node;
    public final int nodesInQueue;
    public final double objectiveIncumbentSolution;

    public CustomProcessNextNodeEvent(Object source, BAPNode node, int nodesInQueue, double objectiveIncumbentSolution) {
        super(source);
        this.node = node;
        this.nodesInQueue = nodesInQueue;
        this.objectiveIncumbentSolution = objectiveIncumbentSolution;
    }
}

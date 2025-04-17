package branchAndPrice;

import java.util.EventObject;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;

public class CustomPruneNodeEvent extends EventObject {
    public final BAPNode node;
    public final double nodeBound;
    public final double bestIntegerSolution;

    public CustomPruneNodeEvent(Object source, BAPNode node, double nodeBound, double bestIntegerSolution) {
        super(source);
        this.node = node;
        this.nodeBound = nodeBound;
        this.bestIntegerSolution = bestIntegerSolution;
    }
}

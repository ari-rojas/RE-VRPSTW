package branchAndPrice;

import java.util.EventObject;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;

public class CustomNodeIsIntegerEvent extends EventObject{
    public final BAPNode node;
    public final double bound;
    public final double objective;

    public CustomNodeIsIntegerEvent(Object source, BAPNode node, double bound, double objective) {
        super(source);
        this.node = node;
        this.bound = bound;
        this.objective = objective;
    }

}

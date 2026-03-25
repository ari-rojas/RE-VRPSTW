package branchAndPrice;

import java.util.EventObject;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;

public class FinishIPRootNodeEvent extends EventObject{

    public final BAPNode node;
    public final double time;

    public FinishIPRootNodeEvent(Object source, BAPNode node, double time) {
        super(source);
        this.node = node;
        this.time = time;
    }

    
}
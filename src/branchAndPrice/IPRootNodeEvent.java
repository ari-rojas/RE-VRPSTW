package branchAndPrice;

import java.util.EventObject;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;

public class IPRootNodeEvent extends EventObject{

    public final BAPNode node;

    public IPRootNodeEvent(Object source, BAPNode node) {
        super(source);
        this.node = node;
    }

    
}

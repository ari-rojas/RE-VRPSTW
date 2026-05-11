package branchAndPrice;

import java.util.EventObject;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;

public class MIPMasterEvent extends EventObject{
    public final BAPNode node;

    public MIPMasterEvent(Object source, BAPNode node) {
        super(source);
        this.node = node;
    }

    
}

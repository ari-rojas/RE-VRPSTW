package branchAndPrice;

import java.util.EventObject;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;

public class CGMasterIsInfeasibleEvent extends EventObject {
    public final BAPNode node;

    public CGMasterIsInfeasibleEvent(Object source, BAPNode node) {
        super(source);
        this.node = node;
    }

}

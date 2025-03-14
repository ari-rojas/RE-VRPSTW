package branchAndPrice;

import java.util.EventObject;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;

public class CGProblemsLBEvent extends EventObject {
    public final BAPNode node;

    public CGProblemsLBEvent(Object source, BAPNode node) {
        super(source);
        this.node = node;
    }
}

package branchAndPrice;

import java.util.EventObject;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;

public class CGFixingByReducedCostEvent extends EventObject {
    public final BAPNode node;
    public final double UB;
    public final double LB;

    public CGFixingByReducedCostEvent(Object source, BAPNode node, double UB, double LB) {
        super(source);
        this.node = node;
        this.UB = UB;
        this.LB = LB;
    }

}

package branchAndPrice;

import java.util.EventObject;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;

public class CGFixingByReducedCostEvent extends EventObject {
    
    public final double UB;
    public final double LB;

    public CGFixingByReducedCostEvent(Object source, double UB, double LB) {
        super(source);
        this.UB = UB;
        this.LB = LB;
    }

}

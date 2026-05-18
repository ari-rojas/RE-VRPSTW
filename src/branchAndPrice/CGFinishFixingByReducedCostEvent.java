package branchAndPrice;

import java.util.List;
import java.util.Map;
import java.util.EventObject;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;

public class CGFinishFixingByReducedCostEvent extends EventObject {
    
    public final Map<Integer, Double> arcs;
    public final double best_rc;

    public CGFinishFixingByReducedCostEvent(Object source, Map<Integer, Double> arcs, double best_rc){
        super(source);
        this.arcs = arcs;
        this.best_rc = best_rc;

    }
}

package branchAndPrice;


import java.util.EventObject;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;

public class FinishLexicographicMasterEvent extends EventObject{

    public final BAPNode node;
    public final double depletion;
    public final double cost;

    public FinishLexicographicMasterEvent(Object source, BAPNode node, double objective, double cost) {
        super(source);
        this.node = node;
        this.depletion = objective;
        this.cost = cost;
    }

    
}

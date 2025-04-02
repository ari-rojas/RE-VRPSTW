package branchAndPrice;


import java.util.EventObject;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;

public class FinishLexicographicMasterEvent extends EventObject{

    public final BAPNode node;
    public final Double objective;

    public FinishLexicographicMasterEvent(Object source, BAPNode node, Double objective) {
        super(source);
        this.node = node;
        this.objective = objective;
    }

    
}

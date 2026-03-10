package branchAndPrice;


import java.util.EventObject;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;

public class FinishLexicographicMasterEvent extends EventObject{

    public final BAPNode node;
    public final boolean found_integer_solution;

    public FinishLexicographicMasterEvent(Object source, BAPNode node, boolean found) {
        super(source);
        this.node = node;
        this.found_integer_solution = found;
    }

    
}

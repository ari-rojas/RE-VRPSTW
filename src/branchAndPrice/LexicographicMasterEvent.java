package branchAndPrice;

import java.util.EventObject;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;

public class LexicographicMasterEvent extends EventObject{
    public final BAPNode node;

    public LexicographicMasterEvent(Object source, BAPNode node) {
        super(source);
        this.node = node;
    }

    
}

package branchAndPrice;

import java.util.List;
import java.util.EventObject;

import columnGeneration.Route;

public class RollbackEvent extends EventObject{

    public final int BL;
    public final int explosion;

    public RollbackEvent(Object source, int BL, int explosion) {
        super(source);
        this.BL = BL;
        this.explosion = explosion;
    }

    
}
package branchAndPrice;

import java.util.EventObject;

public class FinishIPSolutionEvent extends EventObject{

    public final double time;

    public FinishIPSolutionEvent(Object source, double time) {
        super(source);
        this.time = time;
    }

    
}
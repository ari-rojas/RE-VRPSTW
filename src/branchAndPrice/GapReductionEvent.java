package branchAndPrice;

import java.util.EventObject;

public class GapReductionEvent extends EventObject{

    public double gapReduction;

    public GapReductionEvent(Object source, double gR) {
        super(source);
        this.gapReduction = gR;
    }

    
}

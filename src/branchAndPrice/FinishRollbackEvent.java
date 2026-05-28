package branchAndPrice;

import java.util.List;
import java.util.EventObject;

import columnGeneration.Route;

public class FinishRollbackEvent extends EventObject{

    public final List<Route> solution;
    public final double objective;
    public final double ncols;
    public final int nSRCs;
    public final int nVehicleBranches;

    public FinishRollbackEvent(Object source, List<Route> solution, double obj, int n, int nSRCs, int nVB) {
        super(source);
        this.solution = solution;
        this.objective = obj;
        this.ncols = n;
        this.nSRCs = nSRCs;
        this.nVehicleBranches = nVB;
    }

    
}
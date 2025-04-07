package branchAndPrice;

import java.util.EventListener;

public interface ExtendBAPListener extends EventListener {
    
    void CGMasterIsInfeasible(CGMasterIsInfeasibleEvent var1);

    void CGProblemsLB(CGProblemsLBEvent var1);

    void startLexicographicMaster(LexicographicMasterEvent var1);

    void finishLexicographicMaster(FinishLexicographicMasterEvent var1);

}

package branchAndPrice;

import java.util.EventListener;

public interface ExtendBAPListener extends EventListener {
    
    void CGMasterIsInfeasible(CGMasterIsInfeasibleEvent var1);

    void CGProblemsLB(CGProblemsLBEvent var1);

    void startLexicographicMaster(LexicographicMasterEvent var1);

    void finishLexicographicMaster(FinishLexicographicMasterEvent var1);

    void IPRootNode(IPSolutionEvent var1);

    void finishIPRootNode(FinishIPSolutionEvent var1);

    void Rollback(RollbackEvent var1);

    void finishRollback(FinishRollbackEvent var1);

    void gapReduction(GapReductionEvent var1);
}

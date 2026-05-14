package branchAndPrice;

import java.util.EventListener;

public interface ExtendBAPListener extends EventListener {
    
    void CGMasterIsInfeasible(CGMasterIsInfeasibleEvent var1);

    void CGProblemsLB(CGProblemsLBEvent var1);

    void startMIPMaster(MIPMasterEvent var1);

    void finishMIPMaster(FinishMIPMasterEvent var1);

    void IPRootNode(IPRootNodeEvent var1);

    void finishIPRootNode(FinishIPRootNodeEvent var1);

    void Rollback(RollbackEvent var1);

    void finishRollback(FinishRollbackEvent var1);

}

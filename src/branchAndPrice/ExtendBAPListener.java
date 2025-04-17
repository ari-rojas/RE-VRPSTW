package branchAndPrice;

import java.util.EventListener;

public interface ExtendBAPListener extends EventListener {
    
    void CGMasterIsInfeasible(CGMasterIsInfeasibleEvent var1);

    void CGProblemsLB(CGProblemsLBEvent var1);

    void startLexicographicMaster(LexicographicMasterEvent var1);

    void finishLexicographicMaster(FinishLexicographicMasterEvent var1);

    void customNodeIsInteger(CustomNodeIsIntegerEvent var1);

    void customProcessNextNode(CustomProcessNextNodeEvent var1);

    void customStartBAP(CustomStartBAPEvent var1);

    void customPruneNode(CustomPruneNodeEvent var1);

}

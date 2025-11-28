package branchAndPrice;

import java.util.Iterator;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchAndPrice;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;

public class ExtendBAPNotifier{

  private final AbstractBranchAndPrice<?, ?, ?> parent;
  private final List<ExtendBAPListener> customListeners;
    
    public ExtendBAPNotifier(AbstractBranchAndPrice<?, ?, ?> parent) {
      
      this.parent = parent;
      this.customListeners = new ArrayList<ExtendBAPListener>();
      
    }

    public void addExtendBAPListener(ExtendBAPListener listener) {
      this.customListeners.add(listener);
    }

    public void removeExtendBAPListener(ExtendBAPListener listener) {
      this.customListeners.remove(listener);
    }

    public void fireCGMasterIsInfeasibleEvent(BAPNode node) {
      CGMasterIsInfeasibleEvent cgMasterIsInfeasibleEvent = null;

      ExtendBAPListener listener;
      for(Iterator var3 = this.customListeners.iterator(); var3.hasNext(); listener.CGMasterIsInfeasible(cgMasterIsInfeasibleEvent)) {
        listener = (ExtendBAPListener)var3.next();
        if (cgMasterIsInfeasibleEvent == null) {
          cgMasterIsInfeasibleEvent = new CGMasterIsInfeasibleEvent(this.parent, node);
        }
      }

    }

    public void fireCGProblemsLBEvent(BAPNode node) {
      CGProblemsLBEvent cgProblemsLBEvent = null;

      ExtendBAPListener listener;
      for(Iterator var3 = this.customListeners.iterator(); var3.hasNext(); listener.CGProblemsLB(cgProblemsLBEvent)) {
        listener = (ExtendBAPListener)var3.next();
        if (cgProblemsLBEvent == null) {
          cgProblemsLBEvent = new CGProblemsLBEvent(this.parent, node);
        }
      }

    }

    public void fireLexicographicMasterEvent(BAPNode node) {
      LexicographicMasterEvent lexiMasterEvent = null;

      ExtendBAPListener listener;
      for(Iterator var3 = this.customListeners.iterator(); var3.hasNext(); listener.startLexicographicMaster(lexiMasterEvent)) {
        listener = (ExtendBAPListener)var3.next();
        if (lexiMasterEvent == null) {
          lexiMasterEvent = new LexicographicMasterEvent(this.parent, node);
        }
      }
    }

    public void fireFinishLexicographicMasterEvent(BAPNode node, double obj, double cost) {
      FinishLexicographicMasterEvent lexiMasterEvent = null;

      ExtendBAPListener listener;
      for(Iterator var3 = this.customListeners.iterator(); var3.hasNext(); listener.finishLexicographicMaster(lexiMasterEvent)) {
        listener = (ExtendBAPListener)var3.next();
        if (lexiMasterEvent == null) {
          lexiMasterEvent = new FinishLexicographicMasterEvent(this.parent, node, obj, cost);
        }
      }
    }

    public void fireFixingByReducedCostEvent(BAPNode node, double UB, double LB) {
      CGFixingByReducedCostEvent frcEvent = null;

      ExtendBAPListener listener;
      for(Iterator var3 = this.customListeners.iterator(); var3.hasNext(); listener.fixingByReducedCost(frcEvent)) {
        listener = (ExtendBAPListener)var3.next();
        if (frcEvent == null) {
          frcEvent = new CGFixingByReducedCostEvent(this.parent, node, UB, LB);
        }
      }
    }

    public void fireFinishFixingByReducedCostEvent(BAPNode node, Map<Integer, Double> arcs, double best_rc) {
      CGFinishFixingByReducedCostEvent frcEvent = null;

      ExtendBAPListener listener;
      for(Iterator var3 = this.customListeners.iterator(); var3.hasNext(); listener.finishFixingByReducedCost(frcEvent)) {
        listener = (ExtendBAPListener)var3.next();
        if (frcEvent == null) {
          frcEvent = new CGFinishFixingByReducedCostEvent(this.parent, node, arcs, best_rc);
        }
      }
    }


}

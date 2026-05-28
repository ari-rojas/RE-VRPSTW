package branchAndPrice;

import java.util.Iterator;
import java.util.ArrayList;
import java.util.List;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchAndPrice;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;

import columnGeneration.Route;

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

    public void fireFinishLexicographicMasterEvent(BAPNode node, boolean found) {
      FinishLexicographicMasterEvent lexiMasterEvent = null;

      ExtendBAPListener listener;
      for(Iterator var3 = this.customListeners.iterator(); var3.hasNext(); listener.finishLexicographicMaster(lexiMasterEvent)) {
        listener = (ExtendBAPListener)var3.next();
        if (lexiMasterEvent == null) {
          lexiMasterEvent = new FinishLexicographicMasterEvent(this.parent, node, found);
        }
      }
    }

    public void fireRollbackEvent(int BL, int explosion){
      RollbackEvent rbEvent = null;

      ExtendBAPListener listener;
      for(Iterator var3 = this.customListeners.iterator(); var3.hasNext(); listener.Rollback(rbEvent)) {
        listener = (ExtendBAPListener)var3.next();
        if (rbEvent == null) {
          rbEvent = new RollbackEvent(this.parent, BL, explosion);
        }
      }
    }

    public void fireFinishRollbackEvent(List<Route> sol, double obj, int n, int nSRCs, int nVB){
      FinishRollbackEvent rbEvent = null;

      ExtendBAPListener listener;
      for(Iterator var3 = this.customListeners.iterator(); var3.hasNext(); listener.finishRollback(rbEvent)) {
        listener = (ExtendBAPListener)var3.next();
        if (rbEvent == null) {
          rbEvent = new FinishRollbackEvent(this.parent, sol, obj, n, nSRCs, nVB);
        }
      }

    }

    public void fireGapReductionEvent(double gR){
      GapReductionEvent gapEvent = null;

      ExtendBAPListener listener;
      for(Iterator var3 = this.customListeners.iterator(); var3.hasNext(); listener.gapReduction(gapEvent)) {
        listener = (ExtendBAPListener)var3.next();
        if (gapEvent == null) {
          gapEvent = new GapReductionEvent(this.parent, gR);
        }
      }

    }

    public void fireIPSolutionEvent(){
      IPSolutionEvent IPEvent = null;

      ExtendBAPListener listener;
      for(Iterator var3 = this.customListeners.iterator(); var3.hasNext(); listener.IPRootNode(IPEvent)) {
        listener = (ExtendBAPListener)var3.next();
        if (IPEvent == null) {
          IPEvent = new IPSolutionEvent(this.parent);
        }
      }

    }

    public void fireFinishIPSolutionEvent(double time){
      FinishIPSolutionEvent IPEvent = null;

      ExtendBAPListener listener;
      for(Iterator var3 = this.customListeners.iterator(); var3.hasNext(); listener.finishIPRootNode(IPEvent)) {
        listener = (ExtendBAPListener)var3.next();
        if (IPEvent == null) {
          IPEvent = new FinishIPSolutionEvent(this.parent, time);
        }
      }

    }


}

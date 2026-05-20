package branchAndPrice;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecisionListener;
import org.slf4j.LoggerFactory;

import columnGeneration.Route;
import model.EVRPTW;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;
import org.slf4j.Logger;

public class customGraphManipulator{

    protected final Logger logger = LoggerFactory.getLogger(customGraphManipulator.class);
    private BAPNode<EVRPTW, Route> previous_node;
    private final Set<BranchingDecisionListener> listeners;

    public final Map<Integer, ArrayList<Integer>> rootPaths;
	public final Map<Integer, ArrayList<BranchingDecision>> branchingDecisions;

    public customGraphManipulator(BAPNode<EVRPTW, Route> rootNode, Map<Integer, ArrayList<Integer>> rootPaths, Map<Integer, ArrayList<BranchingDecision<EVRPTW,Route>>> branchingDecisions){
        this.previous_node = rootNode;
        this.listeners = new LinkedHashSet<BranchingDecisionListener>();
        this.rootPaths = rootPaths;
        this.branchingDecisions = branchingDecisions;
    }

    public void next(BAPNode<EVRPTW, Route> next_node){
        
        int mutualNodesOnPath = 0;
        ArrayList<Integer> previous_rootPath = rootPaths.get(previous_node.nodeID);
        ArrayList<Integer> next_rootPath = rootPaths.get(next_node.nodeID);

        for(int i = 0; i < Math.min(previous_rootPath.size(), next_rootPath.size()) && previous_rootPath.get(i) == next_rootPath.get(i); i++) {
            mutualNodesOnPath++;
        }

        // Removing the branching decisions of the non-mutual nodes
        for (int n = mutualNodesOnPath; n < previous_rootPath.size(); n++){
            int nodeID = previous_rootPath.get(n);
            for (BranchingDecision bd: branchingDecisions.get(nodeID)) this.rewindBranchingDecision(bd);
        }

        // Adding the remaining branching decisions following the new root path
        for (int n = mutualNodesOnPath; n < next_rootPath.size(); n++){
            int nodeID = next_rootPath.get(n);
            for (BranchingDecision bd: branchingDecisions.get(nodeID)) this.performBranchingDecision(bd);
        }

        this.previous_node = next_node;
    }

    public void restore() {

        for (int nodeID: rootPaths.get(previous_node.nodeID)){
            for (BranchingDecision bd: branchingDecisions.get(nodeID)) this.rewindBranchingDecision(bd);
        }
    }

    protected void addBranchingDecisionListener(BranchingDecisionListener listener) {
        this.listeners.add(listener);
    }

    protected void removeBranchingDecisionListener(BranchingDecisionListener listener) {
        this.listeners.remove(listener);
    }

    private void performBranchingDecision(BranchingDecision bd) {
        for(BranchingDecisionListener listener : this.listeners) {
            listener.branchingDecisionPerformed(bd);
        }

    }

    private void rewindBranchingDecision(BranchingDecision bd) {
        for(BranchingDecisionListener listener : this.listeners) {
            listener.branchingDecisionReversed(bd);
        }

    }
    
}

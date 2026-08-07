package columnGeneration;

import java.util.Comparator;
import java.util.HashSet;
import java.util.Arrays;
import java.util.BitSet;


/**
 * Class that represents a Label for the Labeling Algorithm
 */
public class BackwardLabel{

	public int vertex; 						//current vertex of the label
	public int dominanceVertex;				//vertex id for dominance checks
	public int index; 						//currentLabel index
	public int nextArc; 					//nextVertex (from which it stems this one)
	public int nextLabelIndex; 				//nextLabel index (from which it stems this one)
	public double reducedCost; 				//reduced cost
	public int remainingLoad; 				//remaining load
	public int remainingTime;				//remaining time
	public int[] remainingEnergy; 			//remaining energy
	public int chargingTime; 				//time required to charge
	public boolean[] unreachable; 			//customers that are not reachable by resource limitations
	public boolean[] ng_path; 				//customers that visit them would violate the ng-path cycling restrictions
	public boolean[] eta; 					//number of times modulo 2 that the label has visited customers in S (a triplet in a SRC)
	public HashSet<Integer> srcIndices; 	//SRC indices for which \eta = 1
	public BitSet feasible_Ts;
	public boolean frcHasCandidate;

	/** Creates a new Label.*/
	public BackwardLabel(int nextLabelIndex, double reducedCost, int remainingLoad, int remainingTime, int[] remainingEnergy, int chargingTime, boolean[] unreachable, boolean[] ng_path, boolean[] eta, HashSet<Integer> srcIndices) {
		this.nextLabelIndex = nextLabelIndex;
		this.reducedCost = reducedCost;
		this.remainingLoad = remainingLoad;
		this.remainingTime = remainingTime;
		this.remainingEnergy = remainingEnergy;
		this.chargingTime = chargingTime;
		this.unreachable = unreachable;
		this.ng_path = ng_path;
		this.eta = eta;
		this.srcIndices = srcIndices;
	}

	/** Obtains the string representation of a label. */
	@Override
	public String toString(){
		return "l("+vertex+"): r="+reducedCost+",q="+remainingLoad+",t="+remainingTime+", e="+Arrays.toString(remainingEnergy) + ", b="+chargingTime + ", unreach="+ Arrays.toString(unreachable) + ", ng=" + Arrays.toString(ng_path) + ", eta=" + Arrays.toString(eta);
	}

	public BackwardLabel clone(){
		
		BackwardLabel newLab = new BackwardLabel(this.nextLabelIndex, this.reducedCost, this.remainingLoad, this.remainingTime, this.remainingEnergy, this.chargingTime, this.unreachable, this.ng_path, this.eta, this.srcIndices);
		newLab.vertex = this.vertex;
		newLab.nextArc = this.nextArc;
		return newLab;
	}

	/** @return a negative integer, zero, or a positive integer as this object is less than, equal to, or greater than the specified object. */
	public static class SortLabels implements Comparator<BackwardLabel> {

		public int superDepotID;
		public int T_startID;

		public SortLabels(int superDepot, int T_start){
			this.superDepotID = superDepot;
			this.T_startID = T_start;
		}

		@Override
		public int compare(BackwardLabel L1, BackwardLabel L2) {

			// If the labels are at the Super Depot Node
			if (L1.dominanceVertex == superDepotID) {
				if (L1.chargingTime<L2.chargingTime) return -1;			// Lower charging time gets priority
				if (L1.chargingTime>L2.chargingTime) return 1;
				if (L1.remainingTime>L2.remainingTime) return -1;		// Later departure time gets priority
				if (L1.remainingTime<L2.remainingTime) return 1;
				if (L1.reducedCost<L2.reducedCost) return -1;		// Then lower reduced cost gets priority
				return 1;
			}

			// If the vertex is a charging time period
			if(L1.vertex>T_startID) {
				if (L1.chargingTime<L2.chargingTime) return -1;			// Lower charging time gets priority
				if (L1.chargingTime>L2.chargingTime) return 1;		
				if (L1.reducedCost<L2.reducedCost) return -1;		// Then lower reduced cost gets priority
				return 1;											
			}

			// For non-first customer vertices
			int gamma = L1.remainingEnergy.length-1;
			if (L1.remainingLoad>L2.remainingLoad) return -1;						// Higher remaining load capacity gets priority
			if (L1.remainingLoad<L2.remainingLoad) return 1;
			if (L1.remainingEnergy[gamma]>L2.remainingEnergy[gamma]) return -1;		// Then higher remaining energy capacity
			if (L1.remainingEnergy[gamma]<L2.remainingEnergy[gamma]) return 1;
			if (L1.remainingTime>L2.remainingTime) return -1;						// Then higher remaining time
			if (L1.remainingTime<L2.remainingTime) return 1;
			if (L1.reducedCost<L2.reducedCost) return -1;							// Then lower reduced cost
			if (L1.reducedCost>L2.reducedCost) return 1;
			return 0;																// Finally, is both labels have the same resource consumptions, no priority is set
		}
	}

}
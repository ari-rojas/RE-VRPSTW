package columnGeneration;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;


/**
 * Class that represents a Label for the Labeling Algorithm
 */
public class Label{

	public int vertex; 						//current vertex of the label
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

	/** Creates a new Label.*/
	public Label(int vertex, int nextArc, int nextLabelIndex, double reducedCost, int remainingLoad, int remainingTime, int[] remainingEnergy, int chargingTime, boolean[] unreachable, boolean[] ng_path, boolean[] eta, HashSet<Integer> srcIndices) {
		this.vertex = vertex;
		this.nextArc = nextArc;
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
		return "l("+vertex+"): r="+reducedCost+",q="+remainingLoad+",t="+remainingTime+", e="+Arrays.toString(remainingEnergy) + ", b="+chargingTime + ", ix="+index;
	}

	public Label clone(){
		
		return new Label(this.vertex, this.nextArc, this.nextLabelIndex, this.reducedCost, this.remainingLoad, this.remainingTime, this.remainingEnergy, this.chargingTime, this.unreachable, this.ng_path, this.eta, this.srcIndices);
	}

	/** @return a negative integer, zero, or a positive integer as this object is less than, equal to, or greater than the specified object. */
	public static class SortLabels implements Comparator<Label> {
		@Override
		public int compare(Label L1, Label L2) {

			// For charging vertices
			if(L1.vertex>L1.unreachable.length+1 && L2.vertex>L2.unreachable.length+1) {
				if(L1.chargingTime<L2.chargingTime) return -1;			// Less charging time gets priority
				else if(L1.chargingTime>L2.chargingTime) return 1;		
				else if(L1.reducedCost<L2.reducedCost) return -1;		// Then lower reduced cost gets priority
				else return 1;											
			}

			int gamma = L1.remainingEnergy.length-1;
			// For non-chargin vertices
			if(L1.remainingLoad>L2.remainingLoad) return -1;			// Higher remaining load capacity gets priority
			if(L1.remainingLoad<L2.remainingLoad) return 1;
			if(L1.remainingEnergy[gamma]>L2.remainingEnergy[gamma]) return -1;		// Then higher remaining energy capacity
			if(L1.remainingEnergy[gamma]<L2.remainingEnergy[gamma]) return 1;
			if(L1.remainingTime>L2.remainingTime) return -1;			// Then higher remaining time
			if(L1.remainingTime<L2.remainingTime) return 1;
			if(L1.reducedCost<L2.reducedCost) return -1;				// Then lower reduced cost
			if(L1.reducedCost>L2.reducedCost) return 1;
			return 0;													// Finally, is both labels have the same resource consumptions, no priority is set
		}
	}
}
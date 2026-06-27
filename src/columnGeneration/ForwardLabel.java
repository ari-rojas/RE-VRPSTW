package columnGeneration;

import java.util.Comparator;
import java.util.HashSet;
import java.util.Arrays;


/**
 * Class that represents a Label for the Labeling Algorithm
 */
public class ForwardLabel{

	public int vertex; 						//current vertex of the label
	public int index; 						//currentLabel index
	public int previousArc; 				//nextVertex (from which it stems this one)
	public int previousLabelIndex; 		    //nextLabel index (from which it stems this one)
	public double reducedCost; 				//reduced cost
    public double chargingBound;            //charging reduced cost contribution bound
	public int cumulativeLoad; 				//remaining load
	public int cumulativeTime;				//remaining time
    public int latestDeparture;             //latest departure from the depot
    public int travelTimes;                 //sum of travel and service times along the subpath
	public int[] remainingEnergy; 			//remaining energy
	public int chargingTime; 				//time required to charge
	public boolean[] unreachable; 			//customers that are not reachable by resource limitations
	public boolean[] ng_path; 				//customers that visit them would violate the ng-path cycling restrictions
	public boolean[] eta; 					//number of times modulo 2 that the label has visited customers in S (a triplet in a SRC)
	public HashSet<Integer> srcIndices; 	//SRC indices for which \eta = 1
	public boolean frcHasCandidate;

	/** Creates a new Label.*/
	public ForwardLabel(int nextLabelIndex, double reducedCost, double chargingBound, int cumulativeLoad, int cumulativeTime, int latDeparture, int travelTimes, int[] remainingEnergy, int chargingTime, boolean[] unreachable, boolean[] ng_path, boolean[] eta, HashSet<Integer> srcIndices) {
		this.previousLabelIndex = nextLabelIndex;
		this.reducedCost = reducedCost;
        this.chargingBound = chargingBound;
		this.cumulativeLoad = cumulativeLoad;
		this.cumulativeTime = cumulativeTime;
        this.latestDeparture = latDeparture;
        this.travelTimes = travelTimes;
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
		return "l("+vertex+"): r="+reducedCost+",q="+cumulativeLoad+",t="+cumulativeTime+", e="+Arrays.toString(remainingEnergy) + ", b="+chargingTime + ", unreach="+ Arrays.toString(unreachable) + ", ng=" + Arrays.toString(ng_path) + ", eta=" + Arrays.toString(eta);
	}

	public ForwardLabel clone(){
		
		ForwardLabel newLab = new ForwardLabel(this.previousLabelIndex, this.reducedCost, this.chargingBound, this.cumulativeLoad, this.cumulativeTime, this.latestDeparture, this.travelTimes, this.remainingEnergy, this.chargingTime, this.unreachable, this.ng_path, this.eta, this.srcIndices);
		newLab.vertex = this.vertex;
		newLab.previousArc = this.previousArc;
		return newLab;
	}

	/** @return a negative integer, zero, or a positive integer as this object is less than, equal to, or greater than the specified object. */
	public static class SortForwardLabels implements Comparator<ForwardLabel> {

		public SortForwardLabels(){
		}

		@Override
		public int compare(ForwardLabel L1, ForwardLabel L2) {

			// For non-first customer vertices
			int gamma = L1.remainingEnergy.length-1;
			if (L1.cumulativeLoad<L2.cumulativeLoad) return -1;						            // Higher remaining load capacity gets priority
			if (L1.cumulativeLoad>L2.cumulativeLoad) return 1;
			if (L1.remainingEnergy[gamma]>L2.remainingEnergy[gamma]) return -1;		            // Then higher remaining energy capacity
			if (L1.remainingEnergy[gamma]<L2.remainingEnergy[gamma]) return 1;
			if (L1.cumulativeTime<L2.cumulativeTime) return -1;						            // Then higher remaining time
			if (L1.cumulativeTime>L2.cumulativeTime) return 1;
			if (L1.reducedCost+L1.chargingBound<L2.reducedCost+L2.chargingBound) return -1;		// Then lower reduced cost
			if (L1.reducedCost+L1.chargingBound>L2.reducedCost+L2.chargingBound) return 1;
			return 0;																            // Finally, is both labels have the same resource consumptions, no priority is set
		}
	}

}
package columnGeneration;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import org.jorlib.frameworks.columnGeneration.colgenMain.AbstractColumn;
import model.EVRPTW;


/**
 * Implementation of a column (route) in the mE-VRSPTW
 */
public final class Route extends AbstractColumn<EVRPTW, PricingProblem> {

	public final HashMap<Integer, Integer> route; 	//number of times each customer is visited
	public final int[] routeSequence;				//route sequence
	public final int cost;							//cost of this column in the objective 
	public final int departureTime;					//departure time
	public final int energy;						//energy required by the route
	public final int load;							//load
	public double reducedCost;						//reduced cost (when priced)
	public ArrayList<Integer> arcs;					//used routing arcs
	public ArrayList<Integer> PParcs;				//used PP arcs

	//charging information
	public int lastChargingTime;
	public int chargingTime;

	public int BBnode=-1;						//Node in the BB tree in which it was priced

	/**
	 * Creates a new route (column).
	 * @param creator: description of who created the column (an algorithm)
	 * @param isArtificial: indicated whether the route is artificial (for initialization purposes)
	 */
	public Route(String creator, boolean isArtificial, HashMap<Integer, Integer> route, int[] routeSequence, PricingProblem pricingProblem, int cost, int departureTime, int energy, int load, double reducedCost, ArrayList<Integer> arcs, ArrayList<Integer> PParcs, int initialChargingTime, int chargingTime) {
		super(pricingProblem, isArtificial, creator);
		this.route=route;
		this.routeSequence = routeSequence;
		this.cost= cost;
		this.departureTime = departureTime;
		this.energy = energy;
		this.load = load;
		this.reducedCost = reducedCost;
		this.arcs = arcs;
		this.PParcs = PParcs;
		this.lastChargingTime = initialChargingTime;
		this.chargingTime = chargingTime;
	}

	/** The equals and hashCode methods are important (for the jORlib). **/
	@Override
	public boolean equals(Object o) {
		if(this==o)
			return true;
		if(!(o instanceof Route))
			return false;
		Route other=(Route) o;
		if((this.isArtificialColumn && !other.isArtificialColumn)|| (!this.isArtificialColumn && other.isArtificialColumn)) return false;
		return (this.PParcs.equals(other.PParcs) && this.lastChargingTime == other.lastChargingTime && this.chargingTime == other.chargingTime);
	}

	@Override
	public int hashCode() {
		return route.hashCode()+lastChargingTime+chargingTime;
	}

	/** Returns the string representation of a route. */
	@Override
	public String toString() {
		if (this.value>0)
			return "Value: "+ Math.floor(10000*this.value)/10000+" Cost: "+this.cost +" Energy: " + this.energy + " Route: "+ Arrays.toString(this.routeSequence)+" charging time: [" + (this.lastChargingTime-this.chargingTime+1) + "," + this.lastChargingTime +"]"+ " Departure: "  + this.departureTime +  " Arcs:" + arcs.toString() + " PP Arcs: " + PParcs.toString();
		else 
			return "Reduced Cost: "+ Math.floor(10000*this.reducedCost)/10000+" Cost: "+this.cost +" Energy: " + this.energy +" Route: "+ Arrays.toString(this.routeSequence)+ " charging time: [" + (this.lastChargingTime-this.chargingTime+1) + "," + this.lastChargingTime +"]"+ " Departure: "  + this.departureTime + " Arcs:" + arcs.toString() + " PP Arcs: " + PParcs.toString();
	}

	/** Clones the route. */
	public Route clone() {
		return new Route(this.creator, this.isArtificialColumn, (HashMap<Integer, Integer>) this.route.clone(), (int[]) this.routeSequence.clone(), this.associatedPricingProblem, this.cost, this.departureTime, this.energy, this.load, this.reducedCost, (ArrayList<Integer>) this.arcs.clone(), (ArrayList<Integer>) this.PParcs.clone(), this.lastChargingTime, this.chargingTime);
	}
}
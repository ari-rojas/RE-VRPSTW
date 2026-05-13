package columnGeneration;

import java.util.Comparator;

public class DepotLabel {
    
    public int index = -1;
    public final Label label;
    public final int depot_vertex_id;

    public DepotLabel(Label lab, int depot_vertex_id){

        this.label = lab;
        this.depot_vertex_id = depot_vertex_id;

    }

    @Override
	public String toString(){
		return "Depot Label at i0 "+ depot_vertex_id +"; Label: " + label.toString();
	}

    public static class SortDepotLabels implements Comparator<DepotLabel> {
		@Override
		public int compare(DepotLabel L1, DepotLabel L2) {
			
            Label lab1 = L1.label; Label lab2 = L2.label;

            if(lab1.chargingTime<lab2.chargingTime) return -1;			// Less charging time gets priority
            if(lab1.chargingTime>lab2.chargingTime) return 1;		
            if(lab1.remainingTime>lab2.remainingTime) return -1;        // Later departure time gets priority
            if(lab1.remainingTime<lab2.remainingTime) return 1;
            if(lab1.reducedCost<lab2.reducedCost) return -1;		// Then lower reduced cost gets priority
            
            return 1;
			
		}
	}

}

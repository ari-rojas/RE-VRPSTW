package branchAndPrice;

import java.util.EventObject;

public class CustomStartBAPEvent extends EventObject {
    public final String instanceName;
    public final double objectiveIncumbentSolution;

    public CustomStartBAPEvent(Object source, String instanceName, double objectiveIncumbentSolution) {
        super(source);
        this.instanceName = instanceName;
        this.objectiveIncumbentSolution = objectiveIncumbentSolution;
    }
}

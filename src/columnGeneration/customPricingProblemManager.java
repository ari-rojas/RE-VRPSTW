package columnGeneration;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;
import org.jorlib.frameworks.columnGeneration.pricing.PricingProblemBundle;
import org.jorlib.frameworks.columnGeneration.pricing.PricingProblemManager;

import model.EVRPTW;

public class customPricingProblemManager extends PricingProblemManager{

    public customPricingProblemManager(List<PricingProblem> pricingProblems, Map<Class<? extends AbstractPricingProblemSolver<EVRPTW, Route, PricingProblem>>, PricingProblemBundle<EVRPTW, Route, PricingProblem>> pricingProblemBundles){
        super(pricingProblems, pricingProblemBundles);
    }
    
}

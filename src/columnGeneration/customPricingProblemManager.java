package columnGeneration;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import org.jorlib.frameworks.columnGeneration.colgenMain.AbstractColumn;
import org.jorlib.frameworks.columnGeneration.io.TimeLimitExceededException;
import org.jorlib.frameworks.columnGeneration.util.Configuration;

import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblem;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;
import org.jorlib.frameworks.columnGeneration.pricing.PricingProblemBundle;

public class customPricingProblemManager<T, U extends AbstractColumn<T,V>, V extends AbstractPricingProblem<T>> {
    private static final Configuration config = Configuration.getConfiguration();
    private final Map<Class<? extends AbstractPricingProblemSolver<T, U, V>>, PricingProblemBundle<T, U, V>> pricingProblemBundles;
    private final ExecutorService executor;
    private final List<Future<Void>> futures;

    public customPricingProblemManager(List<V> pricingProblems, Map<Class<? extends AbstractPricingProblemSolver<T, U, V>>, PricingProblemBundle<T, U, V>> pricingProblemBundles) {
        this.pricingProblemBundles = pricingProblemBundles;
        this.executor = Executors.newFixedThreadPool(config.MAXTHREADS);
        this.futures = new ArrayList<Future<Void>>(pricingProblems.size());
    }

    public List<U> solvePricingProblems(Class<? extends AbstractPricingProblemSolver<T, U, V>> solver) throws TimeLimitExceededException {
        PricingProblemBundle<T, U, V> bundle = (PricingProblemBundle<T, U, V>)this.pricingProblemBundles.get(solver);

        for(AbstractPricingProblemSolver<T, U, V> solverInstance : bundle.solverInstances) {
            Future<Void> f = this.executor.submit(solverInstance);
            this.futures.add(f);
        }

        for(Future<Void> f : this.futures) {
            try {
                f.get();
            } catch (ExecutionException e) {
                if (e.getCause() instanceof TimeLimitExceededException) {
                this.shutdownAndAwaitTermination(this.executor);
                throw (TimeLimitExceededException)e.getCause();
                }

                e.printStackTrace();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        List<U> newColumns = new ArrayList<U>();

        for(AbstractPricingProblemSolver<T, U, V> solverInstance : bundle.solverInstances) {
            newColumns.addAll(solverInstance.getColumns());
        }

        return newColumns;
    }

    public void setTimeLimit(long timeLimit) {
        for(PricingProblemBundle<T, U, V> bunddle : this.pricingProblemBundles.values()) {
            for(AbstractPricingProblemSolver<T, U, V> solverInstance : bunddle.solverInstances) {
                solverInstance.setTimeLimit(timeLimit);
            }
        }

    }

    private void shutdownAndAwaitTermination(ExecutorService pool) {
        pool.shutdownNow();

        try {
            if (!pool.awaitTermination(60L, TimeUnit.SECONDS)) {
                System.err.println("Pool did not terminate");
            }
        } catch (InterruptedException var3) {
            pool.shutdownNow();
            Thread.currentThread().interrupt();
        }

    }

    public void close() {
        this.executor.shutdownNow();

        for(PricingProblemBundle<T, U, V> bunddle : this.pricingProblemBundles.values()) {
            for(AbstractPricingProblemSolver<T, U, V> solverInstance : bunddle.solverInstances) {
                solverInstance.close();
            }
        }

    }
}
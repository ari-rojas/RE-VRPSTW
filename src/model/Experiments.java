package model;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.jorlib.frameworks.columnGeneration.util.Configuration;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.io.File;
import java.io.PrintStream;
import java.lang.reflect.Field;

import model.EVRPTW.Arc;

import columnGeneration.Route;
import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

public class Experiments {

    public static void deleteStaticObject(Class<?> clazz, String fieldName) {
        try {
            Field field = clazz.getDeclaredField(fieldName);
            field.setAccessible(true); // Allow access to private or protected fields
            field.set(null, null);     // Set the static field to null
        } catch (NoSuchFieldException | IllegalAccessException e) {
            e.printStackTrace();
            throw new RuntimeException("Unable to delete the static object", e);
        }
    }

    public static List<String> generateNames(String instance_prefix) {

        HashMap<String, Integer> num = new HashMap<>();
        num.put(instance_prefix+"C1", 9); num.put(instance_prefix+"C2", 8); num.put(instance_prefix+"R1", 12);
        num.put(instance_prefix+"R2", 11); num.put(instance_prefix+"RC1", 8); num.put(instance_prefix+"RC2", 8);

        // Create a list to store all names
        List<String> names = new ArrayList<>();

        // Outer loop for numeric values
        for (int n : new int[]{50, 25}) {
            // Loop through instances
            for (String instance : num.keySet()) {
                // Get the max number from the map
                int maxNum = num.get(instance);
                // Loop to generate names
                for (int j = 1; j <= maxNum; j++) {
                    // Format the number to 2 digits
                    String number = String.format("%02d", j);
                    // Create the name
                    String name = instance + number + "-" + n;
                    // Add to the list
                    names.add(name);
                }
            }
        }
        // Return the list of names
        return names;
    }

    public static void determine_number_of_chargers(String name){

        try {
            PrintStream fileOut = new PrintStream("./results/log/B-"+name+".log");
            System.setOut(fileOut);

            if (!name.equals("")) {
                System.out.println(" ========================== "+name+" ========================== ");
                
                boolean same_obj = false;
                Double last_obj = 0.;
        
                int B = 1;
                while (!same_obj){
        
                    EVRPTW evrptw = new EVRPTW(name, 0, B, false, "RE-VRSPTW","tuning");
                    EVRPTWSolver Solver =  new EVRPTWSolver(evrptw, null);
        
                    Double obj = Solver.upperBound;
                    if (obj.doubleValue() == last_obj.doubleValue() && obj.doubleValue() < 100000) same_obj = true;
                    else {
                        B++; last_obj = obj;
                    }
        
                    deleteStaticObject(Configuration.class, "instance");
        
                }
            }

        } catch (Exception ex){
            ex.printStackTrace();
        }
        

    }

    public static void test_number_of_chargers(String name){

        try{
            PrintStream fileOut = new PrintStream("./results/log/KB-"+name+".log");
            System.setOut(fileOut);

            // RUN ALL THE EXPERIMENTS AT ONCE

            File xmlFile = new File("./results/tuning/Num_chargers.xml");
            DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
            DocumentBuilder builder = factory.newDocumentBuilder();
            Document doc = builder.parse(xmlFile);
            Element num_chargers_element = (Element) doc.getElementsByTagName("tuning").item(0);

            System.out.println(" ========================== "+name+" ========================== ");
            
            Element instance_element = (Element) num_chargers_element.getElementsByTagName(name).item(0);
            Element unb_B = (Element) instance_element.getElementsByTagName("unb_B").item(0);
            int B = Integer.parseInt(unb_B.getElementsByTagName("K").item(0).getTextContent());

            EVRPTW evrptw = new EVRPTW(name, 0, B, false, "RE-VRSPTW", "tuning");
            EVRPTWSolver Solver =  new EVRPTWSolver(evrptw, null);

            deleteStaticObject(Configuration.class, "instance");

        } catch (Exception e){
            e.printStackTrace();
        }


    }

    public static void tune_number_of_chargers(){

        try{
            PrintStream fileOut = new PrintStream("./results/log/1. Tune_num_chargers.log");
            System.setOut(fileOut);

            File xmlFile = new File("./data/1. Num_chargers.xml");
            DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
            DocumentBuilder builder = factory.newDocumentBuilder();
            Document doc = builder.parse(xmlFile);
            Element num_chargers_element = (Element) doc.getElementsByTagName("num_chargers").item(0);

            List<String> names = generateNames("");
            for (String name : names){

                try {
                    
                    Element instance_element = (Element) num_chargers_element.getElementsByTagName(name).item(0);
                    boolean needs_tuning = Boolean.parseBoolean(instance_element.getElementsByTagName("needs_tuning").item(0).getTextContent());

                    if (needs_tuning) {
                        System.out.println(" ========================== "+name+" ========================== ");
                        
                        int max_chargers = Integer.parseInt(instance_element.getElementsByTagName("num_chargers").item(0).getTextContent());
                        int min_chargers = Integer.parseInt(instance_element.getElementsByTagName("min_chargers").item(0).getTextContent());

                        Double last_obj = 10000.;
                        for (int B=max_chargers; B >= min_chargers; B--){
                            EVRPTW evrptw = new EVRPTW(name, 0, B, false, "RE-VRSPTW", "tuning");
                            EVRPTWSolver Solver =  new EVRPTWSolver(evrptw, null);

                            Double obj = Solver.upperBound;
                            if (obj.doubleValue() > last_obj.doubleValue()){
                                deleteStaticObject(Configuration.class, "instance");
                                break;
                            } else{
                                last_obj = obj;
                            }
                
                            deleteStaticObject(Configuration.class, "instance");
                        }
                    }

                } catch (Exception ex){
                    ex.printStackTrace();
                    break;
                }
            }

        } catch (Exception e){
            e.printStackTrace();
        }

    }

    public static void run_experiments(String instances_prefix, int gamma, String experiment){

        try{

            // RUN ALL THE EXPERIMENTS AT ONCE
            List<String> names = generateNames(instances_prefix);
            for (String name : names){

                try {
                    if (!name.equals("")) {
                        
                        EVRPTW evrptw = new EVRPTW(name, gamma, 0, true, "RE-VRSPTW", experiment);
                        EVRPTWSolver Solver =  new EVRPTWSolver(evrptw, null);
            
                        deleteStaticObject(Configuration.class, "instance");
                    }

                } catch (Exception ex){
                    ex.printStackTrace();
                    break;
                }
            }

        } catch (Exception e){
            e.printStackTrace();
        }

    }

    public static void run_robustness_experiments(String instance){

        String alg = "ExNewPricing-Long";
        if (instance != ""){
            try {

                int gamma = 0;
                if (instance == "R211-50") gamma = 1;

                while (gamma <= 10){
                    
                    EVRPTW evrptw = new EVRPTW(instance, gamma, 0, true, alg, "Gamma"+gamma);
                    EVRPTWSolver Solver = new EVRPTWSolver(evrptw, null);

                    Solver.solve(32400000L); evrptw.fileOut.close();
                    ArrayList<Route> solution = Solver.close();

                    double obj = Solver.upperBound;
                    boolean isOptimal = Solver.isOptimal;

                    int next_gamma = gamma + 1;
                    if (obj > 1e7 && isOptimal){ // the experiment and all the ones that follow it are infeasible

                        for (int g = gamma+1; g <= 10; g++) printSolutionToLogFile(alg, "Gamma"+g, instance, obj, solution, 0);
                        break;

                    } else if (gamma <= 10) {

                        // Retrieve the solution
                        int nR = solution.size(); int[] departureTimes = new int[nR]; int maxT = 0;

                        int[] nominalEnergy = new int[nR];
                        PriorityQueue<Integer>[] energyDeviations = new PriorityQueue[nR];
                        int[] worstCaseEnergy = new int[nR];
                        for (int r = 0; r < nR; r++){

                            Route route = solution.get(r);

                            int d = route.departureTime; departureTimes[r] = d;
                            if (d > maxT) maxT = d;
                            
                            int nomEnergy = 0;
                            PriorityQueue<Integer> energyDevs = new PriorityQueue<>(Comparator.reverseOrder());
                            for (Integer arcID: route.arcs){
                                Arc arc = evrptw.arcs[arcID];
                                nomEnergy += arc.energy;
                                energyDevs.add(arc.energy_deviation);
                            }

                            nominalEnergy[r] = nomEnergy;

                            int worstEnergy = nomEnergy;
                            for (int g = 1; g <= gamma; g++) {
                                Integer dev = energyDevs.poll();
                                if (dev != null) worstEnergy += dev;
                                else break;
                            }

                            energyDeviations[r] = energyDevs; worstCaseEnergy[r] = worstEnergy;

                        }

                        boolean isRobust = true; 
                        // Assess the robustness of the current solution
                        for (int g = gamma+1; g <= 10; g++){

                            long startTime = System.currentTimeMillis();
                            int[] chargingTimes = new int[nR];
                            
                            for (int r = 0; r < nR; r++){

                                int worstEnergy = worstCaseEnergy[r];
                                Integer nextWorseDev = energyDeviations[r].poll();
                                if (nextWorseDev != null) worstEnergy += nextWorseDev;

                                if (worstEnergy > evrptw.E) { isRobust = false; break; }
                                worstCaseEnergy[r] = worstEnergy;

                                chargingTimes[r] = evrptw.f_inverse[worstEnergy];

                            }

                            if (!isRobust) break;

                            // Solving the Charging Scheduling Model
                            SchedulingResult result = solveChargingScheduling(nR, maxT, chargingTimes, departureTimes, evrptw.B);

                            if (!result.feasible) break;

                            // Retrieving the robust solution
                            int[] startingTimes = result.startingTimes;
                            ArrayList<Route> robust_solution = new ArrayList<>();
                            for (int r = 0; r < nR; r++){

                                Route route = solution.get(r);
                                Route new_route = new Route("initSolution", false, (HashMap<Integer, Integer>) route.route.clone(), (int[]) route.routeSequence.clone(), route.associatedPricingProblem, route.cost, route.departureTime, worstCaseEnergy[r], route.load, route.reducedCost, (ArrayList<Integer>) route.arcs.clone(), startingTimes[r], chargingTimes[r]);
                                new_route.value = 1;

                                robust_solution.add(new_route);
                            }

                            long totalTime = System.currentTimeMillis()-startTime;
                            printSolutionToLogFile(alg, "Gamma"+g, instance, obj, robust_solution, totalTime);

                            next_gamma ++;

                        }

                    }

                    gamma = next_gamma;
                    deleteStaticObject(Configuration.class, "instance");

                }

            } catch (Exception ex) { ex.printStackTrace(); }
        
        }

    }

    public static void printSolutionToLogFile(String algorithm, String experiment, String instance, double obj, List<Route> solution, long time) {

        try {
            PrintStream fileOut = new PrintStream("./results/log/"+algorithm+"/"+experiment+"/"+instance+".log");
            System.setOut(fileOut);

            System.out.println("================ SOLUTION BPC - " + instance +" ================");
            System.out.println("BAP terminated with objective: "+obj);
            System.out.println("Total Number of iterations: "+0);
            System.out.println("Total Number of processed nodes: "+0);
            System.out.println("Total Time spent on master problems (s): "+0+" Total time spent on pricing problems (s): "+0);
            System.out.println("Total running time (s): "+ getTimeInSeconds(time));
            
            System.out.println("Solution is optimal: "+true);
            System.out.println("Columns (only non-zero columns are returned):");

            for (Route column : solution){
                System.out.println(column.toString());
            }

            System.out.println("TIME BRANCHING - Total time is: 0.0");

            fileOut.close();

        } catch (Exception e) { e.printStackTrace();}
	}

    public static SchedulingResult solveChargingScheduling(int n, int maxT, int[] charging_times, int[] departure_times, double B) throws IloException {

		try {
			IloCplex cplex = new IloCplex();
			cplex.setOut(null); 			//disable CPLEX output
			cplex.setParam(IloCplex.Param.Threads, 1);

			// x[t] = x_(t,t+1), t = 0..maxT
			IloNumVar[] x = new IloNumVar[maxT + 1];
			for (int t = 0; t <= maxT; t++) {
				x[t] = cplex.numVar(0.0, Double.MAX_VALUE, "x_" + t + "_" + (t + 1));
			}

			// y[r][t] = y_t^r, r = 0..n-1, b_r <= t < d_r
			IloIntVar[][] y = new IloIntVar[n][maxT + 1];
			for (int r = 0; r < n; r++) {
				for (int t = charging_times[r]; t < departure_times[r]; t++) {
					y[r][t] = cplex.boolVar("y_" + t + "_" + r);
				}
			}

			// -----------------------------
			// Boundary constraints
			// x_(0,1) = 0
			// x_(T,T+1) = 0, with T = maxT
			// -----------------------------
			cplex.addEq(x[0], 0.0, "x_0_1_zero");
			cplex.addEq(x[maxT], 0.0, "x_" + maxT + "_" + (maxT + 1) + "_zero");

			// -----------------------------
			// 1) Charging capacity
			// x_(t,t+1) + sum_{r : b_r <= t < d_r} y_t^r <= B     for all t = 1..maxT
			// -----------------------------
			for (int t = 1; t <= maxT; t++) {
				IloLinearNumExpr lhs = cplex.linearNumExpr();

				lhs.addTerm(1.0, x[t]);
				for (int r = 0; r < n; r++) {
					if (charging_times[r] <= t && t < departure_times[r]) lhs.addTerm(1.0, y[r][t]);
				}

				cplex.addLe(lhs, B, "capacity_" + t);
			}

			// -----------------------------
			// 2) Flow conservation
			// x_(t,t+1) + sum_{r : b_r <= t < d_r} y_t^r
			//   = x_(t-1,t) + sum_{r : t+b_r-1 < d_r} y_{t+b_r-1}^r
			// for all t = 1..maxT
			// -----------------------------
			for (int t = 1; t <= maxT; t++) {
				IloLinearNumExpr lhs = cplex.linearNumExpr();
				IloLinearNumExpr rhs = cplex.linearNumExpr();

				// Left-hand side: x_(t,t+1) + sum_{r : b_r <= t < d_r} y_t^r
				lhs.addTerm(1.0, x[t]);
				for (int r = 0; r < n; r++) {
					if (charging_times[r] <= t && t < departure_times[r])  lhs.addTerm(1.0, y[r][t]);
				}

				// Right-hand side: x_(t-1,t) + sum_{r : t+b_r-1 < d_r} y_{t+b_r-1}^r
				rhs.addTerm(1.0, x[t - 1]);
				for (int r = 0; r < n; r++) {
					int tau = t + charging_times[r] - 1; // tau = t + b_r - 1
					if (tau < departure_times[r]) rhs.addTerm(1.0, y[r][tau]);
				}

				cplex.addEq(lhs, rhs, "flow_" + t);
			}

			// -----------------------------
			// 3) Each EV charges once
			// sum_{t : b_r <= t < d_r} y_t^r = 1     for all r = 0..n-1
			// -----------------------------
			for (int r = 0; r < n; r++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();

				for (int t = charging_times[r]; t < departure_times[r]; t++) expr.addTerm(1.0, y[r][t]);
				cplex.addEq(expr, 1.0, "charge_once_" + r);
			}

			// Solve the model, and if a feasible solution exists, retrieve it
			if (cplex.solve()) {
                // Extract solution
                int[] startingTimes = new int[n];

                for (int r = 0; r < n; r++) {
                    int init_t = 0;
                    for (int t=departure_times[r]-1; t>=charging_times[r]; t--){
						if (cplex.getValue(y[r][t]) > 0.5) {init_t = t-charging_times[r]+1; break;}
					}
                    startingTimes[r] = init_t;
                }

                cplex.end();

                return new SchedulingResult(true, startingTimes);
				
			} else { cplex.end(); return new SchedulingResult(false, null); }


        } catch (IloException e) { e.printStackTrace(); return new SchedulingResult(false, null); }
    }

    /** Returns the real (double) objective (divided by 10). */
	public static double getScaledObjective(double objective) {
		double realCost = objective;
		realCost = realCost*0.1+0.05;
		return Math.floor(realCost*10)/10;
	}

	/** Returns the time in seconds (and considering two decimals). */
	public static double getTimeInSeconds(double time) {
		double realTime = time*0.001;
		realTime = Math.floor(realTime*100)/100; //two decimals
		return realTime;
	}

    public static void run_experiments(int gamma, String experiment){

        run_experiments("", gamma, experiment);

    }

    public static void main(String[] args){

        run_robustness_experiments(args[0]);
    
    }

    public static class SchedulingResult {

        public final boolean feasible;
        // Optional: store solutions (only if feasible)
        public final int[] startingTimes;

        public SchedulingResult(boolean feasible, int[] starts) {
            this.feasible = feasible;
            this.startingTimes = starts;
        }
    }
}

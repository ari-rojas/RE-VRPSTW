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

    public static void determine_number_of_chargers(String instance){

        try {
            PrintStream fileOut = new PrintStream("./results/log/RojasTuneB/" + instance + ".log");
            System.setOut(fileOut);

            if (!instance.equals("")) {
                System.out.println(" ========================== " + instance + " ========================== ");
                
                File xmlFile = new File("./data/" + instance + ".xml");
                DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
                DocumentBuilder builder = factory.newDocumentBuilder();
                Document doc = builder.parse(xmlFile);

                Element infoElement = (Element) doc.getElementsByTagName("info").item(0);
                int B = Integer.parseInt(infoElement.getElementsByTagName("num_chargers").item(0).getTextContent());
        
                boolean feasible = false; long timeLimit = 86400000L; double last_obj = Double.MAX_VALUE;
                while (!feasible){
                    
                    long time = System.currentTimeMillis();
                    EVRPTW evrptw = new EVRPTW(instance, 0, B, false, "RojasTuneB", "","");
                    EVRPTWSolver Solver =  new EVRPTWSolver(evrptw, new ArrayList<>());

                    Solver.solve(timeLimit); ArrayList<Route> solution = Solver.close();
        
                    Double obj = Solver.upperBound;
                    if (obj.doubleValue() < 1e7 && System.currentTimeMillis()-time < timeLimit && Math.abs(obj-last_obj) < 1e-4) feasible = true;
                    else B++;
        
                    deleteStaticObject(Configuration.class, "instance");
                    last_obj = obj+0;
        
                }
            }

        } catch (Exception ex){
            ex.printStackTrace();
        }

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

    public static void main(String[] args){

        determine_number_of_chargers(args[0]);
    
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
package model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.jorlib.frameworks.columnGeneration.util.Configuration;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.io.File;
import java.io.PrintStream;
import java.lang.reflect.Field;

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
        
                    EVRPTW evrptw = new EVRPTW(name, 0, B, false, "tuning");
                    EVRPTWSolver Solver =  new EVRPTWSolver(evrptw);
        
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

            EVRPTW evrptw = new EVRPTW(name, 0, B, false, "tuning");
            EVRPTWSolver Solver =  new EVRPTWSolver(evrptw);

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
                            EVRPTW evrptw = new EVRPTW(name, 0, B, false, "tuning");
                            EVRPTWSolver Solver =  new EVRPTWSolver(evrptw);

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

    public static void run_experiments(String instances_prefix, int gamma){

        try{

            // RUN ALL THE EXPERIMENTS AT ONCE
            List<String> names = generateNames(instances_prefix);
            for (String name : names){

                try {
                    if (!name.equals("")) {
                        
                        EVRPTW evrptw = new EVRPTW(name, gamma, 0, true, "tuning");
                        EVRPTWSolver Solver =  new EVRPTWSolver(evrptw);
            
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

    public static void run_experiments(int gamma){

        run_experiments("", gamma);

    }

    public static void main(String[] args){

        run_experiments(args[0], Integer.parseInt(args[1]));
    
    }
}

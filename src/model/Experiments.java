package model;

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

    public static void determine_number_of_chargers(String instance){

        try {
            PrintStream fileOut = new PrintStream("./results/log/YaminTuneB/" + instance + ".log");
            System.setOut(fileOut);

            if (!instance.equals("")) {
                System.out.println(" ========================== " + instance + " ========================== ");
                
                File xmlFile = new File("./data/New-Yamin24-Tight/" + instance + ".xml");
                DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
                DocumentBuilder builder = factory.newDocumentBuilder();
                Document doc = builder.parse(xmlFile);

                Element infoElement = (Element) doc.getElementsByTagName("info").item(0);
                int B = Integer.parseInt(infoElement.getElementsByTagName("num_chargers").item(0).getTextContent());
        
                boolean feasible = false; long timeLimit = 86400000L;
                while (!feasible){
                    
                    long time = System.currentTimeMillis();
                    EVRPTW evrptw = new EVRPTW(instance, 0, B, false, "New-Yamin24-Tight", "","");
                    EVRPTWSolver Solver =  new EVRPTWSolver(evrptw);

                    Solver.solve(timeLimit);
        
                    Double obj = Solver.upperBound;
                    if (obj.doubleValue() < 1e7 && System.currentTimeMillis()-time < timeLimit) feasible = true;
                    else B++;
        
                    deleteStaticObject(Configuration.class, "instance");
        
                }
            }

        } catch (Exception ex){
            ex.printStackTrace();
        }

    }

    public static void main(String[] args){

        determine_number_of_chargers(args[0]);
    
    }
}

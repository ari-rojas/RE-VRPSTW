package model;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.PriorityQueue;
import java.util.Properties;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.DirectedWeightedMultigraph;
import org.jorlib.frameworks.columnGeneration.model.ModelInterface;
import org.jorlib.frameworks.columnGeneration.util.Configuration;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import columnGeneration.Label;

/**
 * The Electric Vehicle Routing and Overnight Charging Scheduling Problem on a Multigraph
 * @author Daniel Yam�n (Universidad de los Andes)
 */

public final class EVRPTW implements ModelInterface {

	public final String instanceName;						//instance name
	public PrintStream fileOut;
	public final String en_dev;

	public Configuration config;

	//Basic information
	public int C; 											//number of customers
	public int Q; 											//load capacity
	public int V; 											//number of vertices (customer-based graph)
	
	public DirectedWeightedMultigraph<Integer, Arc> graph; 	//routing graph
	public Vertex[] vertices; 								//set of routing vertices
	public Arc[] arcs; 										//set of routing arcs
	public int numArcsRoadNetwork; 							//number of arcs in the road network
	
	public DirectedWeightedMultigraph<Integer, PPArc> PPgraph;	//pricing problem graph
	public PPVertex[] PPvertices;
	public PPArc[] PParcs;
	public int numArcs; 										//number of arcs

	public int gamma;										//uncertainty budget by vehicle (route)

	//Energy information
	public int B; 											//number of chargers
	public int E; 											//energy capacity
	public int T_min; 										//opening time of the depot
	public int last_charging_period; 						//charging end time
	public int[] f_inverse; 								//(inverse) recharging function

	//Acceleration strategies and BPC
	public final int Delta;
	public final int DeltaMax;
	public final double precision = 0.09; 					//precision for the column generation algorithm (it is scaled by 10)
	public long exactPricingTime = 0; 						//time spent on the exact labeling algorithm
	public long heuristicPricingTime = 0; 					//time spent on the heuristic labeling algorithm
	public int columnsRootNode = 0; 						//columns generated at the root node
	public int cutsRootNode = 0; 							//cuts separated at the root node
	public int[] infeasiblePPArcs; 							//infeasible arcs in the pricing problem

	//Log file creation
	public boolean print_log;
	public String algorithm;
	public String experiment;

	public boolean CUTSENABLED;
	public boolean rollbackTrigger;
	public int rollbackBaseLine;
	public int rollbackExplosion;
	public int rollbackFactor = 125;
	public int cut_iterations;
	public int pricingSoftFactor = 5;
	public double gapReductionRequirement = 0.35;

	// Identifiers for the differnt types of vertices in the Pricing Problem Graph
	public static final byte C0 = 0; 	 	// Customer depot nodes
	public static final byte C1 = 1; 		// Non-first customer nodes
	public static final byte Depot = 2;		// Returning depot node
	public static final byte Tt = 3; 		// Charging time period nodes
	public static final byte Source = 4;	// Dummy source node

	public int C0_startID;
	public int C1_startID;
	public int T_startID;
	public int superDepotID;

	public static final String[] VERTEX_TYPE_NAMES = {"i0","i1","0","t","s"};

	// Identifiers for the different types of arcs in the Pricing Problem Graph
	public static final byte AR0 = 0;	// Routing arcs between customer depot nodes and non-first customer nodes
	public static final byte AR1 = 1;	// Routing arcs between non-first customer nodes
	public static final byte AC1 = 2;	// Charging Scheduling arcs to finishing charging times
	public static final byte AC2 = 3;	// Charging Scheduling arcs between consecutive charging time periods
	public static final byte AC3 = 4;	// Charging Scheduling arcs to starting charging times

	public int lenAR0;
	public int lenAR1;

	/**
	 * Constructs a new mE-VRSPTW instance. 
	 * @param instanceName input instance.
	 * @throws IOException Throws IO exception when the instance cannot be found.
	 */
	public EVRPTW(String instanceName, int gamma, int num_chargers, boolean print_log, String algorithm, String experiment, String en_dev) throws IOException {
		
		//Properties
		Properties properties = new Properties();
		properties.setProperty("MAXTHREADS", "1"); //only one thread
		Configuration.readFromFile(properties);
		this.config = Configuration.getConfiguration();

		this.en_dev = en_dev;
		
		this.instanceName = instanceName.trim();
		int start_ix = 0; int end_ix = 2;
		this.Delta = (this.instanceName.substring(start_ix, end_ix).equals("R1") || this.instanceName.substring(start_ix, end_ix).equals("C1") || this.instanceName.substring(start_ix, end_ix+1).equals("RC1")) ? 7 : 12;
		this.DeltaMax = Delta+5;
		this.C = Integer.parseInt(this.instanceName.substring(Math.max(this.instanceName.length() - 2, 0))); // Number of customers must be the last two
		this.V = C+2;
		this.graph = new DirectedWeightedMultigraph<Integer, EVRPTW.Arc>(Arc.class);
		this.PPgraph = new DirectedWeightedMultigraph<Integer, EVRPTW.PPArc>(PPArc.class);
		this.numArcs = 0;
		this.gamma = gamma;
		this.algorithm = algorithm;
		this.experiment = experiment;
		this.print_log = print_log;
		this.B = num_chargers;

		//create a new file output stream.
		//create a new file output stream.
		if (this.print_log){
			this.fileOut = new PrintStream("./results/log/"+this.algorithm+"/"+this.experiment+"/"+this.getName()+".log");
			System.setOut(this.fileOut);
		}

		//read the instance
		readData();
		if (this.print_log){
			System.out.println(" - Number of chargers: " + this.B);
			System.out.println(" - Energy capacity: " + this.E);
			System.out.println(" - Full recharging time: " + this.f_inverse[this.E]);
			System.out.println(" - Charging time periods: "+this.last_charging_period);
			System.out.println(" - Number of Routing Graph arcs: "+this.numArcsRoadNetwork);
			System.out.println(" - Number of arcs in AR0, AR1: "+this.lenAR0+", "+this.lenAR1);
			System.out.println(" - Number of PP arcs: "+this.numArcs);
		}

		/* for (int arcID = 0; arcID < this.numArcs; arcID++){
			System.out.println(PParcs[arcID].toString());
		} */

	}

	/** Name of the current instance */
	@Override
	public String getName() {
		return instanceName;
	}

	/**
	 * Reads an instance.
	 * The file must be stored in ./data/instances and following the guidelines of the VRPREP. 
	 */
	private void readData() {
		try {

			/** Reading input file **/
			if (this.print_log){
				System.out.println(" - ================ LOADING INSTANCE ================");
				System.out.println(" - Instance: " + this.getName());
			}
			
			File xmlFile = new File("./data/" + this.getName() + ".xml");
			
			DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
			DocumentBuilder builder = factory.newDocumentBuilder();
			Document doc = builder.parse(xmlFile);

			// Load initial information
			loadInitialInformation(doc);
			// Retrieve vehicle profile information
			Element vehicleProfileElement = (Element) doc.getElementsByTagName("vehicle_profile").item(0);
			Element customElements = (Element) vehicleProfileElement.getElementsByTagName("custom").item(0);
			this.last_charging_period = Integer.parseInt(customElements.getElementsByTagName("last_charging_period").item(0).getTextContent())+1;
			
			// Build the routing original graph
			build_routing_graph(doc);

			// Build the Pricing Problem graph
			build_pricing_problem_graph();
			
		}
		catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	/** Defines some parameters (according to the information in the .xml file). */
	public void loadInitialInformation(Document doc) {
		Element infoElement = (Element) doc.getElementsByTagName("info").item(0);
		//number of chargers, number of arcs, two alternatives pairs, and average number of alternatives
		if (this.print_log) this.B = Integer.parseInt(infoElement.getElementsByTagName("num_chargers").item(0).getTextContent());
	}

	public void build_pricing_problem_graph(){
		
		/////////////////////////////////////
		/// VERTICES
		/////////////////////////////////////
		
		this.PPvertices = new PPVertex[this.C*2+3+this.last_charging_period];
		int id = 0;

		// Dummy source node
		this.PPvertices[id] = new PPVertex(id, Source, 0); PPgraph.addVertex(id); id ++;

		// Routing customer depot nodes
		this.C0_startID = id-1;
		for (int i = 1; i <= this.C; i++){
			this.PPvertices[id] = new PPVertex(id, C0, this.vertices[i], i); PPgraph.addVertex(id); id ++; }
		// Non-first customer nodes
		this.C1_startID = id-1;
		for (int i = 1; i <= this.C; i++){
			this.PPvertices[id] = new PPVertex(id, C1, this.vertices[i], i); PPgraph.addVertex(id); id ++; }
		// Returning depot node
		this.PPvertices[id] = new PPVertex(id, Depot, this.vertices[this.C+1], 0); PPgraph.addVertex(id); id ++;

		// Charging time vertices
		this.T_startID = id-1;
		for (int t = 1; t <= last_charging_period; t++) {
			this.PPvertices[id] = new PPVertex(id, Tt, t); PPgraph.addVertex(id); id ++; }

		// Dummy Super Depot node to save the Labels of all C0 vertices
		this.PPvertices[id] = new PPVertex(id, C0, this.vertices[0], 0); this.superDepotID = id;
		
		int auxNumArcs = 2*(V*V-V);
		for (int ix = 0; ix <= superDepotID; ix ++) this.PPvertices[ix].unprocessedLabels = new PriorityQueue<Label>(auxNumArcs, new Label.SortLabels(superDepotID, T_startID));
		
		/////////////////////////////////////
		/// ARCS
		////////////////////////////////////
		
		this.PParcs = new PPArc[this.C*(this.C+this.last_charging_period+2)+2*this.last_charging_period];
		id = 0; this.lenAR0 = 0; this.lenAR1 = 0;

		// AR0 Routing arcs from the customer depot nodes to non-first customer nodes
		for (int tail = 1; tail <= this.C; tail++){
			int tail_vertex_id = this.C0_startID + tail;
			
			for (Arc routing_arc: graph.outgoingEdgesOf(tail)){
				int head = routing_arc.head;
				if (head <= this.C && this.graph.getEdge(0,tail).energy + routing_arc.energy + this.graph.getEdge(head,this.C+1).min_energy <= this.E && this.vertices[tail].opening_tw + routing_arc.time <= this.vertices[head].closing_tw){
					int head_vertex_id = this.C1_startID+head;
					PPArc newArc = new PPArc(id, AR0, this.graph.getEdge(tail, head), tail_vertex_id, head_vertex_id); this.PParcs[id] = newArc;
					PPgraph.addEdge(tail_vertex_id, head_vertex_id, newArc); id ++; this.lenAR0++;
				}
			}

			// Routing arc (\in AR2) from the customer depot node i0 to the returning depot
			PPArc newArc = new PPArc(id, AR0, this.graph.getEdge(tail, this.C+1), tail_vertex_id, this.T_startID); this.PParcs[id] = newArc;
			PPgraph.addEdge(tail_vertex_id, this.T_startID, newArc); id ++; this.lenAR0++;
		}

		// AR1 Routing arcs between non-first customer nodes
		for (int tail = 1; tail <= this.C; tail++){
			int tail_vertex_id = this.C1_startID+tail;
			PPVertex vx = PPvertices[tail_vertex_id];

			if (vx.routing_vertex.feasible_nonfirst){
			
				for (Arc routing_arc: graph.outgoingEdgesOf(tail)){
					int head = routing_arc.head;
					if (head <= this.C && this.graph.getEdge(0,tail).min_energy + routing_arc.energy + this.graph.getEdge(head,this.C+1).min_energy <= this.E && this.vertices[tail].open_tw_nonfirst + routing_arc.time <= this.vertices[head].closing_tw){
						int head_vertex_id = this.C1_startID+head;
						PPArc newArc = new PPArc(id, AR1, this.graph.getEdge(tail, head), tail_vertex_id, head_vertex_id); this.PParcs[id] = newArc;
						PPgraph.addEdge(tail_vertex_id, head_vertex_id, newArc); id ++; this.lenAR1++;
					}
				}

				// Routing arc (\in AR2) from the customer depot node i0 to the returning depot
				PPArc newArc = new PPArc(id, AR1, this.graph.getEdge(tail, this.C+1), tail_vertex_id, this.T_startID); this.PParcs[id] = newArc;
				PPgraph.addEdge(tail_vertex_id, this.T_startID, newArc); id ++; this.lenAR1++;
			}

		}

		// AC1 Charging arcs for finishing charging time periods
		for (int head = 1; head <= this.C; head++){
			int head_vertex_id = this.C0_startID+head;
			Vertex vx = vertices[head];

			for (int t = vx.min_chargingTime; t < vx.last_departure; t++){
				int tail_vertex_id = this.T_startID + t;
				PPArc newArc = new PPArc(id, AC1, tail_vertex_id, head_vertex_id); this.PParcs[id] = newArc;
				PPgraph.addEdge(tail_vertex_id, head_vertex_id, newArc); id ++;
			}
		}
		
		for (int t = 1; t < this.last_charging_period; t++){
			int tail_vertex_id = this.T_startID+t;

			// AC2 Charging arcs between consecutive charging time periods
			PPArc newArc = new PPArc(id, AC2, tail_vertex_id, tail_vertex_id+1); this.PParcs[id] = newArc;
			PPgraph.addEdge(tail_vertex_id, tail_vertex_id+1, newArc); id ++;

			// AC3 Charging arcs for starting charging time periods
			newArc = new PPArc(id, AC3, 0, tail_vertex_id); this.PParcs[id] = newArc;
			PPgraph.addEdge(0, tail_vertex_id, newArc); id ++;

		}
		PPArc newArc = new PPArc(id, AC3, 0, this.T_startID+this.last_charging_period); this.PParcs[id] = newArc;
		PPgraph.addEdge(0, this.T_startID+this.last_charging_period, newArc);

		this.numArcs = id+1;

	}

	public void build_routing_graph(Document doc){
			
		// Load vertices
		this.vertices = new Vertex[this.V];
		loadVertices(doc);
		this.T_min = vertices[0].opening_tw;

		//Load arcs
		int auxNumArcs = 2*(V*V-V)+3*last_charging_period;
		arcs = new Arc[auxNumArcs];
		loadArcs(doc);

		// Load fleet (vehicle profile and charging information)
		loadFleet(doc);

		/** Neighborhoods (ng-path). **/
		for (int i = 1; i <= this.C; i++) {

			//Sort arcs
			ArrayList<Arc> incoming = new ArrayList<Arc>(graph.incomingEdgesOf(i));
			Collections.sort(incoming, new SortByCost());
			ArrayList<Arc> outgoing = new ArrayList<Arc>(graph.outgoingEdgesOf(i));
			Collections.sort(outgoing, new SortByCost());

			int outgoingIndex = 0;
			int incomingIndex = 0;
			boolean processOutgoing = true;
			while(true) {
				if((outgoingIndex>=outgoing.size() && incomingIndex>=incoming.size()) || vertices[i].neighbors.size()>=this.Delta) break;
				else if(outgoingIndex>=outgoing.size()) processOutgoing = false;
				else if(incomingIndex>=incoming.size()) processOutgoing = true;
				else if (outgoing.get(outgoingIndex).cost<incoming.get(incomingIndex).cost) {
					processOutgoing=true;
				}else {processOutgoing = false;}
				int nextCustomer = (processOutgoing) ? outgoing.get(outgoingIndex).head: incoming.get(incomingIndex).tail;
				if (nextCustomer!=0 && nextCustomer!=C+1 && !vertices[i].neighbors.contains(nextCustomer)) {
					vertices[i].neighbors.add(nextCustomer);
				}
				if(processOutgoing) outgoingIndex++;
				else incomingIndex++;
			}
			vertices[i].neighbors.add(i); //does not count for the Delta
		}

		/** (A priori) unreachable customers **/
		for (int i = 1; i <= this.C; i++) {
			for (int j = 1; j <= this.C; j++) {
				if(i!=j && !graph.containsEdge(i, j))
					vertices[j].unreachable.add(i); //backward labeling
			}
		}

	}

	/** Defines the vertices (according to the information in the .xml file) */
	public void loadVertices(Document doc) {

		NodeList nodeNodes = doc.getElementsByTagName("node");
		for (int i = 0; i < nodeNodes.getLength(); i++) {

			Node node = nodeNodes.item(i);
			Element nodeElement = (Element) node;

			//id, coordx, coordy, and load
			int id = Integer.parseInt(nodeElement.getAttribute("id"));
			int node_type = Integer.parseInt(nodeElement.getAttribute("type"));
			int coordx = Integer.parseInt(nodeElement.getElementsByTagName("cx").item(0).getTextContent());
			int coordy = Integer.parseInt(nodeElement.getElementsByTagName("cy").item(0).getTextContent());
			int load = Integer.parseInt(nodeElement.getElementsByTagName("load").item(0).getTextContent());

			//time windows
			Element timeWindowsElement = (Element) nodeElement.getElementsByTagName("tw").item(0);
			int opening_tw = Integer.parseInt(timeWindowsElement.getElementsByTagName("start").item(0).getTextContent());
			int closing_tw = Integer.parseInt(timeWindowsElement.getElementsByTagName("end").item(0).getTextContent());

			if (node_type > 0){ // Customer nodes
				Element customElement = (Element) nodeElement.getElementsByTagName("custom").item(0);
				int last_departure = Integer.parseInt(customElement.getElementsByTagName("last_departure").item(0).getTextContent());
				boolean feasible_nonfirst = customElement.getElementsByTagName("non_first_feasible").item(0).getTextContent().equals("true");
				int open_tw_nonfirst = Integer.parseInt(customElement.getElementsByTagName("tw1_start").item(0).getTextContent());
				int min_b = Integer.parseInt(customElement.getElementsByTagName("min_chargingTime").item(0).getTextContent());
				vertices[id] = new Vertex(id, coordx, coordy, load, opening_tw, closing_tw, last_departure, feasible_nonfirst, open_tw_nonfirst, min_b);
			} else  vertices[id] = new Vertex(id, coordx, coordy, load, opening_tw, closing_tw); // Depot nodes

			graph.addVertex(id);
		}
	}

	/** Defines the arcs (according to the information in the .xml file). */
	public void loadArcs(Document doc) {

		NodeList linkNodes = doc.getElementsByTagName("link");
		this.numArcsRoadNetwork = linkNodes.getLength();

		for (int i = 0; i < linkNodes.getLength(); i++) {
			
			Node link = linkNodes.item(i);
			Element linkElement = (Element) link;

			//id, head and tail
			int id = Integer.parseInt(linkElement.getAttribute("id"));
			int tail = Integer.parseInt(linkElement.getAttribute("tail"));
			int head = Integer.parseInt(linkElement.getAttribute("head"));

			//travel cost and travel time
			int cost = Integer.parseInt(linkElement.getElementsByTagName("travel_cost").item(0).getTextContent());
			int time = Integer.parseInt(linkElement.getElementsByTagName("travel_time").item(0).getTextContent());

			//custom elements (energy and minimum values)
			Element customElements = (Element) linkElement.getElementsByTagName("custom").item(0);
			int energy = Integer.parseInt(customElements.getElementsByTagName("energy_consumption").item(0).getTextContent());
			int energy_deviation = Integer.parseInt(customElements.getElementsByTagName("energy_deviation"+this.en_dev).item(0).getTextContent());
			int min_energy = Integer.parseInt(customElements.getElementsByTagName("min_energy").item(0).getTextContent());
			int min_cost = Integer.parseInt(customElements.getElementsByTagName("min_cost").item(0).getTextContent());
			int min_time = Integer.parseInt(customElements.getElementsByTagName("min_time").item(0).getTextContent());
			boolean minCostAlternative = true;
			
			Arc newArc = new Arc(id, tail, head, cost, time, energy, energy_deviation, min_energy, min_cost, min_time, minCostAlternative);
			arcs[id] = newArc;
			graph.addEdge(tail, head, newArc);
			
		}
	}

	/** Defines the arcs (according to the information in the .xml file). */
	public void loadFleet(Document doc) {

		Element vehicleProfileElement = (Element) doc.getElementsByTagName("vehicle_profile").item(0);
		//vehicles' load capacity
		this.Q = Integer.parseInt(vehicleProfileElement.getElementsByTagName("capacity").item(0).getTextContent());

		//custom elements (last charging time and energy capacity)
		Element customElements = (Element) vehicleProfileElement.getElementsByTagName("custom").item(0);
		this.E = Integer.parseInt(customElements.getElementsByTagName("energy_capacity").item(0).getTextContent());

		//inverse recharging function
		this.f_inverse= new int[E+1]; //piecewise-linear function (rounded up)
		Element inverseChargingFunctionElement = (Element) customElements.getElementsByTagName("inverse_recharging_function").item(0);
		NodeList breakpointNodes = inverseChargingFunctionElement.getElementsByTagName("breakpoint");
		for (int i = 0; i < breakpointNodes.getLength(); i++) {
			Node breakpoint = breakpointNodes.item(i);
			Element breakpointElement = (Element) breakpoint;
			int energyLevel = Integer.parseInt(breakpointElement.getElementsByTagName("energy_level").item(0).getTextContent());
			int timesteps = Integer.parseInt(breakpointElement.getElementsByTagName("periods").item(0).getTextContent());
			f_inverse[energyLevel] = timesteps;
		}


	}

	/** Cleans the SRCs Indices of each customer vertex
	 * This method should be called after finishing a CG cycle and after solving a BAP node.
	 */
	public void cleanSRCs(){
		
		int auxNumArcs = 2*(V*V-V);
		for (int i = 1; i <= this.C; i++){
			this.vertices[i].SRCIndices = new ArrayList<>();

			int C0_ix = this.C0_startID + i;
			this.PPvertices[C0_ix].processedLabels = new ArrayList<Label>(auxNumArcs);

			int C1_ix = this.C1_startID + i;
			this.PPvertices[C1_ix].processedLabels = new ArrayList<Label>(auxNumArcs);
		}
	}

	/** Class that represents a vertex. */
	public class Vertex{

		public final int node_id;
		public final int xcoord; 								//x-coordinate
		public final int ycoord; 								//y-coordinate
		public final int load; 								//load of the vertex
		public final int opening_tw; 							//opening time window
		public final int closing_tw; 							//closing time window
		public final int last_departure;
		public final boolean feasible_nonfirst;
		public final int open_tw_nonfirst;
		public final int min_chargingTime;
		
		public final HashSet<Integer> unreachable; 			//(a priori) unreachable customers from this vertex
		public HashSet<Integer> neighbors;
		public ArrayList<Integer> SRCIndices; 			//indices of the SRC containing this vertex

		/**
		 * Creates a new (customer) vertex.
		 * @throws IOException Throws IO exception when the instance cannot be found.
		 */
		public Vertex(int id, int xcoord, int ycoord, int load, int opening_tw, int closing_tw, int last_departure, boolean feas_nonfirst, int tw_nonfirst, int min_b) {
			this.node_id = id;
			this.xcoord = xcoord;
			this.ycoord = ycoord;
			this.load = load;
			this.opening_tw = opening_tw;
			this.closing_tw = closing_tw;
			this.last_departure = last_departure;
			this.feasible_nonfirst = feas_nonfirst;
			this.open_tw_nonfirst = tw_nonfirst;
			this.min_chargingTime = min_b;

			this.unreachable = new HashSet<Integer>(C);
			this.SRCIndices = new ArrayList<>();
			this.neighbors = new HashSet<Integer>(C);
		}

		public Vertex(int id, int xcoord, int ycoord, int load, int opening_tw, int closing_tw) {
			this.node_id = id;
			this.xcoord = xcoord;
			this.ycoord = ycoord;
			this.load = load;
			this.opening_tw = opening_tw;
			this.closing_tw = closing_tw;
			this.last_departure = opening_tw;
			this.feasible_nonfirst = false;
			this.open_tw_nonfirst = opening_tw;
			this.min_chargingTime = 0;

			this.unreachable = null;
			this.SRCIndices = null;
			this.neighbors = null;
		}

		/**
		 * Obtains the string of a vertex
		 */
		@Override
		public String toString(){
			return "" + node_id;
		}
	}

	/** Class that represents an arc. */
	public class Arc{

		public final int tail; 							//tail vertex of the arc
		public final int head; 							//head vertex of the arc
		public final int id; 								//arc id
		public final int cost; 							//cost of the arc
		public final int time; 							//time of the arc
		public final int energy; 							//energy of the arc
		public final int energy_deviation;				//worst-case energy deviation of the arc
		public final int min_energy;
		public final int min_cost;
		public final int min_time;
		public final boolean minCostAlternative;

		/**
		 * Creates a new arc.
		 * @throws IOException Throws IO exception when the instance cannot be found.
		 */
		public Arc(int id, int tail, int head, int cost, int time, int energy, int energy_deviation, int min_energy, int min_cost, int min_time, boolean minCostAlt) {
			this.id = id;
			this.tail = tail;
			this.head = head;
			this.cost = cost;
			this.time = time;
			this.energy = energy;
			this.energy_deviation = energy_deviation;
			this.min_energy = min_energy;
			this.min_cost = min_cost;
			this.min_time = min_time;
			this.minCostAlternative = minCostAlt;
		}

		/** Obtains the string representation of the arc. */
		@Override
		public String toString(){
			return "Routing Arc " + this.id + " ("+tail+","+head+")";
		}

		@Override
		public int hashCode() {
			// TODO Auto-generated method stub
			return id;
		}
	}

	public class PPVertex{

		public final int id; 									//vertex id
		public final byte vertex_type;
		public final Vertex routing_vertex;						//corresponding routing vertex
		
		public final int node_number;

		public ArrayList<Label> processedLabels; 					//labels that have reached the vertex and are non-dominated
		public PriorityQueue<Label> unprocessedLabels; 				//labels that have reached the vertex but have not yet been processed

		/**
		 * Creates a new routing subgraph customer / depot vertex.
		 * @throws IOException Throws IO exception when the instance cannot be found.
		 */
		public PPVertex(int id, byte type, Vertex vertex, int number) {
			this.id = id;
			this.vertex_type = type;
			this.routing_vertex = vertex;
			this.node_number = number;

			int auxNumArcs = 2*(V*V-V);
			this.processedLabels = new ArrayList<Label>(auxNumArcs);
		}

		/**
		 * Creates a new charging subgraph vertex.
		 * @throws IOException Throws IO exception when the instance cannot be found.
		 */
		public PPVertex(int id, byte type, int number) {
			this.id = id;
			this.vertex_type = type;
			this.routing_vertex = null;
			this.node_number = number;

			int auxNumArcs = 2*(V*V-V);
			this.processedLabels = new ArrayList<Label>(auxNumArcs);
		}

		/**
		 * Obtains the string of a vertex
		 */
		@Override
		public String toString(){
			return VERTEX_TYPE_NAMES[this.vertex_type] +" "+ node_number;
		}
	}
	
	public class PPArc extends DefaultWeightedEdge{

		private static final long serialVersionUID = 1L;
		public final int id; 								//arc id
		public final byte arc_type;	
		public final Arc routing_arc;						//routing arc
		
		public final int tail_vertex_id; 						//tail vertex of the arc
		public final int head_vertex_id; 						//head vertex of the arc
		public double modifiedCost; 						//modified cost of the arc

		/**
		 * Creates a new arc belonging to the routing subgraph of the Pricing Problem.
		 * @throws IOException Throws IO exception when the instance cannot be found.
		 */
		public PPArc(int id, byte type, Arc arc, int tail, int head ) {
			this.id = id;
			this.arc_type = type;
			this.routing_arc = arc;
			this.tail_vertex_id = tail;
			this.head_vertex_id = head;
			this.modifiedCost = 0.0;
		}
		
		/**
		 * Creates a new arc belonging to the charging subgraph of the Pricing Problem.
		 * @throws IOException Throws IO exception when the instance cannot be found.
		 */
		public PPArc(int id, byte type, int tail, int head) {
			this.id = id;
			this.arc_type = type;
			this.routing_arc = null;
			this.tail_vertex_id = tail;
			this.head_vertex_id = head;
			this.modifiedCost = 0.0;
		}

		/** Obtains the string representation of the arc. */
		@Override
		public String toString(){
			if (this.arc_type <= AR1) return "PP Arc " + this.id + " ("+PPvertices[tail_vertex_id].toString()+", "+PPvertices[head_vertex_id].toString()+"); " + this.routing_arc.toString();
			else return "PP Arc " + this.id + "("+PPvertices[tail_vertex_id].toString()+", "+PPvertices[head_vertex_id].toString()+")";
		}

		@Override
		public int hashCode() {
			// TODO Auto-generated method stub
			return id;
		}
	}

	/**
	 * @return a negative integer, zero, or a positive integer as this object is less than, equal to, or greater than the specified object.
	 */
	public class SortByCost implements Comparator<Arc> {
		@Override
		public int compare(Arc a1, Arc a2) {
			if(a1.cost<a2.cost) return -1;
			if(a1.cost>a2.cost) return 1;
			return 0;
		}
	}

}
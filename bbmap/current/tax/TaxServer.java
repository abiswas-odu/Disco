package tax;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.net.InetSocketAddress;

import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;

import dna.Parser;
import fileIO.ReadWrite;
import shared.Tools;
import structures.IntList;

import com.sun.net.httpserver.Headers;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

/**
 * @author Shijie Yao, Brian Bushnell
 * @date Dec 13, 2016
 *
 */
public class TaxServer {

	public static void main(String[] args) throws Exception {
		TaxServer ts=new TaxServer(args);
		System.err.println("Ready!");
		//ts.begin();
	}
	
	public TaxServer(String[] args) throws Exception {
		int port_=8321;
		
		//Create a parser object
		Parser parser=new Parser();
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //Strip leading hyphens
			
			
			if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("verbose2")){
				verbose2=Tools.parseBoolean(b);
			}else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				tableFile=b;
				if("auto".equalsIgnoreCase(b)){tableFile=TaxTree.defaultTableFile();}
			}else if(a.equals("tree") || a.equals("taxtree")){
				treeFile=b;
				if("auto".equalsIgnoreCase(b)){treeFile=TaxTree.defaultTreeFile();}
			}else if(a.equals("accession")){
				accessionFile=b;
				if("auto".equalsIgnoreCase(b)){accessionFile=TaxTree.defaultAccessionFile();}
			}else if(a.equals("reverse")){
				reverseOrder=Tools.parseBoolean(b);
			}else if(a.equals("port")){
				port_=Integer.parseInt(b);
			}else{
				throw new RuntimeException(b);
			}
		}
		
		port=port_;
		
		typeMap=makeTypeMap();
		commonMap=makeCommonMap();
		
		
		
		if(tableFile!=null){
			outstream.println("Loading gi table.");
			GiToNcbi.initialize(tableFile);
		}
		if(treeFile!=null){
			outstream.println("Loading tree.");
			tree=ReadWrite.read(TaxTree.class, treeFile, true);
			if(tree.nameMap==null){
				outstream.println("Hashing names.");
				tree.hashNames();
			}
			assert(tree.nameMap!=null);
		}else{
			tree=null;
			throw new RuntimeException("No tree specified.");
		}
		if(accessionFile!=null){
			AccessionToTaxid.tree=tree;
			outstream.println("Loading accession table.");
			AccessionToTaxid.load(accessionFile);
			System.gc();
		}
		
		try {
			httpServer = HttpServer.create(new InetSocketAddress(port), 0);
		} catch (IOException e) {
			throw(e);
		}
		httpServer.createContext("/tax", new TaxHandler());
		httpServer.createContext("/help", new HelpHandler());
		httpServer.setExecutor(null); // creates a default executor
		httpServer.start();
	}

	class HelpHandler implements HttpHandler {
		
		public void handle(HttpExchange t) throws IOException {
			
			{
				Headers h = t.getResponseHeaders();
				String type="text/plain";
				h.add("Content-Type", type);
			}
			
			String response=USAGE;
			if(verbose2){System.out.println("Sending: "+response);}
			
			t.sendResponseHeaders(200, response.length());
			OutputStream os = t.getResponseBody();
			os.write(response.getBytes());
			os.close();
		}
	}

	class TaxHandler implements HttpHandler {
		
		public void handle(HttpExchange t) throws IOException {
			
			//String query = t.getRequestURI().getQuery(); //the KEY=VAL&KEY=VAL params in URL
			String rparam = t.getRequestURI().toString();   //restful style params, KEY/VAL in URL 
			while(rparam.startsWith("/")){
				rparam = rparam.substring(1);
			}
			while(rparam.endsWith("/")){
				rparam = rparam.substring(0, rparam.length()-1);
			}
			if(verbose){System.out.println(rparam);}

			String[] params = rparam.split("/");
			if(verbose2){System.out.println(Arrays.toString(params));}
			
			final String response=toResponse(params);
			
			{
				Headers h = t.getResponseHeaders();
				String type=response.startsWith("{") ? "application/json" : "text/plain";
				h.add("Content-Type", type);
			}
			
			if(verbose2){System.out.println("Sending: "+response);}
			
			t.sendResponseHeaders(200, response.length());
			OutputStream os = t.getResponseBody();
			os.write(response.getBytes());
			os.close();
		}
	}
	

	
	String toResponse(String[] params){
//		System.err.println("a");
		if(params.length<3){
			if(params.length==2 && "advice".equalsIgnoreCase(params[1])){return TAX_ADVICE;}
			return USAGE;
		}
//		System.err.println("b");
		
		if(params.length>4){return USAGE;}
//		System.err.println("c");
		
		final String query=params[params.length-1];
		final String[] names=query.split(",");
//		System.err.println("d");
		final boolean ancestor=(params.length>3 && params[2].equalsIgnoreCase("ancestor"));
//		System.err.println("e");
//		System.err.println(params[2]+", "+ancestor);
//		System.err.println("f");
		final int type;
		{
			String typeS=params[1];
			Integer value=typeMap.get(typeS);
			if(value==null){
				if(typeS.equalsIgnoreCase("advice")){
					return TAX_ADVICE;
				}else{
					return "{\"error\": \"Bad type; should be gi, taxid, or name.\"}";
				}
			}
			type=value.intValue();
		}
		final boolean plaintext=(type>=PT_OFFSET);
		
		if(verbose2){System.out.println("Type: "+type);}
		if(type==NAME || type==PT_NAME){
			for(int i=0; i<names.length; i++){
				if(names[i].contains("%20")){names[i]=names[i].replace("%20", " ");}
			}
			if(verbose2){System.out.println("Revised: "+Arrays.toString(names));}
		}
		
		if(ancestor){
			return toAncestor(type, names, plaintext, query);
		}
		
		if(plaintext){
			return toText(type, names);
		}
		
		if(names.length==1){
			return toJson(type, names[0]).toString();
		}
		ArrayList<JsonObject> list=new ArrayList<JsonObject>();
		for(String name : names){
			list.add(toJson(type, name));
		}
		return JsonObject.toString(list);
	}
	
	String toAncestor(final int type, final String[] names, boolean plaintext, String query){
		IntList ilist=toIntList(type, names);
		int id=FindAncestor.findAncestor(tree, ilist);
		TaxNode tn=(id>-1 ? tree.getNode(id) : null);
		if(tn==null){return new JsonObject(query, "error","Not found.").toString();}
		if(plaintext){return ""+id;}
		
		JsonObject j=new JsonObject(query);
		j.add("name", tn.name);
		j.add("tax_id", ""+tn.id);
		j.add("level", ""+tn.levelString());
		while(tn!=null && tn.level!=TaxTree.LIFE){
			j.add(toJson(tn));
			if(tn.pid==tn.id){break;}
			tn=tree.getNode(tn.pid);
		}
		return j.toString();
	}
	
	IntList toIntList(final int type, final String[] names){
		IntList list=new IntList(names.length);
		if(type==PT_GI || type==GI){
			for(String name : names){
				TaxNode tn=getTaxNodeGi(Integer.parseInt(name));
				if(tn!=null){list.add(tn.id);}
			}
		}else if(type==PT_NAME || type==NAME){
			for(String name : names){
				TaxNode tn=getTaxNodeByName(name);
				if(tn!=null){list.add(tn.id);}
			}
		}else if(type==PT_NCBI || type==NCBI){
			for(String name : names){
				TaxNode tn=getTaxNodeNcbi(Integer.parseInt(name));
				if(tn!=null){list.add(tn.id);}
			}
		}else if(type==PT_ACCESSION || type==ACCESSION){
			for(String name : names){
				int ncbi=AccessionToTaxid.get(name);
				if(ncbi>=0){list.add(ncbi);}
			}
		}else{
			throw new RuntimeException("{\"error\": \"Bad type\"}");
		}
		return list;
	}
	
	String toText(final int type, final String[] names){
		
		StringBuilder sb=new StringBuilder();
		String comma="";
		
		if(type==PT_GI){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeGi(Integer.parseInt(name));
				if(tn==null){sb.append("-1");}
				else{sb.append(tn.id);}
				comma=",";
			}
		}else if(type==PT_NAME){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeByName(name);
				if(tn==null){sb.append("-1");}
				else{sb.append(tn.id);}
				comma=",";
			}
		}else if(type==PT_NCBI){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeNcbi(Integer.parseInt(name));
				if(tn==null){sb.append("-1");}
				else{sb.append(tn.id);}
				comma=",";
			}
		}else if(type==PT_ACCESSION || type==ACCESSION){
			for(String name : names){
				sb.append(comma);
				int ncbi=AccessionToTaxid.get(name);
				sb.append(ncbi);
				comma=",";
			}
//		}else if(type==PT_HEADER){
//			for(String name : names){
//				sb.append(comma);
//				TaxNode tn=getTaxNodeHeader(name);
//				if(tn==null){sb.append("-1");}
//				else{sb.append(tn.id);}
//				comma=",";
//			}
		}else{
			return "Bad type; should be pt_gi or pt_name; e.g. /tax/pt_gi/1234";
		}
		
		return sb.toString();
	}
	
	JsonObject toJson(final int type, final String name){
		TaxNode tn=null;
		
		if(type==GI){
			tn=getTaxNodeGi(Integer.parseInt(name));
		}else if(type==NAME){
			tn=getTaxNodeByName(name);
		}else if(type==NCBI){
			tn=getTaxNodeNcbi(Integer.parseInt(name));
		}else if(type==ACCESSION){
			int ncbi=AccessionToTaxid.get(name);
			tn=(ncbi>=0 ? tree.getNode(ncbi) : null);
//		}else if(type==HEADER){
//			tn=getTaxNodeHeader(name);
		}else{
			return new JsonObject(""+type,"error","Bad type; should be gi, taxid, or name; e.g. /tax/name/homo_sapiens");
		}
		if(verbose2){System.out.println("Got node: "+tn);}
		
		if(tn!=null){
			JsonObject j=new JsonObject(name);
			j.add("name", tn.name);
			j.add("tax_id", ""+tn.id);
			j.add("level", ""+tn.levelString());
			while(tn!=null && tn.level!=TaxTree.LIFE){
				j.add(toJson(tn));
				if(tn.pid==tn.id){break;}
				tn=tree.getNode(tn.pid);
			}
			return j;
		}
		return new JsonObject(name, "error","Not found.");
	}
	
	JsonObject toJson(TaxNode tn){
		JsonObject j=new JsonObject(tn.levelString());
		j.add("name",tn.name);
		j.add("tax_id",""+tn.id);
		return j;
	}
	
	TaxNode getTaxNodeByName(String name){
		if(verbose2){System.out.println("Fetching node for "+name);}
		List<TaxNode> list=tree.getNodesByName(name);
		if(verbose2){System.out.println("Fetched "+list);}
		if(list==null){
			if(verbose2){System.out.println("Fetched in common map "+name);}
			String name2=commonMap.get(name);
			if(verbose2){System.out.println("Fetched "+name2);}
			if(name2!=null){list=tree.getNodesByName(name2);}
		}
		return list==null ? null : list.get(0);
	}
	
	TaxNode getTaxNodeGi(int gi){
		int ncbi=-1;
		try {
			ncbi=GiToNcbi.getID(gi);
		} catch (Throwable e) {
			if(verbose){e.printStackTrace();}
		}
		return ncbi<0 ? null : getTaxNodeNcbi(ncbi);
	}
	
	TaxNode getTaxNodeHeader(String header){
		int ncbi=-1;
		try {
			ncbi=GiToNcbi.getID(header);
		} catch (Throwable e) {
			if(verbose){e.printStackTrace();}
		}
		return ncbi<0 ? null : getTaxNodeNcbi(ncbi);
	}
	
	TaxNode getTaxNodeNcbi(int ncbi){
		TaxNode tn=null;
		try {
			tn=tree.getNode(ncbi);
		} catch (Throwable e) {
			if(verbose){e.printStackTrace();}
		}
		return tn;
	}
	
	/** This is called if the program runs with no parameters */
	private void printOptions(){
		throw new RuntimeException("TODO");
	}

	private HashMap<String, Integer> makeTypeMap() {
		HashMap<String, Integer> map=new HashMap<String, Integer>(63);
		map.put("gi", GI);
		map.put("name", NAME);
		map.put("tax_id", NCBI);
		map.put("ncbi", NCBI);
		map.put("taxid", NCBI);
		map.put("header", HEADER);
		map.put("accession", ACCESSION);
		map.put("pt_gi", PT_GI);
		map.put("pt_name", PT_NAME);
		map.put("pt_tax_id", PT_NCBI);
		map.put("pt_ncbi", PT_NCBI);
		map.put("pt_taxid", PT_NCBI);
		map.put("pt_header", PT_HEADER);
		map.put("pt_header", PT_HEADER);
		map.put("pt_accession", PT_ACCESSION);
		
		return map;
	}
	
	public static HashMap<String, String> makeCommonMap(){
		HashMap<String, String> map=new HashMap<String, String>();
		map.put("human", "homo sapiens");
		map.put("cat", "felis catus");
		map.put("dog", "canis lupus familiaris");
		map.put("mouse", "mus musculus");
		map.put("cow", "bos taurus");
		map.put("bull", "bos taurus");
		map.put("horse", "Equus ferus");
		map.put("pig", "Sus scrofa domesticus");
		map.put("sheep", "Ovis aries");
		map.put("goat", "Capra aegagrus");
		map.put("turkey", "Meleagris gallopavo");
		map.put("fox", "Vulpes vulpes");
		map.put("chicken", "Gallus gallus domesticus");
		map.put("wolf", "canis lupus");
		map.put("fruitfly", "drosophila melanogaster");
		map.put("zebrafish", "Danio rerio");
		map.put("catfish", "Ictalurus punctatus");
		map.put("trout", "Oncorhynchus mykiss");
		map.put("salmon", "Salmo salar");
		map.put("tilapia", "Oreochromis niloticus");
		map.put("e coli", "Escherichia coli");
		map.put("e.coli", "Escherichia coli");

		map.put("lion", "Panthera leo");
		map.put("tiger", "Panthera tigris");
		map.put("bear", "Ursus arctos");
		map.put("deer", "Odocoileus virginianus");
		map.put("coyote", "Canis latrans");

		map.put("corn", "Zea mays subsp. mays");
		map.put("maize", "Zea mays subsp. mays");
		map.put("oat", "Avena sativa");
		map.put("wheat", "Triticum aestivum");
		map.put("rice", "Oryza sativa");
		map.put("potato", "Solanum tuberosum");
		map.put("barley", "Hordeum vulgare");
		map.put("poplar", "Populus alba");
		map.put("lettuce", "Lactuca sativa");
		map.put("beet", "Beta vulgaris");
		map.put("strawberry", "Fragaria x ananassa");
		map.put("orange", "Citrus sinensis");
		map.put("lemon", "Citrus limon");
		map.put("soy", "Glycine max");
		map.put("soybean", "Glycine max");
		map.put("grape", "Vitis vinifera");
		map.put("olive", "Olea europaea");
		map.put("cotton", "Gossypium hirsutum");
		map.put("apple", "Malus pumila");
		map.put("bannana", "Musa acuminata");
		map.put("tomato", "Solanum lycopersicum");
		map.put("sugarcane", "Saccharum officinarum");
		map.put("bean", "Phaseolus vulgaris");
		map.put("onion", "Allium cepa");
		map.put("garlic", "Allium sativum");
		
		map.put("pichu", "mus musculus");
		map.put("pikachu", "mus musculus");
		map.put("vulpix", "Vulpes vulpes");
		map.put("ninetails", "Vulpes vulpes");
		map.put("mareep", "Ovis aries");
		
		return map;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String tableFile=TaxTree.defaultTableFile();
	private String treeFile=TaxTree.defaultTreeFile();
	private String accessionFile=null;
	
	private final TaxTree tree;
	private final HashMap<String, Integer> typeMap;
	private final HashMap<String, String> commonMap;
	
	/** Reverse order for tax lines */
	private boolean reverseOrder=true;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public final int port;
	
	public final HttpServer httpServer;

	public static final int PT_OFFSET=16;
	public static final int UNKNOWN=0, GI=1, NAME=2, NCBI=3, HEADER=4, ACCESSION=5;
	public static final int PT_GI=GI+PT_OFFSET, PT_NAME=NAME+PT_OFFSET, PT_NCBI=NCBI+PT_OFFSET, PT_HEADER=HEADER+PT_OFFSET, PT_ACCESSION=ACCESSION+PT_OFFSET;
	
	public static final String TAX_ADVICE="This site does not give tax advice.";
	public static final String USAGE="Welcome to the JGI taxonomy server!\n"
			+ "This service provides taxonomy information from NCBI taxID numbers, gi numbers, organism names, and accessions.\n"
			+ "The output is formatted as a Json object.\n"
			+ "Usage:\n\n"
			+ "/tax/name/homo_sapiens will give taxonomy information for an organism name.\n"
			+ "/tax/taxid/9606 will give taxonomy information for an NCBI taxID.\n"
			+ "/tax/gi/1234 will give taxonomy information from an NCBI gi number.\n"
			+ "/tax/accession/NZ_AAAA01000057.1 will give taxonomy information from an accession.\n"
//			+ "/tax/header/ will attempt to parse a sequence header such as >gi|7|emb|X51700.1| Bos taurus\n"
			+ "\nComma-delimited lists are accepted, such as tax/gi/1234,7000,42\n"
			+ "For plaintext (non-Json) results of the taxID number alone, use the pt_ prefix.  For example:\n\n"
			+ "/tax/pt_name/homo_sapiens\n"
			+ "/tax/pt_gi/1234\n"
			+ "\nTo find the common ancestor of multiple organisms, add /ancestor/. For example:\n"
			+ "/tax/taxid/ancestor/1234,5678,42\n"
			+ "/tax/name/ancestor/homo_sapiens,canis_lupus,bos_taurus\n";
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false, verbose2=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	
}

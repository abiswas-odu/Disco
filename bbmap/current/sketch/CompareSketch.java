package sketch;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

import dna.Parser;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.Heap;
import tax.PrintTaxonomy;
import tax.TaxNode;
import tax.TaxTree;

/**
 * Compares one or more input sketches to a set of reference sketches.
 * 
 * @author Brian Bushnell
 * @date July 29, 2016
 *
 */
public class CompareSketch extends SketchObject {
	

	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		CompareSketch cs=new CompareSketch(args);
		
		//Run the object
		cs.process(t);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CompareSketch(String[] args){
		
		//Process any config files
		args=Parser.parseConfig(args);
		
		//Detect whether the uses needs help
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		//Print the program name and arguments
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		//Set some shared static variables regarding PIGZ
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		Shared.READ_BUFFER_LENGTH=1;
		
		//Create a parser object
		Parser parser=new Parser();
		parser.out1="stdout.txt";
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //Strip leading hyphens
			
			if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("in")){
				addFiles(b, in);
			}else if(a.equals("ref")){
				addFiles(b, ref);
			}else if(Sketch.parseCoding(a, b)){
				//Do nothing
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Tools.parseKMG(b);
				//Set a variable here
			}else if(a.equals("mincount") || a.equals("minhits")  || a.equals("hits") || a.equals("matches")){
				minCount=Integer.parseInt(b);
			}else if(a.equals("minid") || a.equals("id")){
				minID=Float.parseFloat(b);
				if(minID>1){minID/=100;}
			}else if(a.equals("records")){
				maxRecords=Integer.parseInt(b);
				assert(maxRecords>=1) : "Max records must be at least 1.";
			}else if(a.equals("format")){
				format=Integer.parseInt(b);
			}
			
			
			else if(a.equals("size")){
				size=(int)Tools.parseKMG(b);
			}else if(a.equals("maxfraction")){
				Sketch.maxGenomeFraction=Float.parseFloat(b);
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(a.equals("rcomp")){
				rcomp=Tools.parseBoolean(b);
			}
			
			else if(a.equals("taxtree") || a.equals("tree")){
				taxTreeFile=b;
				if("auto".equalsIgnoreCase(b)){taxTreeFile=TaxTree.defaultTreeFile();}
			}else if(a.equals("level") || a.equals("taxlevel") || a.equals("minlevel")){
				taxLevel=Integer.parseInt(b);
			}else if(a.equals("printtax") || a.equals("printtaxa")){
				printTax=Tools.parseBoolean(b);
			}
			
			else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}
			
			else if(b==null && arg.indexOf('=')<0){//Parse standard flags in the parser
				ref.add(arg);
			}
			
			else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		{//Process parser fields
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			out=parser.out1;
		}
		
		//Ensure there is an input file
		if(in.isEmpty()){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		
		//Ensure there is an ref file
		if(ref.isEmpty()){
			printOptions();
			throw new RuntimeException("Error - at least one reference file is required.");
		}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		ffout=FileFormat.testOutput(out, FileFormat.TEXT, null, false, overwrite, append, false);
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, taxTreeFile) || !Tools.testInputFiles(true, true, in.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read to some input files.\n");
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in.toArray(new String[0])) || !Tools.testForDuplicateFiles(true, ref.toArray(new String[0]))){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		tool=new SketchTool(size, k, minCount, rcomp);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void process(Timer t){
		Timer ttotal=new Timer();
		
		if(taxTreeFile!=null){setTaxtree(taxTreeFile);}
		t.start();
		inSketches=tool.loadSketches_MT(in);
		refSketches=tool.loadSketches_MT(ref);
		final int numLoaded=(inSketches.size()+refSketches.size());
		t.stop();
		outstream.println("Loaded "+numLoaded+" sketches in \t"+t);
		t.start();
		
		final int threads=Shared.threads();
		final AtomicInteger fakeID=new AtomicInteger(minFakeID);
		@SuppressWarnings("unchecked")
		ConcurrentHashMap<Integer, Comparison> maps[]=new ConcurrentHashMap[inSketches.size()];
		for(int i=0; i<maps.length; i++){
			maps[i]=new ConcurrentHashMap<Integer, Comparison>(101);
		}
		
		ArrayList<CompareThread> alct=new ArrayList<CompareThread>(threads);
		for(int i=0; i<threads; i++){
			alct.add(new CompareThread(i, threads, fakeID, maps));
		}
		for(CompareThread ct : alct){ct.start();}
		boolean success=true;
		for(CompareThread ct : alct){
			
			//Wait until this thread has terminated
			while(ct.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					ct.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			
			synchronized(ct){
				success&=ct.success;
			}
		}
		alct=null;
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		
		@SuppressWarnings("unchecked")
		ArrayList<Comparison>[] alca=new ArrayList[inSketches.size()];
		for(int i=0; i<alca.length; i++){
			ConcurrentHashMap<Integer, Comparison> map=maps[i];
			maps[i]=null;
			ArrayList<Comparison> al=alca[i]=new ArrayList<Comparison>(map.size());
			for(Entry<Integer, Comparison> e : map.entrySet()){
				al.add(e.getValue());
			}
			Shared.sort(al);
			Collections.reverse(al);
			while(al.size()>maxRecords){
				al.remove(al.size()-1);
			}
		}
		
		TextStreamWriter tsw=(ffout==null ? null : new TextStreamWriter(ffout));
		if(tsw!=null){tsw.start();}
		if(tsw!=null){tsw.poisonAndWait();}
		errorState&=tsw.errorState;
		
		writeResults(alca);
		
		t.stop();
		outstream.println("\nRan "+(inSketches.size()*refSketches.size())+" comparisons in \t"+t);
		ttotal.stop();
		outstream.println("Total Time: \t"+ttotal);
	}
	
	private void writeResults(ArrayList<Comparison>[] alca){
		if(ffout==null){return;}
		TextStreamWriter tsw=new TextStreamWriter(ffout);
		tsw.start();
		for(int i=0; i<alca.length; i++){
			Sketch s=inSketches.get(i);
			ArrayList<Comparison> al=alca[i];
			writeResults(al, s, tsw);
		}
		tsw.poisonAndWait();
		errorState&=tsw.errorState;
	}
	
	private void writeResults(ArrayList<Comparison> al, Sketch s, TextStreamWriter tsw){
		tsw.println("\nResults for "+s.name()+":\n");
		
		ArrayList<TaxNode> tnl=new ArrayList<TaxNode>();
		StringBuilder sb=new StringBuilder();
		for(Comparison c : al){

			if(format==0){
				tsw.print(String.format("WKID %.2f%%\tKID %.2f%%\tmatches %d\tcompared %d",
						100*c.idWeighted(), 100*c.idMin(), c.hits, c.minIndex)+"\ttaxID "+c.taxID+"\tgSize "+c.genomeSize+"\t"+c.name+"\n");
				if(printTax && taxtree!=null && c.taxID>=0 && c.taxID<minFakeID){
					TaxNode tn=taxtree.getNode(c.taxID);
					if(tn!=null){
						PrintTaxonomy.printTaxonomy(tn, tsw, taxtree, TaxTree.DOMAIN);
					}
				}
				tsw.print("\n");
			}else{
				if(taxtree!=null && c.taxID>=0 && c.taxID<minFakeID){
					TaxNode tn=taxtree.getNode(c.taxID);
					while(tn!=null && tn.pid!=tn.id && tn.level<=TaxTree.DOMAIN){
						tnl.add(tn);
						tn=taxtree.getNode(tn.pid);
					}
				}
				
				sb.append(String.format("WKID %.2f%%\tKID %.2f%%\tmatches %d\tcompared %d\t",
						100*c.idWeighted(), 100*c.idMin(), c.hits, c.minIndex));
				sb.append("\ttaxID ").append(c.taxID).append('\t');
				sb.append(c.name).append('\t');
				
				for(int i=tnl.size()-1; i>=0; i--){
					TaxNode tn=tnl.get(i);
					sb.append(tn.name);
					if(i>0){sb.append(';');}
				}
				sb.append('\n');
				
				tsw.print(sb.toString());
				
				tnl.clear();
				sb.setLength(0);
			}
		}
	}
	
//	private static void compareOneToAll(final Sketch a, final ArrayList<Sketch> others, int minCount, float minID, 
//			TaxTree tree, Heap<Comparison> heap, TextStreamWriter tsw){
//		final int[] buffer=Sketch.makeBuffer();
//		
//		if(tsw!=null){tsw.print("Comparing "+a.name);}
//		
//		for(Sketch b : others){
//			Comparison c=compareOneToOne(a, b, buffer, minCount, minID, heap);
//		}
//		
//		ArrayList<Comparison> list=heap.toList();
//		Collections.reverse(list);
//		if(tsw!=null){
//			for(Comparison c : list){
////				final TaxNode tn=(taxtree==null || c.taxID<1 ? null : taxtree.getNode(c.taxID));
//				tsw.print(String.format("\n%.2f%% WKID, %.2f%% KID, %d matches, %d compared",
//						100*c.idWeighted(), 100*c.idMin(), c.hits, c.minIndex)+"\t"+c.taxID+"\t"+c.name);
//			}
//		}
//		
//		if(tsw!=null){tsw.print("\n");}
//	}
	
	private static Comparison compareOneToOne(final Sketch a, final Sketch b, int[] buffer, int minCount, float minID, Heap<Comparison> heap){
//		assert(heap!=null); //Optional, for testing.
		final int matches=a.countMatches(b, buffer);
		assert(matches==buffer[0]);
		final int div=buffer[1];
		
		if(matches<minCount || matches/(float)div<minID){
//			System.err.print(".");
			return null;
		}
		if(heap!=null && !heap.hasRoom() && heap.peek().hits>matches){return null;}
		
//		System.err.print("*");
		Comparison c=new Comparison(buffer, b);
		if(heap==null || heap.add(c)){return c;}
		return null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private void addFiles(String a, Collection<String> list){
		if(a==null){return;}
		File f=null;
		if(a.indexOf(',')>=0){f=new File(a);}
		if(f==null || f.exists()){
			list.add(a);
		}else{
			for(String s : a.split(",")){list.add(s);}
		}
	}
	
	/** This is called if the program runs with no parameters */
	private void printOptions(){
		throw new RuntimeException("Please run the shellscript with no arguments for usage information."); //TODO
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class CompareThread extends Thread {
		
		CompareThread(final int tid_, final int incr_, final AtomicInteger fake_, ConcurrentHashMap<Integer, Comparison>[] maps_){
			tid=tid_;
			incr=incr_;
			fakeID=fake_;
			maps=maps_;
		}
		
		public void run(){
			success=false;
			final int refLim=refSketches.size();
			final int inLim=inSketches.size();
			
			for(int refNum=tid; refNum<refLim; refNum+=incr){
				final Sketch b=refSketches.get(refNum);
				
				for(int inNum=0; inNum<inLim; inNum++){
					final Sketch a=inSketches.get(inNum);
					processPair(a, b, maps[inNum]);
				}
			}
			synchronized(this){success=true;}
		}
		
		private boolean processPair(Sketch a, Sketch b, ConcurrentHashMap<Integer, Comparison> map){
			Comparison c=compareOneToOne(a, b, buffer, minCount, minID, null);
			if(c==null){return false;}
			if(c.taxID<0){c.taxID=fakeID.getAndIncrement();}
			
			TaxNode tn=(taxtree==null ? null : taxtree.getNode(b.taxID));
			if(tn!=null){
				c.name=tn.name;
				if(tn.level<taxLevel){
					tn=taxtree.getNode(b.taxID, taxLevel);
				}
			}
			Integer key=(tn==null ? c.taxID : tn.id);

			Comparison old=map.get(key);
//			System.err.println("A. Old: "+(old==null ? 0 : old.hits)+", new: "+c.hits);
			if(old!=null && old.compareTo(c)>0){return false;}
			
			old=map.put(key, c);
			while(old!=null && old.compareTo(c)>0){
//				System.err.println("B. Old: "+(old==null ? 0 : old.hits)+", new: "+c.hits);
				c=old;
				old=map.put(key, c);
			}
			return true;
		}
		
		private final int tid, incr;
		private final int[] buffer=Sketch.makeBuffer();
		
		private final AtomicInteger fakeID;
		private ConcurrentHashMap<Integer, Comparison> maps[];
		
		boolean success=false;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private ArrayList<String> in=new ArrayList<String>();
	
	private String out="stdout.txt";
	
	private String taxTreeFile=null;
	
	private ArrayList<String> ref=new ArrayList<String>();

	private int size=10000;
	private int k=31;
	private boolean rcomp=true;
	private final SketchTool tool;
	
	private int maxRecords=100;
	
	private int minCount=3;
	
	private float minID=0.02f;
	
	private int taxLevel=2;
	
	private ArrayList<Sketch> inSketches;
	
	private ArrayList<Sketch> refSketches;
	
	private boolean printTax=true;
	
	private int format=0;
	
	private static final int minFakeID=1900000000;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary output file */
	private final FileFormat ffout;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=false;
	/** Append to existing output files */
	private boolean append=false;
	
}

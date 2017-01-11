package tax;

import java.io.PrintStream;
import java.util.Arrays;

import dna.Parser;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ByteBuilder;
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;
import structures.IntList;

/**
 * @author Brian Bushnell
 * @date May 9, 2016
 *
 */
public class FindAncestor {
	
	public static void main(String[] args){
		Timer t=new Timer();
		FindAncestor fa=new FindAncestor(args);
		fa.process(t);
	}
	
	public FindAncestor(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				tableFile=b;
				if("auto".equalsIgnoreCase(b)){tableFile=TaxTree.defaultTableFile();}
			}else if(a.equals("tree") || a.equals("taxtree")){
				treeFile=b;
				if("auto".equalsIgnoreCase(b)){treeFile=TaxTree.defaultTreeFile();}
			}else if(a.equals("invalid")){
				outInvalid=b;
			}else if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}

		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);
		ffoutInvalid=FileFormat.testOutput(outInvalid, FileFormat.TXT, null, true, overwrite, append, false);
		ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);
		
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
		lifeNode=tree.getNodeByName("life");
	}
	
	void process(Timer t){
		
		ByteFile bf=ByteFile.makeByteFile(ffin1, false);
		ByteStreamWriter bsw=new ByteStreamWriter(ffout1);
		bsw.start();
		
		bsw.print("#Name\tAncestor\tMajority\tTaxonomy...\n".getBytes());
		
		ByteStreamWriter bswInvalid=null;
		if(ffoutInvalid!=null){
			bswInvalid=new ByteStreamWriter(ffoutInvalid);
			bswInvalid.start();
		}
		
//		final HashArray1D counts=countTable ? new HashArray1D(256000, true) : null;
		final IntList giList=new IntList();
		final IntList ncbiList=new IntList();
		final IntList traversal=new IntList();
		
		byte[] line=bf.nextLine();
		ByteBuilder bb=new ByteBuilder();
		
		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=line.length;
				
				giList.clear();
				ncbiList.clear();
				traversal.clear();
				
				final int giCount=getGiNumbers(line, giList, ',');
				final int ncbiCount=getNcbiNumbers(giList, ncbiList);
				
				taxaCounted+=giCount;
				taxaValid+=ncbiCount;
				final boolean valid=(ncbiCount>0);
				
				if(valid){
					linesValid++;
					int ancestor=findAncestor(ncbiList);
					int majority=findMajority(ncbiList);
					
					for(int i=0; i<line.length && line[i]!='\t'; i++){
						bb.append(line[i]);
					}
					bb.append('\t');
					bb.append(ancestor);
					bb.append('\t');
					bb.append(majority);
					bb.append('\t');
					
					fillTraversal(majority, traversal, true);
					writeTraversal(traversal, bb);
					bb.append('\n');
					
					for(int i=0; i<ncbiList.size; i++){
						final int id=ncbiList.get(i);
						fillTraversal(id, traversal, true);
						writeTraversal(traversal, bb);
						bb.append('\n');
					}
					
					bsw.print(bb.toBytes());
					bb.clear();
				}else{
					if(bswInvalid!=null){
						bswInvalid.println(line);
					}
				}
			}
			line=bf.nextLine();
		}
		
		errorState|=bf.close();
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		if(bswInvalid!=null){errorState|=bswInvalid.poisonAndWait();}
		
		t.stop();
		
		double rpnano=linesProcessed/(double)(t.elapsed);
		double bpnano=bytesProcessed/(double)(t.elapsed);

		String rpstring=(linesProcessed<100000 ? ""+linesProcessed : linesProcessed<100000000 ? (linesProcessed/1000)+"k" : (linesProcessed/1000000)+"m");
		String bpstring=(bytesProcessed<100000 ? ""+bytesProcessed : bytesProcessed<100000000 ? (bytesProcessed/1000)+"k" : (bytesProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		outstream.println("Time:                         \t"+t);
		outstream.println("Lines Processed:    "+rpstring+" \t"+String.format("%.2fk lines/sec", rpnano*1000000));
		outstream.println("Bytes Processed:    "+bpstring+" \t"+String.format("%.2fm bytes/sec", bpnano*1000));
		
		outstream.println();
		outstream.println("Valid Lines:       \t"+linesValid);
		outstream.println("Invalid Lines:     \t"+(linesProcessed-linesValid));
//		if(counts!=null){
//			outstream.println("Unique Taxa:       \t"+taxaCounted);
//		}
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private void fillTraversal(int id, IntList traversal, boolean addLife){
		traversal.clear();
		for(TaxNode node=tree.getNode(id); node!=null && node!=lifeNode; node=tree.getNode(node.pid)){
			traversal.add(node.id);
		}
		if(addLife || traversal.size==0){traversal.add(lifeNode.id);}
	}
	
	private void writeTraversal(IntList traversal, ByteBuilder bb){
		for(int i=traversal.size-1; i>=0; i--){
			final int id=traversal.get(i);
			if(id>=0){
				TaxNode tn=tree.getNode(id);
				//			bb.append(tn.level+"_"+tn.name);
				bb.append(/*tn.level+"_"+*/tn.name);
				if(i>0){bb.append('\t');}
			}
		}
	}
	
	private int getGiNumbers(final byte[] line, final IntList list, final char delimiter){
		int i=0;
		
		//Skip name
		while(i<line.length && line[i]!='\t'){i++;}
		
		//Skip whitespaces
		while(i<line.length && Character.isWhitespace(line[i])){i++;}
		
		while(i<line.length){
			while(i<line.length && line[i]==delimiter){i++;}
			int start=i;
			while(i<line.length && line[i]!=delimiter){i++;}
			final int stop=i;
			if(Tools.startsWith(line, prefix, start)){start+=3;}
			assert(start<stop) : "Badly formatted line at "+start+":\n"+new String(line);
//			System.err.println(start+","+stop+",'"+new String(line).substring(start, stop)+"'");
			if(start<stop){
				final int number=Tools.parseInt(line, start, stop);
				list.add(number);
			}
		}
		return list.size;
	}
	
	private int getNcbiNumbers(final IntList giList, final IntList ncbiList){
		final int size=giList.size;
		for(int i=0; i<size; i++){
			final int gi=giList.get(i);
			final int ncbi=GiToNcbi.getID(gi);
//			System.err.println(gi+" -> "+ncbi);
			if(ncbi>=0){ncbiList.add(ncbi);}
		}
//		for(int i=0; i<20; i++){
//			System.err.println(GiToNcbi.getID(i*i*i));
//		}
//		assert(false);
		return ncbiList.size;
	}
	
	private int findAncestor(IntList list){
		return findAncestor(tree, list);
	}
	
	public static int findAncestor(TaxTree tree, IntList list){
		if(list.size<1){
			assert(false);
			return -1;
		}
		int ancestor=list.get(0);
		for(int i=1; i<list.size && ancestor>-1; i++){
			final int id=list.get(i);
//			System.err.println(ancestor+"+"+id+" -> "+tree.commonAncestor(ancestor, id));
			int x=tree.commonAncestor(ancestor, id);
			if(x>-1){
				ancestor=x;
			}
		}
//		System.err.println("Ancestor node: "+tree.getNode(ancestor));
//		System.err.println(list+" -> "+ancestor);
//		if(ancestor<0){ancestor=lifeNode.id;}
		return ancestor;
	}
	
	private int findMajority(IntList list){
		if(list.size<3){return findAncestor(list);}
		final int majority=list.size/2+1;
//		System.err.println("Majority: "+majority);
		
		for(int i=0; i<list.size; i++){
			final int id=list.get(i);
			TaxNode tn=tree.getNode(id);
//			System.err.println("Found node "+tn);
			assert(tn!=null) : "No node for id "+id;
			if(tn!=null){
				tree.percolateUp(tn, 1);
			}
		}
		
		TaxNode best=lifeNode;
		for(int i=0; i<list.size; i++){
			final int id=list.get(i);
			TaxNode tn=tree.getNode(id);
			while(tn!=null && tn!=lifeNode){
				if(tn.countSum>=majority && tn.level<best.level){
					best=tn;
					break;
				}
				tn=tree.getNode(tn.pid);
			}
		}
		
//		System.err.println("Best node: "+best);
		
		for(int i=0; i<list.size; i++){
			final int id=list.get(i);
			TaxNode tn=tree.getNode(id);
			if(tn!=null){
				tree.percolateUp(tn, -1);
			}
		}
		
		return best.id;
	}
	
	/*--------------------------------------------------------------*/
	
	private void printOptions(){
		assert(false) : "printOptions: TODO";
	}
	
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	private String outInvalid=null;

	private String tableFile=null;
	private String treeFile=null;
	
	private final TaxTree tree;
	
	private final TaxNode lifeNode;
	
	/*--------------------------------------------------------------*/
	
	private long taxaCounted=0;
	private long taxaValid=0;
	private long linesProcessed=0;
	private long linesValid=0;
	private long bytesProcessed=0;
	
	private long maxLines=Long.MAX_VALUE;

//	private boolean prefix=false;
	private boolean countTable=true;
//	private boolean keepInvalidSequence=false;
	
	private final String prefix="gi|";
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	private final FileFormat ffoutInvalid;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	
//	private static final byte[] ncbi=">ncbi|".getBytes();
	
}

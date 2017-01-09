package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import dna.Parser;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sort.Shuffle;
import sort.Shuffle.ShuffleThread;
import fileIO.FileFormat;


/**
 * Generates artificial cross-contaminated data by mixing reads.
 * Takes input from multiple files, and writes output to the same number of files.
 * @author Brian Bushnell
 * @date Oct 27, 2014
 *
 */
public class CrossContaminate {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		Timer t=new Timer();
		CrossContaminate cc=new CrossContaminate(args);
		cc.process(t);
	}
	
	public CrossContaminate(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		for(String s : args){if(s.startsWith("out=standardout") || s.startsWith("out=stdout")){outstream=System.err;}}
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		boolean setInterleaved=false; //Whether it was explicitly set.

		
		
		Shared.READ_BUFFER_LENGTH=Tools.min(200, Shared.READ_BUFFER_LENGTH);
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		

		ArrayList<String> inTemp=new ArrayList<String>();
		ArrayList<String> outTemp=new ArrayList<String>();
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens

			if(Parser.isJavaFlag(arg)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(parser.parseCommon(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
//				align2.FastaReadInputStream2.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("in")){
				String[] split2=b.split(",");
				for(String name : split2){
					inNames.add(name);
				}
			}else if(a.equals("out")){
				String[] split2=b.split(",");
				for(String name : split2){
					outNames.add(name);
				}
			}else if(a.equals("innamefile")){
				String[] split2=b.split(",");
				for(String name : split2){
					inTemp.add(name);
				}
			}else if(a.equals("outnamefile")){
				String[] split2=b.split(",");
				for(String name : split2){
					outTemp.add(name);
				}
			}else if(a.equals("shuffle")){
				shuffle=Tools.parseBoolean(b);
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}else if(a.equals("minsinks") || a.equals("ns")){
				minSinks=Integer.parseInt(b);
			}else if(a.equals("maxsinks") || a.equals("xs")){
				maxSinks=Integer.parseInt(b);
			}else if(a.equals("minprob") || a.equals("np")){
				minProb=Double.parseDouble(b);
			}else if(a.equals("maxprob") || a.equals("xp")){
				maxProb=Double.parseDouble(b);
			}else if(a.equals("showspeed")){
				showspeed=Tools.parseBoolean(b);
			}else if(a.equals("shufflethreads")){
				shufflethreads=Integer.parseInt(b);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=parser.overwrite;
			append=parser.append;

			setInterleaved=parser.setInterleaved;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		DecontaminateByNormalization.parseStringsFromFiles(inTemp);
		DecontaminateByNormalization.parseStringsFromFiles(outTemp);
		
		inNames.addAll(inTemp);
		outNames.addAll(outTemp);
		inTemp=outTemp=null;
		
		if(inNames.isEmpty() || inNames.size()!=outNames.size()){
			assert(false) : inNames+"\n"+outNames;
			printOptions();
			throw new RuntimeException("Error - at least one input file is required, and # input files must equal # output files.");
		}
		
		assert(minSinks<=maxSinks);
		minSinks=Tools.max(0, minSinks);
		maxSinks=Tools.min(inNames.size()-1, maxSinks);
		assert(minSinks<=maxSinks) : minSinks+", "+maxSinks;
		
		assert(minProb<=maxProb);
		assert(minProb>=0 && maxProb<=1);

		minProbPow=Math.log(minProb);
		maxProbPow=Math.log(maxProb);
		
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(!Tools.testInputFiles(true, true, inNames.toArray(new String[0]))){
			outstream.println(outNames);
			throw new RuntimeException("Can't find some input files:\n"+inNames+"\n");
		}
		
		if(!Tools.testOutputFiles(overwrite, append, false, outNames.toArray(new String[0]))){
			outstream.println(outNames);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files.\n");
		}
		
		if(seed>0){randy.setSeed(seed);}
		
		vessels=makeVessels(outNames);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	void process(Timer t){
		
		outstream.println("Processing data.");
		for(int i=0; i<inNames.size(); i++){
			try{
				processOneSource(i);
			}catch(Throwable e){
				System.err.println("Failed to open file "+inNames.get(i)+"\nException:"+e+"\n");
				errorState=true;
			}
		}
		
		for(Vessel v : vessels){
			errorState|=v.close();
		}
		
		if(shuffle){
			shuffle(shufflethreads);
		}
		
		t.stop();
		
		if(showspeed){
			double rpnano=readsProcessed/(double)(t.elapsed);
			double bpnano=basesProcessed/(double)(t.elapsed);

			String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
			String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");

			while(rpstring.length()<8){rpstring=" "+rpstring;}
			while(bpstring.length()<8){bpstring=" "+bpstring;}

			outstream.println("Time:                         \t"+t);
			outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
			outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		}
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	void shuffle(final int threads){
		outstream.println("Shuffling output in "+threads+" thread"+(threads==1 ? "." : "s."));
		Shuffle.showSpeed=Shuffle.printClass=false;
		Shuffle.setMaxThreads(threads);
		for(Vessel v : vessels){
			ShuffleThread st=new ShuffleThread(v.fname, null, v.fname, null, Shuffle.SHUFFLE, true);
			st.start();
		}
		Shuffle.waitForFinish();
	}
	
	void processOneSource(int sourceNum){
		String fname=inNames.get(sourceNum);
		
		FileFormat ffin=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, true);
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin, null);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		final boolean paired=cris.paired();
		if(verbose){
			if(!ffin.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		}
		
		ArrayList<Vessel> sinks=assignSinks(vessels, sourceNum);
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			outstream.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((r.mate!=null)==paired);
			}

			while(reads!=null && reads.size()>0){
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					{
						readsProcessed++;
						basesProcessed+=initialLength1;
					}
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=initialLength2;
					}
					
					addRead(r1, sinks);
				}

				cris.returnList(ln.id, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadWrite.closeStream(cris);
	}
	
	private void addRead(Read r, ArrayList<Vessel> list){
		double p=randy.nextDouble();
		for(Vessel v : list){
			if(p>=v.prob){
				v.bsw.println(r, true);
				r=null;
				break;
			}
		}
		assert(r==null) : p+"\n"+list;
	}
	
	private ArrayList<Vessel> makeVessels(ArrayList<String> strings){
		ArrayList<Vessel> list=new ArrayList<Vessel>(strings.size());
		for(String s : strings){
			Vessel v=new Vessel(s, true);
			list.add(v);
		}
		return list;
	}
	
	private ArrayList<Vessel> assignSinks(ArrayList<Vessel> list, int sourceNum){
		int potential=list.size()-1;
		assert(potential>=minSinks && maxSinks<=potential) : potential+", "+minSinks+", "+maxSinks;
		int range=maxSinks-minSinks+1;
		
		int sinks=minSinks+(range>0 ? randy.nextInt(range) : 0);
		assert(sinks>=0);
		
		for(Vessel v : list){v.prob=0;}
		ArrayList<Vessel> sinklist=(ArrayList<Vessel>) list.clone();
		list=null;
		Vessel source=sinklist.remove(sourceNum);
		if(verbose || true){
			System.err.println("Source:   \t"+inNames.get(sourceNum));
			System.err.println("Sinks:    \t"+sinks);
		}
		
		while(sinklist.size()>sinks){
			int x=randy.nextInt(sinklist.size());
			sinklist.set(x, sinklist.get(sinklist.size()-1));
			sinklist.remove(sinklist.size()-1);
		}
//		if(verbose){System.err.println("Sinklist:\n"+sinklist);}
		
		{
			double probRange=maxProbPow-minProbPow;
			
			assert(probRange>=0) : minProb+", "+maxProb+", "+minProbPow+", "+maxProbPow+", "+probRange;
			
			double remaining=1.0;
			for(Vessel v : sinklist){
				double c=Math.pow(Math.E, minProbPow+randy.nextDouble()*probRange)*remaining;
				remaining-=c;
				v.prob=c;
			}
			source.prob=remaining;
			sinklist.add(source);
			if(verbose || true){System.err.println("Sinklist:\t"+sinklist+"\n");}
			double d=0;
			for(Vessel v : sinklist){
				d+=v.prob;
				v.prob=d;
			}
//			if(verbose){System.err.println("Sinklist:\t"+sinklist);}
			d=0;
			for(Vessel v : sinklist){
				double temp=v.prob;
				v.prob=d;
				d=temp;
			}
//			if(verbose){System.err.println("Sinklist:\t"+sinklist);}
		}
		Collections.reverse(sinklist);
		assert(sinklist.get(sinklist.size()-1).prob==0.0) : sinklist;
		
//		if(verbose){
//			System.err.println("Sinklist:\t"+sinklist);
//			System.err.println();
//		}
		if(verbose || true){System.err.println();}
		
		return sinklist;
	}
	
	private void printOptions(){
		assert(false) : "printOptions: TODO";
//		outstream.println("Syntax:\n");
//		outstream.println("java -ea -Xmx512m -cp <path> jgi.ReformatReads in=<infile> in2=<infile2> out=<outfile> out2=<outfile2>");
//		outstream.println("\nin2 and out2 are optional.  \nIf input is paired and there is only one output file, it will be written interleaved.\n");
//		outstream.println("Other parameters and their defaults:\n");
//		outstream.println("overwrite=false  \tOverwrites files that already exist");
//		outstream.println("ziplevel=4       \tSet compression level, 1 (low) to 9 (max)");
//		outstream.println("interleaved=false\tDetermines whether input file is considered interleaved");
//		outstream.println("fastawrap=70     \tLength of lines in fasta output");
//		outstream.println("qin=auto         \tASCII offset for input quality.  May be set to 33 (Sanger), 64 (Illumina), or auto");
//		outstream.println("qout=auto        \tASCII offset for output quality.  May be set to 33 (Sanger), 64 (Illumina), or auto (meaning same as input)");
//		outstream.println("outsingle=<file> \t(outs) Write singleton reads here, when conditionally discarding reads from pairs.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class Vessel{
		
		public Vessel(String fname_, boolean allowSubprocess){
			fname=fname_;
			ff=FileFormat.testOutput(fname, FileFormat.FASTQ, null, allowSubprocess, overwrite, append, false);
			bsw=new ByteStreamWriter(ff);
			bsw.start();
		}
		
		public boolean close(){
			return bsw.poisonAndWait();
		}
		
		public String toString(){
			return fname+", "+String.format("%.6f", prob);
		}
		
		final String fname;
		final FileFormat ff;
		final ByteStreamWriter bsw;
		
		double prob;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/*--------------------------------------------------------------*/

	private ArrayList<String> inNames=new ArrayList<String>();
	private ArrayList<String> outNames=new ArrayList<String>();
	
	private ArrayList<Vessel> vessels;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	private long seed=-1;
	
	private int minSinks=1;
	private int maxSinks=8;
	private double minProb=0.000005;
	private double maxProb=0.025;

	private double minProbPow=Math.log(minProb);
	private double maxProbPow=Math.log(maxProb);
	
//	private double root=3.0;
//
//	private double minProbRoot=Math.pow(minProb, 1/root);
//	private double maxProbRoot=Math.pow(maxProb, 1/root);
	
	private final Random randy=new Random();
	
	long readsProcessed=0;
	long basesProcessed=0;
	
	private int shufflethreads=3;
	
	private boolean shuffle=false;
	private boolean showspeed=true;
	
	/*--------------------------------------------------------------*/
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	
}

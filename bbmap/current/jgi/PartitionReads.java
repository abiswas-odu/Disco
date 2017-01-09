package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.SamLine;
import structures.ListNum;
import dna.Parser;
import fileIO.ByteFile;
import fileIO.ReadWrite;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import fileIO.FileFormat;

/**
 * @author Brian Bushnell
 * @date June 1, 2016
 *
 */
public class PartitionReads {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		PartitionReads as=new PartitionReads(args);
		as.process(t);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public PartitionReads(String[] args){
		
		//Process any config files
		args=Parser.parseConfig(args);
		
		//Detect whether the uses needs help
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		//Print the program name and arguments
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		boolean setInterleaved=false; //Whether interleaved was explicitly set.
		
		//Set some shared static variables regarding PIGZ
		Shared.READ_BUFFER_LENGTH=Tools.max(400, Shared.READ_BUFFER_LENGTH);
		Shared.capBuffers(2);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		SamLine.SET_FROM_OK=true;
		
		//Create a parser object
		Parser parser=new Parser();
		
		String qfout1=null, qfout2=null;
		
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
			}else if(a.equals("ways")){
				ways=Integer.parseInt(b);
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;

			out1=parser.out1;
			out2=parser.out2;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		assert(ways>0) : "Ways must be at least 1.";
		
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		
		//Do qfin # replacement
		if(qfin1!=null && qfin2==null && qfin1.indexOf('#')>-1 && !new File(qfin1).exists()){
			qfin2=qfin1.replace("#", "2");
			qfin1=qfin1.replace("#", "1");
		}
		
		//Do output file # replacement
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		
		//Do output file # replacement
		if(qfout1!=null && qfout2==null && qfout1.indexOf('#')>-1){
			qfout2=qfout1.replace("#", "2");
			qfout1=qfout1.replace("#", "1");
		}
		
		//Adjust interleaved detection based on the number of input files
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		//Ensure there is an input file
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		//Ensure out2 is not set without out1
		if(out1==null){
			if(out2!=null){
				printOptions();
				throw new RuntimeException("Error - cannot define out2 without defining out1.");
			}
		}
		
		//Adjust interleaved settings based on number of output files
		if(!setInterleaved){
			assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}
		
		assert(out1==null || out1.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		assert(out2==null || out2.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		assert(qfout1==null || qfout1.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		assert(qfout2==null || qfout2.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read to some input files.\n");
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		//Create output FileFormat objects
		ffout1=new FileFormat[ways];
		ffout2=new FileFormat[ways];
		qfout1Array=new String[ways];
		qfout2Array=new String[ways];
		for(int i=0; i<ways; i++){
			String temp1=out1==null ? null : out1.replaceFirst("%", ""+i);
			String temp2=out2==null ? null : out2.replaceFirst("%", ""+i);
			ffout1[i]=FileFormat.testOutput(temp1, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
			ffout2[i]=FileFormat.testOutput(temp2, FileFormat.FASTQ, extout, true, overwrite, append, ordered);

			qfout1Array[i]=qfout1==null ? null : qfout1.replaceFirst("%", ""+i);
			qfout2Array[i]=qfout2==null ? null : qfout2.replaceFirst("%", ""+i);
		}
		
		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		useSharedHeader=(ffin1.samOrBam() && ffout1!=null && ffout1.length>0 && ffout1[0]!=null && ffout1[0].samOrBam());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros[]=new ConcurrentReadOutputStream[ways];
		if(ffout1!=null && ffout1.length>0){
			final int buff=1;
			
			for(int i=0; i<ways; i++){
				ros[i]=ConcurrentReadOutputStream.getStream(ffout1[i], ffout2[i], qfout1Array[i], qfout2Array[i], buff, null, useSharedHeader);
				ros[i].start(); //Start the stream
			}
		}
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the read stream
		processInner(cris, ros);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		//Report timing and results
		{
			t.stop();
			
			//Calculate units per nanosecond
			double rpnano=readsProcessed/(double)(t.elapsed);
			double bpnano=basesProcessed/(double)(t.elapsed);
			
			//Add "k" and "m" for large numbers
			String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
			String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");
			
			//Format the strings so they have they are right-justified
			while(rpstring.length()<8){rpstring=" "+rpstring;}
			while(bpstring.length()<8){bpstring=" "+bpstring;}
			
			outstream.println("Time:                         \t"+t);
			outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
			outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		}
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Iterate through the reads */
	void processInner(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros[]){
		
		//Do anything necessary prior to processing
		
		{
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			ArrayList<Read>[] outLists=new ArrayList[ways];
			for(int i=0; i<ways; i++){
				outLists[i]=new ArrayList<Read>();
			}
			
			int nextIndex=0;
			
			//As long as there is a nonempty read list...
			while(reads!=null && reads.size()>0){
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				//Loop through each read in the list
				for(Read r1 : reads){
					final Read r2=r1.mate;
					
					//Track the initial length for statistics
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					//Increment counters
					readsProcessed+=1+r1.mateCount();
					basesProcessed+=initialLength1+initialLength2;
					
					outLists[nextIndex].add(r1);
					nextIndex=(nextIndex+1)%ways;
				}
				
				//Output reads to the output stream
				for(int i=0; i<ways; i++){
					if(ros!=null){ros[i].add(outLists[i], ln.id);}
					outLists[i]=new ArrayList<Read>();
				}
				
				//Notify the input stream that the list was used
				cris.returnList(ln.id, ln.list.isEmpty());
				if(verbose){outstream.println("Returned a list.");}
				
				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		//Do anything necessary after processing
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This is called if the program runs with no parameters */
	private void printOptions(){
		throw new RuntimeException("TODO");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;
	
	private String qfin1=null;
	private String qfin2=null;

	/** Primary output file path */
	private String out1=null;
	/** Secondary output file path */
	private String out2=null;

	private String[] qfout1Array=null;
	private String[] qfout2Array=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	/** Split data into this many output files */
	private int ways=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	/** Secondary input file */
	private final FileFormat ffin2;
	
	/** Primary output file */
	private final FileFormat[] ffout1;
	/** Secondary output file */
	private final FileFormat[] ffout2;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private boolean useSharedHeader=false;
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
	/** This flag has no effect on singlethreaded programs */
	private final boolean ordered=false;
	
}

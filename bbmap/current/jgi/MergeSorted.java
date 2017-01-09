package jgi;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import dna.Parser;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sort.SortByName;
import stream.FASTQ;
import stream.FastaReadInputStream;

public class MergeSorted {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		MergeSorted as=new MergeSorted(args);
		as.process(t);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public MergeSorted(String[] args){

		//Process any config files
		args=Parser.parseConfig(args);
		
		//Print the program name and arguments
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		
		//Set some shared static variables regarding PIGZ
		Shared.READ_BUFFER_LENGTH=Tools.min(200, Shared.READ_BUFFER_LENGTH);
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		String in=null;
		
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
			
			
			if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
//				SortByName.verbose=verbose;
			}else if(a.equals("delete")){
				delete=Tools.parseBoolean(b);
			}else if(a.equals("in")){
				in=b;
			}else if(a.equals("groups") || a.equals("g") || a.equals("sets") || a.equals("ways")){  
				groups=Integer.parseInt(b);
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			parser.processQuality();
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			out1=parser.out1;
			out2=parser.out2;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		if(in==null){
			throw new RuntimeException("Error - at least one input file is required.");
		}
		
		if(groups>0){
			assert(in.contains("%")) : "With a set number of groups, input filename must contain a % symbol where the number is.";
			for(int i=0; i<groups; i++){
				inList.add(in.replaceFirst("%", ""+i));
			}
		}else{
			String[] split=in.split(",");
			for(String s : split){inList.add(s);}
		}
		
		//Ensure there is an input file
		if(inList.isEmpty()){
			throw new RuntimeException("Error - at least one input file is required.");
		}
		
		//Do output file # replacement
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		//Ensure out2 is not set without out1
		if(out1==null){
			if(out2!=null){
				throw new RuntimeException("Error - cannot define out2 without defining out1.");
			}
		}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, inList.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read to some input files.\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	public void process(Timer t){
		
		if(verbose){
			System.err.println("Processing "+inList);
		}
		errorState=SortByName.mergeAndDump(inList, ffout1, ffout2, delete, ffout1.samOrBam());
		
		//Report timing and results
		{
			t.stop();
			
			outstream.println("Time:                         \t"+t);
		}
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary output file path */
	private String out1=null;
	/** Secondary output file path */
	private String out2=null;
	
	private ArrayList<String> inList=new ArrayList<String>();
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;

	private int groups=-1;
	private boolean delete=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary output file */
	private final FileFormat ffout1;
	/** Secondary output file */
	private final FileFormat ffout2;
	
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
	/** This flag has no effect on singlethreaded programs */
	private final boolean ordered=false;
	
}

package jgi;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.KillSwitch;
import stream.Read;
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
 * @date June 20, 2014
 *
 */
public class Shred {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		Shred mb=new Shred(args);
		mb.process(t);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public Shred(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;

		
		
		Shared.READ_BUFFER_LENGTH=Tools.min(100, Shared.READ_BUFFER_LENGTH);
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		long seed=-1;
		Parser parser=new Parser();
		boolean even=false;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens

			if(a.equals("length") || a.equals("len") || a.equals("shredlen") || a.equals("shredlength")){
				shredLength=(int)Tools.parseKMG(b);
			}else if(a.equals("overlap")){
				overlap=(int)Tools.parseKMG(b);
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("even") || a.equals("equal")){
				even=Tools.parseBoolean(b);
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}else if(a.equals("median")){
				median=(int)Tools.parseKMG(b);
			}else if(a.equals("variance")){
				variance=(int)Tools.parseKMG(b);
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		evenLengths=even;
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
			
			extin=parser.extin;
			extout=parser.extout;
			
			minLength=parser.minReadLength;
		}
		
		minLength=Tools.mid(1, minLength, shredLength);
		assert(shredLength>0);
		assert(shredLength>overlap);
		increment=shredLength-overlap;
		incMult=1.0/increment;
		assert(increment>0);
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		
		randy=(seed>=0 ? new Random(seed) : new Random());
		if(median>0 && variance<0){variance=median;}
	}
	
	public boolean parseArgument(String arg, String a, String b){
		if(a.equals("reads") || a.equals("maxreads")){
			maxReads=Tools.parseKMG(b);
			return true;
		}else if(a.equals("some_argument")){
			maxReads=Tools.parseKMG(b);
			return true;
		}
		return false;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			cris.start();
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}

		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=2;
			
			if(cris.paired()){KillSwitch.kill("This program does not support paired reads.");}
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, null, buff, null, false);
			ros.start();
		}else{ros=null;}
		
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the read stream
		processInner(cris, ros);
		
		ReadWrite.closeStreams(cris, ros);
		if(verbose){outstream.println("Finished.");}
		
		errorState|=ReadStats.writeAll();
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		t.stop();
		
		final double rpnano=readsProcessed/(double)(t.elapsed);
		final double bpnano=basesProcessed/(double)(t.elapsed);
		
		outstream.println("Time:                         \t"+t);
		
		{
			String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
			String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");
			
			while(rpstring.length()<8){rpstring=" "+rpstring;}
			while(bpstring.length()<8){bpstring=" "+bpstring;}
			
			outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
			outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		}

		{
			String rpstring=(readsOut<100000 ? ""+readsOut : readsOut<100000000 ? (readsOut/1000)+"k" : (readsOut/1000000)+"m");
			String bpstring=(basesOut<100000 ? ""+basesOut : basesOut<100000000 ? (basesOut/1000)+"k" : (basesOut/1000000)+"m");

			while(rpstring.length()<8){rpstring=" "+rpstring;}
			while(bpstring.length()<8){bpstring=" "+bpstring;}
			
			outstream.println("Reads Out:          "+rpstring);
			outstream.println("Bases Out:          "+bpstring);
		}
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Iterate through the reads */
	void processInner(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){
		
		readsProcessed=0;
		basesProcessed=0;
		
		readsOut=0;
		basesOut=0;
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		if(reads!=null && !reads.isEmpty()){
			Read r=reads.get(0);
			assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
		}

		while(reads!=null && reads.size()>0){
			if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
			
			ArrayList<Read> listOut=new ArrayList<Read>();
			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx);

				final int initialLength1=r1.length();

				readsProcessed++;
				basesProcessed+=initialLength1;
				
				if(median>0){
					processRandomly(r1, listOut);
				}else if(evenLengths){
					processEvenly(r1, listOut);
				}else{
					processUnevenly(r1, listOut);
				}
			}

			if(ros!=null){ros.add(listOut, ln.id);}

			cris.returnList(ln.id, ln.list.isEmpty());
			if(verbose){outstream.println("Returned a list.");}
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	void processRead(final Read r1, final ArrayList<Read> list){
		if(r1.length()<minLength){return;}
		if(r1.length()<=shredLength){
			r1.numericID=readsOut;
			list.add(r1);
			readsOut++;
			basesOut+=r1.length();
			return;
		}
		if(evenLengths){
			processEvenly(r1, list);
		}else{
			processUnevenly(r1, list);
		}
	}
	
	void processEvenly(final Read r1, final ArrayList<Read> list){
		final byte[] bases=r1.bases;
		final byte[] quals=r1.quality;
		final String name=r1.id;
		
		final int chunks=(int)Math.ceil((bases.length-overlap)*incMult);
		assert(chunks>0);
		double inc2=bases.length/(double)chunks;
		
		for(int chunk=0; chunk<chunks; chunk++){
			int a=(int)Math.floor(inc2*chunk);
			int b=(chunk==chunks-1 ? bases.length : overlap+(int)Math.floor(inc2*(chunk+1)));
			b=Tools.min(b, a+shredLength);
			final int length=b-a;
			if(length<minLength){return;}
			final byte[] bases2=KillSwitch.copyOfRange(bases, a, b);
			final byte[] quals2=(quals==null ? null : KillSwitch.copyOfRange(quals, a, b));
			Read shred=new Read(bases2, quals2, readsOut, name+"_"+a+"-"+(b-1));
			readsOut++;
			basesOut+=shred.length();
			list.add(shred);
		}
	}
	
	void processUnevenly(final Read r1, final ArrayList<Read> list){
		final byte[] bases=r1.bases;
		final byte[] quals=r1.quality;
		final String name=r1.id;
		for(int i=0; i<bases.length; i+=increment){
			final int limit=Tools.min(i+shredLength, bases.length);
			final int length=limit-i;
			if(length<minLength){return;}
			final byte[] bases2=KillSwitch.copyOfRange(bases, i, limit);
			final byte[] quals2=(quals==null ? null : KillSwitch.copyOfRange(quals, i, limit));
			Read shred=new Read(bases2, quals2, readsOut, name+"_"+i+"-"+(limit-1));
			readsOut++;
			basesOut+=shred.length();
			list.add(shred);
			if(limit==bases.length){return;}
			assert(limit<bases.length);
		}
	}
	
	void processRandomly(final Read r1, final ArrayList<Read> list){
		final byte[] bases=r1.bases;
		final byte[] quals=r1.quality;
		final String name=r1.id;
		for(int i=0; i<bases.length;){
			int rand=Tools.min(randy.nextInt(2*variance), randy.nextInt(3*variance), 2*variance);
			final int limit=Tools.max(i+minLength, Tools.min(i+rand+median-variance, bases.length));
			final int length=limit-i;
			if(length<minLength || limit>bases.length){return;}
			final byte[] bases2=KillSwitch.copyOfRange(bases, i, limit);
			final byte[] quals2=(quals==null ? null : KillSwitch.copyOfRange(quals, i, limit));
			Read shred=new Read(bases2, quals2, readsOut, name+"_"+i+"-"+(limit-1));
			readsOut++;
			basesOut+=shred.length();
			list.add(shred);
			if(limit==bases.length){return;}
			assert(limit<bases.length);
			i=limit;
		}
	}
	
	/** This is called if the program runs with no parameters */
	private void printOptions(){
		throw new RuntimeException("TODO");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	
	private String extin=null;
	private String extout=null;
	
	/*--------------------------------------------------------------*/

	protected long readsProcessed=0;
	protected long basesProcessed=0;
	protected long readsOut=0;
	protected long basesOut=0;
	
	private long maxReads=-1;
	
	private int median=-1;
	private int variance=-1;
	
	private int shredLength=500;
	private int minLength=1;
	private int overlap=0;
	private final int increment;
	private final double incMult;
	
	private final boolean evenLengths;
	
	private final Random randy;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	
}

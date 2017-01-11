package stream;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.FileFormat;
import shared.Shared;
import shared.Timer;

/**
 * Loads sam files rapidly with multiple threads.
 * 
 * @author Brian Bushnell
 * @date November 4, 2016
 *
 */
public abstract class SamStreamer {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static final void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		int threads=Shared.threads();
		if(args.length>1){threads=Integer.parseInt(args[1]);}
		SamStreamer as=new SamReadStreamer(args[0], threads);
//		SamStreamer as=new SamStreamer(args[0], 1);
		
		//Run the object
		as.start();
		as.test();
		
		t.stop("Time: ");
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public SamStreamer(String fname_, int threads_){
		this(FileFormat.testInput(fname_, FileFormat.SAM, null, true, false), threads_);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public SamStreamer(FileFormat ffin_){
		this(ffin_, DEFAULT_THREADS);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public SamStreamer(FileFormat ffin_, int threads_){
		fname=ffin_.name();
		threads=threads_;
		ffin=ffin_;
		
		inq=new ArrayBlockingQueue<ArrayList<byte[]>>(threads/2+1);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	
	final void test(){
		for(ArrayList<Read> list=nextReads(); list!=null; list=nextReads()){
			if(verbose){outstream.println("Got list of size "+list.size());}
		}
	}
	
	
	/** Create read streams and process all data */
	public final void start(){
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the reads in separate threads
		spawnThreads();
		
		if(verbose){outstream.println("Finished; closing streams.");}
	}

	public final ArrayList<Read> nextList(){return nextReads();}
	public abstract ArrayList<Read> nextReads();
	public abstract ArrayList<SamLine> nextLines();
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/** Spawn process threads */
	abstract void spawnThreads();
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	protected String fname;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Quit after processing this many input reads; -1 means no limit */
	protected long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	final FileFormat ffin;
	
	final ArrayBlockingQueue<ArrayList<byte[]>> inq;
	
	final int threads;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	static final ArrayList<Read> POISON_READS=new ArrayList<Read>(0);
	static final ArrayList<SamLine> POISON_LINES=new ArrayList<SamLine>(0);
	static final ArrayList<byte[]> POISON_BYTES=new ArrayList<byte[]>(0);
	static final int LIST_SIZE=200;
	public static int DEFAULT_THREADS=3;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	protected PrintStream outstream=System.err;
	/** Print verbose messages */
	public static final boolean verbose=false;
	public static final boolean verbose2=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	
}

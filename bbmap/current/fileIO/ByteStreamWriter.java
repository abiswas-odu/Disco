package fileIO;

import java.io.IOException;
import java.io.OutputStream;
import java.util.concurrent.ArrayBlockingQueue;

import kmer.AbstractKmerTable;
import shared.Shared;
import shared.Timer;
import stream.ByteBuilder;
import stream.Read;
import ukmer.AbstractKmerTableU;

import dna.AminoAcid;
import dna.Data;



/**
 * @author Brian Bushnell
 * @date Oct 21, 2014
 *
 */
public class ByteStreamWriter extends Thread {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static void main(String[] args){
		Timer t=new Timer();
		final int alen=1000;
		byte[] array=new byte[alen];
		for(int i=0; i<array.length; i++){
			array[i]=AminoAcid.numberToBase[i&3];
		}
		array[array.length-1]='\n';
		long iters=Long.parseLong(args[1]);
		String fname=args[0];
		ByteStreamWriter bsw=new ByteStreamWriter(fname, true, false, true);
		bsw.start();
		for(long i=0; i<iters; i++){
			bsw.print(array);
		}
		bsw.poisonAndWait();
		t.stop();
		System.err.println("MB/s: \t"+String.format("%.2f", ((alen*iters)/(t.elapsed/1000.0))));
		System.err.println("Time: \t"+t);
	}
	
	public ByteStreamWriter(String fname_, boolean overwrite_, boolean append_, boolean allowSubprocess_){
		this(fname_, overwrite_, append_, allowSubprocess_, 0);
	}
	
	public ByteStreamWriter(String fname_, boolean overwrite_, boolean append_, boolean allowSubprocess_, int format){
		this(FileFormat.testOutput(fname_, FileFormat.TEXT, format, 0, allowSubprocess_, overwrite_, append_, true));
	}
	
	public ByteStreamWriter(FileFormat ff){
		FASTQ=ff.fastq() || ff.text();
		FASTA=ff.fasta();
		BREAD=ff.bread();
		SAM=ff.samOrBam();
		BAM=ff.bam();
		SITES=ff.sites();
		INFO=ff.attachment();
		OTHER=(!FASTQ && !FASTA && !BREAD && !SAM && !BAM && !SITES && !INFO);
		
		
		fname=ff.name();
		overwrite=ff.overwrite();
		append=ff.append();
		allowSubprocess=ff.allowSubprocess();
		assert(!(overwrite&append));
		assert(ff.canWrite()) : "File "+fname+" exists and overwrite=="+overwrite;
		if(append && !(ff.raw() || ff.gzip())){throw new RuntimeException("Can't append to compressed files.");}
		
		if(!BAM || !Data.SAMTOOLS() || !Data.SH()){
			outstream=ReadWrite.getOutputStream(fname, append, true, allowSubprocess);
		}else{
			outstream=ReadWrite.getOutputStreamFromProcess(fname, "samtools view -S -b -h - ", true, append, true, true);
		}
		
		queue=new ArrayBlockingQueue<ByteBuilder>(5);
		buffer=new ByteBuilder(initialLen);
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Primary Method        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	@Override
	public void run() {
		if(verbose){System.err.println("running");}
		assert(open) : fname;
		
		synchronized(this){
			started=true;
			this.notify();
		}
		
		ByteBuilder job=null;

		if(verbose){System.err.println("waiting for jobs");}
		while(job==null){
			try {
				job=queue.take();
//				job.list=queue.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		if(verbose){System.err.println("processing jobs");}
		while(job!=null && job!=POISON2){
			if(job.length()>0){
				try {
					outstream.write(job.array, 0, job.length());
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
			job=null;
			while(job==null){
				try {
					job=queue.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		if(verbose){System.err.println("null/poison job");}
//		assert(false);
		open=false;
		ReadWrite.finishWriting(null, outstream, fname, allowSubprocess);
		if(verbose){System.err.println("finish writing");}
		synchronized(this){notifyAll();}
		if(verbose){System.err.println("done");}
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------      Control and Helpers     ----------------*/
	/*--------------------------------------------------------------*/
	
	
	@Override
	public void start(){
		super.start();
		if(verbose){System.err.println(this.getState());}
		synchronized(this){
			while(!started){
				try {
					this.wait(20);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}

	
	public synchronized void poison(){
		//Don't allow thread to shut down before it has started
		while(!started || this.getState()==Thread.State.NEW){
			try {
				this.wait(20);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		if(!open){return;}
		addJob(buffer);
		buffer=null;
//		System.err.println("Poisoned!");
//		assert(false);
		
//		assert(false) : open+", "+this.getState()+", "+started;
		open=false;
		addJob(POISON2);
	}
	
	public void waitForFinish(){
		while(this.getState()!=Thread.State.TERMINATED){
			try {
				this.join(1000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}
	
	/**
	 * @return true if there was an error, false otherwise
	 */
	public boolean poisonAndWait(){
		poison();
		waitForFinish();
		return errorState;
	}
	
	//TODO Why is this synchronized?
	public synchronized void addJob(ByteBuilder bb){
//		System.err.println("Got job "+(j.list==null ? "null" : j.list.size()));
		
		assert(started) : "Wait for start() to return before using the writer.";
//		while(!started || this.getState()==Thread.State.NEW){
//			try {
//				this.wait(20);
//			} catch (InterruptedException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//		}
		
		boolean success=false;
		while(!success){
			try {
				queue.put(bb);
				success=true;
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				assert(!queue.contains(bb)); //Hopefully it was not added.
			}
		}
	}
	
	/** Called after every write to the buffer */
	private final void flushBuffer(boolean force){
		final int x=buffer.length();
		if(x>=maxLen || (force && x>0)){
			addJob(buffer);
			buffer=new ByteBuilder(initialLen);
		}
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Print             ----------------*/
	/*--------------------------------------------------------------*/
	
	
	@Deprecated
	/** Avoid using this if possible. */
	public void print(CharSequence x){
		if(verbose){System.err.println("Added line '"+x+"'");}
		assert(open) : x;
		buffer.append(x);
		flushBuffer(false);
	}
	
	@Deprecated
	/** Avoid using this if possible. */
	public void print(StringBuilder x){
		if(verbose){System.err.println("Added line '"+x+"'");}
		assert(open) : x;
		buffer.append(x);
		flushBuffer(false);
	}
	
	@Deprecated
	/** Avoid using this if possible. */
	public void print(String x){
		if(verbose){System.err.println("Added line '"+x+"'");}
		assert(open) : x;
		buffer.append(x);
		flushBuffer(false);
	}
	
	public void print(int x){
		if(verbose){System.err.println("Added line '"+(x)+"'");}
		assert(open) : x;
		buffer.append(x);
		flushBuffer(false);
	}
	
	public void print(long x){
		if(verbose){System.err.println("Added line '"+(x)+"'");}
		assert(open) : x;
		buffer.append(x);
		flushBuffer(false);
	}
	
	public void print(float x){
		if(verbose){System.err.println("Added line '"+(x)+"'");}
		assert(open) : x;
		buffer.append(x);
		flushBuffer(false);
	}
	
	public void print(double x){
		if(verbose){System.err.println("Added line '"+(x)+"'");}
		assert(open) : x;
		buffer.append(x);
		flushBuffer(false);
	}
	
	public void print(byte x){
		if(verbose){System.err.println("Added line '"+((char)x)+"'");}
		assert(open) : ((char)x);
		buffer.append(x);
		flushBuffer(false);
	}
	
	public void print(char x){
		if(verbose){System.err.println("Added line '"+(x)+"'");}
		assert(open) : (x);
		buffer.append(x);
		flushBuffer(false);
	}
	
	public void print(byte[] x){
		if(verbose){System.err.println("Added line '"+new String(x)+"'");}
		assert(open) : new String(x);
		buffer.append(x);
		flushBuffer(false);
	}
	
	public void print(char[] x){
		if(verbose){System.err.println("Added line '"+new String(x)+"'");}
		assert(open) : new String(x);
		buffer.append(x);
		flushBuffer(false);
	}
	
	public void print(ByteBuilder x){
		if(verbose){System.err.println("Added line '"+x+"'");}
		assert(open) : x;
		buffer.append(x);
		flushBuffer(false);
	}
	
	public void print(ByteBuilder x, boolean destroy){
		if(!destroy || buffer.length()>0){print(x);}
		else{
			if(verbose){System.err.println("Added line '"+x+"'");}
			assert(open) : x;
			addJob(x);
		}
	}
	
	public void print(Read r){
		assert(!OTHER);
		ByteBuilder x=(FASTQ ? r.toFastq(buffer) : FASTA ? r.toFasta(FASTA_WRAP, buffer) : SAM ? r.toSam(buffer) : 
			SITES ? r.toSitesB(buffer) : INFO ? r.toInfoB(buffer) : r.toText(true, buffer));
		flushBuffer(false);
	}
	
	public void printKmer(long kmer, int count, int k){
		AbstractKmerTable.toBytes(kmer, count, k, buffer);
		flushBuffer(false);
	}
	
	public void printKmer(long kmer, int[] values, int k){
		AbstractKmerTable.toBytes(kmer, values, k, buffer);
		flushBuffer(false);
	}
	
	public void printKmer(long[] array, int count, int k){
		AbstractKmerTableU.toBytes(array, count, k, buffer);
		flushBuffer(false);
	}
	
	public void printKmer(long[] array, int[] values, int k){
		AbstractKmerTableU.toBytes(array, values, k, buffer);
		flushBuffer(false);
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------           Println            ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public void println(){print('\n');}
	public void println(CharSequence x){print(x); print('\n');}
	public void println(String x){print(x); print('\n');}
	public void println(StringBuilder x){print(x); print('\n');}
	public void println(int x){print(x); print('\n');}
	public void println(long x){print(x); print('\n');}
	public void println(float x){print(x); print('\n');}
	public void println(double x){print(x); print('\n');}
	public void println(byte x){print(x); print('\n');}
	public void println(char x){print(x); print('\n');}
	public void println(byte[] x){print(x); print('\n');}
	public void println(char[] x){print(x); print('\n');}
	public void println(ByteBuilder x){print(x); print('\n');}
	public void println(ByteBuilder x, boolean destroy){
		if(destroy){print(x.append('\n'));}else{print(x); print('\n');}
	}
	public void printlnKmer(long kmer, int count, int k){printKmer(kmer, count, k); print('\n');}
	public void printlnKmer(long kmer, int[] values, int k){printKmer(kmer, values, k); print('\n');}
	public void printlnKmer(long[] array, int count, int k){printKmer(array, count, k); print('\n');}
	public void printlnKmer(long[] array, int[] values, int k){printKmer(array, values, k); print('\n');}
	public void println(Read r){print(r); print('\n');}

	
	public void println(Read r, boolean paired){
		println(r);
		if(paired && r.mate!=null){println(r.mate);}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private ByteBuilder buffer;
	
	public int initialLen=36000;
	public int maxLen=32768;
	public final boolean overwrite;
	public final boolean append;
	public final boolean allowSubprocess;
	public final String fname;
	private final OutputStream outstream;
	private final ArrayBlockingQueue<ByteBuilder> queue;
	private boolean open=true;
	private volatile boolean started=false;
	
	/** TODO */
	public boolean errorState=false;
	
	/*--------------------------------------------------------------*/
	
	private final boolean BAM;
	private final boolean SAM;
	private final boolean FASTQ;
	private final boolean FASTA;
	private final boolean BREAD;
	private final boolean SITES;
	private final boolean INFO;
	private final boolean OTHER;
	
	private final int FASTA_WRAP=Shared.FASTA_WRAP;
	
	/*--------------------------------------------------------------*/

//	private static final ByteBuilder POISON=new ByteBuilder("POISON_ByteStreamWriter");
	private static final ByteBuilder POISON2=new ByteBuilder(1);
	
	public static boolean verbose=false;
	
}

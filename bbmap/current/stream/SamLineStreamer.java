package stream;

import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteFile;
import fileIO.FileFormat;

/**
 * Loads sam files rapidly with multiple threads.
 * 
 * @author Brian Bushnell
 * @date November 4, 2016
 *
 */
public class SamLineStreamer extends SamStreamer {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public SamLineStreamer(String fname_, int threads_){
		this(FileFormat.testInput(fname_, FileFormat.SAM, null, true, false), threads_);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public SamLineStreamer(FileFormat ffin_){
		this(ffin_, DEFAULT_THREADS);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public SamLineStreamer(FileFormat ffin_, int threads_){
		super(ffin_, threads_);
		outq=new ArrayBlockingQueue<ArrayList<SamLine>>(threads+1);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public ArrayList<SamLine> nextLines(){
		ArrayList<SamLine> list=null;
		while(list==null){
			try {
				list=outq.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			if(verbose){outstream.println("a. Got list size "+list.size());}
		}
		while(list==POISON_LINES){
			if(verbose){outstream.println("b. Got poison.");}
			try {
				outq.put(list);
				list=null;
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		if(verbose){outstream.println("c. done.");}
		return list;
	}
	
	public ArrayList<Read> nextReads(){
		KillSwitch.kill("Unsupported.");
		return null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/** Spawn process threads */
	void spawnThreads(){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=this.threads+1;
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(i, alpt));
		}
		if(verbose){outstream.println("Spawned threads.");}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		if(verbose){outstream.println("Started threads.");}
		
		//Do anything necessary after processing
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	private class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final int tid_, ArrayList<ProcessThread> alpt_){
			tid=tid_;
			alpt=(tid==0 ? alpt_ : null);
		}
		
		//Called by start()
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			if(tid==0){
				processBytes();
			}else{
				makeReads();
			}
			
			//Indicate successful exit status
			success=true;
		}
		
		void processBytes(){
			if(verbose){outstream.println("tid "+tid+" started processBytes.");}

//			ByteFile.FORCE_MODE_BF1=true;
			ByteFile.FORCE_MODE_BF2=true;
			ByteFile bf=ByteFile.makeByteFile(ffin, false);
			
			ArrayList<byte[]> list=new ArrayList<byte[]>(LIST_SIZE);
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				assert(line!=null);
//				outstream.println("a");
				list.add(line);
				if(list.size()>=LIST_SIZE){
//					outstream.println("b");
//					outstream.println(inq.size()+", "+inq.remainingCapacity());
					putBytes(list);
//					outstream.println("c");
					list=new ArrayList<byte[]>(LIST_SIZE); 
				}
//				outstream.println("d");
			}
			if(verbose){outstream.println("tid "+tid+" ran out of input.");}
			if(list.size()>0){
				putBytes(list);
				list=null;
			}
			if(verbose || verbose2){outstream.println("tid "+tid+" done reading bytes.");}
			putBytes(POISON_BYTES);
			if(verbose || verbose2){outstream.println("tid "+tid+" done poisoning.");}
			bf.close();
			if(verbose || verbose2){outstream.println("tid "+tid+" closed stream.");}
			
			success=true;
			
			//Wait for completion of all threads
			boolean allSuccess=true;
			for(ProcessThread pt : alpt){
				
				//Wait until this thread has terminated
				if(pt!=this){
					if(verbose){outstream.println("Waiting for thread "+pt.tid);}
					while(pt.getState()!=Thread.State.TERMINATED){
						try {
							//Attempt a join operation
							pt.join();
						} catch (InterruptedException e) {
							//Potentially handle this, if it is expected to occur
							e.printStackTrace();
						}
					}

					//Accumulate per-thread statistics
					readsProcessed+=pt.readsProcessedT;
					basesProcessed+=pt.basesProcessedT;
					allSuccess&=pt.success;
				}
			}
			
			putReads(POISON_LINES);
			if(verbose || verbose2){outstream.println("tid "+tid+" done poisoning reads.");}
			
			//Track whether any threads failed
			if(!allSuccess){errorState=true;}
			if(verbose || verbose2){outstream.println("tid "+tid+" finished!");}
			
		}
		
		void putBytes(ArrayList<byte[]> list){
			if(verbose){outstream.println("tid "+tid+" putting blist size "+list.size());}
			while(list!=null){
				try {
					inq.put(list);
					list=null;
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			if(verbose){outstream.println("tid "+tid+" done putting blist");}
		}
		
		void putReads(ArrayList<SamLine> list){
			if(verbose){outstream.println("tid "+tid+" putting rlist size "+list.size());}
			while(list!=null){
				try {
					outq.put(list);
					list=null;
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			if(verbose){outstream.println("tid "+tid+" done putting rlist");}
		}
		
		ArrayList<byte[]> takeBytes(){
			if(verbose){outstream.println("tid "+tid+" taking blist");}
			ArrayList<byte[]> list=null;
			while(list==null){
				try {
					list=inq.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			if(verbose){outstream.println("tid "+tid+" took blist size "+list.size());}
			return list;
		}
		
		/** Iterate through the reads */
		void makeReads(){
			if(verbose){outstream.println("tid "+tid+" started makeReads.");}
			
			ArrayList<byte[]> list=takeBytes();
			while(list!=POISON_BYTES){
				ArrayList<SamLine> reads=new ArrayList<SamLine>(list.size());
				for(byte[] line : list){
					if(line[0]=='@'){
						//ignore;
					}else{
						SamLine sl=new SamLine(line);
						reads.add(sl);
						
						readsProcessedT++;
						basesProcessedT+=(sl.seq==null ? 0 : sl.length());
					}
				}
				if(reads.size()>0){putReads(reads);}
				list=takeBytes();
			}
			if(verbose || verbose2){outstream.println("tid "+tid+" done making reads.");}

			putBytes(POISON_BYTES);
			if(verbose || verbose2){outstream.println("tid "+tid+" done poisoning bytes.");}
			
//			putReads(POISON_LINES);
//			if(verbose || verbose2){outstream.println("tid "+tid+" done poisoning reads.");}
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Thread ID */
		final int tid;
		
		ArrayList<ProcessThread> alpt;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	final ArrayBlockingQueue<ArrayList<SamLine>> outq;
	
}

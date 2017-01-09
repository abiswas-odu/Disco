package fileIO;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.concurrent.ArrayBlockingQueue;

import shared.Timer;


/**
 * Runs a ByteFile1 in a separate thread.  Can speed up disk reading, particularly of compressed files, at cost of slightly more work done.
 * Drop-in compatible with ByteFile1.
 * @author Brian Bushnell
 * @date Sep 23, 2013
 *
 */
public class ByteFile2 extends ByteFile {
	
	
	public static void main(String[] args) throws IOException{
		ByteFile2 tf=new ByteFile2(args.length>0 ? args[0] : "stdin", false, true);
		long first=0, last=100;
		boolean speedtest=false;
		if(args.length>1){
			if(args[1].equalsIgnoreCase("speedtest")){
				speedtest=true;
				first=0;
				last=Long.MAX_VALUE;
			}else{
				first=Integer.parseInt(args[1]);
				last=first+100;
			}
		}
		if(args.length>2){
			last=Integer.parseInt(args[2]);
		}
		speedtest(tf, first, last, !speedtest);
		
		tf.close();
		tf.reset();
		tf.close();
	}
	
	private static void speedtest(ByteFile2 tf, long first, long last, boolean reprint){
		Timer t=new Timer();
		long lines=0;
		long bytes=0;
		for(long i=0; i<first; i++){tf.nextLine();}
		if(reprint){
			for(long i=first; i<last; i++){
				byte[] s=tf.nextLine();
				if(s==null){break;}

				lines++;
				bytes+=s.length;
				System.out.println(new String(s));
			}
			
			System.err.println("\n");
			System.err.println("Lines: "+lines);
			System.err.println("Bytes: "+bytes);
		}else{
			for(long i=first; i<last; i++){
				byte[] s=tf.nextLine();
				if(s==null){break;}
				lines++;
				bytes+=s.length;
			}
		}
		t.stop();
		
		if(!reprint){
			double rpnano=lines/(double)(t.elapsed);
			double bpnano=bytes/(double)(t.elapsed);

			String rpstring=(lines<100000 ? ""+lines : lines<100000000 ? (lines/1000)+"k" : (lines/1000000)+"m");
			String bpstring=(bytes<100000 ? ""+bytes : bytes<100000000 ? (bytes/1000)+"k" : (bytes/1000000)+"m");

			while(rpstring.length()<8){rpstring=" "+rpstring;}
			while(bpstring.length()<8){bpstring=" "+bpstring;}

			System.err.println("Time:                         \t"+t);
			System.err.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk lines/sec", rpnano*1000000));
			System.err.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bytes/sec", bpnano*1000));
		}
	}
	
//	public ByteFile2(String name()){this(name(), false);}
	
	public ByteFile2(String fname, boolean tryAllExtensions, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess_, false), tryAllExtensions);
	}
	
	public ByteFile2(FileFormat ff, boolean tryAllExtensions){
		super(ff, tryAllExtensions);
		if(verbose){System.err.println("ByteFile2("+ff+", "+tryAllExtensions+")");}
		open();
	}
	
	public final void reset(){
		close();
		open();
	}
	
	public synchronized final boolean close(){
		if(verbose){System.err.println("ByteFile2("+name()+").close()");}
		if(isOpen()){
//			errorState|=ReadWrite.killProcess(name());
			thread.shutdown();
			while(thread.getState()!=Thread.State.TERMINATED){
				try {
					thread.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			thread.bf1.close();
		}
		thread=null;
		currentList=null;
		currentLoc=0;
//		assert(numIn==numOut) : numIn+", "+numOut;
		pushBack=null;
		if(verbose){System.err.println("ByteFile2("+name()+").close() returned "+errorState);}
		return errorState;
	}
	
	@Override
	public byte[] nextLine(){
		
		if(pushBack!=null){
			byte[] temp=pushBack;
			pushBack=null;
			return temp;
		}
//		if(verbose){System.err.println("Reading line.");}
//		byte[] r=null;
		
		byte[][] temp=currentList;
		int tempLoc=currentLoc;
		
		if(temp==null || tempLoc>=temp.length || temp[tempLoc]==null){
			boolean b=getBuffer();
			if(!b){
				if(verbose2){System.err.println("nextLine()->getBuffer() returned false.");}
				return null;
			}
			temp=currentList;
			tempLoc=currentLoc;
			if(temp==null || temp==poison || temp[tempLoc]==null){
				return null;
			}
		}
		
		//TODO: This is a race condition; currentList can be changed to null.  A defensive copy could be created.
		//Note that I read the above warning and added "temp" and "temploc" but I'm not sure if that fixed anything.
		assert(temp!=null && temp!=poison);
		assert(tempLoc<temp.length);
		assert(temp[tempLoc]!=null);
		byte[] r=temp[tempLoc];
		assert(r!=null);
		currentLoc++;
//		numOut++;
		return r;
	}
	
	private boolean getBuffer(){
		if(verbose2){System.err.println("Getting new buffer.");}
		currentLoc=0;
		final BF1Thread bft=thread;
		if(bft==null){
			currentList=null;
			if(verbose2){System.err.println("No buffers available.  thread="+thread+", shutdown="+(thread==null ? "X" : ""+thread.shutdown));}
			return false;
		}
		if(currentList==poison){
			if(verbose2){System.err.println("A: Current list is poison.");}
			return false;
		}
		if(currentList!=null){
			Arrays.fill(currentList, null); //MUST be done or lines get recycled at end of file.
			while(currentList!=null){
				try {
					if(verbose2){System.err.println("adding to qEmpty list size "+currentList.length+"\n"+toString(currentList));}
					bft.qEmpty.put(currentList);
					currentList=null;
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		assert(currentList==null);
		while(currentList==null){
			try {
				assert(bft!=null);
				if(verbose2){System.err.println("C: qFull.size()="+bft.qFull.size());}
				currentList=bft.qFull.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		if(verbose2){
			if(currentList==poison){
				System.err.println("B: Current list is poison.");
			}else{
				System.err.println("getBuffer fetched a new buffer of size "+currentList.length);
			}
		}
		return currentList!=poison;
	}
	
	private final synchronized BF1Thread open(){
		if(verbose2){System.err.println("ByteFile2("+name()+").open()");}
		assert(thread==null);
		currentList=null;
		currentLoc=0;
//		numIn=0;
//		numOut=0;
		thread=new BF1Thread(ff);
		thread.start();
		return thread;
	}
	
	private class BF1Thread extends Thread{
		
//		public BF1Thread(String fname){
//			bf1=new ByteFile1(fname, false, allowSubprocess);
//			qFull=new ArrayBlockingQueue<byte[][]>(buffs+2);
//			qEmpty=new ArrayBlockingQueue<byte[][]>(buffs+2);
//			for(int i=0; i<buffs; i++){
//				try {
//					qEmpty.put(new byte[bufflen][]);
//				} catch (InterruptedException e) {
//					// TODO Auto-generated catch block
//					e.printStackTrace();
//				}
//			}
//		}
		
		public BF1Thread(FileFormat ff){
			bf1=new ByteFile1(ff, false);
			qFull=new ArrayBlockingQueue<byte[][]>(buffs+2);
			qEmpty=new ArrayBlockingQueue<byte[][]>(buffs+2);
			for(int i=0; i<buffs; i++){
				try {
					qEmpty.put(new byte[bufflen][]);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		
		@Override
		public void run(){
			if(verbose){System.err.println("ByteFile2("+name()+").run()");}
			byte[] s=null;
			byte[][] list=null;
			while(list==null){
				try {
					list = qEmpty.take();
				} catch (InterruptedException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
			}
			synchronized(this){
				if(list==poison || shutdown){
					shutdown();
					return;
				}
			}
			
			int loc=0;
			long bases=0;
			
			//At this point, list is not null
			for(s=bf1.nextLine(); s!=null; s=bf1.nextLine()){
				bases+=s.length;
				assert(list!=null) : "Somehow the list became null for "+bf1.name()+" at line "+cntr;
				list[loc]=s;
				loc++;
//				numIn++;
//				if(verbose){System.err.println("Added line "+numIn);}
				if(loc>=bufflen || bases>=buffcapacity){
					if(verbose2){System.err.println("Capacity exceeded.");}
					while(list!=null){
						try {
//							synchronized(this){
//								if(!shutdown){
									if(verbose2){
										System.err.println("A: Adding to qFull list of size "+loc);
										System.err.println(ByteFile2.toString(list));
									}
									cntr+=list.length;
									qFull.put(list);
									if(verbose2){System.err.println("A: qFull.size()="+qFull.size());}
//								}
//							}
							list=null;
							loc=0;
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
					//At this point, list is null
					if(shutdown){
						if(verbose2){System.err.println("Break 1");}
						break;
					}
					while(list==null){
						if(verbose2){System.err.println("Taking empty list.");}
						try {
							list = qEmpty.take();
						} catch (InterruptedException e1) {
							// TODO Auto-generated catch block
							e1.printStackTrace();
						}
					}
					//At this point, list is not null
					bases=0;
					if(list==poison){
						if(verbose2){System.err.println("Break 2");}
						break;
					}
					//At this point, list is not null
				}
			}
			if(verbose2){System.err.println("Run loop exit.");}
			
			while(list!=null && loc>0){
				try {
//					synchronized(this){
//						if(!shutdown){
							if(verbose2){System.err.println("B: Adding list of size "+loc);}
							qFull.put(list);
							if(verbose2){System.err.println("B: qFull.size()="+qFull.size());}
//						}
//					}
					list=null;
					loc=0;
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			//At this point, list is null
			shutdown();
			
			if(verbose){System.err.println("ByteFile2("+name()+").run() finished");}
		}
		
		private synchronized void shutdown(){
			if(verbose || verbose2){System.err.println("ByteFile2("+name()+").shutdown()");}
			if(shutdown){return;}
			shutdown=true;
			if(verbose2){System.err.println("Adding poison.");}
			qFull.add(poison);
			qEmpty.add(poison);
			if(verbose2){System.err.println("D: qFull.size()="+qFull.size());}
			if(verbose || verbose2){System.err.println("ByteFile2("+name()+").shutdown() finished");}
		}
		
		private boolean shutdown=false;
		final ByteFile1 bf1;
		final ArrayBlockingQueue<byte[][]> qFull;
		final ArrayBlockingQueue<byte[][]> qEmpty;
		
	}
	
	public boolean isOpen(){
		if(currentList!=null && currentLoc<currentList.length && currentList[currentLoc]!=null){return true;}
		final BF1Thread bft=thread;
		if(bft==null){
			return false;
		}
		return true;
//		synchronized(bft){
//		//NOTE!!!  This cannot be used because qFull.size() will not return a correctly synchronized value.  Poll() may work.
//			assert(bft.bf1.isOpen() || !bft.qFull.isEmpty()) : bft.bf1.isOpen()+", "+bft.qFull.isEmpty()+", "+bft.qFull.size();
//			return (bft.bf1.isOpen() || !bft.qFull.isEmpty());
//		}
	}
	
	/** For debugging */
	private static String toString(byte[][] x){
		StringBuilder sb=new StringBuilder();
		for(byte[] z : x){
			sb.append(z==null ? "null" : new String(z)).append('\n');
		}
		return sb.toString();
	}
	
	public final InputStream is(){return thread==null ? null : thread.bf1.is();}
	
	public final long lineNum(){return thread==null ? -1 : thread.bf1.lineNum();}

	private long cntr;
	private BF1Thread thread=null;
	private byte[][] currentList=null;
	private int currentLoc=0;
//	private int currentSize=0;
	
//	private long numIn=0, numOut=0;
	
	private static final byte[][] poison=new byte[0][];
	public static boolean verbose=false;
	private static final boolean verbose2=false;
	private static final int bufflen=1000;
	private static final int buffs=4;
	private static final int buffcapacity=256000;
	
	private boolean errorState=false;
	
}

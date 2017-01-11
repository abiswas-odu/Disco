package stream;

import java.lang.management.ManagementFactory;
import java.lang.management.OperatingSystemMXBean;
import java.util.Arrays;
//import com.sun.management.OperatingSystemMXBean;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicIntegerArray;

import shared.Shared;

/**
 * Monitors CPU utilization to determine if the program has crashed.
 * @author Brian Bushnell
 * @date Feb 25, 2015
 *
 */
public final class KillSwitch extends Thread {
	
	public static void main(String[] args){
		double seconds=Double.parseDouble(args[0]);
		double load=Double.parseDouble(args[1]);
		launch(seconds, load);
		if(args.length>2){
			
		}
	}
	
	/**
	 * @param seconds
	 * @param load
	 */
	private KillSwitch(double seconds, double load) {
		maxSeconds=seconds;
		minLoad=load;
	}

	public static boolean launch(){
		return launch(600);
	}

	public static boolean launch(double seconds){
		return launch(seconds, 0.002);
	}
	
	public static synchronized boolean launch(double seconds, double load){
		if(count>0){return false;}
		ks=new KillSwitch(seconds, load);
		ks.start();
		return true;
	}
	
	@Override
	public void run(){
		
		boolean success=monitor();
//		System.err.println("success: "+success);
		if(!success || killFlag.get()){
			if(!suppressMessages){
				System.err.println("Process has decided it has crashed, and will abort.\n" +
						"If this decision was incorrect, please re-run with the flag 'monitor=f'"); 
			}
			kill0();
		}
	}
	
	private boolean monitor(){
		
		final OperatingSystemMXBean bean=ManagementFactory.getOperatingSystemMXBean();
		if(bean.getSystemLoadAverage()<0){
			System.err.println("This OS does not support monitor, so monitoring was disabled.");
			return true;
		}
		
		final long start=System.currentTimeMillis();
		final long buffer=(long)(1+maxSeconds*1000);
		long stop=start+buffer;
//		System.err.println("start="+start+", stop="+stop+", buffer="+buffer);
//		System.err.println("shutdownFlag.get()="+shutdownFlag.get());
		while(!shutdownFlag.get()){
			try {
				sleep(500);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			final double load=bean.getSystemLoadAverage();
			final long time=System.currentTimeMillis();
			if(load>minLoad){stop=time+buffer;}
			if(time>stop){return false;}
//			System.err.println("stop-time="+(stop-time)+", load="+load);
		}
//		System.err.println("shutdownFlag.get()="+shutdownFlag.get());
		return true;
	}

	/*--------------------------------------------------------------*/
	
	public static void kill(String s){
		ballast=null;
		Exception e=new Exception(s);
		e.printStackTrace();
		kill0();
	}
	
//	public static void kill(Throwable e){
//		e.printStackTrace();
//		kill0();
//	}
	
	public static void kill(){
		ballast=null;
		Exception e=new Exception("Aborting.");
		e.printStackTrace();
		kill0();
	}
	
	public static void killSilent(){
		ballast=null;
		kill0();
	}
	
	private static void kill0(){
		ballast=null;
		Runtime.getRuntime().halt(1);
	}
	
	public static void shutdown(){
		shutdownFlag.set(true);
	}
	
	public static void setKillFlag(){
		killFlag.set(true);
	}
	
	/*--------------------------------------------------------------*/
	
	public static final void memKill(OutOfMemoryError e){
		ballast=null;
		synchronized(MemKillMessage){
			e.printStackTrace();
			System.err.println(MemKillMessage);
//			Shared.printMemory();
//			killSilent();
			kill0();
		}
	}
	
	
	/*--------------------------------------------------------------*/
	
	public static final AtomicIntegerArray allocAtomicInt(int len){
		AtomicIntegerArray ret=null;
		try {
			ret=new AtomicIntegerArray(len);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
	public static final long[] allocLong1D(int len){
		long[] ret=null;
		try {
			ret=new long[len];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
	public static final int[] allocInt1D(int len){
		int[] ret=null;
		try {
			ret=new int[len];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
	public static final byte[] allocByte1D(int len){
		byte[] ret=null;
		try {
			ret=new byte[len];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
	public static final char[] allocChar1D(int len){
		char[] ret=null;
		try {
			ret=new char[len];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}

	/*--------------------------------------------------------------*/
	
	public static final int[][] allocInt2D(int x){
		int[][] ret=null;
		try {
			ret=new int[x][];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
	public static final int[][] allocInt2D(int x, int y){
		int[][] ret=null;
		try {
			ret=new int[x][y];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}

	/*--------------------------------------------------------------*/
	
	public static int[][][] allocInt3D(int x) {
		int[][][] ret=null;
		try {
			ret=new int[x][][];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
	public static int[][][] allocInt3D(int x, int y) {
		int[][][] ret=null;
		try {
			ret=new int[x][y][];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
	public static int[][][] allocInt3D(int x, int y, int z) {
		int[][][] ret=null;
		try {
			ret=new int[x][y][z];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}

	/*--------------------------------------------------------------*/
	
	public static byte[] copyOf(byte[] buffer, int newLength) {
		byte[] copy=null;
		try {
			copy=Arrays.copyOf(buffer, newLength);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return copy;
	}
	
	public static int[] copyOf(int[] buffer, int newLength) {
		int[] copy=null;
		try {
			copy=Arrays.copyOf(buffer, newLength);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return copy;
	}
	
	public static long[] copyOf(long[] buffer, int newLength) {
		long[] copy=null;
		try {
			copy=Arrays.copyOf(buffer, newLength);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return copy;
	}

	/*--------------------------------------------------------------*/
	
	public static byte[] copyOfRange(byte[] buffer, int start, int limit) {
		byte[] copy=null;
		try {
			copy=Arrays.copyOfRange(buffer, start, limit);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return copy;
	}

	public static int[] copyOfRange(int[] buffer, int start, int limit) {
		int[] copy=null;
		try {
			copy=Arrays.copyOfRange(buffer, start, limit);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return copy;
	}

	public static long[] copyOfRange(long[] buffer, int start, int limit) {
		long[] copy=null;
		try {
			copy=Arrays.copyOfRange(buffer, start, limit);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return copy;
	}
	
	/*--------------------------------------------------------------*/
	
	private final double maxSeconds;
	private final double minLoad;
	
	private static AtomicBoolean shutdownFlag=new AtomicBoolean(false);
	private static AtomicBoolean killFlag=new AtomicBoolean(false);
	private static int count=0;
	private static KillSwitch ks;
	private static boolean suppressMessages=false;
	
	/*--------------------------------------------------------------*/
	
	private final static String MemKillMessage=new String("\nThis program ran out of memory.\n"
			+ "Try increasing the -Xmx flag and using tool-specific memory-related parameters.");

	public static void addBallast() {
		synchronized(KillSwitch.class){
			ballast=new int[20000];
		}
	}
	
	private static int[] ballast;
	
}

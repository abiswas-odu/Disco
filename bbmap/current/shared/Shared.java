package shared;

import java.lang.management.ManagementFactory;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import dna.Data;
import stream.KillSwitch;

public class Shared {
	
	private static int THREADS=setThreads(-1);
	
	public static int READ_BUFFER_LENGTH=200;
	private static int READ_BUFFER_NUM_BUFFERS=setBuffers();
	public static long READ_BUFFER_MAX_DATA=400000;
	
	/** Temporary, for testing; should be made non-global */
	public static boolean AMINO_IN=false;
	
	//TODO:  For some reason, it seems as though GAPBUFFER must equal exactly 1/2 of GAPLEN.  Not good; 1/4 would be far better.
	
	public static final int GAPBUFFER=64; //TODO:  Seems to break less than 64, for some reason
	public static final int GAPBUFFER2=2*GAPBUFFER;
	public static final int GAPLEN=128; //TODO: May break when over 128
	public static final int MINGAP=GAPBUFFER2+GAPLEN;
	public static final int GAPCOST=Tools.max(1, GAPLEN/64);
	public static final byte GAPC='-';
	
	public static String BBMAP_VERSION_STRING="36.84";
	
	public static boolean TRIM_READ_COMMENTS=false;
	
	public static boolean USE_JNI=false;//Data.GENEPOOL;
	public static boolean USE_MPI=false;
	public static boolean MPI_KEEP_ALL=true;
	/** Use ConcurrentReadInputStreamMPI instead of D */
	public static boolean USE_CRISMPI=true;
	public static int MPI_RANK=0;
	public static int MPI_NUM_RANKS=1;
	
	public static int FASTA_WRAP=70;
	public static byte FAKE_QUAL=30;
	
	public static String BBMAP_CLASS=null;
	public static String[] COMMAND_LINE=null;
	public static List<String> JVM_ARGS(){
		return ManagementFactory.getRuntimeMXBean().getInputArguments();
	}
	
	/** Directory in which to write temp files */
	private static String TMPDIR=(System.getenv("TMPDIR")==null ? null : (System.getenv("TMPDIR")+"/").replaceAll("//", "/"));
//	static{assert(false) : "TMPDIR="+TMPDIR;}
	
	public static String tmpdir(){return TMPDIR;}
	
	public static String setTmpdir(String s){
		if(s==null){TMPDIR=null;}
		s=s.replaceAll("\\", "/");
		if(!s.endsWith("/")){s=s+"/";}
		TMPDIR=s.replaceAll("//", "/");
		return TMPDIR;
	}
	
	/** Anomaly probably resolved as of v.20.1 
	 * This variable should be TRUE for normal users and FALSE for me. */
	public static boolean anomaly=!(System.getProperty("user.dir")+"").contains("/bushnell/") && !Data.WINDOWS;
	
	public static final char[] getTLCB(int len){
		char[] buffer=TLCB.get();
		if(buffer==null || buffer.length<len){
			buffer=new char[len];
			if(len<1000000){TLCB.set(buffer);}
		}
		return buffer;
	}
	private static final ThreadLocal<char[]> TLCB=new ThreadLocal<char[]>();
	
	public static int setThreads(String x){
		int y=Data.LOGICAL_PROCESSORS;
		if(x!=null && !x.equalsIgnoreCase("auto")){
			y=Integer.parseInt(x);
		}
		return setThreads(y);
	}
	
	public static int setThreads(int x){
		if(x>0){
			THREADS=x;
		}else{
			THREADS=Tools.max(1, Data.LOGICAL_PROCESSORS);
		}
		setBuffers();
		return THREADS;
	}
	
	public static int threads(){
		assert(THREADS>0);
		return THREADS;
	}
	
	public static int capBuffers(int num){
		return setBuffers(Tools.min(num, READ_BUFFER_NUM_BUFFERS));
	}

	public static int READ_BUFFER_NUM_BUFFERS() {
		return READ_BUFFER_NUM_BUFFERS;
	}
	
	public static int setBuffers(){
		return setBuffersFromThreads(THREADS);
	}
	
	public static int setBuffersFromThreads(int threads){
		return setBuffers(Tools.max(4, (threads*3)/2));
	}
	
	public static int setBuffers(int num){
		num=Tools.max(2, num);
		return READ_BUFFER_NUM_BUFFERS=num;
	}
	
	public static int numBuffers(){
		return READ_BUFFER_NUM_BUFFERS;
	}
	
	public static boolean LOW_MEMORY=false;
	
	/** Ratio of -Xms to -Xmx parameters */
	public static final double xmsRatio(){
		Runtime rt=Runtime.getRuntime();
		return rt.totalMemory()*1.0/rt.maxMemory();
	}
	
	public static long memAvailable(int readThreads){
		long usableMemory;
		{
			long memory=Runtime.getRuntime().maxMemory();
			double xmsRatio=Shared.xmsRatio();
			usableMemory=(long)Tools.max(((memory-48000000-(Tools.max(readThreads, 4)*400000))*(xmsRatio>0.97 ? 0.82 : 0.72)), memory*0.45);
		}
		return usableMemory;
	}
	
	public static long memFree(){
		Runtime rt=Runtime.getRuntime();
		return rt.freeMemory();
	}
	
	public static long memUsed(){
		Runtime rt=Runtime.getRuntime();
		return rt.maxMemory()-rt.freeMemory();
	}
	
	/** Print statistics about current memory use and availability */
	public static final void printMemory(){
		try{
			if(GC_BEFORE_PRINT_MEMORY){
				System.gc();
				System.gc();
			}
			Runtime rt=Runtime.getRuntime();
			long mmemory=rt.maxMemory()/1000000;
			long tmemory=rt.totalMemory()/1000000;
			long fmemory=rt.freeMemory()/1000000;
			long umemory=tmemory-fmemory;
			System.err.println("Memory: "+"max="+mmemory+/*"m, total="+tmemory+*/"m, "+"free="+fmemory+"m, used="+umemory+"m");
		}catch(Throwable t){}
	}
	
	/** Do garbage collection prior to printing memory usage */
	private static final boolean GC_BEFORE_PRINT_MEMORY=false;
	

	
	/*--------------------------------------------------------------*/
	/*----------------            Java 8            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final void sort(int[] array){sort(array, 0, array.length);}
	public static final void sort(int[] array, int from, int to){
		try {
			if(!parallelSort || array.length<=parallelSortLength){
				Arrays.sort(array, from, to); //Supported in pre-Java 8.
				return;
			}
			Arrays.parallelSort(array, from, to);
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
	}

	public static final void sort(long[] array){sort(array, 0, array.length);}
	public static final void sort(long[] array, int from, int to){
		try {
			if(!parallelSort || array.length<=parallelSortLength){
				Arrays.sort(array, from, to); //Supported in pre-Java 8.
				return;
			}
			Arrays.parallelSort(array, from, to);
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
	}

	public static final <T extends Comparable<? super T>> void sort(T[] array){sort(array, 0, array.length);}
	public static final <T extends Comparable<? super T>> void sort(T[] array, int from, int to){
		try {
			if(!parallelSort || array.length<=parallelSortLength || THREADS<2){
				Arrays.sort(array, from, to); //Supported in pre-Java 8.
				return;
			}
			Arrays.parallelSort(array, from, to);
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
	}
	
	public static final <T extends Comparable<? super T>> void sort(ArrayList<T> list){
		try {
			if(!parallelSort || list.size()<=parallelSortLength){
				Collections.sort(list); //Supported in pre-Java 8.
				return;
			}
			
			{//If this block causes compile errors, just replace the whole function body with "Shared.sort(list, comparator);"
				@SuppressWarnings("unchecked")
				T[] array=list.toArray((T[])new Comparable[0]);
				list.clear();
				Arrays.parallelSort(array);
				for(T r : array){list.add(r);}
			}
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
	}
	
	public static final <T> void sort(ArrayList<T> list, Comparator<? super T> comparator){
		try {
			if(!parallelSort){
				Collections.sort(list, comparator); //Supported in pre-Java 8.
				return;
			}
			
			{//If this block causes compile errors, just replace the whole function body with "Shared.sort(list, comparator);"
				if(list.size()<=parallelSortLength || THREADS<2){
					list.sort(comparator);
					return;
				}
				@SuppressWarnings("unchecked")
				T[] array=list.toArray((T[])new Comparable[0]);
				list.clear();
				Arrays.parallelSort(array, comparator);
				for(T r : array){list.add(r);}
			}
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
	}
	
	public static final int parallelSortLength=10000;
	public static boolean parallelSort=testParallelSort();
//	public static boolean listSort=testListSort();
	
	public static double javaVersion=parseJavaVersion();
	
	private static double parseJavaVersion(){
		String s=System.getProperty("java.version");
		if(s==null){return 1.6;}
		int dots=0;
		StringBuilder sb=new StringBuilder();
		for(int i=0; i<s.length() && dots<2; i++){
			char c=s.charAt(i);
			if(c=='.'){dots++;}
			else if(!Character.isDigit(c)){break;}
			if(dots>1){break;}
			sb.append(c);
		}
		return Double.parseDouble(sb.toString());
	}
	
	private static boolean testParallelSort(){
		Method m=null;
		try {
			m=Arrays.class.getMethod("parallelSort", new Class[] {Object[].class, Comparator.class});
		} catch (NoSuchMethodException e) {
			// TODO Auto-generated catch block
//			e.printStackTrace();
		} catch (SecurityException e) {
			// TODO Auto-generated catch block
//			e.printStackTrace();
		} catch (Throwable t) {
			// TODO Auto-generated catch block
//			e.printStackTrace();
		}
		return m!=null;
	}
	
	static{
		KillSwitch.addBallast();
	}
	
}

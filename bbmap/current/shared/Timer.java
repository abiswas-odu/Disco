package shared;

import java.io.PrintStream;

public class Timer {
	
	public Timer(){start();}
	
	public long start(String s){
		outstream.println(s);
		return start();
	}
	
	public long stop(String s){
		long x=stop();
		if(addTab && s!=null && !s.endsWith("\t")){s=s+"\t";}
		outstream.println(s+this);
		return x;
	}
	
	public long start(){
		time1=time2=System.nanoTime();
		elapsed=0;
		return time1;
	}
	
	public long stop(){
		time2=System.nanoTime();
		elapsed=time2-time1;
		return time2;
	}
	
	public String toString(){
		return String.format("%.3f seconds.", elapsed/1000000000d);
	}

	public long time1;
	public long time2;
	/** in nanos */
	public long elapsed;
	
	public PrintStream outstream=System.err;
	public boolean addTab=true;
}

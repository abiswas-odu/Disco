package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Random;

import assemble.Tadpole;
import dna.Parser;
import shared.Shared;
import shared.Tools;

/**
 * Wraps BBNorm and Tadpole to normalize and correct reads.
 * @author Brian Bushnell
 * @date Oct 15, 2015
 *
 */
public class NormAndCorrectWrapper {
	
	public static void main(String[] args){
//		Timer t=new Timer();
//		t.start();
		NormAndCorrectWrapper rr=new NormAndCorrectWrapper(args);
		rr.process();
	}
	
	public NormAndCorrectWrapper(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			System.exit(0);
		}
		
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
//		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens
			
			if(a.equals("in")){
				in=b;
			}else if(a.equals("out")){
				out=b;
			}else if(a.equals("correctfirst")){
				correctFirst=Tools.parseBoolean(b);
			}else if(a.equals("ow") || a.equals("overwrite")){
				out=b;
			}else{
				throw new RuntimeException("Unsupported argument");
			}
		}
		
	}
	
	public void process(){
		
		Random randy=new Random();
		
		String temp;
		if(correctFirst){
			temp=Shared.tmpdir()+"corrected_"+((randy.nextLong()&Long.MAX_VALUE)^in.hashCode())+".fq.gz";
			String[] tadArgs=new String[] {"in="+in, "out="+temp, "mode=correct", "pigz", "unpigz", "ow="+overwrite}; 
			Tadpole.main(tadArgs);
			System.gc();
			String[] normArgs=new String[] {"in="+temp, "out="+out, "bits=32", "min=2", "target=100", "pigz", "unpigz", "ow="+overwrite}; 
			KmerNormalize.main(normArgs);
		}else{
			temp=Shared.tmpdir()+"normalized_"+((randy.nextLong()&Long.MAX_VALUE)^in.hashCode())+".fq.gz";
			String[] normArgs=new String[] {"in="+in, "out="+temp, "bits=32", "min=2", "target=100", "pigz", "unpigz", "ow="+overwrite}; 
			KmerNormalize.main(normArgs);
			System.gc();
			String[] tadArgs=new String[] {"in="+temp, "out="+out, "mode=correct", "pigz", "unpigz", "ow="+overwrite}; 
			Tadpole.main(tadArgs);
		}
		
		File f=new File(temp);
		if(f.exists()){
			f.delete();
		}
	}
	
	public PrintStream outstream=System.err;
	public String in="reads.fq.gz", out="corrected.fq.gz";
	public boolean overwrite=true;
	public boolean correctFirst=false;
	
}

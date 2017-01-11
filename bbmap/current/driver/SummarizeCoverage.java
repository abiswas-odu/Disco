package driver;

import java.io.File;
import java.util.ArrayList;

import dna.Parser;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Apr 29, 2015
 *
 */
public class SummarizeCoverage {
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			assert(false);
			System.exit(0);
		}
		
		//Create a new Seal instance
		SummarizeCoverage sc=new SummarizeCoverage(args);
		
		///And run it
		sc.process();
	}
	
	public SummarizeCoverage(String[] args){
		
		Parser parser=new Parser();
		
		ArrayList<String> names=new ArrayList<String>();
		
		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if("null".equalsIgnoreCase(b)){b=null;}
			while(a.charAt(0)=='-' && (a.indexOf('.')<0 || i>1 || !new File(a).exists())){a=a.substring(1);}
			
			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(!arg.contains("=")){
				String[] x=(new File(arg).exists() ? new String[] {arg} : arg.split(","));
				for(String x2 : x){names.add(x2);}
			}else{
				throw new RuntimeException("Unknown parameter "+arg);
			}
		}
		
		{//Process parser fields
			out=(parser.out1==null ? "stdout" : parser.out1);
			if(parser.in1!=null){
				String[] x=(new File(parser.in1).exists() ? new String[] {parser.in1} : parser.in1.split(","));
				for(String x2 : x){names.add(x2);}
			}
		}

		in=new ArrayList<String>();
		for(String s : names){
			Tools.getFileOrFiles(s, in, false, false, false, true);
		}
	}
	
	public void process(){
		TextStreamWriter tsw=new TextStreamWriter(out, true, false, false);
		tsw.start();
		tsw.print("#File\tPrimary_Name\tPrimary_Count\tOther_Count\tPrimary_MB\tOther_MB\n");
		for(String fname : in){
			String pname=null;
			long pcount=0, ocount=0;
			double pmb=0, omb=0;
			TextFile tf=new TextFile(fname);
			for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
				if(!line.startsWith("#")){
					String[] split=line.split("\t");
					long count=Long.parseLong(split[5]);
					double mb=Double.parseDouble(split[2]);
					if(pcount==0 || mb>pmb || (mb==pmb && count>pcount)){
						pname=split[0];
						ocount+=pcount;
						omb+=pmb;
						pcount=count;
						pmb=mb;
					}else{
						ocount+=count;
						omb+=mb;
					}
				}
			}
			tf.close();
			tsw.print(String.format("%s\t%s\t%d\t%d\t%.5f\t%.5f\n", fname, pname, pcount, ocount, pmb, omb));
		}
		tsw.poisonAndWait();
	}
	
	final ArrayList<String> in;
	final String out;
	
}

package driver;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import dna.Parser;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import jgi.DecontaminateByNormalization;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date June 1, 2016
 *
 */
public class SummarizeCrossblock {
	
	public static void main(String[] args){
		Timer t=new Timer();
		SummarizeCrossblock sc=new SummarizeCrossblock(args);
		sc.process(t);
	}
	
	public SummarizeCrossblock(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		for(String s : args){if(s.startsWith("out=standardout") || s.startsWith("out=stdout")){outstream=System.err;}}
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ReadWrite.verbose=verbose;
			}else if(parser.in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				parser.in1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		String in=null;
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=parser.overwrite;
			append=parser.append;
			
			in=parser.in1;

			out1=parser.out1;
		}
		
		if(in==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		
		if(in.contains(",")){
			for(String s : in.split(",")){
				inList.add(s);
			}
		}else{
			inList.add(in);
			DecontaminateByNormalization.parseStringsFromFiles(inList);
		}
		
		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TEXT, null, true, overwrite, append, false);
	}
	
	void process(Timer t){
		final TextStreamWriter tsw;
		tsw=ffout1!=null ? new TextStreamWriter(ffout1) : null;
		if(tsw!=null){tsw.start();}
		if(tsw!=null){tsw.print("#fname\tcopies\tcontigs\tcontigsDiscarded\tbases\tbasesDiscarded\n");}
		
		int i=1;
		for(String fname : inList){
			ParseCrossblockResults pcr=null;
			try{
				pcr=new ParseCrossblockResults(new String[] {"in="+fname});
				Timer t2=new Timer();
				pcr.process(t2);
				tsw.print(fname+"\t"+i+"\t"+pcr.contigs()+"\t"+pcr.contigsDiscarded()+"\t"+pcr.bases()+"\t"+pcr.basesDiscarded()+"\n");
			}catch(Throwable e){
				System.err.println(e);
				tsw.print(fname+"\tERROR\n");
			}
			i++;
		}
		if(tsw!=null){errorState|=tsw.poisonAndWait();}
	}

	/*--------------------------------------------------------------*/
	
	private void printOptions(){assert(false) : "printOptions: TODO";}
	
	/*--------------------------------------------------------------*/
	
	private ArrayList<String> inList=new ArrayList<String>();
	private String out1=null;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;

	private long basesKept=0;
	private long basesDiscarded=0;
	private long contigsKept=0;
	private long contigsDiscarded=0;

	public long basesKept(){return basesKept;}
	public long basesDiscarded(){return basesDiscarded;}
	public long contigsKept(){return contigsKept;}
	public long contigsDiscarded(){return contigsDiscarded;}
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffout1;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	
}

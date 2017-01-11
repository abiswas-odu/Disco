package driver;

import java.io.File;
import java.io.PrintStream;
import java.util.Arrays;

import dna.Parser;
import fileIO.TextFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date June 9, 2016
 *
 */
public class ParseCrossblockResults {
	
	public static void main(String[] args){
		Timer t=new Timer();
		ParseCrossblockResults pcr=new ParseCrossblockResults(args);
		pcr.process(t);
	}
	
	public ParseCrossblockResults(String[] args){
		
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
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=parser.overwrite;
			append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
		}
		
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TEXT, null, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.TEXT, null, false, false);
	}
	
	void process(Timer t){
		
		final TextFile tf;
		{
			tf=new TextFile(ffin1);
			if(verbose){outstream.println("Started tf");}
		}
		
		long linesProcessed=0;
		long charsProcessed=0;
		
		{
			String line;
			while((maxReads<0 || linesProcessed<maxReads) && (line=tf.nextLine())!=null){
				linesProcessed++;
				charsProcessed+=line.length();
				if(!line.startsWith("#")){
					processLine(line);
				}
			}
		}
		errorState|=tf.close();

		if(ffout1!=null){
			final TextStreamWriter tsw;
			{
				tsw=new TextStreamWriter(ffout1);
				tsw.start();
				if(verbose){outstream.println("Started tsw");}
				errorState|=tsw.poisonAndWait();
			}
		}
		
		t.stop();
		
		double rpnano=linesProcessed/(double)(t.elapsed);
		double bpnano=charsProcessed/(double)(t.elapsed);

		String rpstring=(linesProcessed<100000 ? ""+linesProcessed : linesProcessed<100000000 ? (linesProcessed/1000)+"k" : (linesProcessed/1000000)+"m");
		String bpstring=(charsProcessed<100000 ? ""+charsProcessed : charsProcessed<100000000 ? (charsProcessed/1000)+"k" : (charsProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		outstream.println("Time:                         \t"+t);
		outstream.println("Lines Processed:    "+rpstring+" \t"+String.format("%.2fk lines/sec", rpnano*1000000));
		outstream.println("Chars Processed:    "+bpstring+" \t"+String.format("%.2fm chars/sec", bpnano*1000));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	private void processLine(String line){
		ResultsLine rl=new ResultsLine(line);
		if(rl.removed){
			basesDiscarded+=rl.length;
			contigsDiscarded++;
		}else{
			basesKept+=rl.length;
			contigsKept++;
		}
	}
	
	
	
	/*--------------------------------------------------------------*/
	
	private static class ResultsLine{
		
		public ResultsLine(String s){
			String[] split=s.split("\t");
			length=Integer.parseInt(split[3]);
			removed=Integer.parseInt(split[2])==1;
		}
		
		final int length;
		final boolean removed;
		
	}
	
	/*--------------------------------------------------------------*/
	
	private void printOptions(){assert(false) : "printOptions: TODO";}
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
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
	
	public long contigs(){return contigsKept+contigsDiscarded;}
	public long bases(){return basesKept+basesDiscarded;}
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	
}

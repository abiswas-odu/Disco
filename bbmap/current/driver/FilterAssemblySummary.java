package driver;

import java.io.File;
import java.io.PrintStream;
import java.util.Arrays;

import dna.Parser;
import fileIO.TextFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import tax.TaxFilter;

/**
 * @author Brian Bushnell
 * @date May 1, 2016
 *
 */
public class FilterAssemblySummary {
	
	public static void main(String[] args){
		Timer t=new Timer();
		FilterAssemblySummary mb=new FilterAssemblySummary(args);
		mb.process(t);
	}
	
	public FilterAssemblySummary(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		for(String s : args){if(s.startsWith("out=standardout") || s.startsWith("out=stdout")){outstream=System.err;}}
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");

		
		
		Shared.READ_BUFFER_LENGTH=Tools.min(200, Shared.READ_BUFFER_LENGTH);
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		
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
			}else if(TaxFilter.validArgument(a)){
				//do nothing
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

		ffin1=FileFormat.testInput(in1, FileFormat.TEXT, null, true, true);
		
		//Make the actual filter
		filter=TaxFilter.makeFilter(args);
	}
	
	void process(Timer t){
		
		final TextFile tf;
		{
			tf=new TextFile(ffin1);
			if(verbose){outstream.println("Started tf");}
		}
		
		final TextStreamWriter tsw;
		{
			tsw=new TextStreamWriter(ffout1);
			tsw.start();
			if(verbose){outstream.println("Started tsw");}
		}

		long linesProcessed=0;
		long linesRetained=0;
		long charsProcessed=0;
		
		{
			String line;
			while((maxReads<0 || linesProcessed<maxReads) && (line=tf.nextLine())!=null){
				linesProcessed++;
				charsProcessed+=line.length();
				String result=processLine(line);
				if(result!=null){
					linesRetained++;
					if(tsw!=null){tsw.println(result);}
				}
			}
		}
		
		errorState|=tsw.poisonAndWait();
		errorState|=tf.close();
		
		t.stop();
		
		double rpnano=linesProcessed/(double)(t.elapsed);
		double bpnano=charsProcessed/(double)(t.elapsed);

		String rpstring=(linesProcessed<100000 ? ""+linesProcessed : linesProcessed<100000000 ? (linesProcessed/1000)+"k" : (linesProcessed/1000000)+"m");
		String kpstring=(linesRetained<100000 ? ""+linesRetained : linesRetained<100000000 ? (linesRetained/1000)+"k" : (linesRetained/1000000)+"m");
		String bpstring=(charsProcessed<100000 ? ""+charsProcessed : charsProcessed<100000000 ? (charsProcessed/1000)+"k" : (charsProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(kpstring.length()<8){kpstring=" "+kpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		outstream.println("Time:                         \t"+t);
		outstream.println("Lines Processed:    "+rpstring+" \t"+String.format("%.2fk lines/sec", rpnano*1000000));
		outstream.println("Lines Retained:     "+kpstring);
		outstream.println("Chars Processed:    "+bpstring+" \t"+String.format("%.2fm chars/sec", bpnano*1000));
		
		
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	private String processLine(String line){
//		System.out.println("Processing line "+line);
		if(line.startsWith("#")){return null;}
		String[] split=line.split("\t");
		assert(split.length>6) : split.length+"\n"+"'"+line+"'";
		String id=split[6];
		int number=Integer.parseInt(id);
//		System.out.println("Found number "+number+" from "+id+"; node: "+filter.tree().getNode(number));
		boolean b=filter.passesFilter((int)number);
//		System.out.println("passesFilter? "+b);
		return b ? line : null;
	}
	
	
	/*--------------------------------------------------------------*/
	
	private void printOptions(){assert(false) : "printOptions: TODO";}
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	
	/** The actual filter */
	private final TaxFilter filter;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	
}

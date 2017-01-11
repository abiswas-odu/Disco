package driver;

import java.io.File;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;

import dna.Parser;
import fileIO.TextFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Oct 17, 2014
 *
 */
public class MergeBigelow {
	
	public static void main(String[] args){
		Timer t=new Timer();
		MergeBigelow mb=new MergeBigelow(args);
		mb.process(t);
	}
	
	public MergeBigelow(String[] args){
		
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
			in2=parser.in2;

			out1=parser.out1;
		}
		
		if(in1==null || in2==null){
			printOptions();
			throw new RuntimeException("Error - two input files are required.");
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}

		assert(Tools.testInputFiles(false, true, in1, in2));
		assert(Tools.testForDuplicateFiles(true, in1, in2, out1));
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TEXT, null, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.TEXT, null, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.TEXT, null, true, true);
	}
	
	void process(Timer t){
		
		table=hash(ffin2);
		
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
		long charsProcessed=0;
		
		{
			String line;
			while((line=tf.nextLine())!=null){
//				System.err.println("Processing "+line);
				linesProcessed++;
				charsProcessed+=line.length();
				CharSequence result=processLine(line);
				if(tsw!=null && result!=null){tsw.println(result);}
				if(maxReads>0 && linesProcessed>=maxReads){break;}
			}
		}
		
		errorState|=tsw.poisonAndWait();
		errorState|=tf.close();
		
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
	
	
	private CharSequence processLine(String line){
		String[] split=line.split(delimiter);
		String[] split2=table.get(split[0]);
		if(split2==null){return line;} //Header
		StringBuilder sb=new StringBuilder();
		String tab="";
//		assert(false) : split.length+", "+split2.length;
//		System.err.println(split[1]);
		if(split.length>1){
			if(split[1].contains(" SCGC")){
				split[1]=split[1].substring(0, split[1].indexOf(" SCGC"));
//				System.err.println(split[1]);
			}
			if(split[1].contains(" "+split[0])){
				split[1]=split[1].substring(0, split[1].indexOf(" "+split[0]));
//				System.err.println(split[1]);
			}
			split[1]=split[1].toLowerCase();
//			System.err.println(split[1]);
		}
		for(int i=0; i<split.length; i++){
			sb.append(tab);
			sb.append(split[i].replace(',','_'));
			tab="\t";
		}
		for(int i=1; i<split2.length; i++){
			sb.append(tab);
			sb.append(split2[i].replace(',','_'));
			tab="\t";
		}
		return sb;
	}
	
	private HashMap<String, String[]> hash(FileFormat ff){
		final HashMap<String, String[]> table=new HashMap<String, String[]>();
		final TextFile tf;
		{
			tf=new TextFile(ff);
			if(verbose){outstream.println("Started tf");}
		}
		{
			String line;
			while((line=tf.nextLine())!=null){
				String[] split=line.split(delimiter);
				table.put(split[0], split);
			}
		}
		return table;
	}
	
	
	/*--------------------------------------------------------------*/
	
	private void printOptions(){assert(false) : "printOptions: TODO";}
	
	/*--------------------------------------------------------------*/

	private String in1=null;
	private String in2=null;
	private String out1=null;
	
	private String delimiter="\t";
	private HashMap<String, String[]> table; 
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/

	private final FileFormat ffin1;
	private final FileFormat ffin2;
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	
}

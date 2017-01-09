package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import dna.Parser;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Oct 3, 2014
 *
 */
public class PhylipToFasta {
	


	public static void main(String[] args){
		Timer t=new Timer();
		PhylipToFasta mb=new PhylipToFasta(args);
		mb.process(t);
	}
	
	public PhylipToFasta(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		for(String s : args){if(s.startsWith("out=standardout") || s.startsWith("out=stdout")){outstream=System.err;}}
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		
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
			}else if(parser.out1==null && i==1 && !arg.contains("=")){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
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
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTA, null, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.PHYLIP, ".phylip", true, true);
	}
	
	void process(Timer t){
		
		ArrayList<StringBuilder> data=new ArrayList<StringBuilder>();
		long bases=0;
		
		{
			final TextFile tf=new TextFile(ffin1);
			String s=tf.nextLine(); //first line is some numbers

			for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
				if(s.startsWith("        ")){break;}
				StringBuilder sb=new StringBuilder();
				data.add(sb);
				sb.append('>');
				int pos=0;
				for(pos=0; pos<s.length(); pos++){
					char c=s.charAt(pos);
					if(Character.isWhitespace(c)){break;}
					sb.append(c);
				}
				sb.append('\n');
				while(pos<s.length() && Character.isWhitespace(s.charAt(pos))){pos++;}
				while(pos<s.length()){
					char c=s.charAt(pos);
					if(Character.isLetter(c)){
						sb.append(c);
						bases++;
					}
					pos++;
				}
			}

			final int mod=data.size();
			for(int i=0; s!=null; i++){
				StringBuilder sb=data.get(i%mod);
				for(int pos=0; pos<s.length(); pos++){
					char c=s.charAt(pos);
					if(Character.isLetter(c)){
						sb.append(c);
						bases++;
					}
					pos++;
				}
				s=tf.nextLine();
			}
			errorState|=tf.errorState;
		}
		final long reads=data.size();
		
		if(ffout1!=null){
			TextStreamWriter tsw=new TextStreamWriter(ffout1);
			tsw.start();
			for(int i=0; i<data.size(); i++){
				StringBuilder sb=data.set(i, null);
				sb.append('\n');
				tsw.print(sb);
			}
			tsw.poisonAndWait();
			errorState|=tsw.errorState;
		}
		
		t.stop();
		
		double rpnano=reads/(double)(t.elapsed);
		double bpnano=bases/(double)(t.elapsed);

		String rpstring=(reads<100000 ? ""+reads : reads<100000000 ? (reads/1000)+"k" : (reads/1000000)+"m");
		String bpstring=(bases<100000 ? ""+bases : bases<100000000 ? (bases/1000)+"k" : (bases/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	private void printOptions(){
		assert(false) : "printOptions: TODO";
//		outstream.println("Syntax:\n");
//		outstream.println("java -ea -Xmx512m -cp <path> jgi.ReformatReads in=<infile> in2=<infile2> out=<outfile> out2=<outfile2>");
//		outstream.println("\nin2 and out2 are optional.  \nIf input is paired and there is only one output file, it will be written interleaved.\n");
//		outstream.println("Other parameters and their defaults:\n");
//		outstream.println("overwrite=false  \tOverwrites files that already exist");
//		outstream.println("ziplevel=4       \tSet compression level, 1 (low) to 9 (max)");
//		outstream.println("interleaved=false\tDetermines whether input file is considered interleaved");
//		outstream.println("fastawrap=70     \tLength of lines in fasta output");
//		outstream.println("qin=auto         \tASCII offset for input quality.  May be set to 33 (Sanger), 64 (Illumina), or auto");
//		outstream.println("qout=auto        \tASCII offset for output quality.  May be set to 33 (Sanger), 64 (Illumina), or auto (meaning same as input)");
//		outstream.println("outsingle=<file> \t(outs) Write singleton reads here, when conditionally discarding reads from pairs.");
	}
	
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;

	private String out1=null;
	
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

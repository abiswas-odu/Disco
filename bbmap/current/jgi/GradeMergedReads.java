package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import dna.Parser;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date May 20, 2014
 *
 */
public class GradeMergedReads {


	public static void main(String[] args){
		Timer t=new Timer();
		GradeMergedReads gmr=new GradeMergedReads(args);
		gmr.process(t);
	}
	
	public GradeMergedReads(String[] args){
		
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
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		FASTQ.DETECT_QUALITY=false;
		
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
			}else if(a.equals("raw") || a.equals("raw1")){
				raw1=b;
				if(b.indexOf("#")>=0){
					raw1=b.replace('#', '1');
					raw2=b.replace('#', '2');
				}
			}else if(a.equals("raw2")){
				raw2=b;
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
//				align2.FastaReadInputStream2.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(parser.in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				parser.in1=arg;
			}else{
				System.err.println("Unknown parameter "+i+": "+args[i]);
				assert(false) : "Unknown parameter "+i+": "+args[i];
//					+"\n"+arg+", "+parser.in1+", "+arg.contains("=")+", "+(arg.toLowerCase().startsWith("stdin")+", "+new File(arg).exists()+", "+new File(arg).getAbsolutePath());
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			in=parser.in1;
			extin=parser.extin;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		FASTQ.PARSE_CUSTOM=false;

		ffin=FileFormat.testInput(in, FileFormat.FASTQ, extin, true, true);
	}
	
	void process(Timer t){
		
		long mergeable=0, total=0;
		if(raw1!=null){
			FileFormat ffraw1=FileFormat.testInput(raw1, FileFormat.FASTQ, extin, true, true);
			FileFormat ffraw2=FileFormat.testInput(raw2, FileFormat.FASTQ, extin, true, true);
			ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffraw1, ffraw2, null, null);
			cris.start();
			
			{
				
				ListNum<Read> ln=cris.nextList();
				ArrayList<Read> reads=(ln!=null ? ln.list : null);

				while(reads!=null && reads.size()>0){
					for(int idx=0; idx<reads.size(); idx++){
						final Read r1=reads.get(idx);
						String s=r1.id;
						total++;
						final int insert=parseInsert(r1.id);
						if(insert>0 && insert<r1.length()+r1.mateLength()){
							mergeable++;
						}
					}

					cris.returnList(ln.id, ln.list.isEmpty());
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
				errorState|=ReadWrite.closeStream(cris);
			}
		}
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin, null, null, null);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
		}
		
		long readsProcessed=0;
		long basesProcessed=0;

		long correct=0;
		long tooLong=0;
		long tooShort=0;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			System.err.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin==null || ffin.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(reads!=null && reads.size()>0){
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					
					final int initialLength1=r1.length();
					final int insert=parseInsert(r1.id);
					
					int delta=insert-initialLength1;
					if(delta==0){
						correct++;
					}else if(delta>0){
						tooLong++;
					}else{
						tooShort++;
					}
					
					readsProcessed++;
					basesProcessed+=initialLength1;
				}

				cris.returnList(ln.id, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadWrite.closeStream(cris);
		
		t.stop();
		
		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);

		String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
		String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		long incorrect=tooShort+tooLong;
		double snr=10*Math.log10((correct+incorrect+0.0001)/(incorrect+0.0001));
		
		if(total>0){
			outstream.println("Input Total:            \t"+total+" pairs");
			outstream.println("Input Overlapping:      \t"+String.format("%.5f",mergeable*100.0/total)+"%\t"+mergeable+" reads");
		}
		outstream.println("Correct:                \t"+String.format("%.5f",correct*100.0/readsProcessed)+"%\t"+correct+" reads");
		outstream.println("Incorrect:              \t"+String.format("%.5f",incorrect*100.0/readsProcessed)+"%\t"+incorrect+" reads");
		outstream.println("Too Short:              \t"+String.format("%.5f",tooShort*100.0/readsProcessed)+"%\t"+tooShort+" reads");
		outstream.println("Too Long:               \t"+String.format("%.5f",tooLong*100.0/readsProcessed)+"%\t"+tooLong+" reads");
		outstream.println("SNR:                    \t"+String.format("%.3f",snr));

		outstream.println();
		
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		
		if(errorState){
			throw new RuntimeException("GradeMergedReads terminated in an error state; the output may be corrupt.");
		}
	}
	
	static int parseInsert(String s){
//		int space=s.indexOf(' ');
//		if(space<0){space=s.length();}
		int space=s.length();
		int equals=s.indexOf('=');
		for(int i=equals+1;i<s.length();i++){//For programs that rename my reads!
			if(!Character.isDigit(s.charAt(i))){
				space=i;
				break;
			}
		}
		s=s.substring(equals+1, space);
		int insert=Integer.parseInt(s);
		return insert;
	}
	
	/*--------------------------------------------------------------*/
	
	private void printOptions(){
		assert(false) : "Help: TODO";
		outstream.println("Syntax:\n");
		outstream.println("java -ea -Xmx200m -cp <path> jgi.GradeMergedReads in=<file>");
		outstream.println();
		outstream.println("Other parameters and their defaults:\n");
		outstream.println("overwrite=false  \tOverwrites files that already exist");
		outstream.println("ziplevel=4       \tSet compression level, 1 (low) to 9 (max)");
		outstream.println("interleaved=false\tDetermines whether input file is considered interleaved");
		outstream.println("fastawrap=70     \tLength of lines in fasta output");
		outstream.println("qin=auto         \tASCII offset for input quality.  May be set to 33 (Sanger), 64 (Illumina), or auto");
	}
	
	
	/*--------------------------------------------------------------*/
	
	private String in=null;
	
	private String extin=null;
	
	private String raw1=null;
	private String raw2=null;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	

}

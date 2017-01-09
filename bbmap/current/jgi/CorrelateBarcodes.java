package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ListNum;
import dna.Parser;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ReadWrite;
import fileIO.FileFormat;
import fileIO.TextStreamWriter;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date June 20, 2014
 *
 */
public class CorrelateBarcodes {

	public static void main(String[] args){
		Timer t=new Timer();
		CorrelateBarcodes mb=new CorrelateBarcodes(args);
		mb.process(t);
	}
	
	public CorrelateBarcodes(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		for(String s : args){if(s.startsWith("out=standardout") || s.startsWith("out=stdout")){outstream=System.err;}}
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		boolean setInterleaved=false; //Whether it was explicitly set.

		
		
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
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
//				align2.FastaReadInputStream2.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("outcor") || a.equals("cor")){
				outcor=b;
			}else if(a.equals("bqhist")){
				bqhist=b;
			}else if(a.equals("baqhist")){//average quality
				aqhist=b;
			}else if(a.equals("bmqhist")){//minimum quality
				mqhist=b;
			}else if(a.equals("mmq")){//minimum min quality cutoff
				minBarcodeMinQuality=Integer.parseInt(b);
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
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;

			out1=parser.out1;
			out2=parser.out2;
			
			extin=parser.extin;
			extout=parser.extout;
			
			minBarcodeAverageQuality=parser.minAvgQuality;
		}
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(out1==null){
			if(out2!=null){
				printOptions();
				throw new RuntimeException("Error - cannot define out2 without defining out1.");
			}
		}
		
		if(!setInterleaved){
			assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		if(out2!=null && out2.equalsIgnoreCase("null")){out2=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);

		ffcor=FileFormat.testOutput(outcor, FileFormat.TEXT, extout, true, overwrite, append, false);
		ffaq=FileFormat.testOutput(aqhist, FileFormat.TEXT, extout, true, overwrite, append, false);
		ffmq=FileFormat.testOutput(mqhist, FileFormat.TEXT, extout, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}
	
	void process(Timer t){
		
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
//		if(verbose){
			if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
//		}

		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=4;
			
			if(cris.paired() && out2==null && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}			

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2))) : "out1 and out2 have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, null, null, buff, null, false);
			ros.start();
		}else{ros=null;}
		
		long readsProcessed=0;
		long basesProcessed=0;
		
		long readsTossed=0;
		long basesTossed=0;
		
		ReadStats readstats=null;
		ReadStats.COLLECT_QUALITY_STATS=(bqhist!=null);
		if(ReadStats.COLLECT_QUALITY_STATS){
			ReadStats.QUAL_HIST_FILE=bqhist;
			readstats=new ReadStats();
		}
		
//		assert(false) : bqhist+", "+(ReadStats.COLLECT_QUALITY_STATS)+", "+(readstats==null);
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			outstream.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(reads!=null && reads.size()>0){
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					final byte[] barbases, barquals;
					{
						String[] s=r1.id.split("_");
						barbases=s[0].getBytes();
						barquals=s[1].getBytes();
						for(int i=0; i<barquals.length; i++){
							barquals[i]-=33;
						}
					}
					
					final int qbar=Read.avgQualityByProbability(barbases, barquals, true, 0);
					final int minqbar=Tools.min(barquals);
					aqhistArray[qbar]++;
					mqhistArray[minqbar]++;
					
					if(qbar<minBarcodeAverageQuality || minqbar<minBarcodeMinQuality){
						r1.setDiscarded(true);
						readsTossed++;
						basesTossed+=(initialLength1+initialLength2);
						if(r2!=null){readsTossed++;}
					}
					
//					System.err.println(new String(barquals)+" -> "+qbar);
					
					{
						readsProcessed++;
						basesProcessed+=initialLength1;
						final int q1=r1.avgQualityByProbability(true, 0);
						qualCor1[q1][qbar]++;
					}
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=initialLength2;
						final int q2=r2.avgQualityByProbability(true, 0);
						qualCor2[q2][qbar]++;
					}
					
					if(readstats!=null){
						readstats.addToQualityHistogram(barquals, 0);
					}
					
				}
				
				ArrayList<Read> listOut=reads;
				
				if(ros!=null){
					if(minBarcodeAverageQuality>0){
						listOut=new ArrayList<Read>(reads.size());
						for(Read r : reads){
							if(!r.discarded()){listOut.add(r);}
						}
					}
					ros.add(listOut, ln.id);
				}

				cris.returnList(ln.id, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		if(ffcor!=null){
			TextStreamWriter tsw=new TextStreamWriter(ffcor);
			tsw.start();
			tsw.print("#Read1_Q\tBar_Q\tstdev\tcount\tRead2_Q\tBar_Q\tstdev\tcount\n");
			for(int i=0; i<qualCor1.length; i++){
				long[] array1=qualCor1[i], array2=qualCor2[i];
				long sum1=Tools.sum(array1), sum2=Tools.sum(array2);
				double avg1=Tools.averageHistogram(array1), avg2=Tools.averageHistogram(array2);
				double dev1=Tools.standardDeviationHistogram(array1), dev2=Tools.standardDeviationHistogram(array2);
				tsw.print(String.format("%d\t%.1f\t%.1f\t%d\t%d\t%.1f\t%.1f\t%d\n", i, avg1, dev1, sum1, i, avg2, dev2, sum2));
			}
			tsw.poisonAndWait();
			errorState|=tsw.errorState;
		}
		
		if(aqhist!=null){
			TextStreamWriter tsw=new TextStreamWriter(ffaq);
			tsw.start();
			tsw.print("#Quality\tcount\tfraction\n");
			long sum=Tools.sum(aqhistArray);
			double mult=1.0/Tools.max(1, sum);
			long y=0;
			for(int i=0; i<aqhistArray.length; i++){
				long x=aqhistArray[i];
				tsw.print(String.format("%d\t%d\t%.5f\n", i, x, x*mult));
				y+=x;
				if(y==sum){break;}
			}
			tsw.poisonAndWait();
			errorState|=tsw.errorState;
		}
		
		if(mqhist!=null){
			TextStreamWriter tsw=new TextStreamWriter(ffmq);
			tsw.start();
			tsw.print("#Quality\tcount\tfraction\n");
			long sum=Tools.sum(mqhistArray);
			double mult=1.0/Tools.max(1, sum);
			long y=0;
			for(int i=0; i<mqhistArray.length; i++){
				long x=mqhistArray[i];
				tsw.print(String.format("%d\t%d\t%.5f\n", i, x, x*mult));
				y+=x;
				if(y==sum){break;}
			}
			tsw.poisonAndWait();
			errorState|=tsw.errorState;
		}

		
		if(readstats!=null){
			errorState|=ReadStats.writeAll();
		}
		
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		t.stop();
		
		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);

		String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
		String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		
		if(minBarcodeAverageQuality>0){
			outstream.println();
			outstream.println("Reads Discarded:    "+readsTossed+" \t"+String.format("%.3f%%",readsTossed*100.0/readsProcessed));
			outstream.println("Reads Discarded:    "+basesTossed+" \t"+String.format("%.3f%%",basesTossed*100.0/basesProcessed));
		}
		
		if(errorState){
			throw new RuntimeException("ReformatReads terminated in an error state; the output may be corrupt.");
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
	private String in2=null;
	
	private String qfin1=null;
	private String qfin2=null;

	private String out1=null;
	private String out2=null;
	
	private String extin=null;
	private String extout=null;

	private String outcor=null;
	private String bqhist=null;
	private String aqhist=null;
	private String mqhist=null;

	private int minBarcodeAverageQuality=0;
	private int minBarcodeMinQuality=0;

	private long[][] qualCor1=new long[50][50];
	private long[][] qualCor2=new long[50][50];

	private long[] aqhistArray=new long[100];
	private long[] mqhistArray=new long[100];
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;
	
	private final FileFormat ffcor;
	private final FileFormat ffaq;
	private final FileFormat ffmq;

	private final FileFormat ffout1;
	private final FileFormat ffout2;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	
}

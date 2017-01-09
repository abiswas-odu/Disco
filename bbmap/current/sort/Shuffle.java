package sort;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

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
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import fileIO.FileFormat;


/**
 * Randomizes the order of reads.
 * @author Brian Bushnell
 * @date Oct 27, 2014
 *
 */
public class Shuffle {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		Timer t=new Timer();
		Shuffle sh=new Shuffle(args);
		sh.process(t);
	}
	
	public Shuffle(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		for(String s : args){if(s.startsWith("out=standardout") || s.startsWith("out=stdout")){outstream=System.err;}}
		if(printClass){outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");}
		
		boolean setInterleaved=false; //Whether it was explicitly set.

		
		
		Shared.READ_BUFFER_LENGTH=Tools.min(200, Shared.READ_BUFFER_LENGTH);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		
		int mode_=SHUFFLE;
		
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
			}else if(a.equals("shuffle")){
				mode_=SHUFFLE;
			}else if(a.equals("name")){
				mode_=SORT_NAME;
			}else if(a.equals("coordinate")){
				mode_=SORT_COORD;
			}else if(a.equals("sequence")){
				mode_=SORT_SEQ;
			}else if(a.equals("id")){
				mode_=SORT_ID;
			}else if(a.equals("mode")){
				if(b==null){
					throw new RuntimeException("mode must be shuffle, name, coordinate, sequence, or id.");
				}else if(b.equals("shuffle")){
					mode_=SHUFFLE;
				}else if(b.equals("name")){
					mode_=SORT_NAME;
				}else if(b.equals("coordinate")){
					mode_=SORT_COORD;
				}else if(b.equals("sequence")){
					mode_=SORT_SEQ;
				}else if(b.equals("id")){
					mode_=SORT_ID;
				}else{
					throw new RuntimeException("mode must be shuffle, name, coordinate, sequence, or id.");
				}
			}else if(a.equals("showspeed") || a.equals("ss")){
				showSpeed=Tools.parseBoolean(b);
			}else if(parser.in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				parser.in1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		mode=mode_;
		assert(mode>=1 && mode<=5) : "mode must be shuffle, name, coordinate, sequence, or id.";
		
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
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			
			extin=parser.extin;
			extout=parser.extout;
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

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	void process(Timer t){
		
		ArrayList<Read> bigList=new ArrayList<Read>(65530);
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		long readsProcessed=0;
		long basesProcessed=0;
		
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
					
					{
						readsProcessed++;
						basesProcessed+=initialLength1;
					}
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=initialLength2;
					}
					bigList.add(r1);
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
		errorState|=ReadStats.writeAll();
		
		if(mode==SHUFFLE){
			Collections.shuffle(bigList);
		}else if(mode==SORT_NAME){
			Shared.sort(bigList, ReadComparatorName.comparator);
		}else if(mode==SORT_SEQ){
			Shared.sort(bigList, new ReadComparatorTopological());
		}else if(mode==SORT_COORD){
			Shared.sort(bigList, new ReadComparatorMapping());
		}else if(mode==SORT_ID){
			Shared.sort(bigList, ReadComparatorID.comparator);
		}else{
			assert(false) : "No mode set.";
		}
		
		if(ffout1!=null){
			final ByteStreamWriter bsw1, bsw2;
			if(ffout1!=null){
				bsw1=new ByteStreamWriter(ffout1);
				bsw1.start();
			}else{bsw1=null;}
			if(ffout2!=null){
				bsw2=new ByteStreamWriter(ffout2);
				bsw2.start();
			}else{bsw2=null;}
			final boolean b=(bsw2==null);
			for(int i=0, lim=bigList.size(); i<lim; i++){
				final Read r1=bigList.set(i, null);
				final Read r2=r1.mate;
				bsw1.println(r1, b);
				if(r2!=null && !b){bsw2.println(r2);}
			}
			if(bsw1!=null){errorState|=bsw1.poisonAndWait();}
			if(bsw2!=null){errorState|=bsw2.poisonAndWait();}
		}
		
		t.stop();
		
		if(showSpeed){
			double rpnano=readsProcessed/(double)(t.elapsed);
			double bpnano=basesProcessed/(double)(t.elapsed);

			String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
			String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");
			
			while(rpstring.length()<8){rpstring=" "+rpstring;}
			while(bpstring.length()<8){bpstring=" "+bpstring;}
			
			outstream.println("Time:                         \t"+t);
			outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
			outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		}
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
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
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static class ShuffleThread extends Thread{
		
		public ShuffleThread(String in1_, String in2_, String out1_, String out2_, int mode_, boolean ow_){
			in1=in1_;
			in2=in2_;
			out1=out1_;
			out2=out2_;
			mode=mode_;
			ow=ow_;
		}
		
		@Override
		public void start(){
			addThread(1);
			super.start();
		}
		
		@Override
		public void run(){
			ArrayList<String> list=new ArrayList<String>();
			if(in1!=null){list.add("in1="+in1);}
			if(in2!=null){list.add("in1="+in2);}
			if(out1!=null){list.add("out1="+out1);}
			if(out2!=null){list.add("out2="+out2);}
			list.add("mode="+MODES[mode]);
			list.add("ow="+ow);
			try{
				Shuffle.main(list.toArray(new String[0]));
			}catch(Throwable e){
				System.err.println("Failed to shuffle "+in1+"\nException:"+e+"\n");
			}
			addThread(-1);
		}
		
		final String in1, in2;
		final String out1, out2;
		final int mode;
		final boolean ow;
		
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String in2=null;
	
	private String qfin1=null;
	private String qfin2=null;

	private String out1=null;
	private String out2=null;

	private String qfout1=null;
	private String qfout2=null;
	
	private String extin=null;
	private String extout=null;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
	private final int mode;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;

	private final FileFormat ffout1;
	private final FileFormat ffout2;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private static int maxShuffleThreads=1;
	private static int currentShuffleThreads=0;
	
	public static void setMaxThreads(final int x){
		assert(x>0);
		synchronized(SHUFFLE_LOCK){
			maxShuffleThreads=x;
		}
	}
	
	public static int addThread(final int x){
		synchronized(SHUFFLE_LOCK){
			while(x>0 && currentShuffleThreads>=maxShuffleThreads){
				try {
					SHUFFLE_LOCK.wait(2000);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			currentShuffleThreads+=x;
			if(currentShuffleThreads<maxShuffleThreads){SHUFFLE_LOCK.notify();}
			return currentShuffleThreads;
		}
	}
	
	public static void waitForFinish(){
		synchronized(SHUFFLE_LOCK){
			while(currentShuffleThreads>=maxShuffleThreads){
				try {
					SHUFFLE_LOCK.wait(2000);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	private static String SHUFFLE_LOCK=new String("SHUFFLE_LOCK");
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	public static boolean showSpeed=true;
	public static boolean printClass=true;

	public static final int SHUFFLE=1, SORT_NAME=2, SORT_SEQ=3, SORT_COORD=4, SORT_ID=5;
	public static final String[] MODES={"shuffle", "name", "sequence", "coordinate", "id"};
	
	
}

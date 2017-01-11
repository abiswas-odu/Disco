package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ListNum;
import dna.AminoAcid;
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
import align2.BandedAlignerConcrete;

/**
 * @author Brian Bushnell
 * @date June 20, 2014
 *
 */
public class CountBarcodes {

	public static void main(String[] args){
		Timer t=new Timer();
		CountBarcodes mb=new CountBarcodes(args);
		mb.process(t);
	}
	
	public CountBarcodes(String[] args){
		
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
		

		expectedCodeMap=new HashMap<String,String>(200);
		validCodeMap=new HashMap<String,String>(200);
		
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
			}else if(a.equals("countundefined")){
				countUndefined=Tools.parseBoolean(b);
			}else if(a.equals("printheader") || a.equals("header")){
				printHeader=Tools.parseBoolean(b);
			}else if(a.equals("printrows") || a.equals("rows") || a.equals("maxrows")){
				maxRows=Integer.parseInt(b);
			}else if(a.equals("expected")){
				if(b!=null){
					for(String code : b.split(",")){
						expectedCodeMap.put(code, code);
						validCodeMap.put(code, code);
					}
				}
			}else if(a.equals("valid")){
				if(b!=null){
					for(String code : b.split(",")){
						validCodeMap.put(code, code);
					}
				}
			}else if(a.equals("counts") || a.equals("outc")){
				outCounts=b;
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
//			if(!parser.setOut){
//				out1="stdout";
//			}
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
		
		validCodes=new ArrayList<String>(validCodeMap.size());
		validCodes.addAll(validCodeMap.values());
		expectedCodes=new ArrayList<String>(expectedCodeMap.size());
		expectedCodes.addAll(expectedCodeMap.values());
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffCounts=FileFormat.testOutput(outCounts, FileFormat.TEXT, null, false, overwrite, append, false);

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
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, false);
			ros.start();
		}else{ros=null;}
		
		long readsProcessed=0;
		long basesProcessed=0;
		
		HashMap<String, StringNum> map=new HashMap<String, StringNum>();
		
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
					
					String id=r1.id;
					int colon=id.lastIndexOf(':');
					String barcode=id.substring(colon+1);
					
					if(countUndefined || AminoAcid.isFullyDefined(barcode)){
						StringNum value=map.get(barcode);
						if(value==null){
							value=new StringNum(barcode, 0);
							map.put(barcode, value);
						}
						value.increment();
					}
				}
				
				final ArrayList<Read> listOut=reads;
				
				if(ros!=null){ros.add(listOut, ln.id);}

				cris.returnList(ln.id, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		ArrayList<StringNum> list=new ArrayList<StringNum>(map.size());
		list.addAll(map.values());
		Shared.sort(list);
		Collections.reverse(list);
		
		TextStreamWriter tsw=new TextStreamWriter(ffCounts);
		tsw.start();
		if(printHeader){
			tsw.print("#code\tcount\tHamming_dist\tedit_dist\tvalid\n");
		}
		for(StringNum sn : list){
			if(maxRows==0){break;}
			maxRows--;
			int hdist=calcHdist(sn.s, expectedCodes);
			int edist=hdist;
			if(hdist>1){
				try {
					edist=calcEdist(sn.s, expectedCodes);
				} catch (Exception e) {
					edist=hdist;
				}
			}
			boolean valid=validCodes.contains(sn.s);
			tsw.print(sn+"\t"+hdist+"\t"+edist+"\t"+(valid ? "valid" : "")+"\n");
		}
		tsw.poisonAndWait();
		errorState|=tsw.errorState;
		
		errorState|=ReadStats.writeAll();
		
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
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * @param s
	 * @param expectedCodes
	 * @return
	 */
	static int calcHdist(String s, ArrayList<String> expectedCodes) {
		int min=s.length();
		for(String code : expectedCodes){
			min=Tools.min(min, hdist(s, code));
			if(min<1){break;}
		}
		return min;
	}
	
	static int hdist(String s, String code) {
		final int min=Tools.min(s.length(), code.length());
		int subs=0;
		for(int i=0; i<min; i++){
			if(s.charAt(i)!=code.charAt(i)){subs++;}
		}
		return subs;
	}
	
	/**
	 * @param s
	 * @param expectedCodes
	 * @return
	 */
	int calcEdist(String s, ArrayList<String> expectedCodes) {
		int min=s.length();
		for(String code : expectedCodes){
			min=Tools.min(min, edist(s, code));
			if(min<1){break;}
		}
		return min;
	}
	
	int edist(String s, String code) {
		int x=bandy.alignForward(s.getBytes(), code.getBytes(), 0, 0, s.length(), true);
		return x;
	}
	
	/*--------------------------------------------------------------*/

	private static class StringNum implements Comparable<StringNum>{
		
		public StringNum(String s_, long n_){
			s=s_;
			n=n_;
		}
		
		public long increment(){
			return (n=n+1);
		}
		
		@SuppressWarnings("unused")
		public long increment(long x){
			return (n=n+x);
		}
		
		/* (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(StringNum o) {
			if(n<o.n){return -1;}
			if(n>o.n){return 1;}
			return s.compareTo(o.s);
		}
		
		public String toString(){
			return s+"\t"+n;
		}
		
		@SuppressWarnings("unused")
		public boolean equals(StringNum other){
			if(other==null){return false;}
			if(s==other.s){return true;}
			if(s==null || other.s==null){return false;}
			return compareTo(other)==0;
		}
		
		final String s;
		long n;
		
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
	private String outCounts=null;

	private String qfout1=null;
	private String qfout2=null;
	
	private String extin=null;
	private String extout=null;
	
	/*--------------------------------------------------------------*/

	private boolean reverseComplimentMate=false;
	private boolean reverseCompliment=false;
	private boolean countUndefined=true;
	private boolean printHeader=true;
	private int maxRows=-1;

	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;

	private final FileFormat ffout1;
	private final FileFormat ffout2;
	private final FileFormat ffCounts;

	private final HashMap<String,String> expectedCodeMap;
	private final HashMap<String,String> validCodeMap;
	private final ArrayList<String> expectedCodes;
	private final ArrayList<String> validCodes;
	
	private final BandedAlignerConcrete bandy=new BandedAlignerConcrete(21);
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	
}

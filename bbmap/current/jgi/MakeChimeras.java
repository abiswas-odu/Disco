package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.KillSwitch;
import stream.Read;
import structures.ListNum;
import dna.Parser;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ReadWrite;
import fileIO.FileFormat;
import fileIO.TextStreamWriter;
import shared.Shared;
import shared.Timer;
import shared.Tools;


/**
 * @author Brian Bushnell
 * @date Oct 7, 2014
 *
 */
public class MakeChimeras {

	public static void main(String[] args){
		Timer t=new Timer();
		MakeChimeras mb=new MakeChimeras(args);
		mb.process(t);
	}
	
	public MakeChimeras(String[] args){
		
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
		
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		
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
			}else if(parser.in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				parser.in1=arg;
			}else if(a.equals("forcelength")){
				forceLength=Integer.parseInt(b);
			}else if(a.equals("readsout") || a.equals("chimeras")){
				readsOut=Tools.parseKMG(b);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			readsIn=parser.maxReads;
			
			in1=parser.in1;
			qfin1=parser.qfin1;

			out1=parser.out1;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
	}
	
	void process(Timer t){
		assert(readsOut>0) : "Please set the 'readsout' flag to a positive integer.";
		
		ArrayList<Read> source=new ArrayList<Read>();
		{
			final ConcurrentReadInputStream cris;
			{
				cris=ConcurrentReadInputStream.getReadInputStream(readsIn, false, ffin1, null, qfin1, null);
				if(verbose){outstream.println("Started cris");}
				cris.start(); //4567
			}
			assert(!cris.paired());
			
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
						assert(r1.mate==null);

						final int initialLength1=r1.length();
						
						if(initialLength1>0){
							source.add(r1);
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
			
			outstream.println("Read Time:                    \t"+t);
			outstream.println("Reads In:           "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
			outstream.println("Bases In:           "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		}
		
		
		
		if(readsOut>=0){
			t.start();
			
			final TextStreamWriter tsw;
			if(ffout1==null){
				tsw=null;
			}else{
				tsw=new TextStreamWriter(ffout1);
				tsw.start();
			}
			
			final Random randy=new Random();
			
			long readsProcessed=0;
			long basesProcessed=0;
			final int mod=source.size();
			for(long i=0; i<readsOut; i++){
				Read a=source.get(randy.nextInt(mod));
				Read b=source.get(randy.nextInt(mod));
				Read c=makeChimera(a, b, randy, i);
				if(c==null){
					i--;
				}else{
					if(tsw!=null && c!=null){
						tsw.println(c);
						readsProcessed++;
						basesProcessed+=c.length();
					}
				}
			}
			
			if(tsw!=null){errorState|=tsw.poisonAndWait();}
			
			t.stop();
			
			double rpnano=readsProcessed/(double)(t.elapsed);
			double bpnano=basesProcessed/(double)(t.elapsed);

			String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
			String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");

			while(rpstring.length()<8){rpstring=" "+rpstring;}
			while(bpstring.length()<8){bpstring=" "+bpstring;}

			outstream.println("Write Time:                   \t"+t);
			outstream.println("Reads Out:          "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
			outstream.println("Bases Out:          "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		}

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * @param a
	 * @param b
	 * @param randy
	 * @return
	 */
	private Read makeChimera(Read a, Read b, Random randy, long numericID) {
		final String id=a.id+" ~ "+b.id;

		final Read a2, b2;
		if(forceLength>0){
//			if(a.length()>b.length()){
//				Read c=a;
//				a=b;
//				b=c;
//			}
			a2=getPiece(a, randy, forceLength);
			b2=getPiece(b, randy, b.length()-forceLength);
			if(a2==null || b2==null){return null;}
		}else{
			a2=getPiece(a, randy);
			b2=getPiece(b, randy);
		}
		
		a=b=null;
		
		final byte[] abases=a2.bases, bbases=b2.bases, aquals=a2.quality, bquals=b2.quality;
		final int alen=a2.length(), blen=b2.length();
		final int len=a2.length()+b2.length();
		byte[] bases=new byte[len];
		byte[] quals=(aquals==null || bquals==null) ? null : new byte[len];

		for(int i=0; i<alen; i++){
			bases[i]=abases[i];
			if(quals!=null){quals[i]=aquals[i];}
		}
		for(int i=alen, j=0; j<blen; i++, j++){
			bases[i]=bbases[j];
			if(quals!=null){quals[i]=bquals[j];}
		}
		
		Read r=new Read(bases, -1, -1, -1, id, quals, numericID, 0);
		if(Tools.nextBoolean(randy)){r.reverseComplement();}
		return r;
	}

	/**
	 * @param b
	 * @param randy
	 * @return
	 */
	private Read getPiece(Read a, Random randy) {
		int len=randy.nextInt(a.length())+1;
		
		final int start;
		if(Tools.nextBoolean(randy)){
			if(Tools.nextBoolean(randy)){
				start=0;
			}else{
				start=a.length()-len;
			}
		}else{
			int range=a.length()-len;
			start=randy.nextInt(range+1);
		}
		
		byte[] bases=KillSwitch.copyOfRange(a.bases, start, start+len);
		byte[] quals=a.quality==null ? null : KillSwitch.copyOfRange(a.quality, start, start+len);
		
		Read r=new Read(bases, -1, -1, -1, a.id, quals, a.numericID, 0);
		if(Tools.nextBoolean(randy)){r.reverseComplement();}
		return r;
	}

	/**
	 * @param b
	 * @param randy
	 * @return
	 */
	private Read getPiece(Read a, Random randy, int len) {
		len=Tools.min(len, a.length());
		if(len<1){return null;}
		
		final int start;
		if(Tools.nextBoolean(randy)){
			if(Tools.nextBoolean(randy)){
				start=0;
			}else{
				start=a.length()-len;
			}
		}else{
			int range=a.length()-len;
			start=randy.nextInt(range+1);
		}
		
		byte[] bases=KillSwitch.copyOfRange(a.bases, start, start+len);
		byte[] quals=a.quality==null ? null : KillSwitch.copyOfRange(a.quality, start, start+len);
		
		Read r=new Read(bases, -1, -1, -1, a.id, quals, a.numericID, 0);
		if(Tools.nextBoolean(randy)){r.reverseComplement();}
		return r;
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
	
	private String qfin1=null;

	private String out1=null;
	
	private String extin=null;
	private String extout=null;

	private int forceLength=0;
	
	/*--------------------------------------------------------------*/

	private long readsIn=-1;
	private long readsOut=-1;
	
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

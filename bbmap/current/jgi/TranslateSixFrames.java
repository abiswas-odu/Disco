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
import stream.SamLine;
import structures.ListNum;
import dna.AminoAcid;
import dna.Parser;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ReadWrite;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import fileIO.FileFormat;

/**
 * @author Brian Bushnell
 * @date Sep 11, 2012
 *
 */
public class TranslateSixFrames {

	public static void main(String[] args){
		Timer t=new Timer();
		TranslateSixFrames rr=new TranslateSixFrames(args);
		rr.process(t);
	}
	
	public TranslateSixFrames(String[] args){
		
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
		
		int frames=6;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
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
			}else if(a.equals("tag")){
				addTag=Tools.parseBoolean(b);
			}else if(a.equals("skipquality")){
				skipquality=Tools.parseBoolean(b);
			}else if(a.equals("translatequality")){
				skipquality=!Tools.parseBoolean(b);
			}else if(a.equals("frames")){
				frames=Integer.parseInt(b);
				assert(frames>=0 && frames<=6) : "Frames must be in the range of 0 to 6";
			}else if(a.equals("aain")){
				NT_IN=!Tools.parseBoolean(b);
			}else if(a.equals("ntin")){
				NT_IN=Tools.parseBoolean(b);
			}else if(a.equals("aaout")){
				NT_OUT=!Tools.parseBoolean(b);
			}else if(a.equals("ntout")){
				NT_OUT=Tools.parseBoolean(b);
			}else if(a.equals("frames")){
				frames=Integer.parseInt(b);
				assert(frames>=0 && frames<=6) : "Frames must be in the range of 0 to 6";
			}else if(parser.in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				parser.in1=arg;
			}else if(parser.out1==null && i==1 && !arg.contains("=")){
				parser.out1=arg;
			}else{
				System.err.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		FRAMES=frames;
		Shared.AMINO_IN=!NT_IN;
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;	
			samplerate=parser.samplerate;
			sampleseed=parser.sampleseed;
			
			overwrite=parser.overwrite;
			append=parser.append;

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
			if(FASTQ.FORCE_INTERLEAVED){System.err.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
//		if(maxReads!=-1){ReadWrite.USE_GUNZIP=ReadWrite.USE_UNPIGZ=false;}
		
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
//			if(ReadWrite.isCompressed(in1)){ByteFile.FORCE_MODE_BF2=true;}
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(out1==null){
			if(out2!=null){
				printOptions();
				throw new RuntimeException("Error - cannot define out2 without defining out1.");
			}
			if(!parser.setOut){
				out1="stdout";
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
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTA, extout, true, overwrite, append, false);  
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTA, extout, true, overwrite, append, false);
		
		ffin1=FileFormat.testInput(in1, FileFormat.FASTA, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTA, extin, true, true);
		
		if((ffout1!=null && ffout1.fasta()) || (ffin1!=null && ffin1.fasta())){skipquality=true;}
	}
	
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, useSharedHeader, ffin1, ffin2, qfin1, qfin2);
			cris.setSampleRate(samplerate, sampleseed);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Input is "+(paired ? "paired" : "unpaired"));}

		ConcurrentReadOutputStream ros=null;
		if(out1!=null){
			final int buff=4;
			
			if(cris.paired() && out2==null && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}			

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2))) : "out1 and out2 have same name.";
			
//			System.err.println("Calling ConcurrentReadOutputStream with out1="+out1+", out2="+out2+", qfout1="+qfout1+", qfout2="+qfout2);
			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, useSharedHeader);
			ros.start();
		}
		
		long readsProcessed=0;
		long basesProcessed=0;
		
		long readsOut1=0;
		long readsOut2=0;
		
		long basesOut1=0;
		long basesOut2=0;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			System.err.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(reads!=null && reads.size()>0){
				
				final ArrayList<Read> listOut=new ArrayList<Read>(reads.size()*(NT_IN ? FRAMES : 1));
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					if(NT_IN){//Translate to amino acids
						toFrames(r1, skipquality, addTag, FRAMES, listOut);
					}else{
						listOut.add(r1);
					}
					
					
					{
						readsProcessed++;
						basesProcessed+=initialLength1;
					}
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=initialLength2;
					}
					
					if(addslash){
						if(r1.id==null){r1.id=""+r1.numericID;}
						if(!r1.id.contains(" /1")){r1.id+=" /1";}
						if(r2!=null){
							if(r2.id==null){r2.id=""+r2.numericID;}
							if(!r2.id.contains(" /2")){r2.id+=" /2";}
						}
					}
				}
				
				if(NT_OUT){//Translate to nucleotides
					for(int i=0; i<listOut.size(); i++){
						final Read aa1=listOut.get(i);
						final Read aa2=aa1.mate;
						final Read nt1=aa1.aminoToNucleic();
						if(aa2!=null){
							final Read nt2=aa2.aminoToNucleic();
							nt1.mate=nt2;
							nt2.mate=nt1;
						}
						listOut.set(i, nt1);
					}
				}
				
				for(Read r1 : listOut){
					Read r2=r1.mate;
					readsOut1++;
					basesOut1+=r1.length();

					if(r2!=null){
						readsOut2++;
						basesOut2+=r2.length();
					}
				}
				
				if(ros!=null){ros.add(listOut, ln.id);}

				cris.returnList(ln.id, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		t.stop();
		
		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);

		long readsOut=readsOut1+readsOut2;
		long basesOut=basesOut1+basesOut2;
		
		String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
		String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");
		String rostring=(readsOut<100000 ? ""+readsOut : readsOut<100000000 ? (readsOut/1000)+"k" : (readsOut/1000000)+"m");
		String aastring=(basesOut<100000 ? ""+basesOut : basesOut<100000000 ? (basesOut/1000)+"k" : (basesOut/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}
		while(rostring.length()<8){rostring=" "+rostring;}
		while(aastring.length()<8){aastring=" "+aastring;}
		
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		outstream.println("Reads Out:          "+rostring);
		outstream.println((NT_OUT ? "Bases Out:          " : "Amino Acids Out:    ")+aastring);
		
		if(errorState){
			throw new RuntimeException("TranslateSixFrames terminated in an error state; the output may be corrupt.");
		}
	}
	
	public static final ArrayList<Read> toFrames(Read r1, boolean skipquality, boolean addTag, int frames){
		return toFrames(r1, skipquality, addTag, frames, new ArrayList<Read>(frames));
	}
	
	public static final ArrayList<Read> toFrames(Read r1, boolean skipquality, boolean addTag, int frames, ArrayList<Read> listOut){
		final Read r2=r1.mate;
		final byte[][] bm1=AminoAcid.toAAsSixFrames(r1.bases);
		final byte[][] qm1=(skipquality ? QNULL : AminoAcid.toQualitySixFrames(r1.quality, 0));
		final byte[][] bm2=(r2==null ? null : AminoAcid.toAAsSixFrames(r2.bases));
		final byte[][] qm2=(r2==null ? null : (skipquality ? QNULL : AminoAcid.toQualitySixFrames(r2.quality, 0)));

		for(int i=0; i<frames; i++){
			Read aa1=new Read(bm1[i], r1.chrom, r1.start, r1.stop, (addTag ? r1.id+frametag[i] : r1.id), qm1[i], r1.numericID, r1.flags|Read.AAMASK);
			Read aa2=null;
			if(r2!=null){
				aa2=new Read(bm2[i], r2.chrom, r2.start, r2.stop, (addTag ? r2.id+frametag[i] : r2.id), qm2[i], r2.numericID, r2.flags|Read.AAMASK);
				aa1.mate=aa2;
				aa2.mate=aa1;
			}
			if(aa1.bases!=null || (aa2!=null && aa2.bases!=null)){listOut.add(aa1);}
		}
		return listOut;
	}
	
	/*--------------------------------------------------------------*/
	
	private void printOptions(){
		outstream.println("Syntax:\n");
		outstream.println("java -ea -Xmx512m -cp <path> jgi.TranslateSixFrames in=<infile> in2=<infile2> out=<outfile> out2=<outfile2>");
		outstream.println("\nin2 and out2 are optional.  \nIf input is paired and there is only one output file, it will be written interleaved.\n");
		outstream.println("Other parameters and their defaults:\n");
		outstream.println("overwrite=false  \tOverwrites files that already exist");
		outstream.println("ziplevel=4       \tSet compression level, 1 (low) to 9 (max)");
		outstream.println("interleaved=false\tDetermines whether input file is considered interleaved");
		outstream.println("fastawrap=70     \tLength of lines in fasta output");
		outstream.println("qin=auto         \tASCII offset for input quality.  May be set to 33 (Sanger), 64 (Illumina), or auto");
		outstream.println("qout=auto        \tASCII offset for output quality.  May be set to 33 (Sanger), 64 (Illumina), or auto (meaning same as input)");
	}
	
	
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
	
	/** Add /1 and /2 to paired reads */
	private boolean addslash=false;

	private boolean skipquality=false;
	
	private boolean NT_IN=true;
	private boolean NT_OUT=false;

	private long maxReads=-1;
	private float samplerate=1f;
	private long sampleseed=-1;
	
	private final int FRAMES;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;

	private final FileFormat ffout1;
	private final FileFormat ffout2;
	
	private static final String[] frametag=new String[] {" fr1", " fr2", " fr3", " fr4", " fr5", " fr6"};
	private static final byte[][] QNULL=new byte[6][];
	private boolean addTag=true;
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	private boolean useSharedHeader;
	
}

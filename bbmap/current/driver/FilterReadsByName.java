package driver;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.ReadStreamWriter;
import stream.SamLine;
import structures.ListNum;
import dna.Parser;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ReadWrite;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import fileIO.FileFormat;

/**
 * @author Brian Bushnell
 * @date Oct 8, 2014
 *
 */
public class FilterReadsByName {

	public static void main(String[] args){
		Timer t=new Timer();
		FilterReadsByName mb=new FilterReadsByName(args);
		mb.process(t);
	}
	
	public FilterReadsByName(String[] args){
		
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
		
		SamLine.SET_FROM_OK=true;
		ReadStreamWriter.USE_ATTACHED_SAMLINE=true;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens

			if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
//				align2.FastaReadInputStream2.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("names")){
				if(b!=null){
					String[] x=b.split(",");
					for(String s : x){
						names.add(s);
					}
				}
			}else if(a.equals("substrings") || a.equals("substring")){
				if(b==null){b="t";}
				if(b.equals("header")){
					headerSubstringOfName=true;
				}else if(b.equals("name")){
					nameSubstringOfHeader=true;
				}else{
					nameSubstringOfHeader=headerSubstringOfName=Tools.parseBoolean(b);
				}
			}else if(a.equals("casesensitive") || a.equals("case")){
				ignoreCase=!Tools.parseBoolean(b);
			}else if(a.equals("include") || a.equals("retain")){
				exclude=!Tools.parseBoolean(b);
			}else if(a.equals("exclude") || a.equals("remove")){
				exclude=Tools.parseBoolean(b);
			}else if(a.equals("prefix") || a.equals("prefixmode")){
				prefixmode=Tools.parseBoolean(b);
			}else if(a.equals("minlen") || a.equals("minlength")){
				minLength=(int)Tools.parseKMG(b);
			}else if(a.equals("from")){
				fromPos=(int)Tools.parseKMG(b);
			}else if(a.equals("to")){
				toPos=(int)Tools.parseKMG(b);
			}else if(a.equals("pos") || a.equals("range")){
				String[] split2=b.split("-");
				fromPos=(int)Tools.parseKMG(split2[0]);
				toPos=(int)Tools.parseKMG(split2[1]);
			}else if(a.equals("truncate")){
				truncateWhitespace=truncateHeaderSymbol=Tools.parseBoolean(b);
			}else if(a.equals("truncatewhitespace") || a.equals("tws")){
				truncateWhitespace=Tools.parseBoolean(b);
			}else if(a.equals("truncateheadersymbol") || a.equals("ths")){
				truncateHeaderSymbol=Tools.parseBoolean(b);
			}else if(parser.in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				parser.in1=arg;
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{
			String[] x=names.toArray(new String[names.size()]);
			names.clear();
			for(String s : x){
				Tools.addNames(s, names, true);
			}
		}
		if(ignoreCase){
			String[] x=names.toArray(new String[names.size()]);
			names.clear();
			for(String s : x){
				names.add(s.toLowerCase());
			}
		}
		if(truncateHeaderSymbol || truncateWhitespace){
			String[] x=names.toArray(new String[names.size()]);
			names.clear();
			for(String s : x){
				String s2=s;
				if(truncateHeaderSymbol && s.length()>1 && (s.charAt(0)=='@' || s.charAt(0)=='>')){s2=s.substring(1);}
				if(truncateWhitespace){s2=s.trim();}
				if(s2.length()>0){
					names.add(s2);
				}
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
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		if(ffin1!=null && ffout1!=null && ffin1.samOrBam() && ffout1.samOrBam()){
			useSharedHeader=true;
		}
	}
	
	void process(Timer t){
		
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, useSharedHeader, ffin1, ffin2, qfin1, qfin2);
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
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, useSharedHeader);
			ros.start();
		}else{ros=null;}
		
		long readsProcessed=0;
		long basesProcessed=0;
		
		long readsOut=0;
		long basesOut=0;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			outstream.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			
			
			while(reads!=null && reads.size()>0){
				
				ArrayList<Read> retain=new ArrayList<Read>(reads.size());
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					readsProcessed+=1+r1.mateCount();
					basesProcessed+=initialLength1+initialLength2;
					
					final String header;
					{
						String temp=(ignoreCase ? r1.id.toLowerCase() : r1.id);
						temp=truncateWhitespace ? temp.trim() : temp;
						header=temp;
					}
					String prefix=null;
					for(int x=1; x<header.length(); x++){
						char prev=x>=2 ? header.charAt(x-2) : 'X';
						char c=header.charAt(x-1);
						char next=header.charAt(x);
						if(Character.isWhitespace(c) || (c=='/' && (next=='1' || next=='2'))){
							prefix=header.substring(0, x).trim();
							break;
						}else if(Character.isWhitespace(prev) && (c=='1' || c=='2') && next==':'){
							prefix=header.substring(0, x).trim();
							break;
						}
					}
					
					boolean keepThisRead=(initialLength1>=minLength || initialLength2>=minLength);
					boolean match=false;
					if(keepThisRead){
						match=(names.contains(header) || (prefix!=null && names.contains(prefix)));
						if(!match && (nameSubstringOfHeader || headerSubstringOfName)){
							for(String name : names){
								if((headerSubstringOfName && name.contains(header)) || (nameSubstringOfHeader && header.contains(name))){match=true;}
								else if(prefix!=null && ((headerSubstringOfName && name.contains(prefix)) || (nameSubstringOfHeader && prefix.contains(name)))){match=true;}
							}
						}else if(!match && prefixmode){
							for(String name : names){
								if(header.startsWith(name)){match=true;} //TODO: Fast hashing like in DemuxByName
							}
						}
						keepThisRead=(match!=exclude);
					}
					
//					assert(false) : names.contains(name)+", "+name+", "+prefix+", "+exclude;
					
					if(keepThisRead){
						if(fromPos>=0){
							TrimRead.trimToPosition(r1, fromPos, toPos, 1);
						}
						retain.add(r1);
						readsOut+=1+r1.mateCount();
						basesOut+=r1.length()+r1.mateLength();
					}
				}
				
				final ArrayList<Read> listOut=retain;
				
				if(ros!=null){ros.add(listOut, ln.id);}

				cris.returnList(ln.id, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadStats.writeAll();
		
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		t.stop();
		
		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);

//		String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
//		String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");
//
//		while(rpstring.length()<8){rpstring=" "+rpstring;}
//		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		outstream.println("Time:               "+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+basesProcessed+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		outstream.println("Reads Out:          "+readsOut);
		outstream.println("Bases Out:          "+basesOut);
		
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
	private boolean exclude=true;
	private boolean prefixmode=false;
	private boolean nameSubstringOfHeader=false;
	private boolean headerSubstringOfName=false;
	private boolean ignoreCase=true;
	private boolean truncateHeaderSymbol=false;
	private boolean truncateWhitespace=false;

	private int minLength=0;

	private int fromPos=-1;
	private int toPos=-1;
	
	private LinkedHashSet<String> names=new LinkedHashSet<String>();
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;

	private final FileFormat ffout1;
	private final FileFormat ffout2;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	private boolean useSharedHeader=false;
	
}

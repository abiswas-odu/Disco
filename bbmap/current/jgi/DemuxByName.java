package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;

import stream.ArrayListSet;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.MultiCros;
import stream.Read;
import structures.ListNum;
import dna.Parser;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ReadWrite;
import fileIO.FileFormat;
import fileIO.TextFile;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;


/**
 * @author Brian Bushnell
 * @date Oct 9, 2014
 *
 */
public class DemuxByName {

	public static void main(String[] args){
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		final int oldCap=Shared.READ_BUFFER_NUM_BUFFERS(), oldZipThreads=ReadWrite.MAX_ZIP_THREADS;
		final boolean oldPigz=ReadWrite.USE_PIGZ, oldUnpigz=ReadWrite.USE_UNPIGZ;
		
		Timer t=new Timer();
		DemuxByName mb=new DemuxByName(args);
		mb.process(t);
		
		Shared.setBuffers(oldCap);
		ReadWrite.USE_PIGZ=oldPigz;
		ReadWrite.USE_UNPIGZ=oldUnpigz;
		ReadWrite.MAX_ZIP_THREADS=oldZipThreads;
	}
	
	public DemuxByName(String[] args){
		
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
			}else if(a.equals("names") || a.equals("name") || a.equals("affixes")){
				if(b!=null){
					String[] x=b.split(",");
					for(String s : x){
						names.add(s);
					}
				}
			}else if(a.equalsIgnoreCase("length") || a.equalsIgnoreCase("len") || a.equalsIgnoreCase("affixlength") || a.equalsIgnoreCase("affixlen")){
				fixedAffixLength=Integer.parseInt(b);
			}else if(a.equals("prefixmode") || a.equals("prefix") || a.equals("pm")){
				prefixMode=Tools.parseBoolean(b);
			}else if(a.equals("suffixmode") || a.equals("suffix") || a.equals("sm")){
				prefixMode=!Tools.parseBoolean(b);
			}else if(a.equals("substringmode") || a.equals("substring")){
				substringMode=Tools.parseBoolean(b);
			}else if(a.equals("outu") || a.equals("outu1")){
				outu1=b;
			}else if(a.equals("outu2")){
				outu2=b;
			}else if(a.equals("delimiter")){
				if(b==null){delimiter=null;}
				else if(b.equalsIgnoreCase("space")){
					delimiter=" ";
				}else if(b.equalsIgnoreCase("tab")){
					delimiter="\t";
				}else if(b.equalsIgnoreCase("whitespace")){
					delimiter="\\s+";
				}else{
					delimiter=b;
				}
			}else if(parser.in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				parser.in1=arg;
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
				File f=new File(s);
				if(f.exists() && f.isFile()){
					TextFile tf=new TextFile(s);
					String[] lines=tf.toStringLines();
					for(String s2 : lines){
						names.add(s2);
					}
				}else{
					names.add(s);
				}
			}

//			for(String s : names){
//				assert(affixLen<0 || affixLen==s.length()) : "All names must have the same length.";
//				affixLen=s.length();
//			}
//			assert(affixLen>0) : "Must include at least one non-zero-length affix (name).";
			
			{
				BitSet bs=new BitSet();
				if(fixedAffixLength>0){
					bs.set(fixedAffixLength);
				}
				for(String s : names){
					bs.set(s.length());
				}
				affixLengths=new int[bs.cardinality()];
				for(int i=0, bit=-1; i<affixLengths.length; i++){
					bit=bs.nextSetBit(bit+1);
					affixLengths[i]=bit;
				}
				Arrays.sort(affixLengths);
				Tools.reverseInPlace(affixLengths);
			}
			
			assert((affixLengths.length>0 && affixLengths[0]>0) || delimiter!=null) : "Must include at least one non-zero-length affix (name), or a delimiter.";
			ReadWrite.MAX_ZIP_THREADS=Tools.max(1, Tools.min(ReadWrite.MAX_ZIP_THREADS, (Shared.threads()*2-1)/Tools.max(1, names.size())));
			if(names.size()>8){ReadWrite.USE_PIGZ=false;}
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

		assert(out1==null || out1.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		assert(out2==null || out2.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		assert(qfout1==null || qfout1.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		assert(qfout2==null || qfout2.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		if(outu1!=null && outu2==null && outu1.indexOf('#')>-1){
			outu2=outu1.replace("#", "2");
			outu1=outu1.replace("#", "1");
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
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, outu1, outu2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2+", "+outu1+", "+outu2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+", "+outu1+", "+outu2+"\n");
		}

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		if(ffin1!=null && out1!=null && ffin1.samOrBam()){
			String ext=ReadWrite.rawExtension(out1);
			useSharedHeader=FileFormat.isSamOrBam(ext);
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
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		final MultiCros mcros;
		if(out1!=null){
			final int buff=4;
			
			mcros=(fixedAffixLength>0 || delimiter!=null ? new MultiCros(out1, out2, false, overwrite, append, true, useSharedHeader, FileFormat.FASTQ, buff) : null);
			
			if(paired && out2==null && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}			

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2))) : "out1 and out2 have same name.";
			
			for(String s : names){
				
				String qf1=null, qf2=null;
				if(qfout1!=null){qf1=qfout1.replace("%", s);}
				if(qfout2!=null){qf2=qfout2.replace("%", s);}
				
				FileFormat ffout1=FileFormat.testOutput(out1.replace("%", s), FileFormat.FASTQ, extout, true, overwrite, append, false);
				FileFormat ffout2=(out2==null ? null : FileFormat.testOutput(out2.replace("%", s), FileFormat.FASTQ, extout, true, overwrite, append, false));
				ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qf1, qf2, buff, null, true);
				ros.start();
				nameToStream.put(s, ros);
			}
		}else{
			mcros=null;
		}
		
		final ConcurrentReadOutputStream rosu;
		if(outu1!=null){
			final int buff=4;
			
			FileFormat ffout1=FileFormat.testOutput(outu1, FileFormat.FASTQ, extout, true, overwrite, append, false);
			FileFormat ffout2=(outu2==null ? null : FileFormat.testOutput(outu2, FileFormat.FASTQ, extout, true, overwrite, append, false));
			rosu=ConcurrentReadOutputStream.getStream(ffout1, ffout2, null, null, buff, null, true);
			rosu.start();
		}else{
			rosu=null;
		}
		
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
			
			for(String s : names){
				nameToArray.put(s, new ArrayList<Read>());
			}
			final ArrayListSet als=(fixedAffixLength<1 && delimiter==null ? null : new ArrayListSet(false));
			
			
			while(reads!=null && reads.size()>0){
				
				ArrayList<Read> unmatched=new ArrayList<Read>();
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					final String id=r1.id;
					final int idlen=id.length();
					
					ArrayList<Read> al2=null;
					if(names.size()>0){
						if(substringMode){
							for(String s : names){
								if(id.contains(s)){
									al2=nameToArray.get(s);
									break;
								}
							}
						}else{
							for(int affixLen : affixLengths){
								final String sub=idlen>=affixLen ? prefixMode ? id.substring(0, affixLen) : id.substring(idlen-affixLen) : id;
								al2=nameToArray.get(sub);
								if(al2!=null){break;}
							}
						}
					}
					
					if(al2!=null || als!=null){
						readsOut++;
						basesOut+=initialLength1;
						if(r2!=null){
							readsOut++;
							basesOut+=initialLength2;
						}
						
						if(al2!=null){
							al2.add(r1);
							{
								readsOut++;
								basesOut+=initialLength1;
							}
							if(r2!=null){
								readsOut++;
								basesOut+=initialLength2;
							}
						}else if(als!=null){
							String sub=r1.id;
							if(fixedAffixLength>0){
								sub=(sub.length()<=fixedAffixLength ? sub : prefixMode ? id.substring(0, fixedAffixLength) : id.substring(idlen-fixedAffixLength));
							}else{
								assert(delimiter!=null);
								String[] split=sub.split(delimiter);
								assert(split.length>1) : "Delimiter '"+delimiter+"' was not found in name '"+sub+"'";
								sub=split[prefixMode ? 0 : split.length-1];
							}
							als.add(r1, sub);
						}
						
					}else{
						unmatched.add(r1);
					}
					{
						readsProcessed++;
						basesProcessed+=initialLength1;
					}
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=initialLength2;
					}
				}
				
				for(String s : names){
					ArrayList<Read> listOut=nameToArray.put(s, new ArrayList<Read>());
					ConcurrentReadOutputStream ros=nameToStream.get(s);
					if(ros!=null){ros.add(listOut, ln.id);}
				}
				if(mcros!=null){
					if(als.size()+names.size()>8){ReadWrite.USE_PIGZ=false;}
					mcros.add(als, ln.id);
				}
				if(rosu!=null){rosu.add(unmatched, ln.id);}
				
				cris.returnList(ln.id, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadStats.writeAll();
		
		errorState|=ReadWrite.closeStream(cris);
		
		for(String s : names){
			ConcurrentReadOutputStream ros=nameToStream.get(s);
			errorState|=ReadWrite.closeStream(ros);
		}
		
		if(mcros!=null){
			errorState|=ReadWrite.closeStreams(mcros);
		}
		
		t.stop();
		
		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);
		
		outstream.println("Time:               "+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+basesProcessed+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		outstream.println("Reads Out:    "+readsOut);
		outstream.println("Bases Out:    "+basesOut);
		
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

	private String outu1=null;
	private String outu2=null;

	private String qfout1=null;
	private String qfout2=null;
	
	private String extin=null;
	private String extout=null;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
//	private boolean exclude=true;

	private String delimiter=null;
	private boolean prefixMode=true;
	private boolean substringMode=false;
//	private int affixLen=-1;
	
	private int fixedAffixLength=-1;
	
	private int[] affixLengths;
	
	private HashSet<String> names=new HashSet<String>();
	private HashMap<String, ArrayList<Read>> nameToArray=new HashMap<String, ArrayList<Read>>();
	private HashMap<String, ConcurrentReadOutputStream> nameToStream=new HashMap<String, ConcurrentReadOutputStream>();
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	private boolean useSharedHeader=false;
	
}

package jgi;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import bloom.KCountArray;
import bloom.KmerCount7MT;
import bloom.KmerCountAbstract;


import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.ReadStreamWriter;
import structures.ListNum;
import dna.AminoAcid;
import dna.Parser;
import fileIO.ReadWrite;
import fileIO.FileFormat;
import fileIO.TextFile;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;

/**
 * @author Brian Bushnell
 * @date Aug 14, 2012
 *
 */
public class MateReadsMT {
	
	
	public static void main(String[] args){
		MateReadsMT mr=new MateReadsMT(args);
		mr.process();
	}
	
	public MateReadsMT(String[] args){
		System.err.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		System.err.println("BBMerge version "+version);
		
		Timer ttotal=new Timer();
		ttotal.start();
		
		in1_primary=(args[0].indexOf('=')>0 ? null : args[0]);
		in2_primary=(in1_primary!=null && args.length>1 && args[1].indexOf('=')<0 ? args[1] : null);
		if(in2_primary!=null && "null".equalsIgnoreCase(in2_primary)){in2_primary=null;}
		
		{
			if(in1_primary!=null && !in1_primary.contains(",") && !in1_primary.startsWith("stdin.") && !in1_primary.equals("stdin")){
				File f=new File(in1_primary);
				if(!f.exists() || !f.isFile()){
					in1_primary=null;
//					throw new RuntimeException(in1+" does not exist.");
				}
			}
			if(in2_primary!=null && !in2_primary.contains(",")){
				File f=new File(in2_primary);
				if(!f.exists() || !f.isFile()){
					in2_primary=null;
//					throw new RuntimeException(in2+" does not exist.");
				}else if(in1_primary.equalsIgnoreCase(in2_primary)){
					throw new RuntimeException("Both input files are the same.");
				}
			}
		}
		
		Parser parser=new Parser();
		ReadWrite.MAX_ZIP_THREADS=Shared.threads()-1;
		
		
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=(split.length>1 ? split[1] : "true");

			if(Parser.isJavaFlag(arg)){
				//jvm argument; do nothing
			}else if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQualityAdjust(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(parser.parseTrim(arg, a, b)){
				//do nothing
			}else if(a.equals("in") || a.equals("in1")){
				in1_primary=b;
			}else if(a.equals("in2")){
				in2_primary=b;
			}else if(a.equals("k") || a.equals("kmer")){
				k_G=Integer.parseInt(b);
			}else if(a.equals("minoverlappingbases") || a.equals("minoverlapbases")){
				MIN_OVERLAPPING_BASES=Integer.parseInt(b);
			}else if(a.equals("minoverlap") || a.equals("minoverlappingkmers") || a.equals("minoverlapkmers")){
				MIN_OVERLAPPING_KMERS=Integer.parseInt(b);
			}else if(a.equals("minoverlappingbases0") || a.equals("minoverlapbases0")){
				MIN_OVERLAPPING_BASES_0=Integer.parseInt(b);
			}else if(a.equals("minoverlap0") || a.equals("minoverlappingkmers0") || a.equals("minoverlapkmers0")){
				MIN_OVERLAPPING_KMERS_0=Integer.parseInt(b);
			}else if(a.equals("minoverlapinsert") || a.equals("minoi")){
				MIN_OVERLAP_INSERT=Integer.parseInt(b);
			}else if(a.equals("badlimit")){
				DEFAULT_BADLIMIT=Integer.parseInt(b);
			}else if(a.startsWith("matrixbits")){
				int matrixbits=Integer.parseInt(b);
				assert(matrixbits<63);
				totalcells_G=1L<<matrixbits;
			}else if(a.startsWith("cells")){
				totalcells_G=Tools.parseKMG(b);
			}else if(a.equals("passes")){
				passes_G=Integer.parseInt(b);
			}else if(a.equals("bin")){
				bin=Integer.parseInt(b);
			}else if(a.equals("threads") || a.equals("t")){
				THREADS=Integer.parseInt(b);
			}else if(a.startsWith("minq")){
				MIN_QUALITY=(byte)Integer.parseInt(b);
			}else if(a.startsWith("minqo")){
				MIN_QUALITY_FOR_OVERLAP=(byte)Integer.parseInt(b);
			}else if(a.equals("maxq")){
				Read.MAX_MERGE_QUALITY=(byte)Integer.parseInt(b);
			}else if(a.startsWith("minprob")){
				KmerCountAbstract.minProb=Float.parseFloat(b);
			}else if(a.startsWith("hashes") || a.startsWith("multihash")){
				hashes=Integer.parseInt(b);
				assert(hashes>0 && hashes<25);
			}else if(a.startsWith("cbits") || a.startsWith("cellbits")){
				cbits_G=Integer.parseInt(b);
				int cmax=(1<<cbits_G)-1;
				MAX_HITS_FOR_BAD=Tools.min(MAX_HITS_FOR_BAD, cmax-1);
				MIN_HITS_FOR_GOOD=Tools.min(MIN_HITS_FOR_GOOD, cmax);
			}else if(a.startsWith("minvotes")){
				MIN_VOTES=Integer.parseInt(b);
			}else if(a.endsWith("hitsforbad")){
				MAX_HITS_FOR_BAD=Integer.parseInt(b);
			}else if(a.endsWith("hitsforgood")){
				MIN_HITS_FOR_GOOD=Integer.parseInt(b);
			}else if(a.equals("maxbadbases")){
				DEFAULT_BADLIMIT_FOR_BASE_MATCHING=Integer.parseInt(b);
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads_G=Tools.parseKMG(b);
			}else if(a.equals("tablereads") || a.startsWith("tablereads")){
				tableReads_G=Tools.parseKMG(b);
			}else if(a.startsWith("gap")){
				if(b.equalsIgnoreCase("true") || b.equalsIgnoreCase("null")){
					b="";
					System.err.println("Note - no hash tables will be used as there are no gaps.");
					gap_G=new int[1][0];
				}else{
					String[] h=b.split("[;:]");
					gap_G=new int[h.length][];
					for(int m=0; m<h.length; m++){
						String[] g=h[m].split(",");
						if(g.length==0 || (g.length==1 && (g[0].length()==0 || g[0].equalsIgnoreCase("null")))){
							gap_G[m]=new int[0];
						}else{
							gap_G[m]=new int[g.length];
							for(int j=0; j<g.length; j++){
								gap_G[m][j]=Integer.parseInt(g[j]);
							}
						}
					}
				}
//				for(int m=0; m<gap.length; m++){System.err.println(Arrays.toString(gap[m]));}
//				assert(false);
			}else if(a.startsWith("extra")){
				extra_G=Arrays.asList(b.split(","));
			}else if(a.equals("outgood") || a.startsWith("outpair") || a.equals("outmerged") || a.equals("out")){
				outgood_G=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("temp") || a.equals("tempfile")){
				assert(b.contains("#") && b.contains(".txt"));
				tempfile=b;
			}else if(a.equals("outb") || a.equals("outu") || a.equals("outunmerged") || a.equals("outbad") || a.startsWith("outunpair") || a.startsWith("outsingle")){
				outbad_G=(b==null || b.equals("null") ? null : b);
//				assert(outbad==null || b.contains("#")) : "Unpaired read output filname must contain '#' symbol.";
			}else if(a.startsWith("outinsert") || a.startsWith("outlength")){
				outinsert_G=(b==null || b.equals("null") ? null : b);
			}else if(a.startsWith("outhist3") || a.equals("hist3")){
				outhist3=(b==null || b.equals("null") ? null : b);
			}else if(a.startsWith("outhist2") || a.equals("hist2")){
				outhist2=(b==null || b.equals("null") ? null : b);
			}else if(a.startsWith("outhist") || a.startsWith("hist") || a.equals("ihist")){
				outhist=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outputfailed")){
				OUTPUT_FAILED=Tools.parseBoolean(b);
			}else if(a.equals("mix")){
				MIX_BAD_AND_GOOD=Tools.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("join")){
				join_G=Tools.parseBoolean(b);
			}else if(a.equals("usemapping")){
				USE_MAPPING=Tools.parseBoolean(b);
			}else if(a.equals("useoverlap") || a.equals("usebases") || a.equals("matebyoverlap") || a.equals("matebybases")){
				MATE_BY_OVERLAP=Tools.parseBoolean(b);
			}else if(a.startsWith("skipmated")){
				SKIP_MATED_READS=Tools.parseBoolean(b);
			}else if(a.startsWith("writeintermediatejoined")){
				WRITE_INTERMEDIATE_JOINED=Tools.parseBoolean(b);
			}else if(a.startsWith("fillmiddleinter")){
				FILL_MIDDLE_INTERMEDIATE=Tools.parseBoolean(b);
			}else if(a.startsWith("fillmiddlefinal")){
				FILL_MIDDLE_FINAL=Tools.parseBoolean(b);
			}else if(a.equals("fillmiddle")){
				FILL_MIDDLE_INTERMEDIATE=FILL_MIDDLE_FINAL=Tools.parseBoolean(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Tools.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
			}else if(a.equals("trimonfailure") || a.equals("tof")){
				if(b!=null && Character.isDigit(b.charAt(0))){
					TRIM_ON_OVERLAP_FAILURE=Integer.parseInt(b);
				}else{
					TRIM_ON_OVERLAP_FAILURE=(Tools.parseBoolean(b) ? 1 : 0);
				}
			}else if(a.equals("mi") || a.equals("minins") || a.equals("mininsert")){
				minInsert=Integer.parseInt(b);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			qtrimLeft=parser.qtrimLeft;
			qtrimRight=parser.qtrimRight;
			trimq=parser.trimq;
			qtrim=((qtrimLeft||qtrimRight)&&trimq>=0);
			minReadLength=parser.minReadLength;
			untrim=parser.untrim;
		}
		
		if(in2_primary==null && in1_primary!=null && in1_primary.contains("#") && !new File(in1_primary).exists()){
			in2_primary=in1_primary.replaceFirst("#", "2");
			in1_primary=in1_primary.replaceFirst("#", "1");
		}
		
		if(in2_primary!=null){
			assert(!in1_primary.equalsIgnoreCase(in2_primary));
			FASTQ.TEST_INTERLEAVED=false;
			FASTQ.FORCE_INTERLEAVED=false;
		}else{
			FASTQ.TEST_INTERLEAVED=true;
			FASTQ.FORCE_INTERLEAVED=true;
		}
		
//		assert(false) : MATE_BY_OVERLAP;

		if(FILL_MIDDLE_INTERMEDIATE){
			if(!WRITE_INTERMEDIATE_JOINED){System.err.println("WRITE_INTERMEDIATE_JOINED forced to true.");}
			WRITE_INTERMEDIATE_JOINED=true;
		}
		if(WRITE_INTERMEDIATE_JOINED){
			if(!join_G){System.err.println("Final output forced to be joined reads.");}
			join_G=true;
			//Ultimately I could re-read the initial files, so this is not truly needed. 
		}
	}
	
	void process(){
		Timer ttotal=new Timer();
		ttotal.start();
//		assert(!FASTQ.PARSE_CUSTOM);
		
		final int hwthr=Shared.threads();
		if(THREADS<1){THREADS=hwthr;}
		System.err.println("Detected "+Runtime.getRuntime().availableProcessors()+" hardware threads; using "+THREADS+" for main process.");
		long memory=(Runtime.getRuntime().maxMemory());
		System.err.println("Detected "+(memory/(1L<<20))+" MB available memory.");
//		System.err.println("PARSE_CUSTOM="+FASTQ.PARSE_CUSTOM);
		
		if(gap_G!=null){
			for(int[] g : gap_G){
				maxtables=Tools.max(maxtables, g.length);
				for(int g2 : g){assert(g2>0) : "TODO: Ungapped kmers do not currently work.  Please use gap lengths of >0.";}
			}
		}
		if(maxtables<1 && !USE_MAPPING && !MATE_BY_OVERLAP){
			throw new RuntimeException("No gap sizes have been specified, so there is no work to do.");
		}
		
		if(passes_G>1){totalcells_G*=2;}
		
		if(auto && maxtables>0 && totalcells_G<0){
			final long usable=(long)Tools.max(((memory-(256000000))*.7), memory*0.4);
			long mem=usable;
			totalcells_G=(mem*8)/cbits_G;
			
//			long tablebytes=((1L<<matrixbits)*cbits)/8;
//			if(tablebytes*3<usable){matrixbits++;}
//			System.err.println(tablebytes/1000000+", "+usable/1000000+", "+(tablebytes*3)/1000000);
			
			System.err.println("\nAuto settings:");
			System.err.println("k:          \t"+k_G);
			System.err.println("cbits:      \t"+cbits_G);
//			System.err.println("matrixbits: \t"+matrixbits);
//			System.err.println("matrixbits2:\t"+matrixbits2);
			System.err.println("cells:      \t"+Tools.toKMG(totalcells_G));
			System.err.println("hashes:     \t"+hashes);
			System.err.println();
		}else if(totalcells_G==-1){
			totalcells_G=1L<<34;
		}
		
		String in1=in1_primary, in2=in2_primary;
		
		KCountArray middleTable=null;
		if(FILL_MIDDLE_INTERMEDIATE || FILL_MIDDLE_FINAL){
			maxtables++;
			long cells=totalcells_G/maxtables;
			if(k_G<32 && cells>(1L<<(2*k_G))){cells=(1L<<(2*k_G));}
			middleTable=KmerCount7MT.makeKca(in1, in2, extra_G, MIDDLE_TABLE_K, cbits_G, 0, cells, hashes+1, MIN_QUALITY, true, tableReads_G, 1, 4, 2, 2, null, 0);
			middleTable.shutdown();
			System.err.println("MiddleTable: \tgap = "+middleTable.gap+"   \tmem = "+middleTable.mem()+"   \tused = "+String.format("%.3f%%",middleTable.usedFraction()*100));
		}
		
		final int cmax=(1<<cbits_G)-1;
		assert(MIN_HITS_FOR_GOOD>MAX_HITS_FOR_BAD && MIN_HITS_FOR_GOOD<=cmax && MAX_HITS_FOR_BAD>0);
		
		
		
		final int phases=(gap_G==null ? 1 : gap_G.length);
		
		KmerCountAbstract.PREJOIN=false;
		
		String a1=in1, a2=in2;
		
		int oldzip=ReadWrite.ZIPLEVEL;
		for(int phase=0; phase<phases-1; phase++){
			ReadWrite.ZIPLEVEL=Tools.min(oldzip, 4);
			String temp=tempfile.replaceFirst("#", "#_P"+phase);
			runPhase(gap_G[phase], a1, a2, extra_G, null, temp, null, cbits_G, k_G, totalcells_G, hashes, passes_G, 
					WRITE_INTERMEDIATE_JOINED, Tools.max(maxReads_G, tableReads_G), tableReads_G, true, (FILL_MIDDLE_INTERMEDIATE ? middleTable : null));
			a1=temp.replaceFirst("#", "1");
			a2=temp.replaceFirst("#", "2");
			System.err.println("\nPhase "+(phase+1)+" statistics.");
			System.err.println("Insert range:        \t"+insertMinTotal+" - "+insertMaxTotal);
			System.err.println("90th percentile:     \t"+Tools.percentile(histTotal, .9));
			System.err.println("50th percentile:     \t"+Tools.percentile(histTotal, .5));
			System.err.println("10th percentile:     \t"+Tools.percentile(histTotal, .1));
			KmerCountAbstract.PREJOIN=true;
		}
		ReadWrite.ZIPLEVEL=oldzip;
		if(!FILL_MIDDLE_FINAL){middleTable=null;}
		runPhase((gap_G==null ? null : gap_G[gap_G.length-1]), a1, a2, extra_G, outinsert_G, outgood_G, outbad_G, cbits_G, k_G, totalcells_G, hashes, passes_G, 
				join_G, maxReads_G, tableReads_G, false, middleTable);
		
		if(outhist!=null){
			StringBuilder sb=new StringBuilder();
//			for(int i=0; i<histTotal.length; i++){
//				sb.append(i+"\t"+histTotal[i]+"\n");
//			}
			for(int i=0; i<histTotal.length && i<=insertMaxTotal; i+=bin){
				int x=0;
				int y=0;
				for(int j=i; j<i+bin && j<histTotal.length; j++){
					x+=histTotal[j];
					y++;
				}
				x=(x+bin-1)/y;
				sb.append(i+"\t"+x+"\n");
			}
			ReadWrite.writeStringInThread(sb, outhist);
		}

		if(outhist2!=null){
			StringBuilder sb=new StringBuilder();
			for(int i=0; i<histTotal.length && i<=insertMaxTotal; i+=bin){
				int x=0;
				int y=0;
				for(int j=i; j<i+bin && j<histTotal.length; j++){
					x+=histTotal[j];
					y++;
				}
				x=(x+bin-1)/y;
				sb.append(x+"\n");
			}
			ReadWrite.writeStringInThread(sb, outhist2);
		}

		if(outhist3!=null){
			
			if(!new File(outhist3).exists()){
				StringBuilder sb=new StringBuilder();
				for(int i=0; i<histTotal.length; i+=bin){
					sb.append(i+"\n");
				}
				ReadWrite.writeString(sb, outhist3);
			}
			
			StringBuilder sb=new StringBuilder();
			TextFile tf=new TextFile(outhist3, false, false);
			for(int i=0; i<histTotal.length; i+=bin){
				int x=0;
				int y=0;
				for(int j=i; j<i+bin && j<histTotal.length && i<=insertMaxTotal; j++){
					x+=histTotal[j];
					y++;
				}
				x=(x+bin-1)/y;
				sb.append(tf.readLine()+"\t"+x+"\n");
			}
			tf.close();
			ReadWrite.writeStringInThread(sb, outhist3);
		}
		
		
		ttotal.stop();
		System.err.println("Total time: "+ttotal+"\n");
		
		long sum=correctCountTotal+incorrectCountTotal;
		
		double div=100d/readsProcessedTotal;
		System.err.println("Pairs:       \t"+readsProcessedTotal);
		System.err.println("Joined:      \t"+sum+String.format((sum<10000 ? "       " : "   ")+"\t%.3f%%", sum*div));
		if(FASTQ.PARSE_CUSTOM){
			System.err.println("Correct:     \t"+correctCountTotal+String.format((correctCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", correctCountTotal*div));
			System.err.println("Incorrect:   \t"+incorrectCountTotal+String.format((incorrectCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", incorrectCountTotal*div));
		}
		System.err.println("Ambiguous:   \t"+ambiguousCountTotal+String.format((ambiguousCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", ambiguousCountTotal*div));
		System.err.println("No Solution: \t"+noSolutionCountTotal+String.format((noSolutionCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", noSolutionCountTotal*div));
		if(minInsert>0){System.err.println("Too Short:   \t"+tooShortCountTotal+String.format((tooShortCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", tooShortCountTotal*div));}
		System.err.println("Avg Insert:          \t\t"+String.format("%.1f", (insertSumCorrectTotal+insertSumIncorrectTotal)*1d/(correctCountTotal+incorrectCountTotal)));
		if(FASTQ.PARSE_CUSTOM){
			System.err.println("Avg Insert Correct:  \t\t"+String.format("%.1f", (insertSumCorrectTotal)*1d/(correctCountTotal)));
			System.err.println("Avg Insert Incorrect:\t\t"+String.format("%.1f", (insertSumIncorrectTotal)*1d/(incorrectCountTotal)));
		}
		
		System.err.println("\nPhase "+(phases)+" statistics.");
		System.err.println("Insert range:        \t"+insertMinTotal+" - "+insertMaxTotal);
		System.err.println("90th percentile:     \t"+Tools.percentile(histTotal, .9));
		System.err.println("50th percentile:     \t"+Tools.percentile(histTotal, .5));
		System.err.println("10th percentile:     \t"+Tools.percentile(histTotal, .1));
	}
	
	public static void runPhase(int[] gap, String in1, String in2, List<String> extra, String outinsert, String outgood, String outbad, 
			int cbits, int k, long totalcells, int multihash, int passes, boolean join, long maxReads, long tableReads, boolean perfectonly, KCountArray middleTable){
		
		assert(((USE_MAPPING || MATE_BY_OVERLAP) && MIN_VOTES<2) || 
				(MIN_VOTES>0 && MIN_VOTES<=gap.length)) : "minVotes is set too high.  Should be at most the number of (overlapping) gaps.";
		
		Timer thash=new Timer(), talign=new Timer();
		
		assert(totalcells>1);
		if(middleTable!=null){totalcells=totalcells-middleTable.cells;}
		long cells=totalcells/(gap==null || gap.length==0 ? 1 : gap.length);
		if(k<32 && cells>1L<<(2*k)){cells=1L<<(2*k);}
		
		ConcurrentReadOutputStream rosgood=null;
		ConcurrentReadOutputStream rosbad=null;
		ConcurrentReadOutputStream rosinsert=null;
		
		if(outgood!=null){
			final String out1, out2;
			
//			assert(outgood.contains("#") || sam || fq) : outgood;
			if(outgood.contains("#")){
				out1=outgood.replaceFirst("#", "1");
				out2=outgood.replaceFirst("#", "2");
			}else{
				out1=outgood;
				out2=null;
				if(!join){System.err.println("Writing joinable reads interleaved.");}
				else{System.err.println("Writing joinable reads joined.");}
			}
			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1));
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2)));

			final FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, append, true);
			final FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, append, true);
			assert(!ff1.samOrBam()) : "Sam files need reference info for the header.";
			
			final int buff=Tools.max(16, 2*THREADS);
			rosgood=ConcurrentReadOutputStream.getStream(ff1, ff2, null, null, buff, null, false);
			rosgood.start();
		}
		
		if(outbad!=null){
			final String out1, out2;
			
//			assert(outbad.contains("#") || sam || fq) : outbad;
			if(outbad.contains("#")){
				out1=outbad.replaceFirst("#", "1");
				out2=outbad.replaceFirst("#", "2");
			}else{
				out1=outbad;
				out2=null;
				System.err.println("Writing unjoinable reads interleaved.");
			}
			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1));
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2)));

			final FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, append, true);
			final FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, append, true);
			assert(!ff1.samOrBam()) : "Sam files need reference info for the header.";
			
			final int buff=Tools.max(16, 2*THREADS);
			rosbad=ConcurrentReadOutputStream.getStream(ff1, ff2, null, null, buff, null, false);
			rosbad.start();
		}
		
		if(outinsert!=null){
			final int buff=Tools.max(16, 2*THREADS);
			
			String out1=outinsert.replaceFirst("#", "1");

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1));
			
			ReadStreamWriter.HEADER=header();
			final FileFormat ff=FileFormat.testOutput(out1, FileFormat.ATTACHMENT, ".info", true, overwrite, append, true);
			rosinsert=ConcurrentReadOutputStream.getStream(ff, null, null, null, buff, null, false);
			rosinsert.start();
		}
		
		
		if(rosgood!=null || rosbad!=null || rosinsert!=null){
			System.err.println("Started output threads.");
		}
		
		thash.start();
		
		KCountArray[] kca=new KCountArray[gap==null ? 0 : gap.length];
		for(int i=0; i<kca.length; i++){
			kca[i]=KmerCount7MT.makeKca(in1, in2, extra, k, cbits, gap[i], cells, multihash, MIN_QUALITY, true, tableReads, passes, 4, 2, 2, null, 0);
		}
		
		for(int i=0; i<kca.length; i++){
			kca[i].shutdown();
			System.err.println("Table "+i+":\tgap = "+kca[i].gap+"   \tmem = "+kca[i].mem()+"   \tused = "+String.format("%.3f%%",kca[i].usedFraction()*100));
//			printStatistics(kca[i]);
		}
		
		if(kca!=null && kca.length>0){
			thash.stop();
			System.err.println("Hash time:  "+thash);
		}
		
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
//		assert(paired);//Fails on empty files.
		if(verbose){System.err.println("Paired: "+paired);}
		
		talign.start();
		
		
		MateThread[] pta=new MateThread[THREADS];
		for(int i=0; i<pta.length; i++){
			pta[i]=new MateThread(cris, rosgood, rosbad, rosinsert, k, kca, join, perfectonly, middleTable);
			pta[i].start();
		}

		insertMinTotal=999999999;
		insertMaxTotal=0;
		
		readsProcessedTotal=0;
		matedCountTotal=0;
		correctCountTotal=0;
		ambiguousCountTotal=0;
		tooShortCountTotal=0;
		incorrectCountTotal=0;
		noSolutionCountTotal=0;
		insertSumCorrectTotal=0;
		insertSumIncorrectTotal=0;
		
		Arrays.fill(histTotal, 0);
		
		for(int i=0; i<pta.length; i++){
			MateThread ct=pta[i];
			synchronized(ct){
				while(ct.isAlive()){
					try {
						ct.join(1000);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				
				readsProcessedTotal+=ct.readsProcessed;
				matedCountTotal+=ct.matedCount;
				correctCountTotal+=ct.correctCount;
				ambiguousCountTotal+=ct.ambiguousCount;
				tooShortCountTotal+=ct.tooShortCount;
				incorrectCountTotal+=ct.incorrectCount;
				noSolutionCountTotal+=ct.noSolutionCount;
				insertSumCorrectTotal+=ct.insertSumCorrect;
				insertSumIncorrectTotal+=ct.insertSumIncorrect;
				
				basesTrimmedTotal+=ct.basesTrimmedT;
				readsTrimmedTotal+=ct.readsTrimmedT;

				insertMinTotal=Tools.min(ct.insertMin, insertMinTotal);
				insertMaxTotal=Tools.max(ct.insertMax, insertMaxTotal);
				
//				System.err.println(ct.insertMin+", "+ct.insertMax);
				
				if(ct.hist!=null){
					for(int h=0; h<ct.hist.length; h++){
						histTotal[h]+=ct.hist[h];
					}
				}
			}
		}
		
		System.err.println("Finished reading");
		errorState|=ReadWrite.closeStreams(cris, rosgood, rosbad, rosinsert);
		
		talign.stop();
//		System.err.println("Align time: "+talign);
	}
	


	public static String header(){
		return "#id\tnumericID\tinsert\tstatus\thashHits\thashMisses\tscore\tsum\tvotes\n";
	}
	
	
	/**
	 * @param r
	 * @param mate
	 */
	public static int mateRead(Read a, Read b, int k1, int k2, long mask1, long mask2, KCountArray kca, int[] rvector) {
		assert(false) : "Redo this based off of kca[] version.";
		if(rvector==null){rvector=new int[6];}
		final int width=kca.gap+k1+k2;
		
		int maxInsert=kca.gap+a.length()+b.length()-MIN_OVERLAPPING_KMERS+1;
		int minInsert=width+MIN_OVERLAPPING_KMERS-1;
		
		if(maxInsert<minInsert){
			return -1; //Can't be tested
		}
		
		if(a.obj==null){a.obj=hash(a, k1, mask1, 2*k2);}
		if(b.obj==null){b.obj=hash(b, k2, mask2, 0);}
		long[] half1=(long[])a.obj;
		long[] half2=(long[])b.obj;
		
		int bestInsert=-1;
		int bestScore=-1; //This serves as the threshold for the minimum score to report. 
		int bestGood=-1;
		int bestBad=DEFAULT_BADLIMIT;
		int bestSum=0;
		
		int bestMismatches=DEFAULT_MISMATCHLIMIT;
		int bestMatches=0;
		
		final int pivot=k1+kca.gap+b.length()+(k1-k2);
		boolean ambig=false;
		
		for(int insert=minInsert; insert<=maxInsert; insert++){
			int score=scoreIP(half1, half2, insert, pivot, kca, rvector, bestBad);
			if(score>0){
				
				final int overlap=Tools.min(insert, a.length()+b.length()-insert);
				int matches=0;
				int mismatches=0;
				boolean ok=true;
				if(overlap>=4){
					mismatches=countMismatches(a, b, insert, bestMismatches+10);
					matches=overlap-mismatches;
					ok=(mismatches<3 || mismatches*2<matches);
					bestMismatches=Tools.min(mismatches, bestMismatches);
				}
				
				int good=rvector[1], bad=rvector[2], sum=rvector[3];
				if(ok && good>=MIN_OVERLAPPING_KMERS && bad<=bestBad){
					if(bad<bestBad || good>bestGood || (good==bestGood && sum>bestSum)){
						ambig=(bestBad==0);
						bestScore=score;
						bestInsert=insert;
						bestGood=good;
						bestBad=bad;
						bestSum=sum;
						if(ambig){break;}
					}else if(good==bestGood && sum==bestSum){
						assert(bad==bestBad && sum==bestSum) : bad+"~"+bestBad+", "+good+"~"+bestGood+", "+sum+"~"+bestSum;
						ambig=true;
					}
				}
			}
		}
		
		rvector[0]=bestScore;
		rvector[1]=bestGood;
		rvector[2]=bestBad;
		rvector[3]=bestSum;
		rvector[4]=(ambig ? 1 : 0);
		
		return bestInsert;
	}
	

	public static int mateRead(Read a, Read b, int k1, int k2, long mask1, long mask2, KCountArray kca[], int[] rvector) {
		
//		verbose=a.numericID>145;
		
		if(USE_MAPPING){
			rvector[0]=100;
			rvector[1]=20;
			rvector[2]=0;
			rvector[3]=20; //What is this?
			rvector[4]=0;
			rvector[5]=Tools.max(1, MIN_VOTES);
			return a.insertSizeMapped(ignoreMappingStrand);
		}
		if(rvector==null){rvector=new int[6];}
		
		if(a.obj==null){a.obj=hash(a, k1, mask1, 2*k2);}
		if(b.obj==null){b.obj=hash(b, k2, mask2, 0);}
		long[] half1=(long[])a.obj;
		long[] half2=(long[])b.obj;
		if(half1==null || half2==null){return -1;}
		
		int bestInsert=-1;
		int bestScore=-1; //This serves as the threshold for the minimum score to report. 
		int bestGood=-1;
		int bestBad=DEFAULT_BADLIMIT;
		
		int minGap=kca[0].gap;
		int maxGap=kca[0].gap;

		final int[] pivot=new int[kca.length];
		final int[] minInsert=new int[kca.length];
		final int[] maxInsert=new int[kca.length];
		
		for(int i=0; i<kca.length; i++){
			int gap=kca[i].gap;
			minGap=Tools.min(minGap, gap);
			maxGap=Tools.max(maxGap, gap);
			pivot[i]=k1+gap+b.length()+(k1-k2);
			minInsert[i]=minGap+k1+k2+MIN_OVERLAPPING_KMERS-1;
			maxInsert[i]=maxGap+a.length()+b.length()-MIN_OVERLAPPING_KMERS+1;
		}
		

		int overallMaxInsert=maxGap+a.length()+b.length()-MIN_OVERLAPPING_KMERS+1;
		int overallMinInsert=minGap+k1+k2+MIN_OVERLAPPING_KMERS-1;
		
		if(overallMaxInsert<overallMinInsert){
			return -1; //Can't be tested
		}
		
		int bestMismatches=DEFAULT_MISMATCHLIMIT;
		int bestMatches=0;
		int bestVotes=0;
//		int bestSum=0;
		int bestMatchScore=0;
		
		boolean ambig=false;
		final int minOverlap=MIN_OVERLAPPING_KMERS_0;
		
//		System.err.println("len = "+kca.length);
//		int guess=(kca.length==1 ? 0 : mateRead(a, b, k1, k2, mask1, mask2, new KCountArray[] {kca[1]}, null));
//		boolean vb=(guess==331 || guess==332);
//		assert(vb) : a.insertSize() +", "+ guess+", "+kca.length;
		
//		assert(false) : overallMinInsert+", "+overallMaxInsert+", "+minInsert[0]+", "+maxInsert[0];
		
		for(int insert=overallMinInsert; insert<=overallMaxInsert; insert++){
//			verbose=(insert==174);
			if(verbose){System.err.println("\nTesting read "+a.numericID+", insert "+insert);}
			
			int good=-1;
			int bad=999999;
			int score=-1;
			int votes=0;
//			int sum=0;
			int matchScore=0;
			for(int g=0; g<kca.length; g++){
				if(insert>=minInsert[g] && insert<=maxInsert[g]){
					if(verbose){System.err.println("Testing gap "+kca[g].gap);}
					int x=scoreIP(half1, half2, insert, pivot[g], kca[g], rvector, bestBad);
					final int good0=rvector[1], bad0=rvector[2];
					if(verbose){System.err.println("score="+score+", rvector="+Arrays.toString(rvector));}
					if((good0>MIN_OVERLAPPING_KMERS) && (bad0>bestBad || good0+bad0>=ACCEL_FACTOR)){
//						score=votes==0 ? x : Tools.min(score, x);
//						if(verbose){System.err.println("new score="+score);}
//						good=votes==0 ? good0 : Tools.min(good0, good);
//						bad=votes==0 ? bad0 : Tools.max(bad0, bad);

						score=votes==0 ? x : Tools.max(score, x);
						if(verbose){System.err.println("new score="+score);}
						good=votes==0 ? good0 : Tools.max(good0, good);
						bad=votes==0 ? bad0 : Tools.min(bad0, bad);
						//					sum=votes==0 ? rvector[3] : Tools.max(rvector[3], sum);
						votes++;
						if(bad>bestBad || score<=0){break;}
					}
				}
				
			}
			if(score>0/* && votes>=MIN_VOTES*/){
				
				final int overlap=Tools.min(insert, a.length()+b.length()-insert);
				int matches=0;
				int mismatches=0;
				boolean ok=true;
				if(overlap>=minOverlap){
					mismatches=countMismatches(a, b, insert, bestMismatches+10);
					matches=overlap-mismatches;
					matchScore=matches-mismatches*4;
					ok=(mismatches<3 || mismatches*2<matches);
				}
				
				if(ok && good>=MIN_OVERLAPPING_KMERS && bad<=bestBad){
					

					ambig=true;
					boolean setBest=false;
					boolean quit=(bad==0 && bestBad==0);
					
					if(bad<bestBad){
						setBest=true;
						ambig=false;
					}else if(votes>bestVotes){
						setBest=true;
					}else if(votes==bestVotes && matches>=bestMatches && mismatches<=bestMismatches && (matches>bestMatches || mismatches<bestMismatches)){
						//						if(mismatches>=bestMismatches){ambig=true;}
						setBest=true;
					}else if(matchScore>=bestMatchScore && (good>bestGood /*|| (good==bestGood && sum>bestSum)*/)){
						setBest=true;
					}
					
					if(setBest){
						bestScore=score;
						bestInsert=insert;
						bestGood=good;
						bestBad=bad;
						bestVotes=votes;
						if(overlap>=minOverlap){
							bestMismatches=mismatches;
							bestMatches=matches;
							bestMatchScore=matchScore;
						}
						if(votes<MIN_VOTES){ambig=true;}
					}
					if(quit){break;}
				}
			}
		}

//		if(vb){
//			if(guess==331){
//				if(bestInsert==331){
//					good331++;
//					System.err.println(guess+", "+bestInsert+", id="+a.numericID+", "+Arrays.toString(rvector));
//				}else{
//					bad331++;
//					System.err.println(guess+", "+bestInsert+", id="+a.numericID+", "+Arrays.toString(rvector));
//				}
//			}else if(guess==332){
//				if(bestInsert==332){
//					good332++;
//					System.err.println(guess+", "+bestInsert+", id="+a.numericID+", "+Arrays.toString(rvector));
//				}else{
//					bad332++;
//					System.err.println(guess+", "+bestInsert+", id="+a.numericID+", "+Arrays.toString(rvector));
//				}
//			}
//			if(good331>0 && bad331>0 && good332>0 && bad332>0){assert(false);}
//		}
		
		rvector[0]=bestScore;
		rvector[1]=bestGood;
		rvector[2]=bestBad;
//		rvector[3]=bestSum;
		rvector[4]=(ambig ? 1 : 0);
		rvector[5]=bestVotes;
		
		return bestInsert;
	}
	

	public static int mateByOverlap(Read a, Read b, int[] rvector, final int minOverlap0, final int minOverlap) {
		if(USE_MAPPING){
			rvector[0]=100;
			rvector[1]=20;
			rvector[2]=0;
			rvector[3]=20; //What is this?
			rvector[4]=0;
			rvector[5]=Tools.max(1, MIN_VOTES);
			return a.insertSizeMapped(ignoreMappingStrand);
		}
		if(rvector==null){rvector=new int[6];}
		final byte[] abases=a.bases, bbases=b.bases, aqual=a.quality, bqual=b.quality;
		
		int bestOverlap=-1;
//		int bestScore=-1; //This serves as the threshold for the minimum score to report. 
		int bestGood=-1;
		int bestBad=DEFAULT_BADLIMIT_FOR_BASE_MATCHING;
		final int margin=2;
		
		boolean ambig=false;
		final int maxOverlap=abases.length+bbases.length-Tools.max(minOverlap, MIN_OVERLAP_INSERT);
//		assert(false) : minOverlap+", "+maxOverlap;
//		System.err.print("\nm");
		
		for(int overlap=Tools.max(minOverlap0, 0); overlap<maxOverlap; overlap++){
//			System.err.print("\nn");
//			verbose=(insert==174);
			if(verbose){System.err.println("\nTesting read "+a.numericID+", overlap "+overlap+", insert "+(abases.length+bbases.length-overlap));}
			
			
			int tested=0;
			int good=0, bad=0;
			
			int istart=(overlap<=abases.length ? 0 : overlap-abases.length);
			int jstart=(overlap<=abases.length ? abases.length-overlap : 0);
//			System.err.print("o");
			
			for(int i=istart, j=jstart, badlim=bestBad+margin; i<overlap && i<bbases.length && j<abases.length && bad<badlim; i++, j++){
				assert(j>=0 && j<=abases.length && i>=0 && i<=bbases.length) : "\njstart="+jstart+", j="+j+", istart="+istart+", i="+i+" \n"+
						"overlap="+overlap+", a.length="+a.length()+", b.length="+b.length()+", bad="+bad+", badlim="+badlim+", good="+good+", tested="+tested;
				byte ca=abases[j], cb=bbases[i];
				if(ca=='N' || cb=='N' || (aqual!=null && aqual[j]<MIN_QUALITY_FOR_OVERLAP) || (bqual!=null && bqual[i]<MIN_QUALITY_FOR_OVERLAP)){
					//do nothing
				}else{
					assert(AminoAcid.isFullyDefined(ca) && AminoAcid.isFullyDefined(cb)) : (char)ca+", "+(char)cb;
					tested++;
					if(ca==cb){good++;}
					else{bad++;}
				}
			}
//			System.err.print("p");
			
//			System.err.println(overlap+", "+bestOverlap+", "+bestGood+", "+bestBad+", "+ambig);

//			System.err.print("a");
			if(good>minOverlap){//Candidate
				if(bad<=bestBad){

//					System.err.print("b");
					if(bad<bestBad || (bad==bestBad && good>bestGood)){//Current winner
						if(bad>bestBad-margin){ambig=true;}
						bestOverlap=overlap;
						bestBad=bad;
						bestGood=good;
//						assert(abases.length+bbases.length-bestOverlap<299) : 
//							((abases.length+bbases.length-bestOverlap)+", "+ambig+", "+overlap+", "+good+", "+bad+", "+tested+", "+bestGood+", "+bestBad+", "+a.insertSize());
					}else if(bad==bestBad){
						ambig=true;
					}

//					System.err.print("c");
					if(ambig && bestBad<margin){
//						System.err.print("d");
						rvector[0]=((bestBad==0 ? 8 : 4)*bestGood-6*bestBad);
						rvector[1]=bestGood;
						rvector[2]=bestBad;
//						rvector[3]=bestSum;
						rvector[4]=(ambig ? 1 : 0);
						rvector[5]=0;
//						System.err.print("e");
						return -1;
					}

//					System.err.print("f");
				}
			}else if(bad<margin){
//				System.err.print("g");
				ambig=true;
				rvector[0]=((bestBad==0 ? 8 : 4)*bestGood-6*bestBad);
				rvector[1]=bestGood;
				rvector[2]=bestBad;
//				rvector[3]=bestSum;
				rvector[4]=(ambig ? 1 : 0);
				rvector[5]=0;
				return -1;
			}
//			System.err.print("h");
			
//			if(abases.length+bbases.length-bestOverlap>299){
//				System.err.println((abases.length+bbases.length-bestOverlap)+", "+ambig+", "+rvector[0]+", "+bestGood+", "+bestBad+", "+a.insertSize());
//			}
			
		}
//		System.err.println("i");
		
		rvector[0]=((bestBad==0 ? 8 : 4)*bestGood-6*bestBad);
		rvector[1]=bestGood;
		rvector[2]=bestBad;
//		rvector[3]=bestSum;
		rvector[4]=(ambig ? 1 : 0);
		rvector[5]=0;
		
//		if(abases.length+bbases.length-bestOverlap>299){
//			System.err.println((abases.length+bbases.length-bestOverlap)+", "+ambig+", "+rvector[0]+", "+bestGood+", "+bestBad+", "+a.insertSize());
//		}
		
//		assert(bestOverlap>-1);
		return (bestOverlap<0 ? -1 : abases.length+bbases.length-bestOverlap);
	}
	
	
	public static int countMismatches(Read a, Read b, int insert, int maxMismatches){
		final int lengthSum=a.length()+b.length();
		if(insert>=lengthSum){return 0;}
		final int overlap=Tools.min(insert, lengthSum-insert);
		
		int mismatches=0;
		
		
		int start1=(insert>a.length() ? a.length()-overlap : 0);
		int start2=(insert>=b.length() ? 0 : b.length()-overlap);
//		System.err.println(insert+", "+overlap+", "+start1+", "+start2);
		
		while(start1<0 || start2<0){start1++; start2++;}
		for(int i=start1, j=start2; i<a.length() && j<b.length(); i++, j++){
			final byte ca=a.bases[i], cb=b.bases[j];
			if(ca!=cb){
				final byte qa=a.quality[i], qb=b.quality[j];
				if(ca=='N' || cb=='N' || qa<7 || qb<7){
					//do nothing
				}else{
					mismatches++;
					if(mismatches>maxMismatches){break;}
				}
			}
		}
		return mismatches;
	}
	
	public static int scoreIP(long[] half1, long[] half2, int insert, int pivot, KCountArray kca, int[] rvector, final int badlimit){
		int start1, start2;
		if(insert<=pivot){ //Short mode; start from base 0 of read a
			start1=0;
			start2=pivot-insert;
		}else{ //Long mode; start from base 0 of read b
			start1=insert-pivot;
			start2=0;
		}
		if(verbose){
			System.err.println("ScoreIP.  Insert:  "+insert+", gap="+kca.gap+", badlimit="+badlimit);
		}
		return score(half1, half2, start1, start2, kca, rvector, badlimit);
	}
	
	public static int score(long[] half1, long[] half2, int start1, int start2, KCountArray kca, int[] rvector, final int badlimit){
		int good=0;
		int bad=0;
		int sum=0;
		final int len=Tools.min(half1.length-start1, half2.length-start2);
//		final int incr=Tools.min(len/8, 8); //Accelerates scoring by a factor of 8 for a preview
		final int incr=Tools.min(len/ACCEL_FACTOR, ACCEL_FACTOR);
		
		if(incr>1){
			for(int i=start1, j=start2; i<half1.length && j<half2.length; i+=incr, j+=incr){
				if(half1[i]!=-1 && half2[j]!=-1){
					long key=half1[i]|half2[j];
					int x=kca.read(key);
					sum+=x;
					if(x>=MIN_HITS_FOR_GOOD){good++;}
					else if(x<=MAX_HITS_FOR_BAD){
						bad++;
						if(bad>badlimit){break;}
					}
//					if(verbose){System.err.print("("+Long.toHexString(half1[i])+","+Long.toHexString(half2[j])+","+Long.toHexString(key)+","+x+")");}
				}
			}
			if(verbose){
				System.err.println("\n(incr="+incr+") Good:  "+good+" \tBad:  "+bad);
			}
			if(bad>good || bad>badlimit){
				rvector[0]=incr*((bad==0 ? 8 : 4)*good-6*bad);
				rvector[1]=incr*good;
				rvector[2]=incr*bad;
				rvector[3]=incr*sum;
				return rvector[0];
			}else{
				good=0;
				bad=0;
				sum=0;
			}
		}
		
		for(int i=start1, j=start2; i<half1.length && j<half2.length; i++, j++){
			if(half1[i]!=-1 && half2[j]!=-1){
				long key=half1[i]|half2[j];
				int x=kca.read(key);
				sum+=x;
				if(x>=MIN_HITS_FOR_GOOD){good++;}
				else if(x<=MAX_HITS_FOR_BAD){
					bad++;
					if(bad>badlimit){break;}
				}
//				if(verbose){System.err.print("("+Long.toHexString(half1[i])+","+Long.toHexString(half2[j])+","+Long.toHexString(key)+","+x+")");}
			}
		}
		if(verbose){
			System.err.println("\nGood:  "+good+" \tBad:  "+bad+" \tSum:  "+sum);
		}
		rvector[0]=((bad==0 ? 8 : 4)*good-6*bad);
		rvector[1]=good;
		rvector[2]=bad;
		rvector[3]=sum;
		return rvector[0];
	}
	
	public static void toHex(long[] array){
		for(int i=0; i<array.length; i++){
			System.out.print(Long.toHexString(array[i])+", ");
		}
		System.out.println();
	}
	
	public static void toHex(long[] array1, long[] array2){
		for(int i=0; i<array1.length; i++){
			System.out.print(Long.toHexString(array1[i]|array2[i])+", ");
		}
		System.out.println();
	}
	
	public static long[] hash(Read r, int k, long mask, int offset){
		if(r==null || r.bases==null || r.length()<k){return null;}
		final byte[] bases=r.bases;
		final byte[] quals=r.quality;
		final long[] half=new long[bases.length-k+1];
		Arrays.fill(half, -1);
		
		
//		System.err.println(k+", "+bases.length+", "+half.length+", offset="+offset);
		
		int len=0;
		long kmer=0;
		for(int i=0, j=i-k+1; i<bases.length; i++, j++){
//			System.err.println(len+", "+i+", "+j);
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0 || (quals!=null && quals[i]<MIN_QUALITY)){
				len=0;
				kmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;
				if(len>=k){
					half[j]=(kmer<<offset);
				}
			}
		}
		return half;
	}
	
	
	private static class MateThread extends Thread{
		
		
		public MateThread(ConcurrentReadInputStream cris_, ConcurrentReadOutputStream rosgood_, ConcurrentReadOutputStream rosbad_, ConcurrentReadOutputStream rosinsert_,
				int k_, KCountArray[] kca_, boolean joinReads_, boolean joinperfectonly_, KCountArray middleTable_) {
			cris=cris_;
			rosgood=rosgood_;
			rosbad=rosbad_;
			rosinsert=rosinsert_;
			k=k_;
			kca=kca_;
			joinReads=joinReads_;
			joinperfectonly=joinperfectonly_;
			middleTable=middleTable_;
		}
		
		
		@Override
		public void run(){
			processMate();
		}

		private void processMate() {

			final int k1=((k+1)/2);
			final int k2=k/2;
			//		assert(k1+k2>=1 && k1+k2<20) : k1+", "+k2+", "+(k1+k2);
			assert(USE_MAPPING || MATE_BY_OVERLAP || kca[0].gap>=0);
			final int kbits1=2*k1;
			final int kbits2=2*k2;
			final long mask1=~((-1L)<<(kbits1));
			final long mask2=~((-1L)<<(kbits2));

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(r.mate!=null || WRITE_INTERMEDIATE_JOINED);
			}


			while(reads!=null && reads.size()>0){

				ArrayList<Read> listg=(rosgood==null ? null : new ArrayList<Read>());
				ArrayList<Read> listb=(rosbad==null ? null : new ArrayList<Read>());
				ArrayList<Read> listi=(rosinsert==null ? null : new ArrayList<Read>());

				for(Read r1 : reads){
					final Read r2=r1.mate;
					
					TrimRead tr1=null, tr2=null;
					
					boolean remove=false;
					if(qtrim){
						if(untrim){
							if(r1!=null){
								tr1=TrimRead.trim(r1, qtrimLeft, qtrimRight, trimq, 1);
								int x=(tr1==null ? 0 : tr1.leftTrimmed+tr1.rightTrimmed);
								basesTrimmedT+=x;
								readsTrimmedT+=(x>0 ? 1 : 0);
							}
							if(r2!=null){
								tr2=TrimRead.trim(r2, qtrimLeft, qtrimRight, trimq, 1);
								int x=(tr2==null ? 0 : tr2.leftTrimmed+tr2.rightTrimmed);
								basesTrimmedT+=x;
								readsTrimmedT+=(x>0 ? 1 : 0);
							}
						}else{
							if(r1!=null){
								int x=TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq, 1);
								basesTrimmedT+=x;
								readsTrimmedT+=(x>0 ? 1 : 0);
							}
							if(r2!=null){
								int x=TrimRead.trimFast(r2, qtrimLeft, qtrimRight, trimq, 1);
								basesTrimmedT+=x;
								readsTrimmedT+=(x>0 ? 1 : 0);
							}
						}
					}

					if(minReadLength>0 && !remove){
						int rlen=(r1==null ? 0 : r1.length());
						int rlen2=(r1.mateLength());
						if(rlen<minReadLength && rlen2<minReadLength){
							basesTrimmedT+=(rlen+rlen2);
							remove=true;
						}
					}

					if(!remove){

						//				verbose=(r.numericID==727 || r.numericID==1364);

						if(r2!=null){r2.reverseComplement();}
						readsProcessed++;


						//				System.err.println("True Insert: "+r.insertSize());

						final int[] rvector=new int[6];
						int trueSize=r1.insertSizeMapped(ignoreMappingStrand);

						int bInsert=-1, hInsert=-1;
						int bestInsert;
						int bestScore=-1;
						int bestGood=-1;
						int bestBad=999999, bBad=999999, hBad=999999;
						boolean ambig, tooShort=false;
						boolean bAmbig=true, hAmbig=true;
						int bestVotes=-1;

						boolean didb=false, didh=false;

						//					assert(false) : r+"\n"+(USE_MAPPING)+", "+(r.chrom==r.mate.chrom)+", "+()+", "+()+", "+()+", "+()+", ";

						if(r2==null){
							bestScore=100;
							bestGood=30;
							bestBad=0;
							bestInsert=r1.length();
							assert(r1.length()==r1.insert()) : r1.length()+" != "+r1.insert()+"; actual = "+trueSize;
							//						if(bestInsert!=trueSize){
							//							System.err.println("Bad insert size for pre-joined read "+r.numericID+": len="+r.length()+", insert="+r.insert()+", actual="+trueSize);
							//						}
							bestVotes=1;
							ambig=false;
						}else if(USE_MAPPING && r1.chrom==r2.chrom && r1.start<r1.stop && ((r1.mapped() || r1.synthetic()) && (r2.mapped() || r2.synthetic()))){
							bestScore=100;
							bestGood=30;
							bestBad=0;
							bestInsert=trueSize;
							bestVotes=1;
							ambig=false;
						}else if(SKIP_MATED_READS && r1.insertvalid() && r1.insert()>0){
							bestScore=100;
							bestGood=30;
							bestBad=0;
							bestInsert=r1.insert();
							bestVotes=1;
							ambig=false;
						}else{
							if(MATE_BY_OVERLAP){
								didb=true;
								bInsert=mateByOverlap(r1, r2, rvector, MIN_OVERLAPPING_BASES_0, MIN_OVERLAPPING_BASES);
								bestScore=rvector[0];
								bestGood=rvector[1];
								bBad=rvector[2];
								bAmbig=(rvector[4]==1);
								bestVotes=rvector[5];
								final int len1=r1.length(), len2=r2.length();
								for(int trims=0, q=trimq; trims<TRIM_ON_OVERLAP_FAILURE && !qtrim && bInsert<0 /*&& !bAmbig*/; trims++, q+=8){
//									System.err.println(trims+", "+q);
									Serializable old1=r1.obj;
									Serializable old2=r2.obj;
									tr1=TrimRead.trim(r1, false, true, q, 1+len1*4/10); //r1.length());
									tr2=TrimRead.trim(r2, true, false, q, 1+len2*4/10); //r2.length());
									r1.obj=old1;
									r2.obj=old2;
									if(tr1!=null || tr2!=null){
//										System.err.println(r1.length()+", "+r2.length());
										int x=mateByOverlap(r1, r2, rvector, MIN_OVERLAPPING_BASES_0-1, MIN_OVERLAPPING_BASES);
										if(x>-1){
//											System.err.println(trims);
											bInsert=x;
											bestScore=rvector[0];
											bestGood=rvector[1];
											bBad=rvector[2];
											bAmbig=(rvector[4]==1);
											bestVotes=rvector[5];
											trims=TRIM_ON_OVERLAP_FAILURE;
										}else{
											if(tr1!=null){tr1.untrim();}
											if(tr2!=null){tr2.untrim();}
										}
									}
								}
							}
							if(kca!=null && kca.length>0 && (bAmbig || bInsert<0 || bestBad>0 || bestGood<30)){
								didh=true;
								hInsert=mateRead(r1, r2, k1, k2, mask1, mask2, kca, rvector);
								bestScore=rvector[0];
								bestGood=rvector[1];
								hBad=rvector[2];
								hAmbig=(rvector[4]==1);
								bestVotes=rvector[5];
							}

							if(hInsert==bInsert){
								bestInsert=hInsert;
								bestBad=Tools.min(bBad, hBad);
								ambig=bAmbig && hAmbig;
							}
							//						else if(!didb || bAmbig || bInsert<0){
							//							bestInsert=hInsert;
							//							ambig=hAmbig;
							//						}
							else if(!didh || hAmbig || hInsert<0){
								bestInsert=bInsert;
								bestBad=bBad;
								ambig=bAmbig;
							}else{
								bestInsert=hInsert;
								bestBad=hBad;
								ambig=hAmbig;
							}
						}

						tooShort=(!ambig && bestInsert>0 && bestInsert<minInsert);
						
						if(joinperfectonly && bestBad>0){ambig=true;}

						if(ambig){ambiguousCount++;}
						else if(tooShort){tooShortCount++;}
						else if(bestInsert==-1){noSolutionCount++;}
						else if(bestInsert==trueSize){correctCount++;insertSumCorrect+=bestInsert;}
						else{incorrectCount++;insertSumIncorrect+=bestInsert;}


						if(bestInsert>-1 && !ambig){
							insertMin=Tools.min(bestInsert, insertMin);
							insertMax=Tools.max(bestInsert, insertMax);
							hist[Tools.min(bestInsert, hist.length-1)]++;
						}
						r1.setInsert(ambig ? -1 : bestInsert);

						{//Clear memory.
							if(r1.obj!=null){assert(r1.obj.getClass()==long[].class) : r1.obj.getClass();}
							r1.obj=null;
							if(r2!=null){r2.obj=null;}
						}
						
						if(OUTPUT_FAILED || bestInsert>-1){

							if(untrim && (ambig || bestInsert<0 || !joinReads)){
								if(tr1!=null){tr1.untrim();}
								if(tr2!=null){tr2.untrim();}
							}
							
							if((ambig || bestInsert<0 || tooShort) && (rosbad!=null || !MIX_BAD_AND_GOOD)){
								if(rosbad!=null){
									if(listb!=null){listb.add(r1);}
								}
							}else{
								if(listg!=null){
									Read x=r1;
									if(joinReads && r2!=null){
										x=r1.joinRead(bestInsert);
										//Disabled because ErrorCorrectMT was retired.
//										if(middleTable!=null && x.containsNocalls()){
//											BitSet bs=ErrorCorrectMT.detectNBulk(x);
//											ErrorCorrectMT.correctErrorsBothSides(x, middleTable, MIDDLE_TABLE_K, MIN_HITS_FOR_GOOD, MAX_HITS_FOR_BAD, bs, 9999999);
//										}
									}
									listg.add(x);
								}
							}
							
							if(rosinsert!=null){
								StringBuilder sb=new StringBuilder(40);
								sb.append(r1.id==null ? r1.numericID+"" : r1.id).append('\t');
								sb.append(r1.numericID).append('\t');

								sb.append(bestInsert);
								sb.append('\t');

								if(bestInsert<0){sb.append('F');}//Failed
								else if(ambig){sb.append('A');} //Ambiguous
								else if(tooShort){sb.append('S');} //Short
								else if(bestInsert>0 && bestBad<1){sb.append('P');} //Perfect
								else{sb.append('I');}//Imperfect

								if(bestInsert>0){
									sb.append("\t"+bestGood+"\t"+bestBad+"\t"+bestScore+"\t"+bestVotes);
								}
								r1.obj=sb;
								listi.add(r1);
							}
						}

						//				if(bestInsert!=trueSize && bestInsert>0 && !ambig){
						//					System.err.println("\nIncorrect answer for read "+r.numericID+"; true insert = "+trueSize+", called at "+bestInsert);
						////					verbose=true;
						//					for(int i=0; i<300; i++){
						//						int x=testRead(r, r.mate, k1, k2, mask1, mask2, kca, rvector, i);
						//						if((x>0 && rvector[2]<=bestBad) || i==trueSize || i==bestInsert){
						//							verbose=true;
						//							testRead(r, r.mate, k1, k2, mask1, mask2, kca, rvector, i);
						//							verbose=false;
						//						}
						//					}
						////					verbose=false;
						//				}

						//				assert(r.numericID<200);
						//				assert(false);
						if(r2!=null){r2.reverseComplement();}
					}
				}

				if(rosgood!=null){rosgood.add(listg, ln.id);}
				if(rosbad!=null){rosbad.add(listb, ln.id);}
				if(rosinsert!=null){rosinsert.add(listi, ln.id);}

				//			System.err.println("returning list");
				cris.returnList(ln.id, ln.list.isEmpty());
				//			System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
				//			System.err.println("reads: "+(reads==null ? "null" : reads.size()));
			}
			cris.returnList(ln.id, ln.list.isEmpty());
		}

		int[] hist=new int[1000];

		long readsProcessed=0;
		long matedCount=0;
		long correctCount=0;
		long ambiguousCount=0;
		long tooShortCount=0;
		long incorrectCount=0;
		long noSolutionCount=0;
		long insertSumCorrect=0;
		long insertSumIncorrect=0;
		int insertMax=0;
		int insertMin=999999999;
		
		long basesTrimmedT=0;
		long readsTrimmedT=0;
		
		private final ConcurrentReadInputStream cris;
		private final ConcurrentReadOutputStream rosgood;
		private final ConcurrentReadOutputStream rosbad;
		private final ConcurrentReadOutputStream rosinsert;
		private final int k;
		private final KCountArray[] kca;
		private final boolean joinReads;
		private final boolean joinperfectonly;
		private final KCountArray middleTable;
	}
	
	private String in1_primary;
	private String in2_primary;
	
	private String outgood_G=null;
	private String outbad_G=null;
	private String outinsert_G=null;
	private String outhist=null;
	private String outhist2=null;
	private String outhist3=null;
	
	private long maxReads_G=-1;
	private long tableReads_G=-1;
	private int k_G=K_DEFAULT;
	private int[][] gap_G=null;
	private int cbits_G=2;
	private long totalcells_G=-1;
	private int hashes=2;
	private int passes_G=1;
	private String tempfile="mateReadsTemp#.txt.gz";
	private boolean join_G=true;
	private int maxtables=0;
	private boolean auto=true;
	
	private List<String> extra_G=null;
	
	static boolean errorState=false;
	
	static boolean qtrimRight=false;
	static boolean qtrimLeft=false;
	static boolean untrim=false;
	static byte trimq=6;
	static int minReadLength=0;
	static int minInsert=0;
	static boolean qtrim=false;
	static int TRIM_ON_OVERLAP_FAILURE=1;
	
	
	static int[] histTotal=new int[1000];
	static int bin=1;

	static long readsProcessedTotal=0;
	static long matedCountTotal=0;
	static long correctCountTotal=0;
	static long ambiguousCountTotal=0;
	static long tooShortCountTotal=0;
	static long incorrectCountTotal=0;
	static long noSolutionCountTotal=0;
	static long insertSumCorrectTotal=0;
	static long insertSumIncorrectTotal=0;
	static long basesTrimmedTotal=0;
	static long readsTrimmedTotal=0;
	static int insertMinTotal=999999999;
	static int insertMaxTotal=0;

	public static int MIN_OVERLAPPING_KMERS=10;
	public static int MIN_OVERLAPPING_KMERS_0=4;
	public static int MIN_OVERLAPPING_BASES=12;
	public static int MIN_OVERLAPPING_BASES_0=8;
	public static int MIN_OVERLAP_INSERT=16;
	public static int ACCEL_DIV=10; //Acceleration is actually proportional to inverse of this number.
	public static int ACCEL_FACTOR=ACCEL_DIV; //Max distance between samples
	public static int DEFAULT_BADLIMIT=25;
	public static int DEFAULT_BADLIMIT_FOR_BASE_MATCHING=3;
	public static int DEFAULT_MISMATCHLIMIT=6;
	public static int MIN_HITS_FOR_GOOD=3;
	public static int MAX_HITS_FOR_BAD=1;
	public static int MIN_VOTES=1;
	public static int K_DEFAULT=29;
	public static int MIDDLE_TABLE_K=31;
	public static byte MIN_QUALITY=8;
	public static byte MIN_QUALITY_FOR_OVERLAP=7;
	/** Skip alignment and calculate insert from mapping info */ 
	public static boolean USE_MAPPING=false;
	public static boolean MATE_BY_OVERLAP=true;
	public static boolean SKIP_MATED_READS=false;
	public static boolean OUTPUT_FAILED=true;
	public static boolean MIX_BAD_AND_GOOD=false;
	public static boolean WRITE_INTERMEDIATE_JOINED=false;
	public static boolean FILL_MIDDLE_INTERMEDIATE=false;
	public static boolean FILL_MIDDLE_FINAL=false;
	public static boolean overwrite=true;
	public static boolean append=false;
	public static boolean verbose=false;
	public static boolean ignoreMappingStrand=false;
	
	public static int THREADS=-1;
	public static float version=2.0f;
	
}

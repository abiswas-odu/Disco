package clump;

import java.io.File;
import java.io.PrintStream;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.SynchronousQueue;
import java.util.concurrent.atomic.AtomicInteger;

import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.Read;
import structures.ListNum;
import sun.reflect.Reflection;
import dna.Parser;
import fileIO.ByteFile;
import fileIO.ReadWrite;
import jgi.BBMerge;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sort.ReadComparatorID;
import sort.ReadComparatorName;
import fileIO.FileFormat;
import bloom.KCountArray;

/**
 * @author Brian Bushnell
 * @date June 20, 2014
 *
 */
public class KmerSort3 {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		assert(false);
		main(-1, args);
	}
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(long fileMem, String[] args){
		final boolean pigz=ReadWrite.USE_PIGZ, unpigz=ReadWrite.USE_UNPIGZ;
		final float ztd=ReadWrite.ZIP_THREAD_MULT;
		final int mzt=ReadWrite.MAX_ZIP_THREADS;
		final int oldzl=ReadWrite.ZIPLEVEL;
		Timer t=new Timer();
		KmerSort3 ks=new KmerSort3(fileMem, args);
		ks.process(t);
		ReadWrite.USE_PIGZ=pigz;
		ReadWrite.USE_UNPIGZ=unpigz;
		ReadWrite.ZIP_THREAD_MULT=ztd;
		ReadWrite.MAX_ZIP_THREADS=mzt;
		ReadWrite.ZIPLEVEL=oldzl;
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public KmerSort3(long fileMem_, String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		boolean setInterleaved=false; //Whether it was explicitly set.
		
		Shared.READ_BUFFER_LENGTH=Tools.min(200, Shared.READ_BUFFER_LENGTH);
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
				verbose=KmerComparator.verbose=Tools.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
				assert(k>0 && k<32);
			}else if(a.equals("mincount") || a.equals("mincr")){
				minCount=Integer.parseInt(b);
			}else if(a.equals("rename") || a.equals("addname")){
				addName=Tools.parseBoolean(b);
			}else if(a.equals("shortname") || a.equals("shortnames")){
				if(b!=null && b.equals("shrink")){
					shrinkName=true;
				}else{
					shrinkName=false;
					shortName=Tools.parseBoolean(b);
				}
			}else if(a.equals("rcomp") || a.equals("reversecomplement")){
				rcomp=Tools.parseBoolean(b);
			}else if(a.equals("ecco")){
				ecco=Tools.parseBoolean(b);
			}else if(a.equals("condense") || a.equals("consensus") || a.equals("concensus")){//Note the last one is intentionally misspelled
				condense=Tools.parseBoolean(b);
			}else if(a.equals("correct") || a.equals("ecc")){
				correct=Tools.parseBoolean(b);
			}else if(a.equals("passes")){
				passes=Integer.parseInt(b);
			}
			
			else if(a.equals("dedupe")){
				dedupe=Tools.parseBoolean(b);
			}else if(a.equals("markduplicates")){
				dedupe=ClumpList.markOnly=Tools.parseBoolean(b);
			}else if(a.equals("markall")){
				boolean x=Tools.parseBoolean(b);
				if(x){
					dedupe=ClumpList.markOnly=ClumpList.markAll=true;
				}else{
					ClumpList.markAll=false;
				}
			}else if(a.equals("removeallduplicates") || a.equals("allduplicates")){
				ClumpList.markAll=Tools.parseBoolean(b);
			}else if(a.equals("optical") || a.equals("opticalonly")){
				ClumpList.opticalOnly=Tools.parseBoolean(b);
			}else if(a.equals("dupesubs") || a.equals("duplicatesubs") || a.equals("dsubs") || a.equals("subs")){
				ClumpList.maxSubstitutions=Integer.parseInt(b);
			}else if(a.equals("dupedist") || a.equals("duplicatedistance") || a.equals("ddist") || a.equals("dist") || a.equals("opticaldist")){
				ClumpList.maxOpticalDistance=Float.parseFloat(b);
				ClumpList.opticalOnly=ClumpList.maxOpticalDistance>=0;
			}else if(a.equals("scanlimit") || a.equals("scan")){
				ClumpList.scanLimit=Integer.parseInt(b);
			}
			
			else if(a.equals("prefilter")){
				KmerReduce.prefilter=Tools.parseBoolean(b);
			}else if(a.equals("groups") || a.equals("g") || a.equals("sets") || a.equals("ways")){  
				groups=Integer.parseInt(b);
				splitInput=(groups>1);
			}else if(a.equals("seed")){
				KmerComparator.defaultSeed=Long.parseLong(b);
			}else if(a.equals("hashes")){
				KmerComparator.setHashes(Integer.parseInt(b));
			}else if(a.equals("border")){
				KmerComparator.defaultBorder=Integer.parseInt(b);
			}else if(a.equals("minprob")){
				KmerComparator.minProb=Float.parseFloat(b);
				
			}else if(a.equals("unpair")){
				unpair=Tools.parseBoolean(b);
			}else if(a.equals("repair")){
				repair=Tools.parseBoolean(b);
			}else if(a.equals("namesort") || a.equals("sort")){
				namesort=Tools.parseBoolean(b);
			}else if(a.equals("fetchthreads")){
				fetchThreads=Integer.parseInt(b);
				assert(fetchThreads>0) : fetchThreads+"\nFetch threads must be at least 1.";
			}else if(a.equals("reorder") || a.equals("reorderclumps")){
				//reorder=Tools.parseBoolean(b);
			}else if(a.equals("reorderpaired") || a.equals("reorderclumpspaired")){
//				reorderpaired=Tools.parseBoolean(b);
			}
			
			else if(Clump.parseStatic(arg, a, b)){
				//Do nothing
			}
			
			else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		Clump.rename=condense;
		if(dedupe){KmerComparator.compareSequence=true;}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			setInterleaved=parser.setInterleaved;

			in1=parser.in1;
			in2=parser.in2;

			out1=parser.out1;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
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
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		if(out1==null){ffout=null;}
		else{
			int g=out1.contains("%") ? groups : 1;
			ffout=new FileFormat[g];
			if(g==1){
				ffout[0]=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
			}else{
				ReadWrite.ZIPLEVEL=2;
				ReadWrite.setZipThreadMult(Tools.min(0.5f, 2f/(g+1)));
				for(int i=0; i<g; i++){
					ffout[i]=FileFormat.testOutput(out1.replaceFirst("%", ""+i), FileFormat.FASTQ, extout, (g<10), overwrite, append, false);
				}
			}
		}
		
		if(groups>1 && in1.contains("%") && (splitInput || !new File(in1).exists())){
			ffin1=new FileFormat[groups];
			ffin2=new FileFormat[groups];
			for(int i=0; i<groups; i++){
				ffin1[i]=FileFormat.testInput(in1.replaceFirst("%", ""+i), FileFormat.FASTQ, extin, true, true);
				ffin2[i]=in2==null ? null : FileFormat.testInput(in2.replaceFirst("%", ""+i), FileFormat.FASTQ, extin, true, true);
			}
		}else{
			assert(!in1.contains("%") && groups==1) : "The % symbol must only be present in the input filename if groups>1.";
			ffin1=new FileFormat[1];
			ffin1[0]=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
			ffin2=new FileFormat[1];
			ffin2[0]=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
			groups=1;
		}
		
		long sizeSum=0, expectedMemSum=0;
		for(FileFormat ff : ffin1){
			long x=new File(ff.name()).length();
			sizeSum+=x;
			if(ff.compressed()){x*=40;}else{x*=8;}
			expectedMemSum+=x;
		}
		for(FileFormat ff : ffin2){
			if(ff!=null){
				long x=new File(ff.name()).length();
				sizeSum+=x;
				if(ff.compressed()){x*=40;}else{x*=8;}
				expectedMemSum+=x;
			}
		}

		expectedSizePerGroup=(sizeSum+groups+1)/(Tools.max(groups, 1));
		expectedMemPerGroup=(expectedMemSum+groups+1)/(Tools.max(groups, 1));
		totalMem=Shared.memAvailable(1);
		fileSize=sizeSum;
		fileMem=fileMem_<1 ? 40*fileSize : fileMem_;
		memRatio=fileMem*1.0/Tools.max(1, fileSize);
		
//		if(groups>1){ReadWrite.USE_UNPIGZ=false;} //Not needed since they are not concurrent
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Count kmers */
	void preprocess(){
		if(minCount>1){
			if(groups>1){
				table=ClumpTools.table();
				assert(table!=null);
			}else{
				Timer ctimer=new Timer();
				if(verbose){ctimer.start("Counting pivots.");}
				table=ClumpTools.getTable(in1, in2, k, minCount);
				if(verbose){ctimer.stop("Count time: ");}
			}
		}
	}

	/** Create read streams and process all data */
	void process(Timer t){
		
		preprocess();

		final ConcurrentReadOutputStream rosa[]=(ffout==null ? null : new ConcurrentReadOutputStream[ffout.length]);
		for(int i=0; rosa!=null && i<rosa.length; i++){
			final int buff=1;

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			rosa[i]=ConcurrentReadOutputStream.getStream(ffout[i], null, null, null, buff, null, false);
			rosa[i].start();
		}
		
		readsProcessed=basesProcessed=diskProcessed=memProcessed=0;
		
		//Process the read stream
		processInner(rosa);
		lastMemProcessed=memThisPass;
		
		table=null;
		ClumpTools.clearTable();
		
		errorState|=ReadStats.writeAll();
		
		t.stop();
		
		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);

		String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
		String rpstring2=readsProcessed+"";
		String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");
		
		String cpstring=""+(groups==1 ? clumpsProcessedThisPass : clumpsProcessedTotal);
		String epstring=""+correctionsTotal;
		String dpstring=""+duplicatesTotal;

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}

		while(rpstring2.length()<10){rpstring2=" "+rpstring2;}
		while(cpstring.length()<10){cpstring=" "+cpstring;}
		while(epstring.length()<10){epstring=" "+epstring;}
		while(dpstring.length()<10){dpstring=" "+dpstring;}
		
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		outstream.println();

		outstream.println("Reads In:         "+rpstring2);
		outstream.println("Clumps Formed:    "+cpstring);
		if(correct){
			outstream.println("Errors Corrected: "+epstring);
		}
		if(dedupe){
			outstream.println("Duplicates Found: "+dpstring);
		}
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Collect and sort the reads */
	void processInner(final ConcurrentReadOutputStream rosa[]){
		if(verbose){outstream.println("Making comparator.");}
		KmerComparator kc=new KmerComparator(k, addName, (rcomp || condense || correct));
		
		ClumpList.UNRCOMP=(!rcomp && !condense);
		Timer t=new Timer();
		
		fetchThreads=Tools.min(groups, fetchThreads, ffin1.length);
		assert(fetchThreads>0);
		SynchronousQueue<ArrayList<Read>> listQ=new SynchronousQueue<ArrayList<Read>>();
		ArrayList<FetchThread> alft=fetchReads(kc, fetchThreads, listQ, rosa);
		int poisonCount=0;
		
		for(int group=0; group<groups; group++){
			
			if(verbose){t.start("Fetching reads.");}
			ArrayList<Read> reads=null;
			while(poisonCount<fetchThreads && reads==null){
				try {
					reads=listQ.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(reads==POISON){
					poisonCount++;
					reads=null;
				}
			}
			if(verbose){t.stop("Fetch time: ");}
			
			//Done by FetchThread
//			if(verbose){t.start("Sorting.");}
//			Shared.sort(reads, kc);
//			if(verbose){t.stop("Sort time: ");}
			
			if(verbose){t.start("Making clumps.");}
			readsProcessedThisPass=reads.size();
			ClumpList cl=new ClumpList(reads, k, false);
			clumpsProcessedThisPass=cl.size();
			clumpsProcessedTotal+=clumpsProcessedThisPass;
			if(verbose){t.stop("Clump time: ");}
			
			if(dedupe){
				reads.clear();
				if(verbose){t.start("Deduping.");}
				reads=processClumps(cl, ClumpList.DEDUPE);
				if(verbose){t.stop("Dedupe time: ");}
			}else if(condense){
				reads.clear();
				if(verbose){t.start("Condensing.");}
				reads=processClumps(cl, ClumpList.CONDENSE);
				if(verbose){t.stop("Condense time: ");}
			}else if(correct){
				reads.clear();
				if(verbose){t.start("Correcting.");}
				reads=processClumps(cl, ClumpList.CORRECT);
				if(verbose){t.stop("Correct time: ");}
				
				if(verbose){outstream.println("Seed: "+kc.seed);}
				if(groups>1){
					if(verbose){outstream.println("Reads:        \t"+readsProcessedThisPass);}
					outstream.println("Clumps:       \t"+clumpsProcessedThisPass);
					if(correct){
						outstream.println("Corrections:  \t"+correctionsThisPass);
					}
					outstream.println();
				}
				
				if(passes>1 && groups==1){
					if(verbose){outstream.println("Pass 1.");}
					if(verbose){outstream.println("Reads:        \t"+readsProcessedThisPass);}
					outstream.println("Clumps:       \t"+clumpsProcessedThisPass);
					if(correct){
						outstream.println("Corrections:  \t"+correctionsThisPass);
					}
					outstream.println();
					
					for(int pass=1; pass<passes; pass++){
						kc=new KmerComparator(k, kc.seed<0 ? -1 : kc.seed+1, kc.border-1, kc.hashes, false, kc.rcompReads);
						reads=runOnePass(reads, kc);

						if(verbose){outstream.println("Seed: "+kc.seed);}
						if(verbose){outstream.println("Pass "+(pass+1)+".");}
						if(verbose){outstream.println("Reads:        \t"+readsProcessedThisPass);}
						outstream.println("Clumps:       \t"+clumpsProcessedThisPass);
						if(correct){
							outstream.println("Corrections:  \t"+correctionsThisPass);
						}
						outstream.println();
					}
				}
			}
			
			if(repair || namesort){
				if(groups>1){
					if(verbose){t.start("Name-sorting.");}
					reads=nameSort(reads, false);
					if(verbose){t.stop("Sort time: ");}
				}else{
					if(namesort){
						if(verbose){t.start("Name-sorting.");}
						reads=idSort(reads, repair);
						if(verbose){t.stop("Sort time: ");}
					}else{
						reads=read1Only(reads);
					}
				}
			}
			
			if(doHashAndSplit || groups==0){
				addToRos(rosa, reads, t, kc);
			}else{
				if(group>0){
					ConcurrentReadOutputStream ros=rosa[group-1];
					errorState|=ReadWrite.closeStream(ros);
				}
				rosa[group].add(reads, 0);
			}
			reads=null;
		}
		
		if(verbose){outstream.println("Closing fetch threads.");}
		while(poisonCount<fetchThreads){
			ArrayList<Read> reads=null;
			try {
				reads=listQ.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			if(reads==POISON){
				poisonCount++;
				reads=null;
			}
		}
		long readsThisPass=closeFetchThreads(alft);
		if(verbose){outstream.println("Closed fetch threads.");}
		
		if(rosa!=null){
			if(verbose){outstream.println("Waiting for writing to complete.");}
			for(ConcurrentReadOutputStream ros : rosa){
				errorState=ReadWrite.closeStream(ros)|errorState;
			}
			if(verbose){t.stop("Write time: ");}
		}
		
		if(verbose){outstream.println("Done!");}
	}
	
	private void addToRos(ConcurrentReadOutputStream[] rosa, ArrayList<Read> list, Timer t, KmerComparator old){
		if(rosa==null){return;}
		assert(rosa.length>0);
		if(rosa.length==1){
			if(verbose){t.start("Writing.");}
			rosa[0].add(list, 0);
			return;
		}
		KmerComparator kc=new KmerComparator(old.k, old.seed+1, old.border-1, old.hashes, false, false);
		final int div=rosa.length;
		assert(div==groups);
		@SuppressWarnings("unchecked")
		ArrayList<Read>[] array=new ArrayList[div];
		for(int i=0; i<array.length; i++){
			array[i]=new ArrayList<Read>();
		}
		if(verbose){t.start("Splitting.");}
		hashAndSplit(list, kc, array);
		if(verbose){t.stop("Split time: ");}
		if(verbose){t.start("Writing.");}
		for(int i=0; i<div; i++){
			rosa[i].add(array[i], 0);
			array[i]=null;
		}
		if(verbose){System.err.println("Sent writable reads.");}
	}
	
	private ArrayList<Read> runOnePass(ArrayList<Read> reads, KmerComparator kc){
//		for(Read r : reads){r.obj=null;}//No longer necessary
		
		Timer t=new Timer();
		
		table=null;
		if(minCount>1){
			if(verbose){t.start("Counting pivots.");}
			table=ClumpTools.getTable(reads, k, minCount);
			if(verbose){t.stop("Count time: ");}
		}
		
		if(verbose){t.start("Hashing.");}
		kc.hashThreaded(reads, table, minCount);
		if(verbose){t.stop("Hash time: ");}
		
		if(verbose){t.start("Sorting.");}
		Shared.sort(reads, kc);
		if(verbose){t.stop("Sort time: ");}
		
		if(verbose){t.start("Making clumps.");}
		readsProcessedThisPass=reads.size();
		ClumpList cl=new ClumpList(reads, k, false);
		reads.clear();
		clumpsProcessedThisPass=cl.size();
		clumpsProcessedTotal+=clumpsProcessedThisPass;
		if(verbose){t.stop("Clump time: ");}
		
		if(verbose){t.start("Correcting.");}
		reads=processClumps(cl, ClumpList.CORRECT);
		if(verbose){t.stop("Correct time: ");}
		
		return reads;
	}
	
	public void hashAndSplit(ArrayList<Read> list, KmerComparator kc, ArrayList<Read>[] array){
		int threads=Shared.threads();
		ArrayList<HashSplitThread> alt=new ArrayList<HashSplitThread>(threads);
		for(int i=0; i<threads; i++){alt.add(new HashSplitThread(i, threads, list, kc));}
		for(HashSplitThread ht : alt){ht.start();}
		
		/* Wait for threads to die */
		for(HashSplitThread ht : alt){
			
			/* Wait for a thread to die */
			while(ht.getState()!=Thread.State.TERMINATED){
				try {
					ht.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			for(int i=0; i<groups; i++){
				array[i].addAll(ht.array[i]);
				ht.array[i]=null;
			}
		}
	}
	
	private class HashSplitThread extends Thread{
		
		@SuppressWarnings("unchecked")
		HashSplitThread(int id_, int threads_, ArrayList<Read> list_, KmerComparator kc_){
			id=id_;
			threads=threads_;
			list=list_;
			kc=kc_;
			array=new ArrayList[groups];
			for(int i=0; i<groups; i++){
				array[i]=new ArrayList<Read>();
			}
		}
		
		@Override
		public void run(){
			for(int i=id; i<list.size(); i+=threads){
				Read r=list.get(i);
				kc.hash(r, null, 0, true);
				ReadKey key=(ReadKey)r.obj;
				array[(int)(kc.hash(key.kmer)%groups)].add(r);
			}
		}
		
		final int id;
		final int threads;
		final ArrayList<Read> list;
		final KmerComparator kc;
		final ArrayList<Read>[] array;
	}
	
	private ArrayList<Read> nameSort(ArrayList<Read> list, boolean pair){
		Shared.sort(list, ReadComparatorName.comparator);
		if(!pair){return list;}
		
		ArrayList<Read> list2=new ArrayList<Read>(1+list.size()/2);
		
		Read prev=null;
		for(Read r : list){
			if(prev==null){
				prev=r;
				assert(prev.mate==null);
			}else{
				if(prev.id.equals(r.id) || FASTQ.testPairNames(prev.id, r.id, true)){
					prev.mate=r;
					r.mate=prev;
					prev.setPairnum(0);
					r.setPairnum(1);
					list2.add(prev);
					prev=null;
				}else{
					list2.add(prev);
					prev=r;
				}
			}
		}
		return list2;
	}
	
	private ArrayList<Read> idSort(ArrayList<Read> list, boolean pair){
		Shared.sort(list, ReadComparatorID.comparator);
		if(!pair){return list;}
		
		ArrayList<Read> list2=new ArrayList<Read>(1+list.size()/2);
		
		Read prev=null;
		for(Read r : list){
			if(prev==null){
				prev=r;
				assert(prev.mate==null);
			}else{
				if(prev.numericID==r.numericID){
					assert(prev.pairnum()==0 && r.pairnum()==1) : prev.id+"\n"+r.id;
					prev.mate=r;
					r.mate=prev;
					prev.setPairnum(0);
					r.setPairnum(1);
					list2.add(prev);
					prev=null;
				}else{
					list2.add(prev);
					prev=r;
				}
			}
		}
		return list2;
	}
	
	private ArrayList<Read> read1Only(ArrayList<Read> list){
		ArrayList<Read> list2=new ArrayList<Read>(1+list.size()/2);
		for(Read r : list){
			assert(r.mate!=null);
			if(r.pairnum()==0){list2.add(r);}
		}
		return list2;
	}
	
	@Deprecated
	//No longer needed
	public int countClumps(ArrayList<Read> list){
		int count=0;
		long currentKmer=-1;
		for(final Read r : list){
			final ReadKey key=(ReadKey)r.obj;
			if(key.kmer!=currentKmer){
				currentKmer=key.kmer;
				count++;
			}
		}
		return count;
	}
	
	public ArrayList<FetchThread> fetchReads(final KmerComparator kc, final int fetchThreads, SynchronousQueue<ArrayList<Read>> listQ, ConcurrentReadOutputStream[] rosa){
		AtomicInteger nextGroup=new AtomicInteger(0);
		if(verbose){outstream.println("Making "+fetchThreads+" fetch thread"+(fetchThreads==1 ? "." : "s."));}
		ArrayList<FetchThread> alft=new ArrayList<FetchThread>(fetchThreads);
		for(int i=0; i<fetchThreads; i++){alft.add(new FetchThread(kc, listQ, nextGroup, rosa));}
		
		if(verbose){outstream.println("Starting threads.");}
		for(FetchThread ft : alft){ft.start();}
		
		return alft;
	}
	

	private long closeFetchThreads(ArrayList<FetchThread> alft){
		readsThisPass=0;
		memThisPass=0;
		/* Wait for threads to die */
		for(FetchThread ft : alft){

			/* Wait for a thread to die */
			while(ft.getState()!=Thread.State.TERMINATED){
				try {
					ft.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			readsThisPass+=ft.readsProcessedT;
			basesProcessed+=ft.basesProcessedT;
			diskProcessed+=ft.diskProcessedT;
			memThisPass+=ft.memProcessedT;
			
			errorState|=ft.errorState;
		}
		readsProcessed+=readsThisPass;
		memProcessed+=memThisPass;
		
		ecco=false;
		return readsThisPass;
	}
	
	public ArrayList<Read> processClumps(ClumpList cl, int mode){
		long[] rvector=new long[2];
		ArrayList<Read> out=cl.process(Shared.threads(), mode, rvector);
		correctionsThisPass=rvector[0];
		correctionsTotal+=correctionsThisPass;
		duplicatesThisPass=rvector[1];
		duplicatesTotal+=duplicatesThisPass;
		cl.clear();
		return out;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This is called if the program runs with no parameters */
	private void printOptions(){
		throw new RuntimeException("TODO");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class FetchThread extends Thread{
		
		FetchThread(final KmerComparator kc_, SynchronousQueue<ArrayList<Read>> listQ_, AtomicInteger nextGroup_, ConcurrentReadOutputStream[] rosa_){
			kc=kc_;
			listQ=listQ_;
			nextGroup=nextGroup_;
			rosa=rosa_;
		}
		
		@Override
		public void run(){
			for(ArrayList<Read> reads=fetchNext(); reads!=null; reads=fetchNext()){
				boolean success=false;
				while(!success){
					try {
						listQ.put(reads);
						success=true;
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			boolean success=false;
			while(!success){
				try {
					listQ.put(POISON);
					success=true;
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		
		private ArrayList<Read> fetchNext(){
			final int group=nextGroup.getAndIncrement();
			if(group>=groups){return null;}
			
			final ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ffin1[group], ffin2[group], null, null);
			cris.start();
			
			//Check for file size imbalance
			if(!Clump.forceProcess){
				final long size;
				double expectedMem;
				{
					size=new File(ffin1[group].name()).length()+(ffin2[group]==null ? 0 : new File(ffin2[group].name()).length());
//					expectedMem=size*(ffin1[group].compressed() ? 40 : 8);
					expectedMem=(fileMem*(double)size)/fileSize;
				}
				
				if(size-1000>expectedSizePerGroup*3 || expectedMem*3>totalMem){
					
//					outstream.println("size="+size+", expectedSizePerGroup="+expectedSizePerGroup+", expectedMem="+expectedMem+", expectedSizePerGroup="+expectedMemPerGroup);
//					outstream.println("totalMem="+totalMem+", Shared.memFree()="+Shared.memFree()+", Shared.memUsed()="+Shared.memUsed());
//					outstream.println("fileSize="+fileSize+", "+"fileMem="+fileMem+", "+"memRatio="+memRatio+", ");
//					outstream.println((size-1000>expectedSizePerGroup*3)+", "+(expectedMem*3>totalMem));
					
					//TODO: This could also be based on the number of FetchThreads.
					
					if(expectedMem>0.3*totalMem && (fileMem<1 || fileMem>0.8*totalMem)){
						outstream.println("\n***  Warning  ***\n"
								+ "A temp file may be too big to store in memory, due to uneven splitting:");

						outstream.println("expectedMem="+((long)expectedMem)+", fileMem="+fileMem+", available="+totalMem);

						if(repair || namesort){
							outstream.println("It cannot be streamed to output unaltered because "+(namesort ? " namesort=t" : "repair=t"));
							outstream.println("If this causes the program to crash, please re-run with more memory or groups.\n");
						}else{
							outstream.println(
									"It will be streamed to output unaltered.\n"
											+ "To avoid this behavior, increase memory or increase groups.\n"
											+ "Set the flag forceprocess=t to disable this check.\n");
//							Timer t=new Timer();
							ArrayList<Read> list=streamNext_inner(cris);
//							t.stop("Stream time: ");
							return list;
						}
					}
				}
			}
			
			return fetchNext_inner(cris);
		}
		
		private ArrayList<Read> streamNext_inner(ConcurrentReadInputStream cris){
			StreamToOutput sto=new StreamToOutput(cris, rosa, kc, (repair || namesort));
			errorState|=sto.process();
			return new ArrayList<Read>();
		}
		
		private ArrayList<Read> fetchNext_inner(ConcurrentReadInputStream cris){
			
//			Timer t=new Timer();
//			if(verbose){t.start("Making hash threads.");}
			final int subthreads=Tools.mid(1, Shared.threads()/2, 8);
			ArrayList<FetchSubThread> alht=new ArrayList<FetchSubThread>(subthreads);
			long readsThisGroup=0, memThisGroup=0;
			for(int i=0; i<subthreads; i++){alht.add(new FetchSubThread(i, cris, kc, unpair));}
			
//			if(verbose){outstream.println("Starting threads.");}
			for(FetchSubThread ht : alht){ht.start();}
			
//			if(verbose){outstream.println("Waiting for threads.");}
			/* Wait for threads to die */
			for(FetchSubThread ht : alht){
				
				/* Wait for a thread to die */
				while(ht.getState()!=Thread.State.TERMINATED){
					try {
						ht.join();
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}
				readsThisGroup+=ht.readsProcessedST;
				basesProcessedT+=ht.basesProcessedST;
				diskProcessedT+=ht.diskProcessedST;
				memThisGroup+=ht.memProcessedST;
			}
			readsProcessedT+=readsThisGroup;
			memProcessedT+=memThisGroup;

//			if(verbose){t.stop("Hash time: ");}
//			if(verbose){System.err.println("Closing input stream.");}
			errorState=ReadWrite.closeStream(cris)|errorState;
			
//			if(verbose){t.start("Combining thread output.");}
			assert(readsProcessedT<=Integer.MAX_VALUE && readsProcessedT>=0);
			ArrayList<Read> list=new ArrayList<Read>((int)readsThisGroup);
			for(int i=0; i<subthreads; i++){
				FetchSubThread ht=alht.set(i, null);
				list.addAll(ht.storage);
			}
			assert(list.size()==readsThisGroup || (cris.paired() && list.size()*2==readsThisGroup)) : list.size()+", "+readsThisGroup+", "+cris.paired();
//			if(verbose){t.stop("Combine time: ");}
			
//			if(verbose){t.start("Sorting.");}
			Shared.sort(list, kc);
			
//			if(verbose){t.stop("Sort time: ");}
			return list;
		}
		
		final SynchronousQueue<ArrayList<Read>> listQ;
		final AtomicInteger nextGroup;
		final KmerComparator kc;
		final ConcurrentReadOutputStream[] rosa;
		
		protected long readsProcessedT=0;
		protected long basesProcessedT=0;
		protected long diskProcessedT=0;
		protected long memProcessedT=0;
		protected boolean errorState=false;
		
		
		private class FetchSubThread extends Thread{
			
			FetchSubThread(int id_, ConcurrentReadInputStream cris_, KmerComparator kc_, boolean unpair_){
				id=id_;
				cris=cris_;
				kc=kc_;
				storage=new ArrayList<Read>();
				unpair=unpair_;
			}
			
			@Override
			public void run(){
				ListNum<Read> ln=cris.nextList();
				final boolean paired=cris.paired();
				ArrayList<Read> reads=(ln!=null ? ln.list : null);
				
				while(reads!=null && reads.size()>0){
					
					for(Read r : reads){
						if(!r.validated()){
							r.validate(true);
							if(r.mate!=null){r.mate.validate(true);}
						}
						readsProcessedST+=1+r.mateCount();
						basesProcessedST+=r.length()+r.mateLength();
						diskProcessedST+=r.countFastqBytes()+r.countMateFastqBytes();
						memProcessedST+=r.countBytes()+r.countMateBytes()+ReadKey.overhead;
						if(shrinkName){
							Clumpify.shrinkName(r);
							Clumpify.shrinkName(r.mate);
						}else if(shortName){
							Clumpify.shortName(r);
							Clumpify.shortName(r.mate);
						}
					}
					
					if(ecco){
						for(Read r : reads){
							Read r2=r.mate;
							assert(r.obj==null) : "TODO: Pivot should not have been generated yet, though it may be OK.";
							assert(r2!=null) : "ecco requires paired reads.";
							if(r2!=null){
								int x=BBMerge.findOverlapStrict(r, r2, true);
								if(x>=0){
									r.obj=null;
									r2.obj=null;
								}
							}
						}
					}
					
					ArrayList<Read> hashList=reads;
					if(paired && unpair){
						hashList=new ArrayList<Read>(reads.size()*2);
						for(Read r1 : reads){
							Read r2=r1.mate;
							assert(r2!=null);
							hashList.add(r1);
							hashList.add(r2);
							if(groups>1 || !repair || namesort){
								r1.mate=null;
								r2.mate=null;
							}
						}
					}
					
					kc.hash(hashList, table, minCount, true);
					storage.addAll(hashList);
					cris.returnList(ln.id, ln.list.isEmpty());
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
				if(ln!=null){
					cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
				}
				
				//Optimization for TimSort
				if(parallelSort){
					storage.sort(kc);
//					Shared.sort(storage, kc); //Already threaded; this is not needed.
				}else{
					Collections.sort(storage, kc);
				}
			}

			final int id;
			final ConcurrentReadInputStream cris;
			final KmerComparator kc;
			final ArrayList<Read> storage;
			final boolean unpair;

			protected long readsProcessedST=0;
			protected long basesProcessedST=0;
			protected long diskProcessedST=0;
			protected long memProcessedST=0;
		}
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private int k=31;
	private int minCount=0;
	
	private int groups=1;
	
	KCountArray table=null;
	
	/*--------------------------------------------------------------*/
	/*----------------          I/O Fields          ----------------*/
	/*--------------------------------------------------------------*/

	private String in1=null;
	private String in2=null;

	private String out1=null;
	
	private String extin=null;
	private String extout=null;
	
	/*--------------------------------------------------------------*/
	
	protected long readsProcessed=0;
	protected long basesProcessed=0;
	protected long diskProcessed=0;
	protected long memProcessed=0;

	protected long readsThisPass=0;
	protected long memThisPass=0;
	
	protected static long lastMemProcessed=0;
	
	protected long readsProcessedThisPass=0;
	protected long clumpsProcessedThisPass=0;
	protected long correctionsThisPass=0;
	
	protected long duplicatesThisPass=0;
	protected static long duplicatesTotal=0;
	
	protected long clumpsProcessedTotal=0;
	protected static long correctionsTotal=0;

	protected int passes=1;
	
	private long maxReads=-1;
	private boolean addName=false;
	private boolean shortName=false;
	private boolean shrinkName=false;
	private boolean rcomp=false;
	private boolean condense=false;
	private boolean correct=false;
	private boolean dedupe=false;
	private boolean splitInput=false;
	private boolean ecco=false;
	private boolean unpair=false;
	private boolean repair=false;
	private boolean namesort=false;
	
	private final boolean parallelSort=Shared.parallelSort;
	
	private final long expectedSizePerGroup;
	private final long expectedMemPerGroup;
	private final long totalMem;
	private final long fileMem;
	private final long fileSize;
	
	private final double memRatio;
	
	/*--------------------------------------------------------------*/

	private static final ArrayList<Read> POISON=new ArrayList<Read>();
	protected static int fetchThreads=2;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private final FileFormat ffin1[];
	private final FileFormat ffin2[];

	private final FileFormat ffout[];
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=true;
	public static boolean doHashAndSplit=true;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	
}

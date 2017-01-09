package sort;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.PriorityQueue;
import java.util.concurrent.atomic.AtomicLong;

import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.KillSwitch;
import stream.ConcurrentReadOutputStream;
import stream.CrisContainer;
import stream.Read;
import stream.SamLine;
import structures.ListNum;
import var2.ScafMap;
import dna.Parser;
import fileIO.ByteFile;
import fileIO.ReadWrite;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import fileIO.FileFormat;

/**
 * Sorts reads by name, potentially from multiple input files.
 * 
 * @author Brian Bushnell
 * @date September 21, 2016
 *
 */
public class SortByName {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		final boolean oldFI=FASTQ.FORCE_INTERLEAVED, oldTI=FASTQ.TEST_INTERLEAVED;
		SortByName as=new SortByName(args);
		as.process(t);
		FASTQ.FORCE_INTERLEAVED=oldFI;
		FASTQ.TEST_INTERLEAVED=oldTI;
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public SortByName(String[] args){
		
		//Process any config files
		args=Parser.parseConfig(args);
		
		//Detect whether the uses needs help
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		//Print the program name and arguments
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		boolean setInterleaved=false; //Whether interleaved was explicitly set.
		
		//Set some shared static variables regarding PIGZ
		Shared.READ_BUFFER_LENGTH=Tools.min(200, Shared.READ_BUFFER_LENGTH);
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		//Create a parser object
		Parser parser=new Parser();
		boolean ascending=true;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //Strip leading hyphens
			
			
			if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("verbose2")){
				assert(false) : "Verbose2 is disabled.";
//				verbose2=Tools.parseBoolean(b);
			}else if(a.equals("delete")){
				delete=Tools.parseBoolean(b);
			}else if(a.equals("ascending")){
				ascending=Tools.parseBoolean(b);
			}else if(a.equals("descending")){
				ascending=!Tools.parseBoolean(b);
			}else if(a.equals("length")){
				comparator=ReadLengthComparator.comparator;
			}else if(a.equals("name")){
				comparator=ReadComparatorName.comparator;
			}else if(a.equals("quality")){
				comparator=ReadQualityComparator.comparator;
			}else if(a.equals("position")){
				comparator=ReadComparatorPosition.comparator;
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		comparator.setAscending(ascending);
		SamLine.SET_FROM_OK=true;
		
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
		}
		
//		assert(false) : setInterleaved+", "+FASTQ.FORCE_INTERLEAVED+", "+FASTQ.TEST_INTERLEAVED;
//		assert(!FASTQ.FORCE_INTERLEAVED) : setInterleaved+", "+FASTQ.FORCE_INTERLEAVED+", "+FASTQ.TEST_INTERLEAVED;
//		assert(!FASTQ.TEST_INTERLEAVED) : setInterleaved+", "+FASTQ.FORCE_INTERLEAVED+", "+FASTQ.TEST_INTERLEAVED;
		
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		
		//Do output file # replacement
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		
		//Adjust interleaved detection based on the number of input files
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		//Ensure there is an input file
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		//Ensure out2 is not set without out1
		if(out1==null){
			if(out2!=null){
				printOptions();
				throw new RuntimeException("Error - cannot define out2 without defining out1.");
			}
		}
		
		//Adjust interleaved settings based on number of output files
		if(!setInterleaved){
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read to some input files.\n");
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		//Create output FileFormat objects
//		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
//		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		
		tempExt=".fq.gz";
		if(extout==null){
			if(ffout1!=null){
				tempExt=ffout1.fasta() ? ".fa.gz" : ffout1.samOrBam() ? ".sam" : ".fq.gz";
			}
		}else{
			tempExt=extout;
		}
		
		if(comparator==ReadComparatorPosition.comparator){
			if(ReadComparatorPosition.scafMap==null){
				ReadComparatorPosition.scafMap=ScafMap.loadHeader(in1);
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			useSharedHeader=(ffin1.samOrBam() && ffout1!=null && ffout1.samOrBam());
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, useSharedHeader, ffin1, ffin2, qfin1, qfin2);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
//		//Optionally create a read output stream
//		final ConcurrentReadOutputStream ros;
//		if(ffout1!=null){
//			final int buff=4;
//			
//			if(cris.paired() && out2==null && (in1!=null && !ffin1.samOrBam() && !ffout1.samOrBam())){
//				outstream.println("Writing interleaved.");
//			}
//			
//			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, false);
//			ros.start(); //Start the stream
//		}else{ros=null;}
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the read stream
		processInner(cris);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris);
		
		//Report timing and results
		{
			t.stop();
			
			//Calculate units per nanosecond
			double rpnano=readsProcessed/(double)(t.elapsed);
			double bpnano=basesProcessed/(double)(t.elapsed);
			
			//Add "k" and "m" for large numbers
			String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
			String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");
			
			//Format the strings so they have they are right-justified
			while(rpstring.length()<8){rpstring=" "+rpstring;}
			while(bpstring.length()<8){bpstring=" "+bpstring;}
			
			outstream.println("Time:                         \t"+t);
			outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
			outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		}
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Iterate through the reads */
	public void processInner(final ConcurrentReadInputStream cris){
		//Do anything necessary prior to processing
		
		ArrayList<Read> storage=new ArrayList<Read>();
		
		final long maxMem=Shared.memAvailable(1);
		final long memLimit=(long)(maxMem*.75);
		final long currentLimit=(long)(maxMem*.28);
		long currentMem=0;
		long dumped=0;
		AtomicLong outstandingMem=new AtomicLong();
		
		if(verbose){System.err.println("maxMem="+maxMem+", memLimit="+memLimit+", currentLimit="+currentLimit+", currentLimit="+currentLimit);}
		
		{
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			//As long as there is a nonempty read list...
			while(reads!=null && reads.size()>0){
				if(verbose2){outstream.println("Fetched "+reads.size()+" reads.");}
				
				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					//Track the initial length for statistics
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					//Increment counters
					readsProcessed+=1+r1.mateCount();
					basesProcessed+=initialLength1+initialLength2;
					
					currentMem+=r1.countBytes()+(r2==null ? 0 : r2.countBytes());
					
					storage.add(r1);
				}
				
				if(currentMem>=currentLimit){
					if(verbose){System.err.println("currentMem: "+currentMem+" >= "+currentLimit+", dumping. ");}
					outstandingMem.addAndGet(currentMem);
					sortAndDump(storage, currentMem, outstandingMem, null, false);
					storage=new ArrayList<Read>();
					dumped+=currentMem;
					currentMem=0;
					if(verbose){System.err.println("Waiting on memory; outstandingMem="+outstandingMem);}
					waitOnMemory(outstandingMem, memLimit);
					if(verbose){System.err.println("Done waiting; outstandingMem="+outstandingMem);}
				}
				
				//Notify the input stream that the list was used
				cris.returnList(ln.id, ln.list.isEmpty());
//				if(verbose){outstream.println("Returned a list.");}
				
				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}

		outstandingMem.addAndGet(currentMem);
		if(dumped==0){
			if(out1!=null){
				sortAndDump(storage, currentMem, outstandingMem, out1, useSharedHeader);
				storage=null;
				waitOnMemory(outstandingMem, 0);
			}
		}else{
			sortAndDump(storage, currentMem, outstandingMem, null, false);
			storage=null;
			waitOnMemory(outstandingMem, 0);
			mergeAndDump(outTemp, useSharedHeader);
		}
		
	}
	
	private void waitOnMemory(AtomicLong outstandingMem, long target){
		if(outstandingMem.get()>target){
			if(verbose){outstream.println("Syncing; outstandingMem="+outstandingMem);}
			while(outstandingMem.get()>target){
				try {
					synchronized(outstandingMem){
						outstandingMem.wait(2000);
					}
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private boolean mergeAndDump(ArrayList<String> fnames, boolean useHeader) {
		return mergeAndDump(fnames, ffout1, ffout2, delete, useHeader);
	}
	
	public static boolean mergeAndDump(ArrayList<String> fnames, FileFormat ffout1, FileFormat ffout2, boolean delete, boolean useHeader) {
		
		System.err.println("Merging "+fnames);
//		System.err.println("zl="+ReadWrite.ZIPLEVEL);
//		System.err.println("ztd="+ReadWrite.ZIP_THREAD_DIVISOR);
//		System.err.println("mzt="+ReadWrite.MAX_ZIP_THREADS);
//		System.err.println("pigz="+ReadWrite.USE_PIGZ);
		
		
		boolean errorState=false;
		final ConcurrentReadOutputStream ros;
		if(ffout1!=null){
			final int buff=1;
			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, null, null, buff, null, useHeader);
			ros.start(); //Start the stream
		}else{ros=null;}
		
		ArrayList<CrisContainer> cclist=new ArrayList<CrisContainer>(fnames.size());
		PriorityQueue<CrisContainer> q=new PriorityQueue<CrisContainer>(fnames.size());
		for(String fname : fnames){
			CrisContainer cc=new CrisContainer(fname, comparator);
			if(cc.peek()!=null){
				cclist.add(cc);
				q.add(cc);
			}
		}
		
		mergeAndDump(q, ros);
		if(verbose){
			System.err.println("Finished processing "+fnames);
		}
		
		for(CrisContainer cc : cclist){
			errorState|=cc.close();
		}
		if(delete){
			for(String fname : fnames){
				new File(fname).delete();
			}
		}
		if(ros!=null){errorState|=ReadWrite.closeStream(ros);}
		return errorState;
	}
	
	private static void mergeAndDump(final PriorityQueue<CrisContainer> q, final ConcurrentReadOutputStream ros) {
		
		for(CrisContainer cc : q){
			assert(!cc.cris().paired()) : FASTQ.TEST_INTERLEAVED+", "+FASTQ.FORCE_INTERLEAVED;
		}
		
		final int limit=100000;
		ArrayList<Read> buffer=new ArrayList<Read>(2*limit);
		while(!q.isEmpty()){
			if(verbose2){System.err.println("q size: "+q.size());}
			for(int i=0; !q.isEmpty() && (buffer.size()<limit || i<1); i++){//For loop to force it to run at least once
				CrisContainer cc=q.poll();
				if(verbose2){System.err.println("Polled a cc.");}
				ArrayList<Read> list=cc.fetch();
				if(verbose2){System.err.println("Grabbed "+list.size()+" reads.");}
				buffer.addAll(list);
				if(verbose2){System.err.println("Buffer size: "+buffer.size());}
				if(cc.hasMore()){
					if(verbose2){System.err.println("Returned cc.");}
					q.add(cc);
				}else{
					if(verbose2){System.err.println("Discarded cc.");}
				}
			}
			Shared.sort(buffer, comparator);
			if(verbose2){System.err.println("Sorted buffer.");}
//			if(repair){
//				assert(false) : "TODO";
//			}
			
			for(int i=1; i<buffer.size(); i++){
				Read a=buffer.get(i-1);
				Read b=buffer.get(i);
				assert(comparator.compare(a, b)<=0) : "\n"+a+"\n"+b;
				assert(a.mate==null) : "\n"+a+"\n"+b;
			}
			
			final Read peek=(q.isEmpty() ? null : q.peek().peek());
			final int maxIndex=(peek==null ? buffer.size() : indexOfLowestAbovePivot(buffer, peek));
			ArrayList<Read> list=new ArrayList<Read>(maxIndex);
			
			for(int index=0; index<maxIndex; index++){
				Read r=buffer.get(index);
				assert(peek==null || comparator.compare(peek, r)>0) : "\n"+peek+"\n"+r;
				list.add(r);
			}
			if(ros!=null){ros.add(list, 0);}
			
			ArrayList<Read> oldbuffer=buffer;
			buffer=new ArrayList<Read>(2*limit);
			for(int i=maxIndex, size=oldbuffer.size(); i<size; i++){
				buffer.add(oldbuffer.get(i));
			}
			if(verbose2){System.err.println("Buffer contains "+buffer.size()+" reads.");}
		}
		
		assert(buffer.isEmpty());
	}

//	
//	final Read peek=(q.isEmpty() ? null : q.peek().peek());
//	
//	final int maxIndex=(peek==null ? buffer.size() : indexOfLowestAbovePivot(buffer, peek));
//	
//	for(int index=0; index<maxIndex; index++){
//		Read r=buffer.get(index);
//		list.add(r);
//		if(list.size()>=200){
//			if(ros!=null){ros.add(list, 0);}
//			list=new ArrayList<Read>(200);
//		}
//	}
//	
//	if(!list.isEmpty()){
//		if(ros!=null){ros.add(list, 0);}
//		list=null;
//	}
////	if(verbose2){System.err.println("Stopped at index "+index+".");}
//	
//	ArrayList<Read> oldbuffer=buffer;
//	buffer=new ArrayList<Read>(2*limit);
//	for(int i=maxIndex, size=oldbuffer.size(); i<size; i++){
//		buffer.add(oldbuffer.get(i));
//	}
//	if(verbose2){System.err.println("Buffer contains "+buffer.size()+" reads.");}
//}
	
	private static final int indexOfLowestAbovePivot(final ArrayList<Read> list, final Read pivot){
		final int size=list.size();
		final int maxIndex=binarySearch(list, pivot);
		if(maxIndex<0){return 0;}
		if(maxIndex>=size){return size;}
		Read r=list.get(maxIndex);
		final int x=comparator.compare(pivot, r);
		assert(x<=0) : x+"\n"+pivot.id+"\n"+r.id;
		final int ret=(x==0 ? maxIndex+1 : maxIndex);
		assert(ret==size || comparator.compare(pivot, list.get(ret))<0);
		assert(ret==0 || comparator.compare(pivot, list.get(ret-1))>=0);
		return ret;
	}
	
	/** Returns highest index of an equal value, or the insert position */
	private static final int binarySearch(final ArrayList<Read> list, final Read pivot){
		int a=0, b=list.size()-1;
		while(b>a){
			final int mid=(a+b)/2;
			final Read r=list.get(mid);
			final int x=comparator.compare(pivot, r);
			if(x<0){b=mid;}
			else if(x>0){a=mid+1;}
			else{return mid;}
		}
		
//		for(int i=1; i<list.size(); i++){
//			assert(comparator.compare(list.get(i-1), list.get(i))<=0) : i+"\n"+list.get(i-1)+"\n"+ list.get(i)+"\n";
//		}
//		
//		assert(a>=list.size() || comparator.compare(pivot, list.get(a))<=0) : 
//			comparator.compare(pivot, list.get(a))+"\n"+
//					a+", "+b+", "+list.size()+"\n"+pivot.id+"\n"+list.get(a).id;
//			comparator.compare(pivot, list.get(a))+", "+a+", "+b+", "+list.size()+"\n"+pivot.id+"\n"+list.get(a-1).id+"\n"+list.get(a).id+"\n"+list.get(a+1).id;
		while(a>=0 && a<list.size()-1){
			final int x=comparator.compare(pivot, list.get(a+1));
			assert(x<=0) : x+", "+a+", "+b+", "+list.size()+"\n"+pivot.id+"\n"+list.get(a).id;
			if(x<0){break;}
			a++;
		}
		
		if(a>=0 && a==list.size()-1){//Corner case; insert after last element was not handled
			if(comparator.compare(pivot, list.get(a))>0){a++;}
		}
		
		assert(a>=list.size() || comparator.compare(pivot, list.get(a))<=0) : 
			comparator.compare(pivot, list.get(a))+"\n"+
					a+", "+b+", "+list.size()+"\n"+pivot.id+"\n"+list.get(a).id;
		
		return a;
	}
	
	private void sortAndDump(final ArrayList<Read> storage, final long currentMem, final AtomicLong outstandingMem, String fname, boolean useHeader) {
		String temp=fname;
		if(temp==null){
			synchronized(outTemp){
				if(verbose2){System.err.println("Synced to outTemp to make a temp file.");}
				File tmpfile=new File(".");//(Shared.tmpdir()==null ? null : new File(Shared.tmpdir()));
				if(tmpfile!=null && !tmpfile.exists()){tmpfile.mkdirs();}
				try {
					temp=File.createTempFile("sort_temp_", tempExt, tmpfile).toString();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
					KillSwitch.kill(e.getMessage());
				}
				if(verbose2){System.err.println("Temp file: "+temp);}
				outTemp.add(temp);
			}
		}

		if(verbose || true){System.err.println("Created a WriteThread for "+temp);}
		WriteThread wt=new WriteThread(storage, currentMem, outstandingMem, temp, useHeader);
		wt.start();
	}

	/** This is called if the program runs with no parameters */
	private void printOptions(){
		throw new RuntimeException("TODO");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static class WriteThread extends Thread{
		
		public WriteThread(final ArrayList<Read> storage_, final long currentMem_, final AtomicLong outstandingMem_, String fname_, boolean useHeader_){
			storage=storage_;
			currentMem=currentMem_;
			outstandingMem=outstandingMem_;
			fname=fname_;
			useHeader=useHeader_;
		}
		
		public void run(){
			
			if(verbose){System.err.println("Started a WriteThread.");}
			final FileFormat ffout=FileFormat.testOutput(fname, FileFormat.FASTQ, null, true, false, false, false);
			final ConcurrentReadOutputStream ros;
			if(ffout!=null){
				final int buff=4;
				ros=ConcurrentReadOutputStream.getStream(ffout, null, null, null, buff, null, useHeader);
				ros.start(); //Start the stream
			}else{ros=null;}

			if(verbose){System.err.println("Started a ros.");}
			Shared.sort(storage, comparator);
			
			if(verbose){System.err.println("Sorted reads.");}
			
			ArrayList<Read> buffer=new ArrayList<Read>(200);
			long id=0;
			for(int i=0, lim=storage.size(); i<lim; i++){
				Read r=storage.set(i, null);
				buffer.add(r);
				if(buffer.size()>=200){
					if(ros!=null){ros.add(buffer, id);}
					id++;
					buffer=new ArrayList<Read>(200);
				}
			}
			if(ros!=null && buffer.size()>0){ros.add(buffer, id);}
			errorState|=ReadWrite.closeStream(ros);
			if(verbose){System.err.println("Closed ros.");}
			
			synchronized(outstandingMem){
				outstandingMem.addAndGet(-currentMem);
				if(verbose){System.err.println("Decremented outstandingMem: "+outstandingMem);}
				outstandingMem.notify();
				if(verbose){System.err.println("Notified outstandingMem.");}
			}
		}
		
		final ArrayList<Read> storage;
		final long currentMem;
		final AtomicLong outstandingMem;
		final String fname;
		final boolean useHeader;
		boolean errorState=false;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;
	
	private String qfin1=null;
	private String qfin2=null;

	/** Primary output file path */
	private String out1=null;
	/** Secondary output file path */
	private String out2=null;
	
	private ArrayList<String> outTemp=new ArrayList<String>();
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	private String tempExt=null;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	private boolean delete=true;
	
	private boolean useSharedHeader=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	/** Secondary input file */
	private final FileFormat ffin2;
	
	/** Primary output file */
	private final FileFormat ffout1;
	/** Secondary output file */
	private final FileFormat ffout2;
	
	private static ReadComparator comparator=ReadComparatorName.comparator;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** Print verbose messages */
	public static final boolean verbose2=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=false;
	/** Append to existing output files */
	private boolean append=false;
	/** This flag has no effect on singlethreaded programs */
	private final boolean ordered=false;
	
}

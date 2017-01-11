package sketch;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

import dna.AminoAcid;
import dna.Parser;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import structures.LongList;
import tax.GiToNcbi;
import tax.ImgRecord;
import tax.TaxNode;
import tax.TaxTree;

/**
 * Creates MinHashSketches rapidly.
 * 
 * @author Brian Bushnell
 * @date July 6, 2016
 *
 */
public class SketchMaker extends SketchObject {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		SketchMaker sm=new SketchMaker(args);
		
		//Run the object
		sm.process(t);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public SketchMaker(String[] args){
		
		//Process any config files
		args=Parser.parseConfig(args);
		
		//Detect whether the uses needs help
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		//Print the program name and arguments
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		//Set some shared static variables regarding PIGZ
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		Shared.READ_BUFFER_LENGTH=1;
		
		//Create a parser object
		Parser parser=new Parser();
		
		int size_=10000, minsize_=100;
		int k_=31;
		int files_=1;
		float maxGenomeFraction_=0.01f;
		boolean rcomp_=true;
		int mode_=ONE_SKETCH;
		String imgDump=null;
		
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
			}else if(a.equals("files")){
				files_=Integer.parseInt(b);
			}else if(a.equals("size")){
				size_=(int)Tools.parseKMG(b);
			}else if(a.equals("minsize")){
				minsize_=(int)Tools.parseKMG(b);
			}else if(a.equals("prefilter")){
				prefilter=Tools.parseBoolean(b);
			}else if(a.equals("maxfraction")){
				Sketch.maxGenomeFraction=Float.parseFloat(b);
			}else if(a.equals("k")){
				k_=Integer.parseInt(b);
			}else if(a.equals("rcomp")){
				rcomp_=Tools.parseBoolean(b);
			}else if(a.equals("name")){
				outName=b;
			}else if(a.equals("taxid")){
				outTaxid=Integer.parseInt(b);
			}else if(Sketch.parseCoding(a, b)){
				//Do nothing
			}else if(a.equals("mode")){
				if(b.equalsIgnoreCase("single")){
					mode_=ONE_SKETCH;
				}else if(b.equalsIgnoreCase("taxa")){
					mode_=PER_TAXA;
				}else if(b.equalsIgnoreCase("sequence")){
					mode_=PER_SEQUENCE;
				}else if(b.equalsIgnoreCase("img")){
					mode_=IMG;
				}
			}else if(b==null && a.equals("single")){
				mode_=ONE_SKETCH;
			}else if(b==null && a.equals("taxa")){
				mode_=PER_TAXA;
			}else if(b==null && a.equals("sequence")){
				mode_=PER_SEQUENCE;
			}else if(a.equals("img")){
				mode_=IMG;
				if(b!=null){
					imgDump=b;
					if(imgDump.equals("auto")){imgDump=ImgRecord.DefaultDumpFile;}
				}
			}else if(a.equals("imgdump")){
				imgDump=b;
				if(imgDump.equals("auto")){imgDump=ImgRecord.DefaultDumpFile;}
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Tools.parseKMG(b);
				//Set a variable here
			}
			
			else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				giTableFile=b;
				if("auto".equalsIgnoreCase(b)){giTableFile=TaxTree.defaultTableFile();}
			}else if(a.equals("taxtree") || a.equals("tree")){
				taxTreeFile=b;
				if("auto".equalsIgnoreCase(b)){taxTreeFile=TaxTree.defaultTreeFile();}
			}
			
			else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			} 
			
			else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		size=size_;
		minsize=minsize_;
		k=k_;
		rcomp=rcomp_;
		mode=mode_;
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;
			in2=parser.in2;

			out1=parser.out1;
			
			extin=parser.extin;
		}
		files=(out1==null ? 0 : files_);
		
		assert(mode!=ONE_SKETCH || files<2) : "Multiple output files are not allowed in single-sketch mode.";
		
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
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
		
		ffout=makeFFArray(out1, files, overwrite, append);
		
//		//Ensure output files can be written
//		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
//			outstream.println((out1==null)+", "+out1);
//			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
//		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2, taxTreeFile, giTableFile, imgDump)){
			throw new RuntimeException("\nCan't read to some input files.\n");
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, taxTreeFile, giTableFile)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		tool=new SketchTool(size, k, 1, rcomp);
		
		if(taxTreeFile!=null){setTaxtree(taxTreeFile);}
		
		if(giTableFile!=null){
			loadGiToNcbi();
		}
		if(imgDump!=null){
			ImgRecord.imgMap=ImgRecord.toMap(imgDump);
		}
	}
	
	private static FileFormat[] makeFFArray(String fname0, int files, boolean overwrite, boolean append){
		if(files<1 || fname0==null){return null;}
		String[] fnames=new String[files];
		FileFormat[] ff=new FileFormat[files];
		for(int i=0; i<files; i++){
			String fname=fname0;
			if(files>1){
				assert(fname.indexOf('#')>-1) : "Output name requires # symbol for multiple files.";
				fname=fname.replaceFirst("#", ""+i);
			}
			fnames[i]=fname;
			ff[i]=FileFormat.testOutput(fname, FileFormat.TEXT, null, true, overwrite, append, false);
		}
		
		if(!Tools.testOutputFiles(overwrite, append, false, fnames)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+Arrays.toString(fnames)+"\n");
		}
		
		return ff;
	}
	
//	private static ByteStreamWriter[] makeTSWArray(FileFormat[] ff){
//		if(ff==null || ff.length==0){return null;}
//		ByteStreamWriter[] tsw=new ByteStreamWriter[ff.length];
//		for(int i=0; i<ff.length; i++){
//			tsw[i]=new ByteStreamWriter(ff[i]);
//			tsw[i].start();
//		}
//		return tsw;
//	}
	
	private static ByteStreamWriter[] makeTSWArray(FileFormat[] ff){
		if(ff==null || ff.length==0){return null;}
		ByteStreamWriter[] tsw=new ByteStreamWriter[ff.length];
		for(int i=0; i<ff.length; i++){
			tsw[i]=new ByteStreamWriter(ff[i]);
			tsw[i].start();
		}
		return tsw;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private LongList sizeList(){
		LongList sizes=new LongList();
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		
		//Grab the first ListNum of reads
		ListNum<Read> ln=cris.nextList();
		//Grab the actual read list from the ListNum
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		//As long as there is a nonempty read list...
		while(reads!=null && reads.size()>0){
			//					if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access

			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx);

				int taxID=-1;
				TaxNode tn=null;
				if(taxtree!=null){
					tn=taxtree.getNode(r1.id);
					while(tn!=null && tn.pid!=tn.id && tn.level<taxLevel){
						TaxNode temp=taxtree.getNode(tn.pid);
						if(temp==null || temp.level>=TaxTree.LIFE){break;}
						tn=temp;
					}
					if(tn!=null){taxID=tn.id;}
				}
				
				if(taxID>0){
					sizes.increment(taxID, r1.length()+r1.mateLength());
				}
			}

			//Notify the input stream that the list was used
			cris.returnList(ln.id, ln.list.isEmpty());
			//					if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access

			//Fetch a new list
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}

		//Notify the input stream that the final list was used
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
		
		errorState|=ReadWrite.closeStream(cris);
		outstream.println("Created prefilter.");
		
		return sizes;
	}

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		if(prefilter){
			sizeList=sizeList();
		}
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the reads in separate threads
		spawnThreads(cris);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStream(cris);
		
		//TODO: Write sketch
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
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
			
			outstream.println("Wrote "+sketchesWritten+" of "+sketchesMade+" sketches.");
			outstream.println("Time:                         \t"+t);
			outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
			outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		}
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Tools.mid(1, Shared.threads(), 8);
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		
		if(mode==PER_TAXA || mode==IMG){
			intMap=new ConcurrentHashMap<Long, SketchHeap>();
		}else if(mode==PER_SEQUENCE){
			tsw=makeTSWArray(ffout);
		}
		
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, i));
		}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		
		//Wait for completion of all threads
		boolean success=true;
		SketchHeap heap=null;
		for(ProcessThread pt : alpt){
			
			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			
			//Accumulate per-thread statistics
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			kmersProcessed+=pt.kmersProcessedT;
			
			sketchesMade+=pt.sketchesMadeT;
			sketchesWritten+=pt.sketchesWrittenT;
			
//			System.err.println(pt.sketchesMadeT+" "+pt.sketchesWrittenT);
			
//			System.err.println("pt.readsProcessedT="+pt.readsProcessedT);
			if(mode==ONE_SKETCH){
				SketchHeap temp=pt.heap;
				
				if(temp==null){
					//do nothing
				}else if(heap==null){heap=pt.heap;}
				else{heap.add(pt.heap);}
				
				if(heap!=null){
					if(outTaxid>=0){heap.taxID=outTaxid;}
					heap.setTaxName(outName);
					heap.genomeSize=kmersProcessed;
				}
			}
			success&=pt.success;
		}
		
		if(heap!=null && heap.name0()==null){
			heap.setName0(ReadWrite.stripToCore(in1));
		}
		
		if(ffout!=null){
			if(mode==PER_TAXA || mode==IMG){
				tsw=makeTSWArray(ffout);
				success&=writeMap(intMap);
			}else if(mode==ONE_SKETCH){
				tool.write(new Sketch(heap), ffout[0]);
				sketchesMade++;
				sketchesWritten++;
			}
		}
		
		if(tsw!=null){
			for(int i=0; i<tsw.length; i++){tsw[i].poisonAndWait();}
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		
		//Do anything necessary after processing
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          I/O Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	private boolean writeMap(ConcurrentHashMap<Long, SketchHeap> map){
		
		//Determine how many threads may be used
		final int threads=files;

		//Fill a list with WriteThreads
		ArrayList<WriteThread> alwt=new ArrayList<WriteThread>(threads);
		
		@SuppressWarnings("unchecked")
		ArrayDeque<SketchHeap>[] heaps=new ArrayDeque[threads];
		for(int i=0; i<threads; i++){
			heaps[i]=new ArrayDeque<SketchHeap>();
			WriteThread wt=new WriteThread(i, heaps[i]);
			alwt.add(wt);
		}
		
//		Set<Entry<Long, SketchHeap>> set=intMap.entrySet();
		for(Entry<Long, SketchHeap> entry : intMap.entrySet()){
//			set.remove(entry);  This will probably not work
			SketchHeap heap=entry.getValue();
			sketchesMade++;
			if(heap.size()>0 && heap.genomeSize>=minsize){
				heaps[(entry.hashCode()&Integer.MAX_VALUE)%threads].add(heap);
			}
		}
//		intMap.clear();
		intMap=null;
		
		//Start the threads
		for(WriteThread wt : alwt){wt.start();}
		
		//Wait for completion of all threads
		boolean success=true;
		for(WriteThread wt : alwt){

			//Wait until this thread has terminated
			while(wt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					wt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
//			sketchesMade+=wt.sketchesMadeT;
			sketchesWritten+=wt.sketchesWrittenT;
			success&=wt.success;
		}
		return success;
	}
	
	private class WriteThread extends Thread{
		
		WriteThread(int tnum_, ArrayDeque<SketchHeap> queue_){
			tnum=tnum_;
			queue=queue_;
		}
		
		public void run(){
			success=false;
			for(SketchHeap heap=queue.poll(); heap!=null; heap=queue.poll()){
				Sketch s=new Sketch(heap);
				tool.write(s, tsw[tnum]);
				sketchesWrittenT++;
			}
			success=true;
			queue=null;
		}
		
		ArrayDeque<SketchHeap> queue;
		final int tnum;
//		long sketchesMadeT=0;
		long sketchesWrittenT=0;
		boolean success=false;
	}
	
//	private void writeOutput(ConcurrentHashMap<Integer, SketchHeap> map){
//		ByteStreamWriter tsw=new ByteStreamWriter(ffout);
//		tsw.start();
//		KeySetView<Integer, SketchHeap> y=map.keySet();
//		for(Integer x : map.keySet()){
//			SketchHeap heap=map.get(x);
//			Sketch s=tool.toSketch(heap);
//			tool.write(s, tsw);
//		}
//		tsw.poisonAndWait();
//	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Tax Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	private void loadGiToNcbi(){
		Timer t=new Timer();
		outstream.println("Loading gi to taxa translation table.");
		GiToNcbi.initialize(giTableFile);
		t.stop();
		if(true){
			outstream.println("Time: \t"+t);
			Shared.printMemory();
			outstream.println();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This is called if the program runs with no parameters */
	private void printOptions(){
		throw new RuntimeException("TODO"); //TODO
	}
	
	private final long toValue(long kmer, long rkmer){
		long value=(rcomp ? Tools.max(kmer, rkmer) : kmer);
		return value;
	}
	
	public static final long parseImgId(String name){
		int idx=name.indexOf('.');
		if(idx<1){
			System.err.println("Could not parse number from "+name);
			return -1;
		}
		long id=-1;
		try {
			id=Long.parseLong(name.substring(0, idx));
		} catch (NumberFormatException e) {
			System.err.println("Could not parse number from "+name);
		}
		return id;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final int tid_){
			cris=cris_;
			tid=tid_;
			
			shift=2*k;
			shift2=shift-2;
			mask=~((-1L)<<shift);
			
			heap=new SketchHeap(size);
		}
		
		//Called by start()
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			processInner();
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Iterate through the reads */
		void processInner(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired()); //Disabled due to non-static access
			}
			
//			long cntr1=0, cntr2=0, cntr3=0, cntr4=0;

			//As long as there is a nonempty read list...
			while(reads!=null && reads.size()>0){
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access

				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;

					processReadPair(r1, r2);
				}

				//Notify the input stream that the list was used
				cris.returnList(ln.id, ln.list.isEmpty());
//				if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access

				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
			
//			System.out.println(cntr1+", "+cntr2+", "+cntr3+", "+cntr4);
		}
		
		void processReadPair(Read r1, Read r2){
			//Track the initial length for statistics
			final int initialLength1=r1.length();
			final int initialLength2=r1.mateLength();

			//Increment counters
			readsProcessedT+=1+r1.mateCount();
			basesProcessedT+=initialLength1+initialLength2;

			processRead(r1);
			if(r2!=null){processRead(r2);}
			
			{
				int taxID=-1;
				TaxNode tn=null;
				if(taxtree!=null && (mode==PER_TAXA || mode==PER_SEQUENCE || (mode==ONE_SKETCH && heap.taxID<0))){
					tn=taxtree.getNode(r1.id);
					while(tn!=null && tn.pid!=tn.id && tn.level<taxLevel){
						TaxNode temp=taxtree.getNode(tn.pid);
						if(temp==null || temp.level>=TaxTree.LIFE){break;}
						tn=temp;
					}
					if(tn!=null){taxID=tn.id;}
//					System.err.println("Node: "+r1.id+"\n->\n"+tn);
				}
				
				if(mode==PER_TAXA && tn!=null){
					if(taxID==0 || (tn.level>taxLevel && tn.level>=TaxTree.PHYLUM)){return;}
					TaxNode parent=taxtree.getNode(tn.pid);
					if(parent.pid==parent.id){return;}
					if(prefilter && sizeList.get(taxID)<minsize){return;}
				}
				
				if(heap.taxID<0){heap.taxID=taxID;}
				if(heap.name0()==null){
					heap.setName0(r1.id);
				}
				if(heap.taxName()==null && tn!=null){
					heap.setTaxName(tn.name);
				}
				assert(heap.taxID<0 || heap.taxName()!=null) : heap.taxID+", "+heap.taxName()+", "+heap.name()+", "+tn;
				
//				System.err.println("heap: "+heap);
				
				if(mode==PER_SEQUENCE){
					sketchesMadeT++;
					if(heap.size()>0 && heap.genomeSize>=minsize){
						Sketch sketch=new Sketch(heap);
						if(tsw!=null){
							final int choice=(sketch.hashCode()&Integer.MAX_VALUE)%files;
							synchronized(tsw[choice]){
								tool.write(sketch, tsw[choice]);
								sketchesWrittenT++;
							}
						}
					}
					heap.clear();
				}else if(mode==PER_TAXA){
					if(heap.size()>0){
						if(taxID<0 && heap.genomeSize<minsize){
							heap.clear();
						}else{
							final Long key=taxID>-1 ? taxID : new Long(nextUnknown.getAndIncrement());

							final SketchHeap old;
							synchronized(intMap){
								old=intMap.get(key);
								if(old==null){
									intMap.put(key, heap);
								}
							}
							if(old==null){
								heap=new SketchHeap(size);
							}else{
								synchronized(old){
									old.add(heap);
								}
								heap.clear();
							}
						}
					}
				}else if(mode==IMG){
//					cntr1++;
					if(heap.size()>0){
//						cntr2++;
						long imgID=parseImgId(r1.id);
						assert(imgID>0) : r1.id;
						if(imgID>=0){
//							cntr3++;
							Long key=new Long(imgID);
							
							final SketchHeap old;
							synchronized(intMap){
								old=intMap.get(key);
								if(old==null){
									intMap.put(key, heap);
								}
							}
							if(old==null){
								heap=new SketchHeap(size);
							}else{
								synchronized(old){
									old.add(heap);
								}
								heap.clear();
							}
						}
					}
				}
			}
		}
		
		void processRead(final Read r){
			final byte[] bases=r.bases;
			long kmer=0;
			long rkmer=0;
			int len=0;
			
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				long x=AminoAcid.baseToNumber[b];
				long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=(rkmer>>>2)|(x2<<shift2);
				if(x<0){len=0;}else{len++;}
				if(len>=k){
					kmersProcessedT++;
					heap.genomeSize++;
					long z=toValue(kmer, rkmer);
					long hash=SketchTool.hash(z);
					if(hash>0){heap.add(hash);}
				}
			}
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		/** Number of kmers processed by this thread */
		protected long kmersProcessedT=0;
		
		long sketchesMadeT=0;
		long sketchesWrittenT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Thread ID */
		final int tid;

		SketchHeap heap;

		final int shift;
		final int shift2;
		final long mask;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Hashing            ----------------*/
	/*--------------------------------------------------------------*/
	
//	private static long[][] makeCodes(int symbols, int modes){
//		Random randy=new Random(12345);
//		long[][] r=new long[symbols][modes];
//		for(int i=0; i<symbols; i++){
//			for(int j=0; j<modes; j++){
//				r[i][j]=randy.nextLong();
//			}
//		}
//		return r;
//	}
//	
//	private static final long hash(long kmer){
//		long code=kmer;
//		for(int i=0; i<8; i++){
//			int x=(int)(kmer&0xFF);
//			code^=codes[i][x];
//		}
//		return code;
//	}
//	
//	private static final long[][] codes=makeCodes(8, 256);
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;

	/** Primary output file path */
	private String out1=null;
	
	/** Override input file extension */
	private String extin=null;
	
	private String giTableFile=null;
	private String taxTreeFile=null;
	
	private String outName=null;
	private int outTaxid=-1;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;
	/** Number of bases processed */
	protected long kmersProcessed=0;
	/** Number of sketches started */
	protected long sketchesMade=0;
	/** Number of sketches written */
	protected long sketchesWritten=0;
	
	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	private LongList sizeList=null;

	private ConcurrentHashMap<Long, SketchHeap> intMap;
	private ByteStreamWriter tsw[];
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	/** Secondary input file */
	private final FileFormat ffin2;
	
	/** Primary output files */
	private final FileFormat ffout[];
	/** Number of output files */
	private final int files;
	
	private final boolean rcomp;
	private final int mode;
	
	private final SketchTool tool;
	private final int size;
	private final int minsize;
	private final int k;
	
	private int taxLevel=1;
	private boolean prefilter=false;
	
	private final AtomicInteger nextUnknown=new AtomicInteger(2000000000);
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=false;
	/** Append to existing output files */
	private boolean append=false;
	
}

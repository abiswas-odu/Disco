package hiseq;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ListNum;
import dna.Parser;
import fileIO.ByteFile;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import jgi.Dedupe;
import kmer.AbstractKmerTable;
import kmer.HashArray1D;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import fileIO.FileFormat;

/**
 * Analyzes a flow cell for low-quality areas.
 * Removes reads in the low-quality areas.
 * 
 * @author Brian Bushnell
 * @date August 31, 2016
 *
 */
public class AnalyzeFlowCell {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		AnalyzeFlowCell as=new AnalyzeFlowCell(args);
		as.process(t);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public AnalyzeFlowCell(String[] args){
		
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
//		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		//Create a parser object
		Parser parser=new Parser();
		parser.qtrimRight=trimRight;
		parser.trimq=trimq;
		parser.minReadLength=minlen;
		
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
			}else if(a.equals("divisor") || a.equals("size")){
				Tile.xSize=Tile.ySize=Integer.parseInt(b);
			}else if(a.equals("xdivisor") || a.equals("xsize")){
				Tile.xSize=Integer.parseInt(b);
			}else if(a.equals("ydivisor") || a.equals("ysize")){
				Tile.ySize=Integer.parseInt(b);
			}else if(a.equals("target")){
				targetAverageReads=Integer.parseInt(b);
			}else if(a.equals("dump")){
				dump=b;
			}else if(a.equals("indump") || a.equals("ind") || a.equals("dumpin")){
				dumpIn=b;
			}else if(a.equals("pound")){
				pound=Tools.parseBoolean(b);
			}else if(a.equals("loadkmers") || a.equals("usekmers")){
				loadKmers=Tools.parseBoolean(b);
			}else if(a.equals("lqo") || a.equals("lowqualityonly")){
				discardOnlyLowQuality=Tools.parseBoolean(b);
			}else if(a.equals("dl") || a.equals("discardlevel")){
				discardLevel=Integer.parseInt(b);
			}
			
			else if(a.equals("deviations") || a.equals("d")){
				qDeviations=uDeviations=eDeviations=Float.parseFloat(b);
			}else if(a.equals("qdeviations") || a.equals("qd") || a.equals("dq")){
				qDeviations=Float.parseFloat(b);
			}else if(a.equals("udeviations") || a.equals("ud") || a.equals("du")){
				uDeviations=Float.parseFloat(b);
			}else if(a.equals("edeviations") || a.equals("ed") || a.equals("de")){
				eDeviations=Float.parseFloat(b);
			}else if(a.equals("qfraction") || a.equals("qf")){
				qualFraction=Float.parseFloat(b);
			}else if(a.equals("efraction") || a.equals("uf")){
				uniqueFraction=Float.parseFloat(b);
			}else if(a.equals("efraction") || a.equals("ef")){
				errorFreeFraction=Float.parseFloat(b);
			}else if(a.equals("qabsolute") || a.equals("qa")){
				qualAbs=Float.parseFloat(b);
			}else if(a.equals("uabsolute") || a.equals("ua")){
				uniqueAbs=Float.parseFloat(b);
			}else if(a.equals("eabsolute") || a.equals("ea")){
				errorFreeAbs=Float.parseFloat(b);
			}
			
			else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
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
			

			trimq=parser.trimq;
			minlen=parser.minReadLength;
			trimLeft=parser.qtrimLeft;
			trimRight=parser.qtrimRight;
		}
		
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
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, dump)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read to some input files.\n");
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2, dump)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		if(dumpIn==null){
			flowcell=new FlowCell();
			if(loadKmers){loadKmers();}
			fillTiles();
			keySets=null;
		}else{
			flowcell=new FlowCell(dumpIn);
			
			if(flowcell.avgReads<targetAverageReads){
				flowcell=flowcell.widen(targetAverageReads);
			}
			
			avgQuality=flowcell.avgQuality;
			avgUnique=flowcell.avgUnique;
			avgErrorFree=flowcell.avgErrorFree;
			stdQuality=flowcell.stdQuality;
			stdUnique=flowcell.stdUnique;
			stdErrorFree=flowcell.stdErrorFree;

			long readsToDiscard=markTiles(flowcell.toList(), flowcell.avgReads);
		}
		processReads(t);
	}

	/** Create read streams and process all data */
	void loadKmers(){
		Timer t2=new Timer();
		outstream.print("Loading kmers:  \t");
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
//		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		//Process the read stream
		loadKmersInner(cris);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris);
		
		t2.stop();
		outstream.println(t2);
	}

	/** Create read streams and process all data */
	void fillTiles(){
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
//		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		//Process the read stream
		fillTilesInner(cris);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris);
	}

	/** Create read streams and process all data */
	void processReads(Timer t){
		
		if(ffout1!=null){
			Timer t2=new Timer();
			outstream.print("Filtering reads:\t");

			//Create a read input stream
			final ConcurrentReadInputStream cris;
			{
				cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
				cris.start(); //Start the stream
				if(verbose){outstream.println("Started cris");}
			}
			boolean paired=cris.paired();
			//		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}

			//Optionally create a read output stream
			final ConcurrentReadOutputStream ros;
			if(ffout1!=null){
				final int buff=4;

				ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, false);
				ros.start(); //Start the stream
			}else{ros=null;}

			//Process the read stream
			processInner(cris, ros);

			if(verbose){outstream.println("Finished; closing streams.");}

			//Close the read streams
			errorState|=ReadWrite.closeStreams(cris, ros);

			t2.stop();
			outstream.println(t2);
		}
		
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
			
			outstream.println();
			outstream.println("Time:                         \t"+t);
			outstream.println();
			outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
			outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
			
			if(ffout1!=null){

				rpstring=(readsDiscarded<100000 ? ""+readsDiscarded : readsDiscarded<100000000 ? (readsDiscarded/1000)+"k" : (readsDiscarded/1000000)+"m");
				bpstring=(basesDiscarded<100000 ? ""+basesDiscarded : basesDiscarded<100000000 ? (basesDiscarded/1000)+"k" : (basesDiscarded/1000000)+"m");
				while(rpstring.length()<8){rpstring=" "+rpstring;}
				while(bpstring.length()<8){bpstring=" "+bpstring;}
				outstream.println();
				outstream.println("Reads Discarded:    "+rpstring+" \t"+String.format("%.3f%%", readsDiscarded*100.0/readsProcessed));
				outstream.println("Bases Discarded:    "+bpstring+" \t"+String.format("%.3f%%", basesDiscarded*100.0/basesProcessed));
				outstream.println();

			}
		}
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Iterate through the reads */
	void processInner(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){
		
		//Do anything necessary prior to processing					
		readsProcessed=0;
		basesProcessed=0;
		
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
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
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
					
					boolean keep=processReadPair(r1, r2);
					if(!keep){
						reads.set(idx, null);
						readsDiscarded+=1+r1.mateCount();
						basesDiscarded+=initialLength1+initialLength2;
					}
				}
				
				//Output reads to the output stream
				if(ros!=null){ros.add(reads, ln.id);}
				
				//Notify the input stream that the list was used
				cris.returnList(ln.id, ln.list.isEmpty());
				if(verbose){outstream.println("Returned a list.");}
				
				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		//Do anything necessary after processing
		
	}
	
	/** Iterate through the reads */
	public void loadKmersInner(final ConcurrentReadInputStream cris){
		
		keySets=new AbstractKmerTable[WAYS];

		//Initialize tables
		for(int j=0; j<WAYS; j++){
			keySets[j]=new HashArray1D(512000, true);
		}
		//Do anything necessary prior to processing
		
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
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					//Track the initial length for statistics
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					if(initialLength1>=k && randy.nextBoolean()){
						final long kmer=toKmer(r1.bases, r1.quality, randy.nextInt(initialLength1-k2), k);
						if(kmer>=0){
							AbstractKmerTable table=keySets[(int)(kmer%WAYS)];
							table.increment(kmer);
						}
					}
					
					if(initialLength2>=k && randy.nextBoolean()){
						final long kmer=toKmer(r2.bases, r2.quality, randy.nextInt(initialLength2-k2), k);
						if(kmer>=0){
							AbstractKmerTable table=keySets[(int)(kmer%WAYS)];
							table.increment(kmer);
						}
					}
				}
				
				//Notify the input stream that the list was used
				cris.returnList(ln.id, ln.list.isEmpty());
				if(verbose){outstream.println("Returned a list.");}
				
				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		//Do anything necessary after processing
	}
	
	/** Iterate through the reads */
	public void fillTilesInner(final ConcurrentReadInputStream cris){


		Timer t2=new Timer();
		outstream.print("Filling tiles:  \t");
		
		
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
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
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
					
					final MicroTile mt=flowcell.getMicroTile(r1.id);
					
					if(loadKmers){
						if(initialLength1>=k){
							final long kmer=toKmer(r1.bases, r1.quality, randy.nextInt(initialLength1-k2), k);
							if(kmer>=0){
								AbstractKmerTable table=keySets[(int)(kmer%WAYS)];
								if(table.getValue(kmer)>0){mt.hits++;}
								else{mt.misses++;}
							}else{mt.misses++;}
						}
						
						if(initialLength1>=k){
							final long kmer=toKmer(r1.bases, r1.quality, randy.nextInt(initialLength1-k2), k);
							if(kmer>=0){
								AbstractKmerTable table=keySets[(int)(kmer%WAYS)];
								if(table.getValue(kmer)>0){mt.hits++;}
								else{mt.misses++;}
							}else{mt.misses++;}
						}
					}
					
					if(initialLength1>1){
						mt.qualityCount++;
						mt.qualitySum+=r1.avgQualityByProbabilityDouble(true, initialLength1);
						mt.errorFreeSum+=100*r1.probabilityErrorFree(true, initialLength1);
					}
					
					if(initialLength2>=1){
						mt.qualityCount++;
						mt.qualitySum+=r2.avgQualityByProbabilityDouble(true, initialLength2);
						mt.errorFreeSum+=100*r2.probabilityErrorFree(true, initialLength2);
					}
				}
				
				//Notify the input stream that the list was used
				cris.returnList(ln.id, ln.list.isEmpty());
				if(verbose){outstream.println("Returned a list.");}
				
				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		t2.stop();
		outstream.println(t2);
		
		ArrayList<MicroTile> mtList=flowcell.calcStats();
		if(flowcell.avgReads<targetAverageReads){
			flowcell=flowcell.widen(targetAverageReads);
			mtList=flowcell.toList();
		}
		avgQuality=flowcell.avgQuality;
		avgUnique=flowcell.avgUnique;
		avgErrorFree=flowcell.avgErrorFree;
		
		stdQuality=flowcell.stdQuality;
		stdUnique=flowcell.stdUnique;
		stdErrorFree=flowcell.stdErrorFree;
		
		long readsToDiscard=markTiles(mtList, flowcell.avgReads);
		
		if(dump!=null){
			TextStreamWriter tsw=new TextStreamWriter(dump, overwrite, append, false);
			tsw.start();

			tsw.println("#xSize\t"+Tile.xSize);
			tsw.println("#ySize\t"+Tile.ySize);
			tsw.println("#avgReads\t"+String.format("%.1f", flowcell.avgReads));
			
			tsw.println("#avgQuality\t"+String.format("%.3f", avgQuality));
			tsw.println("#avgUnique\t"+String.format("%.3f", avgUnique));
			tsw.println("#avgErrorFree\t"+String.format("%.3f", avgErrorFree));
			
			tsw.println("#stdQuality\t"+String.format("%.5f", stdQuality));
			tsw.println("#stdUnique\t"+String.format("%.5f", stdUnique));
			tsw.println("#stdErrorFree\t"+String.format("%.5f", stdErrorFree));
			
			tsw.println((pound ? "#" : "") + "lane\ttile\tx1\tx2\ty1\ty2\treads\tunique\tquality\tprobErrorFree\tdiscard");
			
			for(Lane lane : flowcell.lanes){
				if(lane!=null){
					for(Tile tile : lane.tiles){
						if(tile!=null){
							tsw.print(tile.toString());
						}
					}
				}
			}
			tsw.poisonAndWait();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	boolean processReadPair(final Read r1, final Read r2){
		boolean passes=processReadPair_inner(r1, r2);
		if(passes){return true;}
		if(trimq>0){
			TrimRead.trimFast(r1, trimLeft, trimRight, trimq, 0);
			if(r2!=null){TrimRead.trimFast(r2, trimLeft, trimRight, trimq, 0);}
			return r1.length()>=minlen && (r2==null || r2.length()>=minlen);
		}else{
			return false;
		}
	}
	
	/**
	 * Process a single read pair.
	 * @param r1 Read 1
	 * @param r2 Read 2 (may be null)
	 * @return True if the reads should be kept, false if they should be discarded.
	 */
	boolean processReadPair_inner(final Read r1, final Read r2){
		
		MicroTile mt=flowcell.getMicroTile(r1.id);
		if(mt==null){
			if(!warned){
				outstream.println("\nWarning - a read was found with no corresponding MicroTile:\n"+r1.id);
				warned=true;
			}
			return true;
		}
		if(mt.discard<discardLevel){return true;}
		if(!discardOnlyLowQuality){return false;}
		
		int len1=r1.length(), len2=r1.mateLength();
		if(len1>0){
			double qual=r1.avgQualityByProbabilityDouble(true, len1);
			double prob=100*r1.probabilityErrorFree(true, len1);
			if(qual<avgQuality-qDeviations*stdQuality){return false;}
			if(prob<avgErrorFree-eDeviations*stdErrorFree){return false;}
		}
		if(len2>0){
			double qual=r2.avgQualityByProbabilityDouble(true, len2);
			double prob=100*r2.probabilityErrorFree(true, len2);
			if(qual<avgQuality-qDeviations*stdQuality){return false;}
			if(prob<avgErrorFree-eDeviations*stdErrorFree){return false;}
		}
		return true;
	}
	
	/** This is called if the program runs with no parameters */
	private void printOptions(){
		throw new RuntimeException("TODO");
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Generate a kmer from specified start location
	 * @param bases
	 * @param start
	 * @param klen kmer length
	 * @return kmer
	 */
	private final long toKmer(final byte[] bases, byte[] quals, final int start, final int klen){
//		if(minprob>0 && quals!=null){
//			float prob=toProb(quals, start, klen);
//			if(prob<minprob){return -1;}
//		}
		final int stop=start+klen;
		assert(stop<=bases.length) : klen+", "+bases.length;
		long kmer=0;
		
		for(int i=start; i<stop; i++){
			final byte b=bases[i];
			final long x=Dedupe.baseToNumber[b];
			kmer=((kmer<<2)|x);
		}
		return kmer;
	}
	
	private final long markTiles(ArrayList<MicroTile> mtList, double avgReads){
		for(MicroTile mt : mtList){
			mt.discard=0;
		}
		long readsToDiscard=0;
		
		cDiscards=qDiscards=eDiscards=uDiscards=mtDiscards=mtRetained=0;
		
		for(MicroTile mt : mtList){
			double q=mt.averageQuality();
			double e=mt.percentErrorFree();
			double u=mt.uniquePercent();

			double dq=avgQuality-q;
			double de=avgErrorFree-e;
			double du=u-avgUnique;
			
			if(mt.qualityCount<10 && mt.qualityCount<0.02f*avgReads){
				mt.discard++;
				cDiscards++;
			}
			if(dq>qDeviations*stdQuality && dq>avgQuality*qualFraction && dq>qualAbs){
				mt.discard++;
				qDiscards++;
			}
			if(de>eDeviations*stdErrorFree && de>avgErrorFree*errorFreeFraction && de>errorFreeAbs){
				mt.discard++;
				eDiscards++;
			}
			if(avgUnique>2 && avgUnique<98){
				if(du>uDeviations*stdUnique && du>avgUnique*uniqueFraction && du>uniqueAbs){
					mt.discard++;
					uDiscards++;
				}
			}
			if(mt.discard>0){
				mtDiscards++;
				readsToDiscard+=mt.qualityCount;
			}
			else{mtRetained++;}
		}
		
		outstream.println();
		outstream.println("Flagged "+mtDiscards+" of "+(mtDiscards+mtRetained)+" micro-tiles, containing "+readsToDiscard+" reads:");
		outstream.println(uDiscards+" exceeded uniqueness thresholds.");
		outstream.println(qDiscards+" exceeded quality thresholds.");
		outstream.println(eDiscards+" exceeded error-free probability thresholds.");
		outstream.println(cDiscards+" had too few reads to calculate statistics.");
		outstream.println();
		
		return readsToDiscard;
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

	private String qfout1=null;
	private String qfout2=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	private boolean pound=true;
	private String dump=null;
	private String dumpIn=null;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Number of reads discarded */
	protected long readsDiscarded=0;
	/** Number of bases discarded */
	protected long basesDiscarded=0;

	protected long cDiscards=0;
	protected long uDiscards=0;
	protected long qDiscards=0;
	protected long eDiscards=0;
	protected long mtDiscards=0;
	protected long mtRetained=0;
	
	
	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	/** Hold kmers.  A kmer X such that X%WAYS=Y will be stored in keySets[Y] */
	private AbstractKmerTable[] keySets;
	
	private int targetAverageReads=800;

	private float minprob=0;
	private static final int WAYS=31;
	private static final int k=31, k2=30;
	
	private final java.util.concurrent.ThreadLocalRandom randy=java.util.concurrent.ThreadLocalRandom.current();
	private FlowCell flowcell;
	
//	private int[] avgQualityArray;
//	private int[] avgUniqueArray;
	
	private long minCountToUse=0;

	private float qDeviations=2f;
	private float uDeviations=1.5f;
	private float eDeviations=2f;
	
	private float qualFraction=0.01f;
	private float uniqueFraction=0.01f;
	private float errorFreeFraction=0.01f;
	
	private float qualAbs=1f;
	private float uniqueAbs=1f;
	private float errorFreeAbs=1f;
	
	private double avgQuality;
	private double avgUnique;
	private double avgErrorFree;
	
	private double stdQuality;
	private double stdUnique;
	private double stdErrorFree;

	private boolean loadKmers=true;
	private boolean discardOnlyLowQuality=true;
	private int discardLevel=1;
	
	private int minlen=30;
	private byte trimq=-1;
	private boolean trimLeft=false;
	private boolean trimRight=true;
	
	private boolean warned=false;
	
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
	/** This flag has no effect on singlethreaded programs */
	private final boolean ordered=false;
	
}

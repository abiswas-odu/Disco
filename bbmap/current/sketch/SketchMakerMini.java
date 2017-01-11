package sketch;

import java.io.PrintStream;
import java.util.ArrayList;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ListNum;
import tax.TaxNode;

/**
 * Creates MinHashSketches rapidly.
 * 
 * @author Brian Bushnell
 * @date July 6, 2016
 *
 */
public class SketchMakerMini extends SketchObject {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public SketchMakerMini(SketchTool tool_, boolean rcomp_){
		
		tool=tool_;
		rcomp=rcomp_;
		k=tool.k;
		
		shift=2*k;
		shift2=shift-2;
		mask=~((-1L)<<shift);
		
		heap=new SketchHeap(tool.size);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	public Sketch toSketch(String fname){
		heap.clear();
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			FileFormat ffin=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(-1, true, ffin, null, null, null);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		processInner(cris);
		
		errorState|=ReadWrite.closeStream(cris);
		Sketch sketch=new Sketch(heap);
		sketchesMade++;
		return sketch;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private final long toValue(long kmer, long rkmer){
		long value=(rcomp ? Tools.max(kmer, rkmer) : kmer);
		return value;
	}
	
	
	/** Iterate through the reads */
	void processInner(ConcurrentReadInputStream cris){

		//Grab the first ListNum of reads
		ListNum<Read> ln=cris.nextList();
		//Grab the actual read list from the ListNum
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

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
	}

	void processReadPair(Read r1, Read r2){
		//Track the initial length for statistics
		final int initialLength1=r1.length();
		final int initialLength2=r1.mateLength();

		//Increment counters
		readsProcessed+=1+r1.mateCount();
		basesProcessed+=initialLength1+initialLength2;

		processRead(r1);
		if(r2!=null){processRead(r2);}

		if(heap.name0()==null){
			heap.setName0(r1.id);
		}

		TaxNode tn=null;
		if(taxtree!=null && heap.taxID<0){
			try {
				tn=taxtree.getNode(r1.id);
			} catch (Throwable e) {}
			if(tn!=null){
				heap.taxID=tn.id;
				if(heap.taxName()==null){
					heap.setTaxName(tn.name);
				}
			}
		}
		assert(heap.taxID<0 || heap.taxName()!=null) : heap.taxID+", "+heap.taxName()+", "+heap.name()+", "+tn;
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
				kmersProcessed++;
				heap.genomeSize++;
				long z=toValue(kmer, rkmer);
				long hash=SketchTool.hash(z);
				if(hash>0){heap.add(hash);}
			}
		}
	}

	/** True only if this thread has completed successfully */
	boolean success=false;

	SketchHeap heap;

	final int shift;
	final int shift2;
	final long mask;
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;
	/** Number of bases processed */
	protected long kmersProcessed=0;
	/** Number of sketches started */
	protected long sketchesMade=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final boolean rcomp;
	
	private final SketchTool tool;
	private final int k;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	
}

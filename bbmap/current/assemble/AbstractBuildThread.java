package assemble;

import java.util.ArrayList;

import stream.ByteBuilder;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.LongList;

/**
 * @author Brian Bushnell
 * @date Jul 18, 2015
 *
 */
abstract class AbstractBuildThread extends Thread {
	
	public AbstractBuildThread(int id_, int mode_, ConcurrentReadInputStream[] crisa_){
		id=id_;
		crisa=crisa_;
		mode=mode_;
	}
	
	/** Input read stream */
	final ConcurrentReadInputStream[] crisa;
	
	final int mode;
	int minCountSeedCurrent;

	final int[] leftCounts=new int[4];
	final int[] rightCounts=new int[4];
	final ByteBuilder builderT=new ByteBuilder();
	
	final LongList insertSizes=new LongList();
	
	ArrayList<Read> contigs=new ArrayList<Read>();
	
	long readsInT=0;
	long basesInT=0;
	long lowqReadsT=0;
	long lowqBasesT=0;
	final int id;
	
}

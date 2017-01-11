package assemble;

import stream.ByteBuilder;
import ukmer.Kmer;

/**
 * Searches for dead ends.
 * @author Brian Bushnell
 * @date Jul 20, 2015
 *
 */
abstract class AbstractExploreThread extends ShaveObject implements Runnable {

	/**
	 * Constructor
	 */
	public AbstractExploreThread(int id_, int kbig_){
		id=id_;
		myKmer=new Kmer(kbig_);
		thread=new Thread(this);
	}

	@Override
	public final void run(){
		//TODO:

		//With processNextVictims enabled, the number of dead ends found drops from the first pass to the next, then stabilizes.
		//So, they are not being reset correctly.

		//Also, the number found - even with one thread - is nondeterministic if both are enabled.
		//Unstable whether or not processNextVictims is disabled.  But that's probably to be expected as the count is not exact.
		//What should be exact is the number of kmers removed for being dead ends.

		//The number is lower than expected.  65k for 600k reads with errors.  Most are bubbles, but 40% should be dead ends, or 240k.

		while(processNextTable(myKmer)){}
		while(processNextVictims(myKmer)){}
		
		for(int i=0; i<removeMatrixT.length; i++){
			for(int j=0; j<removeMatrixT.length; j++){
				if((i==FORWARD_BRANCH || i==BACKWARD_BRANCH) && (j==FORWARD_BRANCH || j==BACKWARD_BRANCH)){
					bubblesFoundT+=removeMatrixT[i][j];
				}
			}
		}
	}

	boolean processNextTable(){return processNextTable(myKmer);}
	abstract boolean processNextTable(final Kmer kmer);

	boolean processNextVictims(){return processNextVictims(myKmer);}
	abstract boolean processNextVictims(final Kmer kmer);

	/*--------------------------------------------------------------*/	

	public final void start(){thread.start();}
	public final Thread.State getState(){return thread.getState();}
	public final void join() throws InterruptedException{thread.join();}

	/*--------------------------------------------------------------*/	
	
	long kmersTestedT=0;
	long deadEndsFoundT=0;
	long bubblesFoundT=0;
	
	final int id;
	final Kmer myKmer;

	final int[] leftCounts=new int[4];
	final int[] rightCounts=new int[4];
	final ByteBuilder builderT=new ByteBuilder();

	long[][] countMatrixT=new long[8][8];
	long[][] removeMatrixT=new long[8][8];
	
	public final Thread thread;
	
}

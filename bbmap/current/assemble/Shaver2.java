package assemble;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicIntegerArray;

import kmer.AbstractKmerTableSet;
import shared.Tools;
import stream.ByteBuilder;
import ukmer.AbstractKmerTableU;
import ukmer.HashArrayU1D;
import ukmer.HashForestU;
import ukmer.Kmer;
import ukmer.KmerNodeU;
import ukmer.KmerTableSetU;
import dna.AminoAcid;

/**
 * Designed for removal of dead ends (aka hairs).
 * @author Brian Bushnell
 * @date Jun 26, 2015
 *
 */
public class Shaver2 extends Shaver {
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public Shaver2(KmerTableSetU tables_, int threads_){
		this(tables_, threads_, 1, 1, 1, 100, 100, true, true);
	}
	
	public Shaver2(KmerTableSetU tables_, int threads_, 
			int minCount_, int maxCount_, int minSeed_, int maxLengthToDiscard_, int maxDistanceToExplore_, 
			boolean removeHair_, boolean removeBubbles_){
		super(tables_, threads_, minCount_, maxCount_, minSeed_, maxLengthToDiscard_, maxDistanceToExplore_, removeHair_, removeBubbles_);
		tables=tables_;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	final AbstractExploreThread makeExploreThread(int id_){return new ExploreThread(id_);}
	@Override
	final AbstractShaveThread makeShaveThread(int id_){return new ShaveThread(id_);}
	
	
	/*--------------------------------------------------------------*/
	/*----------------       Dead-End Removal       ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public boolean exploreAndMark(Kmer kmer, ByteBuilder bb, int[] leftCounts, int[] rightCounts, int minCount, int maxCount, 
			int maxLengthToDiscard, int maxDistanceToExplore, boolean prune
			, long[][] countMatrixT, long[][] removeMatrixT
			){
		bb.clear();
		assert(kmer.len>=kmer.kbig);
		if(findOwner(kmer)>STATUS_UNEXPLORED){return false;}
		
		bb.appendKmer(kmer);
		final int a=explore(kmer, bb, leftCounts, rightCounts, minCount, maxCount, maxDistanceToExplore);
		
		bb.reverseComplementInPlace();
		kmer=tables.rightmostKmer(bb, kmer);
		final int b=explore(kmer, bb, leftCounts, rightCounts, minCount, maxCount, maxDistanceToExplore);

		final int min=Tools.min(a, b);
		final int max=Tools.max(a, b);
		
		countMatrixT[min][max]++;
		
		if(a==TOO_LONG || a==TOO_DEEP || a==LOOP || a==FORWARD_BRANCH){
			claim(bb, STATUS_EXPLORED, false, kmer);
			return false;
		}
		
		if(b==TOO_LONG || b==TOO_DEEP || b==LOOP || b==FORWARD_BRANCH){
			claim(bb, STATUS_EXPLORED, false, kmer);
			return false;
		}
		
		if(bb.length()-kbig>maxLengthToDiscard){
			claim(bb, STATUS_EXPLORED, false, kmer);
			return false;
		}
		
		if(removeHair && min==DEAD_END){
			if(max==DEAD_END || max==BACKWARD_BRANCH){
				removeMatrixT[min][max]++;
				boolean success=claim(bb, STATUS_REMOVE, false, kmer);
				if(verbose || verbose2){System.err.println("Claiming ("+a+","+b+") length "+bb.length()+": "+bb);}
				assert(success);
				return true;
			}
		}
		
		if(removeBubbles){
			if(a==BACKWARD_BRANCH && b==BACKWARD_BRANCH){
				removeMatrixT[min][max]++;
				boolean success=claim(bb, STATUS_REMOVE, false, kmer);
				if(verbose || verbose2){System.err.println("Claiming ("+a+","+b+") length "+bb.length()+": "+bb);}
				assert(success);
				return true;
			}
		}
		
		claim(bb, STATUS_EXPLORED, false, kmer);
		return false;
	}
	
	/** Explores a single unbranching path in the forward direction.
	 * Returns reason for ending in this direction:
	 *  DEAD_END, TOO_LONG, TOO_DEEP, FORWARD_BRANCH, BACKWARD_BRANCH */
	public int explore(Kmer kmer, ByteBuilder bb, int[] leftCounts, int[] rightCounts, int minCount, int maxCount, int maxLength0){
		if(verbose){outstream.println("Entering explore with bb.length()="+bb.length());}
		assert(bb.length()>=kmer.kbig && kmer.len>=kmer.kbig);
		if(bb.length()==0){bb.appendKmer(kmer);}
		
		final int initialLength=bb.length();
		final int maxLength=maxLength0+kbig;
		
		final long firstKey=kmer.xor();
		HashArrayU1D table=tables.getTable(kmer);
		int count=table.getValue(kmer);
		assert(count>=minCount && count<=maxCount);
		
		int nextRightMaxPos=fillRightCounts(kmer, rightCounts);
		int nextRightMax=rightCounts[nextRightMaxPos];
		if(nextRightMax<minCount){
			if(verbose){outstream.println("Returning DEAD_END: rightMax="+nextRightMax);}
			return DEAD_END;
		}
		
		while(bb.length()<=maxLength){
			
			final int rightMaxPos=nextRightMaxPos;
			final int rightMax=rightCounts[rightMaxPos];
			final int rightSecondPos=Tools.secondHighestPosition(rightCounts);
			final int rightSecond=rightCounts[rightSecondPos];
			
			if(verbose){
				outstream.println("kmer: "+toText(kmer));
				outstream.println("Right counts: "+count+", "+Arrays.toString(rightCounts));
				outstream.println("rightMaxPos="+rightMaxPos);
				outstream.println("rightMax="+rightMax);
				outstream.println("rightSecondPos="+rightSecondPos);
				outstream.println("rightSecond="+rightSecond);
			}
			
			final int prevCount=count;
			
			//Generate the new base
			final byte b=AminoAcid.numberToBase[rightMaxPos];
			final long x=rightMaxPos;
			long evicted=kmer.addRightNumeric(x);
			
			//Now consider the next kmer
			if(kmer.xor()==firstKey){
				if(verbose){outstream.println("Returning LOOP");}
				return LOOP;
			}
			table=tables.getTable(kmer);
			
			assert(table.getValue(kmer)==rightMax);
			count=rightMax;
			
			
			{//Fill right and look for dead end
				nextRightMaxPos=fillRightCounts(kmer, rightCounts);
				nextRightMax=rightCounts[nextRightMaxPos];
				if(nextRightMax<minCount){
					if(verbose){outstream.println("Returning DEAD_END: rightMax="+rightMax);}
					return DEAD_END;
				}
			}
			
			
			{//Look left
				final int leftMaxPos=fillLeftCounts(kmer, leftCounts);
				final int leftMax=leftCounts[leftMaxPos];
				final int leftSecondPos=Tools.secondHighestPosition(leftCounts);
				final int leftSecond=leftCounts[leftSecondPos];
				
//				assert(leftMax==1 || leftMax==0) : prevCount+" -> "+Arrays.toString(leftCounts)+", "+count+", "+Arrays.toString(rightCounts);
				
				if(verbose){
					outstream.println("Left counts: "+count+", "+Arrays.toString(leftCounts));
					outstream.println("leftMaxPos="+leftMaxPos);
					outstream.println("leftMax="+leftMax);
					outstream.println("leftSecondPos="+leftSecondPos);
					outstream.println("leftSecond="+leftSecond);
				}
				
				if(leftSecond>=minCount || leftMax>prevCount){//Backward branch
//					assert(leftSecond==1) : prevCount+" -> "+Arrays.toString(leftCounts)+", "+count+", "+Arrays.toString(rightCounts);
					if(leftMax>prevCount){
						if(verbose){outstream.println("Returning BACKWARD_BRANCH_LOWER: " +
								"count="+count+", prevCount="+prevCount+", leftMax="+leftMax+", leftSecond="+leftSecond);}
						return BACKWARD_BRANCH;
					}else{
						assert(leftMax==prevCount);
						if(leftMax>=2*leftSecond){//This constant is adjustable
//							assert(false) : prevCount+" -> "+Arrays.toString(leftCounts)+", "+count+", "+Arrays.toString(rightCounts);
							//keep going
						}else{
							if(verbose){outstream.println("Returning BACKWARD_BRANCH_SIMILAR: " +
								"count="+count+", prevCount="+prevCount+", leftMax="+leftMax+", leftSecond="+leftSecond);}
							return BACKWARD_BRANCH;
						}
					}
				}
				
			}
			
			//Look right
			if(rightSecond>=minCount){
				if(verbose){outstream.println("Returning FORWARD_BRANCH: rightSecond="+rightSecond);}
				return FORWARD_BRANCH;
			}
			
			if(count>maxCount){
				if(verbose){outstream.println("Returning TOO_DEEP: rightMax="+rightMax);}
				return TOO_DEEP;
			}
			
			bb.append(b);
			if(verbose){outstream.println("Added base "+(char)b);}
		}
		
		assert(bb.length()>maxLength);
		if(verbose){outstream.println("Returning TOO_LONG: length="+bb.length());}
		return TOO_LONG;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         ExploreThread        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Searches for dead ends. 
	 */
	class ExploreThread extends AbstractExploreThread{
		
		/**
		 * Constructor
		 */
		public ExploreThread(int id_){
			super(id_, kbig);
		}
		
		@Override
		boolean processNextTable(final Kmer kmer){
			final int tnum=nextTable.getAndAdd(1);
			if(tnum>=tables.ways){return false;}
			final HashArrayU1D table=tables.getTable(tnum);
			final int max=table.arrayLength();
			for(int cell=0; cell<max; cell++){
				int x=processCell(table, cell, kmer);
				deadEndsFoundT+=x;
			}
			return true;
		}
		
		@Override
		boolean processNextVictims(final Kmer kmer){
			final int tnum=nextVictims.getAndAdd(1);
			if(tnum>=tables.ways){return false;}
			final HashArrayU1D table=tables.getTable(tnum);
			final HashForestU forest=table.victims();
			final int max=forest.arrayLength();
			for(int cell=0; cell<max; cell++){
				KmerNodeU kn=forest.getNode(cell);
				int x=traverseKmerNodeU(kn, kmer);
				deadEndsFoundT+=x;
			}
			return true;
		}
		
		private int processCell(HashArrayU1D table, int cell, Kmer kmer0){
			int count=table.readCellValue(cell);
			if(count<minSeed || count>maxCount){return 0;}
			int owner=table.getCellOwner(cell);
			if(owner>STATUS_UNEXPLORED){return 0;}
			Kmer kmer=table.fillKmer(cell, kmer0);
			if(kmer==null){return 0;}
			if(verbose){outstream.println("id="+id+" processing cell "+cell+"; \tkmer="+kmer);}
			
			return processKmer(kmer);
		}
		
		private int traverseKmerNodeU(KmerNodeU kn, Kmer kmer){
			int sum=0;
			if(kn!=null){
				sum+=processKmerNodeU(kn, kmer);
				if(kn.left()!=null){
					sum+=traverseKmerNodeU(kn.left(), kmer);
				}
				if(kn.right()!=null){
					sum+=traverseKmerNodeU(kn.right(), kmer);
				}
			}
			return sum;
		}
		
		private int processKmerNodeU(KmerNodeU kn, Kmer kmer){
			kmer.setFrom(kn.pivot());
			final int count=kn.getValue(kmer);
			if(count<minSeed || count>maxCount){return 0;}
			int owner=kn.getOwner(kmer);
			if(owner>STATUS_UNEXPLORED){return 0;}
			
			return processKmer(kmer);
		}
		
		private int processKmer(Kmer kmer){
			kmersTestedT++;
			boolean b=exploreAndMark(kmer, builderT, leftCounts, rightCounts, minCount, maxCount, maxLengthToDiscard, maxDistanceToExplore, true
					, countMatrixT, removeMatrixT
				);
			return b ? 1 : 0;
		}
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          ShaveThread         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Removes dead-end kmers. 
	 */
	private class ShaveThread extends AbstractShaveThread{

		/**
		 * Constructor
		 */
		public ShaveThread(int id_){
			super(id_);
		}
		
		@Override
		boolean processNextTable(){
			final int tnum=nextTable.getAndAdd(1);
			if(tnum>=tables.ways){return false;}
//			long x=0;
			final HashArrayU1D table=tables.getTable(tnum);
			final AtomicIntegerArray owners=table.owners();
			final int[] values=table.values();
			final int max=table.arrayLength();
			for(int cell=0; cell<max; cell++){
				if(owners.get(cell)==STATUS_REMOVE){
//					x++;
					values[cell]=0;
				}
			}
			for(KmerNodeU kn : table.victims().array()){
				if(kn!=null){traverseKmerNodeU(kn);}
			}
			
			table.clearOwnership();
			kmersRemovedT+=table.regenerate(0);
//			outstream.println(x);
			return true;
		}
		
		private void traverseKmerNodeU(KmerNodeU kn){
			if(kn==null){return;}
			if(kn.owner()==STATUS_REMOVE){kn.set(0);}
			traverseKmerNodeU(kn.left());
			traverseKmerNodeU(kn.right());
		}
		
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Recall Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	int getCount(Kmer kmer){return tables.getCount(kmer);}
	boolean claim(Kmer kmer, int id){return tables.claim(kmer, id);}
	boolean doubleClaim(ByteBuilder bb, int id/*, long rid*/, Kmer kmer){return tables.doubleClaim(bb, id, kmer/*, rid*/);}
	boolean claim(ByteBuilder bb, int id, /*long rid, */boolean earlyExit, Kmer kmer){return tables.claim(bb, id/*, rid*/, earlyExit, kmer);}
	boolean claim(byte[] array, int len, int id, /*long rid, */boolean earlyExit, Kmer kmer){return tables.claim(array, len, id/*, rid*/, earlyExit, kmer);}
	int findOwner(Kmer kmer){return tables.findOwner(kmer);}
	int findOwner(ByteBuilder bb, int id, Kmer kmer){return tables.findOwner(bb, id, kmer);}
	int findOwner(byte[] array, int len, int id, Kmer kmer){return tables.findOwner(array, len, id, kmer);}
	void release(ByteBuilder bb, int id, Kmer kmer){tables.release(bb, id, kmer);}
	void release(byte[] array, int len, int id, Kmer kmer){tables.release(array, len, id, kmer);}
	int fillRightCounts(Kmer kmer, int[] counts){return tables.fillRightCounts(kmer, counts);}
	int fillLeftCounts(Kmer kmer, int[] counts){return tables.fillLeftCounts(kmer, counts);}
	static StringBuilder toText(Kmer kmer){return AbstractKmerTableU.toText(kmer);}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	AbstractKmerTableSet tables(){return tables;}
	
	private final KmerTableSetU tables;
	
}

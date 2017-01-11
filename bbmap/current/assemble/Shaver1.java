package assemble;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicIntegerArray;

import kmer.AbstractKmerTable;
import kmer.AbstractKmerTableSet;
import kmer.HashArray1D;
import kmer.HashForest;
import kmer.KmerNode;
import kmer.KmerTableSet;
import shared.Tools;
import stream.ByteBuilder;
import ukmer.Kmer;
import dna.AminoAcid;

/**
 * Designed for removal of dead ends (aka hairs).
 * @author Brian Bushnell
 * @date Jun 26, 2015
 *
 */
public class Shaver1 extends Shaver {
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/

	
	public Shaver1(KmerTableSet tables_, int threads_){
		this(tables_, threads_, 1, 1, 1, 100, 100, true, true);
	}
	
	public Shaver1(KmerTableSet tables_, int threads_, 
			int minCount_, int maxCount_, int minSeed_, int maxLengthToDiscard_, int maxDistanceToExplore_, 
			boolean removeHair_, boolean removeBubbles_){
		super(tables_, threads_, minCount_, maxCount_, minSeed_, maxLengthToDiscard_, maxDistanceToExplore_, removeHair_, removeBubbles_);
		tables=tables_;
		k=tables.k;
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
	
	
	public boolean exploreAndMark(long kmer, ByteBuilder bb, int[] leftCounts, int[] rightCounts, int minCount, int maxCount, 
			int maxLengthToDiscard, int maxDistanceToExplore, boolean prune
			, long[][] countMatrixT, long[][] removeMatrixT
			){
		bb.clear();
		if(findOwner(kmer)>STATUS_UNEXPLORED){return false;}
		
		bb.appendKmer(kmer, k);
		final int a=explore(kmer, bb, leftCounts, rightCounts, minCount, maxCount, maxDistanceToExplore);
		
		bb.reverseComplementInPlace();
		kmer=tables.rightmostKmer(bb);
		final int b=explore(kmer, bb, leftCounts, rightCounts, minCount, maxCount, maxDistanceToExplore);

		final int min=Tools.min(a, b);
		final int max=Tools.max(a, b);
		
		countMatrixT[min][max]++;
		
		if(a==TOO_LONG || a==TOO_DEEP || a==LOOP || a==FORWARD_BRANCH){
			claim(bb, STATUS_EXPLORED, false);
			return false;
		}
		
		if(b==TOO_LONG || b==TOO_DEEP || b==LOOP || b==FORWARD_BRANCH){
			claim(bb, STATUS_EXPLORED, false);
			return false;
		}
		
		if(bb.length()-k>maxLengthToDiscard){
			claim(bb, STATUS_EXPLORED, false);
			return false;
		}
		
		if(removeHair && min==DEAD_END){
			if(max==DEAD_END || max==BACKWARD_BRANCH){
				removeMatrixT[min][max]++;
				boolean success=claim(bb, STATUS_REMOVE, false);
				if(verbose || verbose2){System.err.println("Claiming ("+a+","+b+") length "+bb.length()+": "+bb);}
				assert(success);
				return true;
			}
		}
		
		if(removeBubbles){
			if(a==BACKWARD_BRANCH && b==BACKWARD_BRANCH){
				removeMatrixT[min][max]++;
				boolean success=claim(bb, STATUS_REMOVE, false);
				if(verbose || verbose2){System.err.println("Claiming ("+a+","+b+") length "+bb.length()+": "+bb);}
				assert(success);
				return true;
			}
		}
		
		claim(bb, STATUS_EXPLORED, false);
		return false;
	}
	
	/** Explores a single unbranching path in the forward direction.
	 * Returns reason for ending in this direction:
	 *  DEAD_END, TOO_LONG, TOO_DEEP, FORWARD_BRANCH, BACKWARD_BRANCH */
	public int explore(long kmer, ByteBuilder bb, int[] leftCounts, int[] rightCounts, int minCount, int maxCount, int maxLength0){
		if(verbose){outstream.println("Entering explore with bb.length()="+bb.length());}
		assert(bb.length()==0 || tables.rightmostKmer(bb)==kmer);
		if(bb.length()==0){bb.appendKmer(kmer, k);}
		
		final int initialLength=bb.length();
		final int maxLength=maxLength0+k;
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=~((-1L)<<shift);
		
		long rkmer=AminoAcid.reverseComplementBinaryFast(kmer, k);
		
		long key=toValue(kmer, rkmer);
		final long firstKey=key;
		HashArray1D table=tables.getTableForKey(key);
		int count=table.getValue(key);
		assert(count>=minCount && count<=maxCount);
		
		int nextRightMaxPos=fillRightCounts(kmer, rkmer, rightCounts, mask, shift2);
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
				outstream.println("kmer: "+toText(kmer, k)+", "+toText(rkmer, k));
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
			final long x2=AminoAcid.numberToComplement[rightMaxPos];
			kmer=((kmer<<2)|(long)x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			
			//Now consider the next kmer
			key=toValue(kmer, rkmer);
			if(key==firstKey){
				if(verbose){outstream.println("Returning LOOP");}
				return LOOP;
			}
			table=tables.getTableForKey(key);
			
			assert(table.getValue(key)==rightMax);
			count=rightMax;
			
			{//Fill right and look for dead end
				nextRightMaxPos=fillRightCounts(kmer, rkmer, rightCounts, mask, shift2);
				nextRightMax=rightCounts[nextRightMaxPos];
				if(nextRightMax<minCount){
					if(verbose){outstream.println("Returning DEAD_END: rightMax="+rightMax);}
					return DEAD_END;
				}
			}
			
			
			{//Look left
				final int leftMaxPos=fillLeftCounts(kmer, rkmer, leftCounts, mask, shift2);
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
	
	/** Explores a single unbranching path in the forward direction.
	 * Returns reason for ending in this direction:
	 *  DEAD_END, TOO_LONG, TOO_DEEP, FORWARD_BRANCH, BACKWARD_BRANCH */
	public int explore2(long kmer, ByteBuilder bb, int[] leftCounts, int[] rightCounts, int minCount, int maxCount, int maxLength0){
		if(verbose){outstream.println("Entering explore with bb.length()="+bb.length());}
		assert(bb.length()==0 || tables.rightmostKmer(bb)==kmer);
		if(bb.length()==0){bb.appendKmer(kmer, k);}
		
		final int maxLength=maxLength0+k;
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=~((-1L)<<shift);
		
		long rkmer=AminoAcid.reverseComplementBinaryFast(kmer, k);
		
		long key=toValue(kmer, rkmer);
		final long firstKey=key;
		HashArray1D table=tables.getTableForKey(key);
		int count=table.getValue(key);
		assert(count>=minCount && count<=maxCount);

		int rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts, mask, shift2);
		int rightMax=rightCounts[rightMaxPos];
		int rightSecondPos=Tools.secondHighestPosition(rightCounts);
		int rightSecond=rightCounts[rightSecondPos];

		if(verbose){
			outstream.println("kmer: "+toText(kmer, k)+", "+toText(rkmer, k));
			outstream.println("Right counts: "+count+", "+Arrays.toString(rightCounts));
			outstream.println("rightMaxPos="+rightMaxPos);
			outstream.println("rightMax="+rightMax);
			outstream.println("rightSecondPos="+rightSecondPos);
			outstream.println("rightSecond="+rightSecond);
		}

		if(rightMax<minCount){
			if(verbose){outstream.println("Returning DEAD_END: rightMax="+rightMax);}
			return DEAD_END;
		}else if(rightSecond>=minCount){
			if(verbose){outstream.println("Returning FORWARD_BRANCH: rightSecond="+rightSecond);}
			return FORWARD_BRANCH;
		}else if(rightMax>maxCount){
			if(verbose){outstream.println("Returning TOO_DEEP: rightMax="+rightMax);}
			return TOO_DEEP;
		}
		
		while(bb.length()<=maxLength){
			final int prevCount=count;
			
			//Generate the new base
			final byte b=AminoAcid.numberToBase[rightMaxPos];
			final long x=rightMaxPos;
			final long x2=AminoAcid.numberToComplement[rightMaxPos];
			kmer=((kmer<<2)|(long)x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			
			//Now consider the next kmer
			key=toValue(kmer, rkmer);
			if(key==firstKey){
				if(verbose){outstream.println("Returning LOOP");}
				return LOOP;
			}
			table=tables.getTableForKey(key);
			
			assert(table.getValue(key)==rightMax);
			count=rightMax;
			
			{//Look right
				rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts, mask, shift2);
				rightMax=rightCounts[rightMaxPos];
				rightSecondPos=Tools.secondHighestPosition(rightCounts);
				rightSecond=rightCounts[rightSecondPos];

				if(verbose){
					outstream.println("kmer: "+toText(kmer, k)+", "+toText(rkmer, k));
					outstream.println("Right counts: "+count+", "+Arrays.toString(rightCounts));
					outstream.println("rightMaxPos="+rightMaxPos);
					outstream.println("rightMax="+rightMax);
					outstream.println("rightSecondPos="+rightSecondPos);
					outstream.println("rightSecond="+rightSecond);
				}

				if(rightMax<minCount){
					if(verbose){outstream.println("Returning DEAD_END: rightMax="+rightMax);}
					return DEAD_END;
				}else if(rightSecond>=minCount){
					if(verbose){outstream.println("Returning FORWARD_BRANCH: rightSecond="+rightSecond);}
					return FORWARD_BRANCH;
				}
			}
			
			{//Look left
				int leftMaxPos=fillLeftCounts(kmer, rkmer, leftCounts, mask, shift2);
				int leftMax=leftCounts[leftMaxPos];
				int leftSecondPos=Tools.secondHighestPosition(leftCounts);
				int leftSecond=leftCounts[leftSecondPos];
				
//				assert(leftMax==1 || leftMax==0) : prevCount+" -> "+Arrays.toString(leftCounts)+", "+count+", "+Arrays.toString(rightCounts);
				
				if(verbose){
					outstream.println("Left counts: "+count+", "+Arrays.toString(leftCounts));
					outstream.println("leftMaxPos="+leftMaxPos);
					outstream.println("leftMax="+leftMax);
					outstream.println("leftSecondPos="+leftSecondPos);
					outstream.println("leftSecond="+leftSecond);
				}
				
				if(leftSecond>=minCount){//Backward branch
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
			
			if(count>maxCount){
				if(verbose){outstream.println("Returning TOO_DEEP: rightMax="+rightMax);}
				return TOO_DEEP;
			}
			
			bb.append(b);
			if(verbose){outstream.println("Added base "+(char)b);}
		}
		
		assert(bb.length()>maxLength);
		if(verbose){outstream.println("Returning TOO_LONG: length="+bb.length()+", rightMax="+rightMax);}
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
			super(id_, k);
		}
		
		@Override
		boolean processNextTable(Kmer kmer){
			final int tnum=nextTable.getAndAdd(1);
			if(tnum>=tables.ways){return false;}
			final HashArray1D table=tables.getTable(tnum);
			final int max=table.arrayLength();
			for(int cell=0; cell<max; cell++){
				int x=processCell(table, cell);
				deadEndsFoundT+=x;
			}
			return true;
		}
		
		@Override
		boolean processNextVictims(Kmer kmer){
			final int tnum=nextVictims.getAndAdd(1);
			if(tnum>=tables.ways){return false;}
			final HashArray1D table=tables.getTable(tnum);
			final HashForest forest=table.victims();
			final int max=forest.arrayLength();
			for(int cell=0; cell<max; cell++){
				KmerNode kn=forest.getNode(cell);
				int x=traverseKmerNode(kn);
				deadEndsFoundT+=x;
			}
			return true;
		}
		
		private int processCell(HashArray1D table, int cell){
			int count=table.readCellValue(cell);
			if(count<minSeed || count>maxCount){return 0;}
			int owner=table.getCellOwner(cell);
			if(owner>STATUS_UNEXPLORED){return 0;}
			long key=table.getKmer(cell);
			if(verbose){outstream.println("id="+id+" processing cell "+cell+"; \tkmer="+key+"\t"+AminoAcid.kmerToString(key, k));}
			
			return processKmer(key);
		}
		
		private int traverseKmerNode(KmerNode kn){
			int sum=0;
			if(kn!=null){
				sum+=processKmerNode(kn);
				if(kn.left()!=null){
					sum+=traverseKmerNode(kn.left());
				}
				if(kn.right()!=null){
					sum+=traverseKmerNode(kn.right());
				}
			}
			return sum;
		}
		
		private int processKmerNode(KmerNode kn){
			final long key=kn.pivot();
			final int count=kn.getValue(key);
			if(count<minSeed || count>maxCount){return 0;}
			int owner=kn.getOwner(key);
			if(owner>STATUS_UNEXPLORED){return 0;}
			
			return processKmer(key);
		}
		
		private int processKmer(long key){
			kmersTestedT++;
			boolean b=exploreAndMark(key, builderT, leftCounts, rightCounts, minCount, maxCount, maxLengthToDiscard, maxDistanceToExplore, true
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
	class ShaveThread extends AbstractShaveThread{

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
			final HashArray1D table=tables.getTable(tnum);
			final AtomicIntegerArray owners=table.owners();
			final int[] values=table.values();
			final int max=table.arrayLength();
			for(int cell=0; cell<max; cell++){
				if(owners.get(cell)==STATUS_REMOVE){
//					x++;
					values[cell]=0;
				}
			}
			for(KmerNode kn : table.victims().array()){
				if(kn!=null){traverseKmerNode(kn);}
			}
			
			table.clearOwnership();
			kmersRemovedT+=table.regenerate(0);
//			outstream.println(x);
			return true;
		}
		
		private void traverseKmerNode(KmerNode kn){
			if(kn==null){return;}
			if(kn.owner()==STATUS_REMOVE){kn.set(0);}
			traverseKmerNode(kn.left());
			traverseKmerNode(kn.right());
		}
		
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Recall Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	private final long toValue(long kmer, long rkmer){return tables.toValue(kmer, rkmer);}
	int getCount(long kmer, long rkmer){return tables.getCount(kmer, rkmer);}
	boolean claim(long kmer, int id){return claim(kmer, AminoAcid.reverseComplementBinaryFast(kmer, k), id);}
	boolean claim(long kmer, long rkmer, int id){return tables.claim(kmer, rkmer, id);}
	boolean doubleClaim(ByteBuilder bb, int id/*, long rid*/){return tables.doubleClaim(bb, id/*, rid*/);}
	boolean claim(ByteBuilder bb, int id, /*long rid, */boolean earlyExit){return tables.claim(bb, id/*, rid*/, earlyExit);}
	boolean claim(byte[] array, int len, int id, /*long rid, */boolean earlyExit){return tables.claim(array, len, id/*, rid*/, earlyExit);}
	int findOwner(long kmer){return tables.findOwner(kmer);}
	int findOwner(ByteBuilder bb, int id){return tables.findOwner(bb, id);}
	int findOwner(byte[] array, int len, int id){return tables.findOwner(array, len, id);}
	void release(ByteBuilder bb, int id){tables.release(bb, id);}
	void release(byte[] array, int len, int id){tables.release(array, len, id);}
	int fillRightCounts(long kmer, long rkmer, int[] counts, long mask, int shift2){return tables.fillRightCounts(kmer, rkmer, counts, mask, shift2);}
	int fillLeftCounts(long kmer, long rkmer, int[] counts, long mask, int shift2){return tables.fillLeftCounts(kmer, rkmer, counts, mask, shift2);}
	static StringBuilder toText(long kmer, int k){return AbstractKmerTable.toText(kmer, k);}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	AbstractKmerTableSet tables(){return tables;}
	
	private final KmerTableSet tables;
	private final int k;
	
}

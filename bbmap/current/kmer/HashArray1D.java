package kmer;

import java.util.ArrayList;
import java.util.Arrays;

import shared.Primes;
import shared.Tools;

/**
 * Stores kmers in a long[] and counts in an int[], with a victim cache.
 * @author Brian Bushnell
 * @date Oct 25, 2013
 *
 */
public final class HashArray1D extends HashArray {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public HashArray1D(int initialSize, boolean autoResize_){
		super(initialSize, autoResize_, false);
		values=allocInt1D(prime+extra);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final int increment(final long kmer){
		int cell=(int)(kmer%prime);
		
		for(final int max=cell+extra; cell<max; cell++){
			long n=array[cell];
			if(n==kmer){
				values[cell]++;
				if(values[cell]<0){values[cell]=Integer.MAX_VALUE;}
				return values[cell];
			}else if(n==NOT_PRESENT){
				array[cell]=kmer;
				size++;
				values[cell]=1;
				if(autoResize && size+victims.size>sizeLimit){resize();}
				return 1;
			}
		}
		int x=victims.increment(kmer);
		if(autoResize && size+victims.size>sizeLimit){resize();}
		return x;
	}
	
	@Override
	public final int incrementAndReturnNumCreated(final long kmer){
		int cell=(int)(kmer%prime);
		
		for(final int max=cell+extra; cell<max; cell++){
			long n=array[cell];
			if(n==kmer){
				values[cell]++;
				if(values[cell]<0){values[cell]=Integer.MAX_VALUE;}
				return 0;
			}else if(n==NOT_PRESENT){
				array[cell]=kmer;
				size++;
				values[cell]=1;
				if(autoResize && size+victims.size>sizeLimit){resize();}
				return 1;
			}
		}
		return victims.incrementAndReturnNumCreated(kmer);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final int readCellValue(int cell) {
		return values[cell];
	}
	
	@Override
	protected final int[] readCellValues(int cell, int[] singleton) {
		singleton[0]=values[cell];
		return singleton;
	}
	
	@Override
	protected final void insertValue(long kmer, int v, int cell) {
		assert(array[cell]==kmer);
		values[cell]=v;
	}
	
	@Override
	protected final void insertValue(long kmer, int[] vals, int cell) {
		assert(array[cell]==kmer);
		assert(vals.length==1);
		values[cell]=vals[0];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final boolean canRebalance() {return false;}
	
	@Override
	protected synchronized void resize(){
//		assert(false);
//		System.err.println("Resizing from "+prime+"; load="+(size*1f/prime));
		if(prime>=maxPrime){
			sizeLimit=0xFFFFFFFFFFFFL;
			return;
		}
		
		final long oldSize=size, oldVSize=victims.size;
		final long totalSize=oldSize+oldVSize;
		
		final long maxAllowedByLoadFactor=(long)(totalSize*minLoadMult);
		final long minAllowedByLoadFactor=(long)(totalSize*maxLoadMult);

//		sizeLimit=Tools.min((long)(maxLoadFactor*prime), maxPrime);
		
		assert(maxAllowedByLoadFactor>=minAllowedByLoadFactor);
		if(maxAllowedByLoadFactor<prime){
			sizeLimit=(long)(maxLoadFactor*prime);
			return;
		}
		
		long x=10+(long)(prime*resizeMult);
		x=Tools.max(x, minAllowedByLoadFactor);
		x=Tools.min(x, maxAllowedByLoadFactor);
		
		int prime2=(int)Tools.min(maxPrime, Primes.primeAtLeast(x));
		
		if(prime2<=prime){
			sizeLimit=(long)(maxLoadFactor*prime);
			assert(prime2==prime) : "Resizing to smaller array? "+totalSize+", "+prime+", "+x;
			return;
		}
		
		prime=prime2;
//		System.err.println("Resized to "+prime+"; load="+(size*1f/prime));
		long[] oldk=array;
		int[] oldc=values;
		KmerNode[] oldv=victims.array;
		array=allocLong1D(prime2+extra);
		Arrays.fill(array, NOT_PRESENT);
		values=allocInt1D(prime2+extra);
		ArrayList<KmerNode> list=victims.toList();
		Arrays.fill(oldv, null);
		victims.size=0;
		size=0;
		sizeLimit=Long.MAX_VALUE;
		
		for(int i=0; i<oldk.length; i++){
			if(oldk[i]>NOT_PRESENT){set(oldk[i], oldc[i]);}
		}

		for(KmerNode n : list){
			if(n.pivot>NOT_PRESENT){set(n.pivot, n.value());}
		}
		
		assert(oldSize+oldVSize==size+victims.size) : oldSize+", "+oldVSize+" -> "+size+", "+victims.size;
		
		sizeLimit=(long)(maxLoadFactor*prime);
	}
	
	@Deprecated
	@Override
	public void rebalance(){
		throw new RuntimeException("Unimplemented.");
	}
	
	@Override
	public long regenerate(final int limit){
		long sum=0;
		assert(owners==null) : "Clear ownership before regeneration.";
		for(int pos=0; pos<values.length; pos++){
			final long key=array[pos];
			if(key>=0){
				final int value=values[pos];
				values[pos]=NOT_PRESENT;
				array[pos]=NOT_PRESENT;
				size--;
				if(value>limit){
					set(key, value);
				}else{
					sum++;
				}
			}
		}
		
		ArrayList<KmerNode> nodes=victims.toList();
		victims.clear();
		for(KmerNode node : nodes){
			int value=node.value();
			if(value<=limit){
				sum++;
			}else{
				set(node.pivot, node.value());
			}
		}
		
		return sum;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private int[] values;
	
	public int[] values(){return values;}
	

	
}

package ukmer;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import stream.ByteBuilder;
import shared.Primes;
import shared.Tools;
import fileIO.ByteStreamWriter;
import fileIO.TextStreamWriter;

/**
 * Stores kmers in a long[] and values in an int[][], with a victim cache.
 * @author Brian Bushnell
 * @date Nov 7, 2014
 *
 */
public abstract class HashArrayU extends AbstractKmerTableU {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	HashArrayU(int initialSize, int kbig_, boolean autoResize_, boolean twod){
		if(initialSize>1){
			initialSize=(int)Tools.min(maxPrime, Primes.primeAtLeast(initialSize));
		}else{
			initialSize=1;
		}
		prime=initialSize;
		sizeLimit=(long)(sizeLimit=(long)(maxLoadFactor*prime));
		kbig=kbig_;
		mult=Kmer.getMult(kbig);
		arrays=new long[mult][];
		for(int i=0; i<mult; i++){
			arrays[i]=allocLong1D(prime+extra);
			Arrays.fill(arrays[i], NOT_PRESENT);
		}
		victims=new HashForestU(Tools.max(10, initialSize/8), autoResize_, twod);
		autoResize=autoResize_;
		TWOD=twod;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
//	public final int set_Test(final long kmer, final int v){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD){
//			int[] old=getValues(kmer, new int[1]);
//			assert(old==null || contains(kmer, old));
//			if(verbose){System.err.println("Fetched "+Arrays.toString(old));}
//			x=set0(kmer, v);
//			assert(old==null || contains(kmer, old)) : "old="+Arrays.toString(old)+", v="+v+", kmer="+kmer+
//				", get(kmer)="+(Arrays.toString(getValues(kmer, new int[1])));
//			assert(contains(kmer, v));
//		}else{
//			int old=getValue(kmer);
//			assert(old==0 || old==-1 || contains(kmer, old));
//			x=set0(kmer, v);
//			assert(contains(kmer, v)) : "old="+old+", v="+v+", kmer="+kmer+", get(kmer)="+getValue(kmer);
//			assert(v==old || !contains(kmer, old));
//		}
//		return x;
//	}
//	
//	public final int set_Test(final long kmer, final int v[]){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD){
//			final int[] singleton=new int[1];
//			int[] old=getValues(kmer, singleton);
//			assert(old==null || contains(kmer, old));
//			if(verbose){System.err.println("Before: old="+Arrays.toString(old)+", v="+Arrays.toString(v));}
//			x=set0(kmer, v);
//			if(verbose){System.err.println("After:  old="+Arrays.toString(old)+", v="+Arrays.toString(v)+", get()="+Arrays.toString(getValues(kmer, singleton)));}
//			assert(old==null || contains(kmer, old)) : "old="+Arrays.toString(old)+", v="+Arrays.toString(v)+", kmer="+kmer+
//				", get(kmer)="+(Arrays.toString(getValues(kmer, new int[1])));
//			assert(contains(kmer, v)) : "old="+Arrays.toString(old)+", v="+Arrays.toString(v)+", kmer="+kmer+
//				", get(kmer)="+(Arrays.toString(getValues(kmer, new int[1])));
//		}else{
//			int old=getValue(kmer);
//			assert(old==0 || old==-1 || contains(kmer, old));
//			x=set0(kmer, v);
//			assert(contains(kmer, v)) : "old="+old+", v="+v+", kmer="+kmer+", get(kmer)="+getValue(kmer);
//			assert(v[0]==old || !contains(kmer, old));
//		}
//		return x;
//	}
//	
//	public final int setIfNotPresent_Test(long kmer, int v){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD){
////			int[] vals=getValues(kmer, null);
////			assert(vals==null || contains(kmer, vals));
////			x=setIfNotPresent(kmer, v);
////			assert(contains(kmer, vals));
////			assert(contains(kmer, v));
//			x=0;
//			assert(false);
//		}else{
//			int old=getValue(kmer);
//			assert(old==0 || old==-1 || contains(kmer, old));
//			x=setIfNotPresent0(kmer, v);
//			assert((old<1 && contains(kmer, v)) || (old>0 && contains(kmer, old))) : kmer+", "+old+", "+v;
//		}
//		return x;
//	}
	
	@Override
	public final int set(Kmer kmer, final int[] v){
		final int cell=findKmerOrEmpty(kmer);
		
		if(cell==HASH_COLLISION){
			if(verbose){System.err.println("C2: Adding "+kmer+", "+v+", "+cell);}
			final int x=victims.set(kmer, v);
			if(autoResize && size+victims.size>sizeLimit){resize();}
			if(verbose){System.err.println("C2: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}
			return x;
		}
		final long[] key=kmer.key();
		
		assert(cell>=0);
		
		final boolean notpresent=(arrays[0][cell]==NOT_PRESENT);
		if(notpresent){
			if(verbose){System.err.println("B2: Setting cell "+cell+" to kmer "+kmer);}
			setKmer(kmer.key(), cell);
		}
		
		if(verbose){System.err.println("A2: Adding "+kmer+", "+Arrays.toString(v)+", "+cell);}
		insertValue(key, v, cell);
		if(verbose){System.err.println("A2: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}
		
		if(notpresent){
			size++;
			if(autoResize && size+victims.size>sizeLimit){resize();}
			return 1;
		}else{
			return 0;
		}
	}
	
	public final void setKmer(long[] key, int cell){
		if(verbose){System.err.println();}
		for(int i=0; i<mult; i++){
			arrays[i][cell]=key[i];
		}
	}
	
	@Override
	public final int set(final Kmer kmer, final int v){
		assert(kmer.mult==mult && kmer.len>=kmer.kbig);
		final int cell=findKmerOrEmpty(kmer);
//		assert(kmer.verify(false)); //123
		
		if(cell==HASH_COLLISION){
			if(verbose){System.err.println("C2: Adding "+kmer+", "+v+", "+cell);}
			final int x=victims.set(kmer, v);
			if(autoResize && size+victims.size>sizeLimit){resize();}
			if(verbose){System.err.println("C2: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}
			return x;
		}
		assert(cell>=0);
		final long[] key=kmer.key();
		
		final boolean notpresent=(arrays[0][cell]==NOT_PRESENT);
		if(notpresent){
			if(verbose){System.err.println("B2: Setting cell "+cell+" to kmer "+kmer);}
			setKmer(key, cell);
		}
		
		if(verbose){System.err.println("A2: Adding "+kmer+", "+v+", "+cell);}
		insertValue(key, v, cell);
		if(verbose){System.err.println("A2: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}
		
		if(notpresent){
			size++;
			if(autoResize && size+victims.size>sizeLimit){resize();}
			return 1;
		}else{
			return 0;
		}
	}


//	protected LongList ll=new LongList(); //123
//	protected IntList il=new IntList();
	
	@Override
	public final int setIfNotPresent(Kmer kmer, int value){
		final int cell=findKmerOrEmpty(kmer);
		
		if(cell==HASH_COLLISION){
			int x=victims.setIfNotPresent(kmer, value);
			if(autoResize && size+victims.size>sizeLimit){resize();}
			return x;
		}
		assert(cell>=0);
		final long[] key=kmer.key();
		
		if(cell==NOT_PRESENT){
			setKmer(key, cell);
			insertValue(key, value, cell);
			size++;
			if(autoResize && size+victims.size>sizeLimit){resize();}
			return 1;
		}else{
			return 0;
		}
	}
	
	@Override
	public final int getValue(Kmer kmer){
		int cell=findKmer(kmer);
		if(cell==NOT_PRESENT){return NOT_PRESENT;}
		if(cell==HASH_COLLISION){return victims.getValue(kmer);}
		return readCellValue(cell);
	}

	/* (non-Javadoc)
	 * @see ukmer.AbstractKmerTableU#getValue(long[], long)
	 */
	@Override
	public int getValue(long[] key, long xor) {
		throw new RuntimeException("Unimplemented");
	}
	
	protected final long[] fillKey(int cell, long[] temp) {
		return fillKey(cell, temp, arrays);
	}
	
	public final Kmer fillKmer(int cell, Kmer kmer) {
		return fillKmer(cell, kmer, arrays); 
	}
	
	public final Kmer fillKmer(int cell, Kmer kmer, long[][] matrix) {
		long[] x=fillKey(cell, kmer.array1(), matrix);
//		assert(false) : x+"\ngetKmer("+cell+", kmer, matrix)"; //123
		if(x==null){return null;}
		kmer.fillArray2();
		if(verbose){System.err.println("Filled kmer "+kmer+": a1="+Arrays.toString(kmer.array1())+", a2="+Arrays.toString(kmer.array2())+", key="+Arrays.toString(kmer.key()));}
		return kmer;
	}
	
	protected final long[] fillKey(int cell, long[] temp, long[][] matrix) {
		assert(temp.length==mult);
		if(matrix[0][cell]<0){
//			assert(false) : matrix[0][cell]+"\ngetKmer("+cell+", kmer, matrix)\n"+Arrays.toString(matrix[0]); //123
			return null;
		}
		for(int i=0; i<temp.length; i++){
			temp[i]=matrix[i][cell];
		}
		if(verbose){System.err.println("cell="+cell+", matrix[0][cell]="+matrix[0][cell]+", temp="+Arrays.toString(temp)+"\nmatrix[0]="+Arrays.toString(matrix[0]));}
		return temp;
	}
	
	@Override
	public final int[] getValues(Kmer kmer, int[] singleton){
		int cell=findKmer(kmer);
		if(cell==NOT_PRESENT){
			singleton[0]=NOT_PRESENT;
			return singleton;
		}
		if(cell==HASH_COLLISION){return victims.getValues(kmer, singleton);}
		return readCellValues(cell, singleton);
	}
	
	@Override
	public final boolean contains(Kmer kmer){
		int cell=findKmer(kmer);
		if(cell==NOT_PRESENT){return false;}
		if(cell==HASH_COLLISION){return victims.contains(kmer);}
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Ownership           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final void initializeOwnership(){
		assert(owners==null);
		owners=allocAtomicInt(arrays[0].length);
		for(int i=0; i<arrays[0].length; i++){
			owners.set(i, NO_OWNER);
		}
		victims.initializeOwnership();
	}
	
	@Override
	public final void clearOwnership(){
		owners=null;
		victims.clearOwnership();
	}
	
	@Override
	public final int setOwner(final Kmer kmer, final int newOwner){
		final int cell=findKmer(kmer);
		assert(cell!=NOT_PRESENT);
		if(cell==HASH_COLLISION){return victims.setOwner(kmer, newOwner);}
		return setOwner(kmer, newOwner, cell);
	}
	
	public final int setOwner(final Kmer kmer, final int newOwner, final int cell){
//		kmer.verify(true);
		assert(matches(kmer.key(), cell)) : "cell="+cell+", key="+Arrays.toString(kmer.key())+", row="+Arrays.toString(cellToArray(cell))+"\n" +
				"kmer="+kmer+", array1="+Arrays.toString(kmer.array1())+", array2="+Arrays.toString(kmer.array2())+", row="+AbstractKmerTableU.toText(cellToArray(cell), kmer.k);
		final int original=owners.get(cell);
		int current=original;
		while(current<newOwner){
			boolean success=owners.compareAndSet(cell, current, newOwner);
			if(!success){current=owners.get(cell);}
			else{current=newOwner;}
		}
		assert(current>=original) : "original="+original+", current="+current+", newOwner="+newOwner+", re-read="+owners.get(cell);
		return current;
	}
	
	@Override
	public final boolean clearOwner(final Kmer kmer, final int owner){
		final int cell=findKmer(kmer);
		assert(cell!=NOT_PRESENT);
		if(cell==HASH_COLLISION){return victims.clearOwner(kmer, owner);}
		return clearOwner(kmer, owner, cell);
	}
	
	public final boolean clearOwner(final Kmer kmer, final int owner, final int cell){
		assert(matches(kmer.key(), cell));
		boolean success=owners.compareAndSet(cell, owner, NO_OWNER);
		return success;
	}
	
	@Override
	public final int getOwner(final Kmer kmer){
		final int cell=findKmer(kmer);
		assert(cell!=NOT_PRESENT);
		if(cell==HASH_COLLISION){return victims.getOwner(kmer);}
		return getCellOwner(cell);
	}
	
	public final int getCellOwner(final int cell){
		return owners.get(cell);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	protected abstract void insertValue(final long[] kmer, final int v, final int cell);
	
	protected abstract void insertValue(final long[] kmer, final int[] vals, final int cell);

	protected abstract int readCellValue(int cell);
	protected abstract int[] readCellValues(int cell, int[] singleton);

	@Override
	final Object get(long[] kmer){
		throw new RuntimeException("Unimplemented.");
	}
		
	final int findKmer(Kmer kmer){
		return findKmer(kmer.key(), kmer.xor());
	}

	final int findKmer(long[] key, long xor){
		int cell=(int)(xor%prime);
		for(final int max=cell+extra; cell<max; cell++){
			final long n=arrays[0][cell];
			if(n==key[0]){
				boolean success=true;
				for(int i=1; i<mult && success; i++){
					if(key[i]!=arrays[i][cell]){success=false;}
				}
				if(success){return cell;}
			}else if(n==NOT_PRESENT){return NOT_PRESENT;}
		}
		return HASH_COLLISION;
	}
	
	final int findKmerOrEmpty(Kmer kmer){
		int cell=(int)(kmer.xor()%prime);
		if(verbose){System.err.println("Started at cell "+cell+" for "+kmer);}
		
		long[] key=kmer.key();
		for(final int max=cell+extra; cell<max; cell++){
			final long n=arrays[0][cell];
			if(n==NOT_PRESENT){
				if(verbose){System.err.println("Chose empty cell "+cell+" for "+kmer);}
				return cell;
			}
			boolean success=true;
			for(int i=0; i<mult && success; i++){
				if(key[i]!=arrays[i][cell]){success=false;}
			}
			if(success){
				if(verbose){System.err.println("Found cell "+cell+" containing "+kmer);}
				return cell;
			}
		}
		return HASH_COLLISION;
	}
	
	final boolean matches(long[] key, int cell){
		assert(cell>=0);
		boolean success=true;
		for(int i=0; i<mult && success; i++){
			if(key[i]!=arrays[i][cell]){success=false;}
		}
		return success;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	final boolean canResize() {return true;}
	
	@Override
	final public long size() {return size;}
	
	@Override
	final public int arrayLength() {return arrays[0].length;}
	
	@Override
	protected abstract void resize();
	
	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/
	
	protected long[] cellToArray(int cell){throw new RuntimeException("Unimplemented");}
	
	@Override
	public final boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount){
		final long[] key=new long[mult];
		final int alen=arrays[0].length;
		if(TWOD){
			final int[] singleton=new int[1];
			for(int i=0; i<alen; i++){
				long[] temp=fillKey(i, key);
				if(temp!=null){
					tsw.print(toText(temp, readCellValues(i, singleton), k).append('\n'));
				}
			}
		}else{
			for(int i=0; i<alen; i++){
				long[] temp=fillKey(i, key);
				if(temp!=null && readCellValue(i)>=mincount){
					tsw.print(toText(temp, readCellValue(i), k).append('\n'));
				}
			}
		}
		if(victims!=null){
			victims.dumpKmersAsText(tsw, k, mincount);
		}
		return true;
	}
	
	@Override
	public final boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount){
		final long[] key=new long[mult];
		final int alen=arrays[0].length;
		if(TWOD){
			final int[] singleton=new int[1];
			for(int i=0; i<alen; i++){
				long[] temp=fillKey(i, key);
				if(temp!=null){
					bsw.printlnKmer(temp, readCellValues(i, singleton), k);
				}
			}
		}else{
			for(int i=0; i<alen; i++){
				long[] temp=fillKey(i, key);
				if(temp!=null && readCellValue(i)>=mincount){
					bsw.printlnKmer(temp, readCellValue(i), k);
				}
			}
		}
		if(victims!=null){
			victims.dumpKmersAsBytes(bsw, k, mincount);
		}
		return true;
	}
	
	@Override
	public final boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount){
		final long[] key=new long[mult];
		final int alen=arrays[0].length;
		if(TWOD){
			final int[] singleton=new int[1];
			for(int i=0; i<alen; i++){
				long[] temp=fillKey(i, key);
				if(temp!=null){
					toBytes(temp, readCellValues(i, singleton), k, bb);
					bb.append('\n');
					if(bb.length()>=16000){
						ByteBuilder bb2=new ByteBuilder(bb);
						synchronized(bsw){bsw.addJob(bb2);}
						bb.clear();
					}
				}
			}
		}else{
			for(int i=0; i<alen; i++){
				long[] temp=fillKey(i, key);
				if(temp!=null && readCellValue(i)>=mincount){
					toBytes(temp, readCellValue(i), k, bb);
					bb.append('\n');
					if(bb.length()>=16000){
						ByteBuilder bb2=new ByteBuilder(bb);
						synchronized(bsw){bsw.addJob(bb2);}
						bb.clear();
					}
				}
			}
		}
		if(victims!=null){
			victims.dumpKmersAsBytes_MT(bsw, bb, k, mincount);
		}
		return true;
	}
	
	@Override
	public final void fillHistogram(long[] ca, int max){
		final int alen=arrays[0].length;
		for(int i=0; i<alen; i++){
			long kmer=arrays[0][i];
			if(kmer!=NOT_PRESENT){
				int count=Tools.min(readCellValue(i), max);
				ca[count]++;
			}
		}
		if(victims!=null){
			victims.fillHistogram(ca, max);
		}
	}
	
	@Override
	public final void countGC(long[] gcCounts, int max){
		final int alen0=arrays.length;
		final int alen=arrays[0].length;
		for(int i=0; i<alen; i++){
			long kmer0=arrays[0][i];
			if(kmer0!=NOT_PRESENT){
				int count=Tools.min(readCellValue(i), max);
				for(int j=0; j<alen0; j++){
					gcCounts[count]+=gc(arrays[j][i]);
				}
			}
		}
		if(victims!=null){
			victims.countGC(gcCounts, max);
		}
	}
	
	public HashForestU victims(){
		return victims;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	AtomicIntegerArray owners;
	long[][] arrays;
	int prime;
	long size=0;
	long sizeLimit;
	final HashForestU victims;
	final boolean autoResize;
	final int kbig;
	final int mult;//Length of Kmer arrays.
	public final boolean TWOD;
	private final Lock lock=new ReentrantLock();
	
	public AtomicIntegerArray owners() {return owners;}
	
	@Override
	final Lock getLock(){return lock;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	final static int extra=21;
	final static int maxPrime=(int)Primes.primeAtMost(Integer.MAX_VALUE-extra);
	final static float resizeMult=2f; //Resize by a minimum of this much
	final static float minLoadFactor=0.58f; //Resize by enough to get the load above this factor
	final static float maxLoadFactor=0.905f; //Reaching this load triggers resizing
	final static float minLoadMult=1/minLoadFactor;
	final static float maxLoadMult=1/maxLoadFactor;
	
}

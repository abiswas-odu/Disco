package clump;

import java.io.Serializable;

import stream.KillSwitch;
import stream.Read;
import structures.IntList;

class ReadKey implements Serializable, Comparable<ReadKey> {
	
//	public static ReadKey makeKeyIfNeeded(Read r){
//		if(r.obj==null){
//			return makeKey(r, true);
//		}
//		return (ReadKey)r.obj;
//	}
	
	public static ReadKey makeKey(Read r, boolean setObject){
		assert(r.obj==null);
		try {
			ReadKey rk=new ReadKey(r);
			if(setObject){r.obj=rk;}
			return rk;
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
			throw new RuntimeException(e);
		}
	}
	
	private ReadKey(Read r){
		this(r, 0, 0, true);
	}
	
	private ReadKey(Read r, long kmer_, int position_, boolean plus_){
		kmer=kmer_;
		position=position_;
		clump=null;
		assert(!r.swapped());
		flipped=false;
		kmerMinusStrand=false;
	}
	
	protected ReadKey(){}
	
	public void set(long kmer_, int position_, boolean minus_){
		setKmer(kmer_);
		setPosition(position_);
//		setClump(null);
		kmerMinusStrand=minus_;
	}
	
	private long setKmer(long kmer_){
		return kmer=kmer_;
	}
	
	private int setPosition(int position_){
		return position=position_;
	}
	
	public Clump setClump(Clump clump_){
		return clump=clump_;
	}
	
	private boolean setFlipped(boolean flipped_){
		assert(flipped!=flipped_);
		return flipped=flipped_;
	}
	
	public void clear(){
		setKmer(0);
		setPosition(0);
		setClump(null);
		kmerMinusStrand=false;
	}
	
	public void flip(Read r, int k){
		assert(r.swapped()==flipped);
		r.reverseComplement();
		r.setSwapped(!r.swapped());
		setFlipped(!flipped);
		setPosition(r.length()-position+k-2);
		assert(r.swapped()==flipped);
	}
	
	@Override
	public int compareTo(ReadKey b){
		if(kmer!=b.kmer){return kmer>b.kmer ? 1 : -1;}
		if(kmerMinusStrand!=b.kmerMinusStrand){return kmerMinusStrand ? 1 : -1;}
		if(position!=b.position){return position<b.position ? 1 : -1;}
		return 0;
	}
	
	@Override
	public boolean equals(Object b){
		return equals((ReadKey)b);
	}
	
	public boolean equals(ReadKey b){
		return b==null ? false : compareTo(b)==0;
	}
	
	public String toString(){
		return position+","+(kmerMinusStrand ? ",t" : ",f")+","+kmer;
	}
	
	public long kmer;
	public int position;
	public boolean flipped;
	public boolean kmerMinusStrand;
	public Clump clump;
	public IntList vars;
	
	public static final int overhead=overhead();
	private static int overhead(){
		return 16+ //self
				1*(8)+ //kmer
				1*(4)+ //int fields
				2*(4)+ //booleans
				2*(8); //pointers
	}
	
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
}

package structures;

import java.util.Arrays;

import shared.Shared;
import shared.Tools;
import stream.KillSwitch;



public final class IntList{
	
	public IntList(){this(256);}
	
	public IntList(int initial){
		assert(initial>0) : initial+"\n"+this;
		array=new int[initial];
	}
	
	public IntList copy() {
		IntList copy=new IntList(size);
		copy.addAll(this);
		return copy;
	}
	
	public void clear(){size=0;}
	
	public final void set(int loc, int value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]=value;
		size=max(size, loc+1);
	}
	
	public final void increment(int loc, int value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]+=value;
		size=max(size, loc+1);
	}
	
	public final int get(int loc){
		return(loc>=size ? 0 : array[loc]);
	}
	
	public final void add(int x){
		if(size>=array.length){
			resize(size*2L+1);
		}
		array[size]=x;
		size++;
	}
	
	public void addAll(IntList counts) {
		final int[] array2=counts.array;
		final int size2=counts.size;
		for(int i=0; i<size2; i++){add(array2[i]);}
	}
	
	public boolean contains(int x) {
		for(int i=0; i<size; i++){
			if(array[i]==x){return true;}
		}
		return false;
	}
	
	public final void setSize(final int size2) {
		if(size2<array.length){resize(size2);}
		size=size2;
	}
	
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;
		final int size3=(int)Tools.min(Integer.MAX_VALUE, size2);
		assert(size3>size) : "Overflow: "+size+", "+size2+" -> "+size3;
		array=KillSwitch.copyOf(array, size3);
	}
	
	public final void shrink(){
		if(size==array.length){return;}
		array=KillSwitch.copyOf(array, size);
	}
	
	public final void shrinkToUnique(){
		//Assumes sorted.
		if(size<=0){
			shrink();
			return;
		}
		
		int unique=1;
		
		for(int i=1; i<size; i++){
			assert(array[i]>=array[i-1]);
			if(array[i]!=array[i-1]){unique++;}
		}
		if(unique==array.length){return;}
		int[] alt=new int[unique];
		
		alt[0]=array[0];
		for(int i=1, j=1; j<unique; i++){
			if(array[i]!=array[i-1]){
				alt[j]=array[i];
				j++;
			}
		}
		
		array=alt;
		size=alt.length;
	}
	
	public int[] toArray(){
		return KillSwitch.copyOf(array, size);
	}
	
	public String toString(){
		return toStringListView();
	}
	
	public String toStringSetView(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<size; i++){
			if(array[i]!=0){
				sb.append(comma+"("+i+", "+array[i]+")");
				comma=", ";
			}
		}
		sb.append(']');
		return sb.toString();
	}
	
	public String toStringListView(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<size; i++){
				sb.append(comma+array[i]);
				comma=", ";
		}
		sb.append(']');
		return sb.toString();
	}
	
	public void sort() {
		if(size>1){Shared.sort(array, 0, size);}
	}
	
	public void reverse() {
		Tools.reverseInPlace(array, 0, size);
	}
	
	/** Assumes this is sorted.
	 * Reduces the list to a set of unique values;
	 * stores their counts in a second list. */
	public void getUniqueCounts(IntList counts) {
		counts.size=0;
		if(size<=0){return;}

		int unique=1;
		int count=1;
		
		for(int i=1; i<size; i++){
			assert(array[i]>=array[i-1]);
			if(array[i]==array[i-1]){
				count++;
			}else{
				array[unique]=array[i];
				unique++;
				counts.add(count);
				count=1;
			}
		}
		if(count>0){
			counts.add(count);
		}
		size=unique;
		assert(counts.size==size);
	}
	
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	public int[] array;
	public int size=0;
	
}

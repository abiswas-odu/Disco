package structures;

import java.util.Arrays;



public final class IntList2{
	
	public IntList2(){this(256);}
	
	public IntList2(int initial){
		assert(initial>0);
		array=new int[initial][];
	}
	
	public final void set(int loc, int[] value){
		if(loc>=array.length){
			resize((loc+1)*2);
		}
		array[loc]=value;
		size=max(size, loc+1);
	}
	
	public final void increment(int loc, int value){
		throw new RuntimeException("Unsupported");
	}
	
	public final int[] get(int loc){
		return(loc>=size ? null : array[loc]);
	}
	
	@Deprecated
	public final void add(int x){
		throw new RuntimeException("Unsupported");
	}
	
	public final void add(int[] x){
		if(size>=array.length){
			resize(max(size*2, 1));
		}
		array[size]=x;
		size++;
	}
	
	public final void resize(int size2){
		assert(size2>size);
		array=Arrays.copyOf(array, size2);
	}
	
	public final void shrink(){
		if(size==array.length){return;}
		array=Arrays.copyOf(array, size);
	}
	
	public final void shrinkToUnique(){
		throw new RuntimeException("Unsupported");
	}
	
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<size; i++){
			if(array[i]!=null){
				sb.append(comma+"("+i+", "+Arrays.toString(array[i])+")");
				comma=", ";
			}
		}
		sb.append(']');
		return sb.toString();
	}
	
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	public int[][] array;
	public int size=0;
	
}

package structures;

/**
 * Maintains a heap of unique values.
 * @author Brian Bushnell
 * @date July 6, 2016
 *
 */
public class LongHeapSet {
	
	public LongHeapSet(int limit_){
		limit=limit_;
		heap=new LongHeap(limit, true);
		set=new LongHashSet(limit*2);
	}
	
	public boolean add(long value){
		if(heap.hasRoom()){
			if(set.add(value)){
				heap.add(value);
				return true;
			}
			return false;
		}
		
		final long bottom=heap.peek();
		if(value>bottom){
			if(set.add(value)){
				set.remove(bottom);
				assert(set.size()<=limit);
				heap.add(value);
				assert(heap.size()<=limit);
				return true;
			}
		}
		return false;
	}
	
	public void clear(){
		heap.clear();
		set.clear();
	}
	
	public void add(LongHeapSet b){
		final long[] array=b.heap.array();
		final int blen=b.heap.size();
		for(int i=1; i<=blen; i++){
			add(array[i]);
		}
	}
	
	public int size(){return heap.size();}
	
	final int limit;
	public LongHeap heap;
	public LongHashSet set;
	
}

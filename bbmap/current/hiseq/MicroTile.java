package hiseq;

public class MicroTile {

	public MicroTile(){}

	public MicroTile(int lane_, int tile_, int x1_, int x2_, int y1_, int y2_){
		lane=lane_;
		tile=tile_;
		x1=x1_;
		x2=x2_;
		y1=y1_;
		y2=y2_;
	}
	
	public boolean contains(int x, int y){
		return x>=x1 && x<=x2 && y>=y1 && y<=y2;
	}
	
	public String toString(){
		return lane+", "+tile+", "+x1+", "+x2+", "+y1+", "+y2;
	}
	
	public double averageQuality(){
		return qualityCount==0 ? 0 : qualitySum/qualityCount;
	}
	
	public double percentErrorFree(){
		return qualityCount==0 ? 0 : errorFreeSum/qualityCount;
	}
	
	public double hitPercent(){
		long count=kmerCount();
		return count==0 ? 0 : hits*100.0/count;
	}
	
	public double uniquePercent(){
		long count=kmerCount();
		return count==0 ? 0 : misses*100.0/count;
	}

	public long kmerCount(){return hits+misses;}
	
	public void add(MicroTile mt) {
		hits+=mt.hits;
		misses+=mt.misses;
		qualityCount+=mt.qualityCount;
		qualitySum+=mt.qualitySum;
		errorFreeSum+=mt.errorFreeSum;
	}
	
	public long hits;
	public long misses;
	public long qualityCount;
	public double qualitySum;
	public double errorFreeSum;
	
	public int discard=0;
	
	public int lane;
	public int tile;
	public int x1, x2;
	public int y1, y2;
	
}

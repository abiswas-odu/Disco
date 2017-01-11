package sketch;

public final class Comparison extends SketchObject implements Comparable<Comparison> {
	
	public Comparison(){}
	
	public Comparison(int[] buffer){
		this(buffer, null);
	}
	
	public Comparison(Sketch b){
		this(null, b);
	}
	
	public Comparison(int[] buffer, Sketch b){
		if(buffer!=null){
			hits=buffer[0];
			minIndex=buffer[1];
			minLen=buffer[2];
			maxLen=buffer[3];
		}
		
		if(b!=null){
			name=b.name();
			taxID=b.taxID;
			imgID=b.imgID;
			genomeSize=b.genomeSize;
		}
	}

	@Override
	public int compareTo(Comparison b) {
		int x=hits-b.hits;
		if(x!=0){return x;}
		x=b.minIndex-minIndex;
		if(x!=0){return x;}
		x=b.minLen-minLen;
		if(x!=0){return x;}
		x=b.maxLen-maxLen;
		if(x!=0){return x;}
		x=taxID-b.taxID;
		if(x!=0){return x;}
		if(name!=null && b.name!=null){
			return name.compareTo(b.name);
		}
		return 0;
	}
	
	public boolean equals(Object b){
		if(b==null || b.getClass()!=this.getClass()){return false;}
		return compareTo((Comparison)b)==0;
	}
	
	public float idWeighted(){
		return hits/(float)minIndex;
	}
	
	public float idMin(){
		return hits/(float)minLen;
	}
	
	public float idMax(){
		return hits/(float)maxLen;
	}
	
	public String name;
	public int taxID;
	public long imgID;
	public long genomeSize;
	
	public int hits;
	public int minIndex;
	public int minLen;
	public int maxLen;
	
}

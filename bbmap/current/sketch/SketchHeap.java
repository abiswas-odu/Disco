package sketch;

import structures.LongHeapSet;

public class SketchHeap extends LongHeapSet {
	
	public SketchHeap(int limit_){
		super(limit_);
	}
	
	public void clear(){
		taxID=-1;
		imgID=-1;
		genomeSize=0;
		taxName=null;
		name0=null;
		super.clear();
	}
	
	public void add(SketchHeap b){
		if(taxID<0){taxID=b.taxID;}
		if(imgID<0){imgID=b.imgID;}
		if(taxName==null){taxName=b.taxName;}
		if(name0==null){name0=b.name0;}
		genomeSize+=b.genomeSize;
		super.add(b);
	}
	
	public StringBuilder toHeader(){
		StringBuilder sb=new StringBuilder();
		sb.append("#SZ:"+size());
		if(genomeSize>0){sb.append("\tGS:"+genomeSize);}
		if(taxID>=0){sb.append("\tID:"+taxID);}
		if(imgID>=0){sb.append("\tIMG:"+imgID);}
		if(taxName!=null){sb.append("\tNM:"+taxName);}
		if(name0!=null){sb.append("\tNM0:"+name0);}
		return sb;
	}
	
	public String toString(){return toHeader().toString();}

	public String name(){return taxName==null ? name0 : taxName;}
	public String taxName(){return taxName;}
	public String name0(){return name0;}
	public void setTaxName(String s){taxName=s;}
	public void setName0(String s){name0=s;}
	
	private String taxName;
	private String name0;
	public long taxID=-1;
	public long imgID=-1;
	public long genomeSize=0;
	
}

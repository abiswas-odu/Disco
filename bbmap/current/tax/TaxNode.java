package tax;

import java.io.Serializable;
import java.util.Comparator;

import shared.Tools;

/**
 * Represents a taxonomic identifier, such as a specific genus.
 * Includes the name, NCBI numeric id, parent id, and taxonomic level.
 * @author Brian Bushnell
 * @date Mar 6, 2015
 *
 */
public class TaxNode implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -926484721933977347L;
	
	public TaxNode(int id_, String name_){
		this(id_, -1, -1, name_);
	}
	
	public TaxNode(int id_, int parent_, int level_, String name_){
		id=id_;
		pid=parent_;
		level=level_;
		name=name_;
	}

	/**
	 * @param split
	 * @param i
	 * @return
	 */
	public boolean matchesName(String[] split, int idx, TaxTree tree) {
		if(idx<0){return true;}
		if(!split[idx].equalsIgnoreCase(name)){return false;}
		return tree.getNode(pid).matchesName(split, idx-1, tree);
	}
	
	public String toString(){
		return "("+id+","+pid+","+countRaw+","+countSum+",'"+(level<0 ? "?" : TaxTree.levelToString(level))+",'"+(canonical ? "T" : "F")+",'"+name+"')";
	}
	
	public boolean equals(TaxNode b){
		if(id!=b.id || pid!=b.pid || level!=b.level || canonical!=b.canonical){return false;}
		if(name==b.name){return true;}
		if((name==null) != (b.name==null)){return false;}
		return name.equals(b.name);
	}
	
	public long incrementRaw(long amt){
		if(amt==0){return countRaw;}
		if(verbose){System.err.println("incrementRaw("+amt+") node: "+this);}
		countRaw+=amt;
		assert(countRaw>=0) : "Overflow! "+countRaw+", "+amt;
		return countRaw;
	}
	
	public long incrementSum(long amt){
		if(amt==0){return countSum;}
		if(verbose){System.err.println("incrementSum("+amt+") node: "+this);}
		countSum+=amt;
		assert(countSum>=0) : "Overflow! "+countSum+", "+amt;
		return countSum;
	}
	
	public String levelString(){return level<0 ? "unknown" : TaxTree.levelToString(level);}

	public String levelToStringShort() {return level<0 ? "x" : TaxTree.levelToStringShort(level);}
	
	public static class CountComparator implements Comparator<TaxNode>{
		
		@Override
		public int compare(TaxNode a, TaxNode b) {
			long x=b.countSum-a.countSum;
//			System.err.println("x="+x+" -> "+Tools.longToInt(x));
			if(x!=0){return Tools.longToInt(x);}
			return a.level==b.level ? a.id-b.id : a.level-b.level;
		}
		
	}
	
	@Override
	public final int hashCode(){return id;}

	public final int id;
	public final String name;
	public int pid;
	public int level;
	public boolean canonical=true;

	public long countRaw=0;
	public long countSum=0;
	
	public static final boolean verbose=false;
	public static final CountComparator countComparator=new CountComparator();
	
	
}

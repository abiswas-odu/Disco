package sketch;

import java.util.Arrays;

import tax.TaxTree;

public class SketchObject {
	
	static synchronized void setTaxtree(String taxTreeFile){
		if(taxTreeFile==null){
			taxtree=null;
			return;
		}
		if(treefile!=null){
			assert(!treefile.equals(taxTreeFile));
			if(treefile.equals(taxTreeFile)){return;}
			treefile=taxTreeFile;
		}
		taxtree=TaxTree.loadTaxTree(taxTreeFile, System.err, false);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final int RAW=0, HEX=1, A48=2;
	public static final char[] codingArray={'R', 'H', 'A'};
	
	public static int CODING=A48;
	public static boolean delta=true;
	public static float maxGenomeFraction=0.02f;
	
	static final byte[] hexTable=new byte[128];
	static {
		Arrays.fill(hexTable, (byte)-1);
		for(int i='0'; i<='9'; i++){
			hexTable[i]=(byte)(i-'0');
		}
		for(int i='A'; i<='F'; i++){
			hexTable[i]=hexTable[i+'a'-'A']=(byte)(i-'A'+10);
		}
		hexTable['x']=hexTable['X']=hexTable['-']=hexTable['+']=0;
	}
	
	static TaxTree taxtree=null;
	private static String treefile=null;
	
	static final int ONE_SKETCH=1, PER_SEQUENCE=2, PER_TAXA=3, IMG=4;
	
}

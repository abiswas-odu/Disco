package var2;

import dna.Gene;
import shared.Tools;
import stream.SamLine;
import structures.CoverageArray;
import structures.CoverageArray2;
import structures.CoverageArray3;

public class Scaffold {
	
	/** Assumes SAM format.
	 * e.g.<br> @SQ	SN:scaffold_0	LN:1785514	AS:build 9 */
	public Scaffold(byte[] line, int scafnum){
		assert(Tools.startsWith(line, "@SQ\t")) : new String(line);
		number=scafnum;
		int a=0, b=0;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(line);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 1: "+new String(line);
		assert(Tools.startsWith(line, "SN:", a));
		name=new String(line, a+3, b-a-3);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 2: "+new String(line);
		assert(Tools.startsWith(line, "LN:", a));
		length=Tools.parseInt(line, a+3, b);
		b++;
		a=b;
	}
	
	public Scaffold(String name_, int scafnum_, int len_){
		name=name_;
		number=scafnum_;
		length=len_;
	}
	
	public void add(SamLine sl){
		int start=sl.pos-1;
		int stop=sl.stop(start, false, false);
		increment(start, stop, sl.strand());
	}
	
	public synchronized void increment(int from, int to, int strand){
		if(ca==null){
			ca=useCA3 ? new CoverageArray3(number, length) : new CoverageArray2(number, length);
		}
		ca.incrementRange(from, to, 1);
		
		if(trackStrand && strand==Gene.MINUS){
			if(caMinus==null){
				caMinus=useCA3 ? new CoverageArray3(number, length) : new CoverageArray2(number, length);
			}
			caMinus.incrementRange(from, to, 1);
		}
	}
	
	public int calcCoverage(Var v){
		return calcCoverage(v, ca);
	}
	
	public int minusCoverage(Var v){
		assert(trackStrand);
		return calcCoverage(v, caMinus);
	}
	
	public int calcCoverage(Var v, CoverageArray ca){
		final int a=v.start;
		final int b=v.stop;
//		assert(false) : ca.maxIndex+", "+a;
		if(ca==null || ca.maxIndex<a){return 0;}
		final int type=v.type();
		final int avg;
		final int rlen=v.reflen();
		long sum=0;
		if(type==Var.SUB || type==Var.NOCALL || type==Var.DEL){
			for(int i=a; i<b; i++){
				sum+=ca.get(i);
			}
			avg=Math.round(sum/(float)rlen);
		}else if(type==Var.INS){
			assert(rlen==0 && a==b);
//			if(a<=0){sum=2*ca.get(0);}
//			else if(b>ca.maxIndex)
			if(b>=ca.maxIndex){
				sum=2*ca.get(ca.maxIndex);
				avg=(int)(sum/2);
			}else{
				sum=ca.get(a)+ca.get(b);
				avg=(int)Math.ceil(sum/2);
			}
		}else{
			throw new RuntimeException("Unknown type "+type);
		}
		return avg;
	}
	
	public String toString(){
		return "@SQ\tSN:"+name+"\tLN:"+length+"\tID:"+number;
	}
	
	public void clearCoverage(){
		ca=null;
		caMinus=null;
	}
	
	public final String name;
	public final int number;
	public final int length;
	public CoverageArray ca;
	public CoverageArray caMinus;
	public byte[] bases;

	public static boolean useCA3=false;
	public static boolean trackStrand=true;
	
}

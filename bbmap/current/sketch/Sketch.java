package sketch;

import shared.Tools;
import stream.ByteBuilder;
import structures.LongList;
import tax.ImgRecord;

/**
 * @author Brian Bushnell
 * @date July 7, 2016
 *
 */
public class Sketch extends SketchObject implements Comparable<Sketch> {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	public Sketch(long[] array_){
		this(array_, -1, -1, -1, null, null);
	}
	
	public Sketch(SketchHeap heap){
		this(SketchTool.toSketchArray(heap.heap, (int)(2+maxGenomeFraction*heap.genomeSize)), 
				(int)heap.taxID, heap.imgID, heap.genomeSize, heap.taxName(), heap.name0());
	}
	
	public Sketch(long[] array_, int taxID_, long imgID_, long gSize_, String taxName_, String name0_){
		array=array_;
		taxID=taxID_;
		imgID=imgID_;
		genomeSize=gSize_;
		
		taxName=taxName_;
		name0=name0_;
		
		if(ImgRecord.imgMap!=null && imgID>=0 && taxID<0){
			ImgRecord record=ImgRecord.imgMap.get(imgID);
			if(record!=null){
				if(record.name!=null && taxName==null){taxName=record.name;}
				taxID=record.taxID;
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	public int countMatches(Sketch other, int[] buffer){
		return countMatches(array, other.array, buffer);
	}
	
	public void toBinary(final int bits){
		assert(binaryArray==null);
		binaryCardinality=0;
		final int len=(bits+63)/64;
		binaryArray=new long[len];
		for(long x : array){
			final int bitIndex=(int)(x%bits);
			final int index=bitIndex/64;
			final int shift=bitIndex-64*index;
			binaryArray[index]|=(1L<<shift);
		}
		for(long x : binaryArray){binaryCardinality+=Long.bitCount(x);}
	}
	
	public void add(Sketch other, int maxlen){
		final long[] a=array;
		final long[] b=other.array;
		if(maxlen<1){
			assert(false);
			maxlen=1000000;
		}
		LongList list=new LongList(Tools.min(maxlen, a.length+b.length));
		binaryArray=null;
		
		for(int i=0, j=0; i<a.length && j<b.length; ){
			final long ka=a[i], kb=b[j];
			if(ka==kb){//match
				list.add(ka);
				i++;
				j++;
			}else if(ka<kb){
				list.add(ka);
				i++;
			}else{
				list.add(kb);
				j++;
			}
			if(list.size()>=maxlen){break;}
		}
		
		if(array.length==list.size()){
			for(int i=0; i<list.size; i++){
				array[i]=list.array[i];
			}
		}else{
			array=list.toArray();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Comparison          ----------------*/
	/*--------------------------------------------------------------*/

//	public static final float[] makeBuffer(){return new float[5];}
	

	public static final int[] makeBuffer(){return new int[4];}
	
	/** Buffer: {matches, weightedDivisor, minLen, maxLen} */
	public static final int countMatches(long[] a, long[] b, int[] buffer){
		assert(a.length>0 && b.length>0);
		int matches=0;
		int i=0, j=0;
		for(; i<a.length && j<b.length; ){
			final long ka=a[i], kb=b[j];
			if(ka==kb){
				matches++;
				i++;
				j++;
			}else if(ka<kb){
				i++;
			}else{
				j++;
			}
		}
		if(buffer!=null){
			final int divisorWeighted=Tools.max(1, Tools.min(i, j));
			buffer[0]=matches;
			buffer[1]=divisorWeighted;
			buffer[2]=Tools.min(a.length, b.length);
			buffer[3]=Tools.max(a.length, b.length);
		}
		return matches;
	}
	
	public static long countMatchesBinary(long[] a, long[] b){
		long matches=0;
		for(int i=0, j=0; i<a.length && j<b.length; ){
			final long ka=a[i], kb=b[j];
			matches+=Long.bitCount(ka&kb);
		}
		return matches;
	}
	
	public float identityBinary(Sketch b){
		long matches=countMatchesBinary(binaryArray, b.binaryArray);
		return matches/(float)(Tools.max(1, Tools.min(binaryCardinality, b.binaryCardinality)));
	}
	
//	public float identity(Sketch b, float[] ret){
//		if(ret!=null){Arrays.fill(ret, 0);}
//		return identityWeighted(array, b.array, ret);
//	}
//	
//	public static float identity(long[] a, long[] b){
//		int matches=countMatches(a, b);
//		return matches/(float)(Tools.max(1, Tools.min(a.length, b.length)));
//	}
	
	@Override
	public int hashCode(){
		return taxID>=0 ? taxID : name()!=null ? name().hashCode() : Integer.rotateRight(super.hashCode(), 4);
	}
	
	@Override
	public int compareTo(Sketch b){
		if(this==b){return 0;}
		if(taxID>-1 && b.taxID>-1){return taxID-b.taxID;}
		int x=taxName.compareTo(b.taxName);
		if(x!=0){return x;}
		if(name0!=null && b.name0!=null){return name0.compareTo(b.name0);}
		return name0!=null ? 1 : b.name0!=null ? -1 : 0;
	}
	
	@Override
	public boolean equals(Object b){
		if(this==b){return true;}
		if(b==null || this.getClass()!=b.getClass()){return false;}
		return equals((Sketch)b);
	}
	
	public boolean equals(Sketch b){
		return compareTo(b)==0;
	}
	
	public ByteBuilder toHeader(){
		ByteBuilder sb=new ByteBuilder();
		sb.append("#SZ:").append(array.length);
		sb.append("\tCD:");
		sb.append(codingArray[CODING]);
		if(delta){sb.append('D');}
		if(genomeSize>0){sb.append("\tGS:").append(genomeSize);}
		if(taxID>=0){sb.append("\tID:").append(taxID);}
		if(imgID>=0){sb.append("\tIMG:").append(imgID);}
		if(taxName!=null){sb.append("\tNM:").append(taxName);}
		if(name0!=null){sb.append("\tNM0:").append(name0);}
		return sb;
	}
	
	public ByteBuilder toBytes(){
		long prev=0;
		ByteBuilder sb=toHeader();
		sb.append("\n");
		byte[] temp=null;
		if(CODING==A48){temp=new byte[12];}
		for(int i=0; i<array.length; i++){
			long key=array[i];
			long x=key-prev;
			if(CODING==A48){
				appendA48(x, sb, temp);
				sb.append('\n');
			}else if(CODING==HEX){
				sb.append(Long.toHexString(x)).append('\n');
			}else if(CODING==RAW){
				sb.append(x).append('\n');
			}else{
				assert(false);
			}
			//if(delta){prev=key;}
		}
		return sb;
	}
	
	public static final void appendA48(long value, ByteBuilder sb, byte[] temp){
		int i=0;
		while(value!=0){
			byte b=(byte)(value&0x3F);
			temp[i]=b;
			value=value>>6;
			i++;
		}
		if(i==0){
			sb.append((byte)'0');
		}else{
			for(i--;i>=0;i--){
				sb.append((char)(temp[i]+48));
			}
		}
	}
	
	public String toString(){
		return toBytes().toString();
	}
	
	public static long parseA48(byte[] line){
		if(line.length==0){return 0;}
		long x=0;
		for(byte b : line){
			x<<=6;
			x|=(((long)b)-48);
		}
		return x;
	}
	
	public static long parseHex(byte[] line){
		if(line.length==0){return 0;}
		long x=0;
		for(byte b : line){
			x<<=4;
			x|=hexTable[b];
		}
		if(line[0]=='-'){x*=-1;}
		return x;
	}
	
	public static boolean parseCoding(String a, String b){
		if(a.equals("delta")){
			delta=Tools.parseBoolean(b);
		}else if(a.equals("a33") || a.equals("a48")){
			boolean x=Tools.parseBoolean(b);
			if(x){CODING=A48;}
			else if(CODING==A48){CODING=HEX;}
		}else if(a.equals("hex")){
			boolean x=Tools.parseBoolean(b);
			if(x){CODING=HEX;}
			else if(CODING==HEX){CODING=A48;}
		}else if(a.equals("raw")){
			boolean x=Tools.parseBoolean(b);
			if(x){CODING=RAW;}
			else if(CODING==RAW){CODING=A48;}
		}else{
			return false;
		}
		return true;
	}

	public String name(){return taxName==null ? name0 : taxName;}
	public String taxName(){return taxName;}
	public String name0(){return name0;}
	public void setTaxName(String s){taxName=s;}
	public void setName0(String s){name0=s;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public long[] array;
	public int taxID;
	public long imgID;
	public long genomeSize;
	private String taxName;
	private String name0;
	
	public long[] binaryArray;
	public long binaryCardinality;
}

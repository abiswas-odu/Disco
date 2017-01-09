package dna;
import java.util.Arrays;
import java.util.HashMap;

import align2.QualityTools;
import shared.Tools;
import stream.Read;


/**
 * @author Brian Bushnell
 * @date July 1, 2010
 *
 */
public final class AminoAcid {
	
	
	public static void main(String[] args){
//		for(String s : stringToAA.keySet()){
//			System.out.println(s+"\t->\t"+stringToAA.get(s));
//		}
		
		String bases="atctgatTGGcgcgatatatcg";
		String acids=stringToAAs(bases);
		
		System.out.println(bases+" -> "+acids);
		
	}
	
	
	private AminoAcid(){
		this(null);
		assert(false);
		System.exit(0);
	}
	
	private AminoAcid(String line){
		String[] s2=line.split(", ");
		String[] s3=new String[s2.length-3];
		for(int i=3; i<s2.length; i++){
			s3[i-3]=s2[i];
		}
		
		name=s2[0];
		symbol=s2[1];
		letter=s2[2].charAt(0);
		codeStrings=s3;
	}
	
	private AminoAcid(String n, String c3, String c1, String[] bases){
		name=n;
		symbol=c3;
		letter=c1.charAt(0);
		codeStrings=bases;
	}
	
	public String toString(){
		return name+", "+symbol+", "+letter+", "+Arrays.toString(codeStrings);
	}
	
	public static String kmerToString(long kmer, int k){
		StringBuilder sb=new StringBuilder(k);
		for(int i=0; i<k; i++){
			int x=(int)(kmer&3);
			sb.append((char)numberToBase[x]);
			kmer>>=2;
		}
		return sb.reverse().toString();
	}
	
	
	public final String name;
	public final String symbol;
	public final char letter;
	public final String[] codeStrings;
	
	
	//a=1
	//c=2
	//g=4
	//t=8
	
//	R 	G A (puRine)
//	Y 	T C (pYrimidine)
//	K 	G T (Ketone)
//	M 	A C (aMino group)
//	S 	G C (Strong interaction)
//	W 	A T (Weak interaction)
//	B 	G T C (not A) (B comes after A)
//	D 	G A T (not C) (D comes after C)
//	H 	A C T (not G) (H comes after G)
//	V 	G C A (not T, not U) (V comes after U)
//	N 	A G C T (aNy)
//	X 	masked
//	- 	gap of indeterminate length
	
	public static final byte[] numberToBase={
		'A','C','G','T','N'
	};
	
	public static final byte[] numberToComplementaryBase={
		'T','G','C','A','N'
	};
	
	public static final byte[] numberToComplement={
		3,2,1,0,4
	};
	
	public static final byte[] numberToBaseExtended={
		' ','A','C','M','G','R','S','V', //0-7
		'T','W','Y','H','K','D','B','N', //8-15
		'X',' ',' ',' ',' ',' ',' ',' ', //16-23
	};
	
	/** Has 'N' in position 0.  Mainly for translating compressed arrays containing zeroes to bases. */
	public static final byte[] numberToBaseExtended2={
		'N','A','C','M','G','R','S','V', //0-7
		'T','W','Y','H','K','D','B','N', //8-15
		'X',' ',' ',' ',' ',' ',' ',' ', //16-23
	};
	
	public static final byte[] degenerateBases={
		' ',' ',' ','M',' ','R','S','V', //0-7
		' ','W','Y','H','K','D','B',' ', //8-15
		' ',' ',' ',' ',' ',' ',' ',' ', //16-23
	};
	
	public static final byte[] numberToComplementaryBaseExtended={
		' ','T','G','K','C','Y','W','B', //0-7
		'A','S','R','D','M','H','V','N', //8-15
		'X',' ',' ',' ',' ',' ',' ',' ', //16-23
	};
	
	/** Element i is: 0 for 'A', 1 for 'C', 2 for 'G', 3 for 'T', -1 otherwise */
	public static final byte[] baseToNumber=new byte[128];
	
	/** Element i is: 3 for 'A', 2 for 'C', 1 for 'G', 0 for 'T', -1 otherwise */
	public static final byte[] baseToComplementNumber=new byte[128];
	
	/** Element i is: 0 for 'A', 1 for 'C', 2 for 'G', 3 for 'T', 4 for 'N', -1 otherwise */
	public static final byte[] baseToNumberACGTN=new byte[128];
	
	/** Element i is: 0 for 'A', 1 for 'C', 2 for 'G', 3 for 'T', 0 for 'N', -1 otherwise */
	public static final byte[] baseToNumberACGTN2=new byte[128];
	
	/** Element i is: 0 for 'A', 1 for 'C', 2 for 'G', 3 for 'T', 4 otherwise */
	public static final byte[] baseToNumberACGTother=new byte[128];
	
	/** A>A, C>C, G>G, T/U>T, other>N */
	public static final byte[] baseToACGTN=new byte[128];
	
	public static final byte[] baseToComplementExtended=new byte[128];
	
	/** Uracil to Thymine, everything else unchanged */
	public static final byte[] uToT=new byte[256];
	
	/** Element i is the bitwise OR of constituent IUPAC base numbers in baseToNumber.<br> 
	 * For example, baseToNumberExtended['M'] = (baseToNumber['A'] | baseToNumber['C']) = (1 | 2) = 3 <br>
	 * Invalid characters are -1 */
	public static final byte[] baseToNumberExtended=new byte[128];
	public static final AminoAcid[] AlphabeticalAAs=new AminoAcid[21];
	public static final AminoAcid[] codeToAA=new AminoAcid[66];
	public static final char[] codeToChar=new char[66];
	public static final byte[] codeToByte=new byte[66];
	public static final byte[] aminoToCode=new byte[128];
	public static final HashMap<String, AminoAcid> stringToAA=new HashMap<String, AminoAcid>(512);
	
	public static final AminoAcid Alanine=new AminoAcid("Alanine, Ala, A, GCU, GCC, GCA, GCG");
	public static final AminoAcid Arginine=new AminoAcid("Arginine, Arg, R, CGU, CGC, CGA, CGG, AGA, AGG");
	public static final AminoAcid Asparagine=new AminoAcid("Asparagine, Asn, N, AAU, AAC");
	public static final AminoAcid AsparticAcid=new AminoAcid("AsparticAcid, Asp, D, GAU, GAC");
	public static final AminoAcid Cysteine=new AminoAcid("Cysteine, Cys, C, UGU, UGC");
	public static final AminoAcid GlutamicAcid=new AminoAcid("GlutamicAcid, Glu, E, GAA, GAG");
	public static final AminoAcid Glutamine=new AminoAcid("Glutamine, Gln, Q, CAA, CAG");
	public static final AminoAcid Glycine=new AminoAcid("Glycine, Gly, G, GGU, GGC, GGA, GGG");
	public static final AminoAcid Histidine=new AminoAcid("Histidine, His, H, CAU, CAC");
	public static final AminoAcid Isoleucine=new AminoAcid("Isoleucine, Ile, I, AUU, AUC, AUA");
	public static final AminoAcid Leucine=new AminoAcid("Leucine, Leu, L, UUA, UUG, CUU, CUC, CUA, CUG");
	public static final AminoAcid Lysine=new AminoAcid("Lysine, Lys, K, AAA, AAG");
	public static final AminoAcid Methionine=new AminoAcid("Methionine, Met, M, AUG");
	public static final AminoAcid Phenylalanine=new AminoAcid("Phenylalanine, Phe, F, UUU, UUC");
	public static final AminoAcid Proline=new AminoAcid("Proline, Pro, P, CCU, CCC, CCA, CCG");
	public static final AminoAcid Serine=new AminoAcid("Serine, Ser, S, UCU, UCC, UCA, UCG, AGU, AGC");
	public static final AminoAcid Threonine=new AminoAcid("Threonine, Thr, T, ACU, ACC, ACA, ACG");
	public static final AminoAcid Tryptophan=new AminoAcid("Tryptophan, Trp, W, UGG");
	public static final AminoAcid Tyrosine=new AminoAcid("Tyrosine, Tyr, Y, UAU, UAC");
	public static final AminoAcid Valine=new AminoAcid("Valine, Val, V, GUU, GUC, GUA, GUG");
	public static final AminoAcid END=new AminoAcid("End, End, *, UAA, UGA, UAG");
	public static final AminoAcid ANY=new AminoAcid("Any, Any, X, XXX");
	
	


	public static final byte[][] COLORS=new byte[][] {
		{0, 1, 2, 3},
		{1, 0, 3, 2},
		{2, 3, 0, 1},
		{3, 2, 1, 0}
	};
	
	/** Returns a new reverse-complemented array in ASCII coding*/
	public static final byte[] reverseComplementBases(final byte[] in){
		byte[] out=new byte[in.length];
		final int last=in.length-1;
		for(int i=0; i<in.length; i++){
			out[i]=baseToComplementExtended[in[last-i]];
		}
		return out;
	}
	

	public static final void reverseComplementBasesInPlace(final byte[] in){
		if(in!=null){reverseComplementBasesInPlace(in, in.length);}
	}
	
	public static final void reverseComplementBasesInPlace(final byte[] in, final int length){
		if(in==null){return;}
		final int last=length-1;
		final int max=length/2;
		for(int i=0; i<max; i++){
			byte a=in[i];
			byte b=in[last-i];
//			assert(b>0 && b<baseToComplementExtended.length) : ((int)b)+"\t"+((char)b)+"\t"+Arrays.toString(in);
//			System.out.println((char)a+", "+(char)b+", "+i+", "+last);
			in[i]=baseToComplementExtended[b];
			in[last-i]=baseToComplementExtended[a];
		}
		if((length&1)==1){//Odd length; process middle
			in[max]=baseToComplementExtended[in[max]];
		}
	}
	
	public static final String reverseComplementBases(String in){
		return in==null ? null : new String(reverseComplementBases(in.getBytes()));
	}
	
	public static final int reverseComplementBinary(int kmer, int k){
		int out=0;
		kmer=~kmer;
		for(int i=0; i<k; i++){
			out=((out<<2)|(kmer&3));
			kmer>>=2;
		}
		return out;
	}
	
	public static final long reverseComplementBinary(long kmer, int k){
		long out=0;
		kmer=~kmer;
		for(int i=0; i<k; i++){
			out=((out<<2)|(kmer&3L));
			kmer>>=2;
		}
		return out;
	}
	
	public static final int reverseComplementBinaryFast(int kmer, int k){
		int out=0;
		int extra=k&3;
		for(int i=0; i<extra; i++){
			out=((out<<2)|((~kmer)&3));
			kmer>>=2;
		}
		k-=extra;
		for(int i=0; i<k; i+=4){
			out=((out<<8)|(rcompBinaryTable[kmer&0xFF]));
			kmer>>=8;
		}
		return out;
	}
	
	public static final long reverseComplementBinaryFast(long kmer, int k){
		long out=0;
		int extra=k&3;
		for(int i=0; i<extra; i++){
			out=((out<<2)|((~kmer)&3L));
			kmer>>=2;
		}
		k-=extra;
		for(int i=0; i<k; i+=4){
			out=((out<<8)|(rcompBinaryTable[(int)(kmer&0xFFL)]));
			kmer>>=8;
		}
		return out;
	}
	
	public static final byte baseToColor(byte base1, byte base2){
		byte a=baseToNumber[base1];
		byte b=baseToNumber[base2];
		if(a<0 && b<0){return 'N';}
		if(a<0){a=3;}
		if(b<0){b=3;}
		return COLORS[a][b];
	}
	
	public static final byte colorToBase(byte base1, byte color){
		if(!isFullyDefined(base1) || color<0 || color>3){
			return (byte)'N';
		}
		byte a=baseToNumber[base1];
		
		return numberToBase[COLORS[a][color]];
	}
	
	public static final byte toNumber(String code){
		return toNumber(code.charAt(0), code.charAt(1), code.charAt(2));
	}
	
	public static final AminoAcid toAA(String code){
		return toAA(code.charAt(0), code.charAt(1), code.charAt(2));
	}
	
	public static final char toChar(String code){
		return toChar(code.charAt(0), code.charAt(1), code.charAt(2));
	}
	
	public static final char[] splitBase(char c){
		byte b=baseToNumberExtended[c];
		int len=Integer.bitCount(b);
		char[] out=new char[len];
		
		int index=0;
		for(int i=0; i<4; i++){
			if(((1<<i)&b)!=0){
				out[index]=(char)numberToBase[i];
				index++;
			}
		}
		return out;
	}
	

	
	
	public static final byte[] numberToBases(int code, int n){
		
		byte[] bytes=new byte[n];
		
		for(int i=n-1; i>=0; i--){
			int temp=code&3;
			code>>=2;
			bytes[i]=numberToBase[temp];
		}
		
		return bytes;
	}
	
	public static final int baseTupleToNumber(byte[] tuple){
		
		int r=0;
		for(int i=0; i<tuple.length; i++){
			int temp=baseToNumberACGTN[tuple[i]];
			if(temp<0 || temp>3){return -1;}
			r=((r<<2)|temp);
		}
		
		return r;
	}
	
	public static boolean isFullyDefined(char base){
		return baseToNumber[base]>=0;
	}
	
	public static boolean isFullyDefined(byte base){
		return base>=0 && baseToNumber[base]>=0;
	}
	
	public static boolean isACGTN(char base){
		return baseToNumberACGTN[base]>=0;
	}
	
	public static boolean isACGTN(byte base){
		return base>=0 && baseToNumberACGTN[base]>=0;
	}
	
	public static boolean containsOnlyACGTN(String s){
		if(s==null || s.length()==0){return true;}
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(baseToNumberACGTN[c]<0){return false;}
		}
		return true;
	}
	
	public static boolean containsOnlyACGTNQ(String s){
		if(s==null || s.length()==0){return true;}
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(c!='?' && baseToNumberACGTN[c]<0){return false;}
		}
		return true;
	}
	
	public static boolean containsOnlyACGTN(byte[] array){
		if(array==null || array.length==0){return true;}
		for(int i=0; i<array.length; i++){
			byte b=array[i];
			if(b<0 || baseToNumberACGTN[b]<0){return false;}
		}
		return true;
	}
	
	public static boolean isFullyDefined(String s){
		for(int i=0; i<s.length(); i++){
			if(!isFullyDefined(s.charAt(i))){return false;}
		}
		return true;
	}
	
	public static boolean isFullyDefined(byte[] s){
		for(int i=0; i<s.length; i++){
			if(!isFullyDefined(s[i])){return false;}
		}
		return true;
	}
	
	public static int countUndefined(byte[] s){
		int x=0;
		for(int i=0; i<s.length; i++){
			if(!isFullyDefined(s[i])){x++;}
		}
		return x;
	}
	
	
	public static final byte toNumber(char c1, char c2, char c3){
		assert(baseToNumberACGTN2[c1]>=0 && baseToNumberACGTN2[c2]>=0 && baseToNumberACGTN2[c3]>=0);
		int x=(baseToNumberACGTN2[c1]<<4)|(baseToNumberACGTN2[c2]<<2)|(baseToNumberACGTN2[c3]);
		return (byte)x;
	}
	
	public static final AminoAcid toAA(char c1, char c2, char c3){
		assert(baseToNumberACGTN2[c1]>=0 && baseToNumberACGTN2[c2]>=0 && baseToNumberACGTN2[c3]>=0);
		int x=(baseToNumberACGTN2[c1]<<4)|(baseToNumberACGTN2[c2]<<2)|(baseToNumberACGTN2[c3]);
		return codeToAA[x];
	}
	
	public static final char toChar(char c1, char c2, char c3){
		assert(baseToNumberACGTN2[c1]>=0 && baseToNumberACGTN2[c2]>=0 && baseToNumberACGTN2[c3]>=0);
		int x=(baseToNumberACGTN2[c1]<<4)|(baseToNumberACGTN2[c2]<<2)|(baseToNumberACGTN2[c3]);
		return codeToChar[x];
	}
	
	public static final byte toByte(byte c1, byte c2, byte c3){
		int a=baseToNumber[c1], b=baseToNumber[c2], c=baseToNumber[c3];
		if(a<0 || b<0 || c<0){return (byte)'.';}
		int x=((a<<4)|(b<<2)|c);
		return codeToByte[x];
	}
	
	public static final char toChar(byte c1, byte c2, byte c3){
		assert(baseToNumberACGTN2[c1]>=0 && baseToNumberACGTN2[c2]>=0 && baseToNumberACGTN2[c3]>=0);
		byte n1=baseToNumberACGTN2[c1], n2=baseToNumberACGTN2[c2], n3=baseToNumberACGTN2[c3];
		if(n1>3 || n2>3 || n3>3){return '?';}
		int x=(n1<<4)|(n2<<2)|(n3);
//		return (x<codeToChar.length ? codeToChar[x] : '?');
		return codeToChar[x];
	}
	
	public static final String stringToAAs(String bases){
		StringBuilder sb=new StringBuilder(bases.length()/3);
		for(int i=2; i<bases.length(); i+=3){
			char a=toAA(bases.charAt(i-2), bases.charAt(i-1), bases.charAt(i)).letter;
			sb.append(a);
		}
		return sb.toString();
	}
	
	public static final byte[][] toAAsSixFrames(byte[] bases){
		byte[][] out=new byte[6][];
		if(bases!=null && bases.length>2){
			for(int i=0; i<3; i++){
				out[i]=toAAs(bases, i);
			}
			byte[] rcomp=reverseComplementBases(bases);
			for(int i=0; i<3; i++){
				out[i+3]=toAAs(rcomp, i);
			}
		}
		return out;
	}
	
	public static final byte[][] toQualitySixFrames(byte[] quals, int offset){
		byte[][] out=new byte[6][];
		if(quals!=null && quals.length>2){
			for(int i=0; i<3; i++){
				out[i]=toAAQuality(quals, i);
			}
			Tools.reverseInPlace(quals);
			for(int i=0; i<3; i++){
				out[i+3]=toAAQuality(quals, i);
			}
			Tools.reverseInPlace(quals);
		}
		
		if(offset!=0){
			for(byte[] array : out){
				if(array!=null){
					for(int i=0; i<array.length; i++){
						array[i]+=offset;
					}
				}
			}
		}
		
		return out;
	}
	
	public static final byte[] toAAs(byte[] bases, int frame){
		assert(frame>=0 && frame<3);
		if(bases==null){return null;}
		int blen=bases.length-frame;
		if(blen<3){return null;}
		blen=blen-(blen%3);
		final int stop=frame+blen;
		final int alen=blen/3;
		
		byte[] out=new byte[alen];
		for(int i=2+frame, j=0; i<stop; i+=3, j++){
			byte a=toByte(bases[i-2], bases[i-1], bases[i]);
			out[j]=a;
		}
		return out;
	}
	
	public static final byte[] toAAQuality(byte[] quals, int frame){
		assert(frame>=0 && frame<3);
		int blen=quals.length-frame;
		if(blen<3){return null;}
		blen=blen-(blen%3);
		final int stop=frame+blen;
		final int alen=blen/3;
		
		byte[] out=new byte[alen];
		for(int i=2+frame, j=0; i<stop; i+=3, j++){
			byte qa=quals[i-2], qb=quals[i-1], qc=quals[i];
			float pa=QualityTools.PROB_CORRECT[qa], pb=QualityTools.PROB_CORRECT[qb], pc=QualityTools.PROB_CORRECT[qc];
			float p=pa*pb*pc;
			byte q=QualityTools.probCorrectToPhred(p);
			out[j]=q;

//			System.out.println();
//			System.out.println(qa+", "+qb+", "+qc+" -> "+q);
//			System.out.println(pa+", "+pb+", "+pc+" -> "+p);
			
		}
//		System.out.println(Arrays.toString(out));
		return out;
	}
	
	public static final byte[] toNTs(final byte[] aminos){
		if(aminos==null){return null;}
		final int alen=aminos.length;
		final int blen=alen*3;
		
		final byte[] out=new byte[blen];
		for(int i=0, j=0; i<alen; i++, j+=3){
			int code=aminoToCode[aminos[i]];
			out[j+2]=numberToBase[(code&3)];
			out[j+1]=numberToBase[((code>>2)&3)];
			out[j]=numberToBase[((code>>4)&3)];
		}
		return out;
	}
	
	public static final short[] rcompBinaryTable=makeBinaryRcompTable(4);
	
	private static final short[] makeBinaryRcompTable(int k){
		int bits=2*k;
		short[] r=new short[1<<bits];
		for(int i=0; i<r.length; i++){
			r[i]=(short)reverseComplementBinary(i, k);
		}
		return r;
	}
	
	static {
		
		for(int i=0; i<uToT.length; i++){uToT[i]=(byte)i;}
		uToT['u']='t';
		uToT['U']='T';
		
		Arrays.fill(baseToACGTN, (byte)'N');
		
		Arrays.fill(baseToNumberExtended, (byte)-1);
		for(int i=0; i<numberToBaseExtended.length; i++){
			char x=(char)numberToBaseExtended[i];
			if(!Character.isWhitespace(x)){
				baseToNumberExtended[x]=(byte)i;
				baseToNumberExtended[Character.toLowerCase(x)]=(byte)i;
			}
		}
		baseToNumberExtended['U']=8;
		baseToNumberExtended['u']=8;
		
		Arrays.fill(baseToNumberACGTN, (byte)-1);
		Arrays.fill(baseToNumberACGTother, (byte)4);
		for(int i=0; i<numberToBase.length; i++){
			char x=(char)numberToBase[i];
			if(!Character.isWhitespace(x)){
				baseToNumberACGTN[x]=baseToNumberACGTother[x]=(byte)i;
				baseToNumberACGTN[Character.toLowerCase(x)]=baseToNumberACGTother[Character.toLowerCase(x)]=(byte)i;
				baseToACGTN[x]=baseToACGTN[Character.toLowerCase(x)]=(byte)x;
			}
		}
		baseToNumberACGTN['U']=baseToNumberACGTN['u']=3;
		baseToNumberACGTother['U']=baseToNumberACGTother['u']=3;
		baseToACGTN['U']=baseToACGTN['u']=(byte)'T';
		
		for(int i=0; i<baseToNumberACGTN.length; i++){baseToNumberACGTN2[i]=baseToNumberACGTN[i];}
		baseToNumberACGTN2['N']=0;
		baseToNumberACGTN2['n']=0;
		
		Arrays.fill(baseToNumber, (byte)-1);
		for(int i=0; i<numberToBase.length; i++){
			char x=(char)numberToBase[i];
			if(x=='A' || x=='C' || x=='G' || x=='T'){
				baseToNumber[x]=(byte)i;
				baseToNumber[Character.toLowerCase(x)]=(byte)i;
			}
		}
		baseToNumber['U']=3;
		baseToNumber['u']=3;
		
		Arrays.fill(baseToComplementNumber, (byte)-1);
		baseToComplementNumber['A']=baseToComplementNumber['a']=3;
		baseToComplementNumber['C']=baseToComplementNumber['c']=2;
		baseToComplementNumber['G']=baseToComplementNumber['g']=1;
		baseToComplementNumber['T']=baseToComplementNumber['t']=0;
		baseToComplementNumber['U']=baseToComplementNumber['u']=0;
		
		Arrays.fill(baseToComplementExtended, (byte)-1);
		for(int i=0; i<numberToBaseExtended.length; i++){
			char x=(char)numberToBaseExtended[i];
			char x2=(char)numberToComplementaryBaseExtended[i];
			baseToComplementExtended[x]=(byte)x2;
			baseToComplementExtended[Character.toLowerCase(x)]=(byte)Character.toLowerCase(x2);
		}
		baseToComplementExtended['U']=(byte)'A';
		baseToComplementExtended['u']=(byte)'a';
		baseToComplementExtended['?']=(byte)'?';
		baseToComplementExtended[' ']=(byte)' ';
		baseToComplementExtended['-']=(byte)'-';
		baseToComplementExtended['*']=(byte)'*';
		baseToComplementExtended['.']=(byte)'.';
		
		
		AlphabeticalAAs[0]=Alanine;
		AlphabeticalAAs[1]=Arginine;
		AlphabeticalAAs[2]=Asparagine;
		AlphabeticalAAs[3]=AsparticAcid;
		AlphabeticalAAs[4]=Cysteine;
		AlphabeticalAAs[5]=GlutamicAcid;
		AlphabeticalAAs[6]=Glutamine;
		AlphabeticalAAs[7]=Glycine;
		AlphabeticalAAs[8]=Histidine;
		AlphabeticalAAs[9]=Isoleucine;
		AlphabeticalAAs[10]=Leucine;
		AlphabeticalAAs[11]=Lysine;
		AlphabeticalAAs[12]=Methionine;
		AlphabeticalAAs[13]=Phenylalanine;
		AlphabeticalAAs[14]=Proline;
		AlphabeticalAAs[15]=Serine;
		AlphabeticalAAs[16]=Threonine;
		AlphabeticalAAs[17]=Tryptophan;
		AlphabeticalAAs[18]=Tyrosine;
		AlphabeticalAAs[19]=Valine;
		AlphabeticalAAs[20]=END;
//		AlphabeticalAAs[21]=ANY;
		
		Arrays.fill(aminoToCode, (byte)-1);
		for(int i=0; i<AlphabeticalAAs.length; i++){
			AminoAcid aa=AlphabeticalAAs[i];
			
			stringToAA.put(aa.name, aa);
			stringToAA.put(aa.symbol, aa);
			stringToAA.put(aa.letter+"", aa);
			for(int j=0; j<aa.codeStrings.length; j++){
				String s=aa.codeStrings[j];
				stringToAA.put(s, aa);
				aa.codeStrings[j]=s.replace('U', 'T');
				stringToAA.put(aa.codeStrings[j], aa);
				
				int x=toNumber(s);
//				System.out.println("x="+x+", aa="+aa);
				codeToAA[x]=aa;
				codeToChar[x]=aa.letter;
				codeToByte[x]=(byte)(aa.letter);
				if(j==0){
					aminoToCode[aa.letter]=(byte)x;
					aminoToCode[Character.toLowerCase(aa.letter)]=(byte)x;
				}
			}
		}
		aminoToCode['X']=aminoToCode['x']=65;
		codeToAA[65]=ANY;
		codeToChar[65]='X';
		codeToByte[65]='X';
		
		stringToAA.put("X", ANY);
		stringToAA.put("Start", Methionine);
		stringToAA.put("Begin", Methionine);
		stringToAA.put("Stop", END);
		stringToAA.put("Aspartic Acid", AsparticAcid);
		stringToAA.put("Glutamic Acid", GlutamicAcid);
		
		String[] temp=stringToAA.keySet().toArray(new String[0]);
		
		for(String s : temp){
			AminoAcid aa=stringToAA.get(s);
			assert(aa!=null);
			stringToAA.put(s.toLowerCase(), aa);
		}
		
	}
	
}

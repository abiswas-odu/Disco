package shared;

import java.io.Serializable;
import java.util.Arrays;

import align2.QualityTools;
import stream.ByteBuilder;
import stream.KillSwitch;
import stream.Read;
import stream.SamLine;
import stream.SiteScore;

import dna.AminoAcid;
import dna.Gene;

/**
 * Helper class for processes that do inline quality trimming.
 * @author Brian Bushnell
 * @date Mar 15, 2013
 *
 */
public final class TrimRead implements Serializable {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 8791743639124592480L;
	
	public static void main(String[] args){
		byte[] bases=args[0].getBytes();
		byte[] quals=(args.length<2 ? null : args[1].getBytes());
		if(quals!=null){
			for(int i=0; i<quals.length; i++){quals[i]-=32;}
		}
		byte[] match=(args.length<3 ? null : args[2].getBytes());
		int minq=(args.length<4 ? 5 : Integer.parseInt(args[3]));
		Read r=new Read(bases, quals, 1);
		r.match=match;
		System.out.println("Before trim:\n"+r.toFastq()+(r.match==null ? "" : "\n"+new String(r.match)));
		System.out.println(Arrays.toString(r.quality));
		TrimRead tr=trim(r, true, true, minq, 1);
		System.out.println("\nAfter trim:\n"+r.toFastq()+(r.match==null ? "" : "\n"+new String(r.match)));
		if(r.match==null){
			r.match=new byte[r.length()];
			for(int i=0; i<r.length(); i++){r.match[i]='m';}
		}
		tr.untrim();
		System.out.println("\nAfter untrim:\n"+r.toFastq()+(r.match==null ? "" : "\n"+new String(r.match)));
	}
	
	public static TrimRead trim(Read r, boolean trimLeft, boolean trimRight, int trimq, int minlen){
		if(r==null || r.bases==null){return null;}
		
		final int a, b;
		if(optimalMode){
			long packed=testOptimal(r.bases, r.quality, QualityTools.PROB_ERROR[trimq]);
			a=trimLeft ? (int)((packed>>32)&0xFFFFFFFFL) : 0;
			b=trimRight ? (int)((packed)&0xFFFFFFFFL) : 0;
		}else if(windowMode){
			a=0;
			b=(trimRight ? testRightWindow(r.bases, r.quality, (byte)trimq, windowLength) : 0);
		}else{
			a=(trimLeft ? testLeft(r.bases, r.quality, (byte)trimq) : 0);
			b=(trimRight ? testRight(r.bases, r.quality, (byte)trimq) : 0);
		}
		return (a+b==0 ? null : new TrimRead(r, a, b, trimq, minlen));
	}
	
	/**
	 * @param r Read to trim
	 * @param trimLeft Trim left side
	 * @param trimRight Trim right side
	 * @param trimq Maximum quality to trim
	 * @param minResult Ensure trimmed read is at least this long
	 * @return Number of bases trimmed
	 */
	public static int trimFast(Read r, boolean trimLeft, boolean trimRight, int trimq, int minResult){
		return trimFast(r, trimLeft, trimRight, trimq, minResult, 0);
	}
	
	/**
	 * This method allows a "discardUnder" parameter, so that qtrim=r will still discard
	 * reads that if left-trimmed, would be shorter than the desired minimum length.
	 * It is not presently used.
	 * @param r Read to trim
	 * @param trimLeft Trim left side
	 * @param trimRight Trim right side
	 * @param trimq Maximum quality to trim
	 * @param minResult Ensure trimmed read is at least this long
	 * @param discardUnder Resulting reads shorter than this are not wanted 
	 * @return Number of bases trimmed
	 */
	public static int trimFast(Read r, boolean trimLeft, boolean trimRight, int trimq, int minResult, int discardUnder){
		final byte[] bases=r.bases, qual=r.quality;
		if(bases==null || bases.length<1){return 0;}
		
//		assert(false) : trimLeft+", "+trimRight+", "+trimq+", "+minResult+", "+discardUnder;
		
		final int len=r.length();
		final int a0, b0;
		final int a, b;
		if(optimalMode){
			long packed=testOptimal(bases, qual, QualityTools.PROB_ERROR[trimq]);
			a0=(int)((packed>>32)&0xFFFFFFFFL);
			b0=(int)((packed)&0xFFFFFFFFL);
			if(trimLeft!=trimRight && discardUnder>0 && len-a0-b0<discardUnder){
				a=trimLeft ? len : 0;
				b=trimRight ? len : 0;
			}else{
				a=trimLeft ? a0 : 0;
				b=trimRight ? b0 : 0;
			}
//			assert(false) : a0+", "+b0+" -> "+a+", "+b;
		}else if(windowMode){
			a=0;
			b=(trimRight ? testRightWindow(bases, qual, (byte)trimq, windowLength) : 0);
		}else{
			a=(trimLeft ? testLeft(bases, qual, (byte)trimq) : 0);
			b=(trimRight ? testRight(bases, qual, (byte)trimq) : 0);
		}
		return trimByAmount(r, a, b, minResult);
	}
	
	public static boolean untrim(Read r){
		if(r==null || r.obj==null){return false;}
		if(r.obj.getClass()!=TrimRead.class){return false;}
		TrimRead tr=(TrimRead)r.obj;
		return tr.untrim();
	}
	
//	public TrimRead(Read r_, boolean trimLeft, boolean trimRight, int trimq_, int minlen_){
//		this(r_, (trimLeft ? testLeft(r_.bases, r_.quality, (byte)trimq_) : 0), (trimRight ? testRight(r_.bases, r_.quality, (byte)trimq_) : 0), trimq_, minlen_);
//	}
	
	public TrimRead(Read r_, int trimLeft, int trimRight, int trimq_, int minlen_){
		minlen_=Tools.max(minlen_, 0);
		r=r_;
		bases1=r.bases;
		qual1=r.quality;
		trimq=(byte)trimq_;
		assert(bases1!=null || qual1==null) : "\n"+new String(bases1)+"\n"+new String(qual1)+"\n";
		assert(bases1==null || qual1==null || bases1.length==qual1.length) : "\n"+new String(bases1)+"\n"+new String(qual1)+"\n";
		int trimmed=trim(trimLeft, trimRight, minlen_);
		if(trimmed>0){
			assert(bases2==null || bases2.length>=minlen_ || bases1.length<minlen_) : bases1.length+", "+bases2.length+", "+minlen_+", "+trimLeft+", "+trimRight;
			r.bases=bases2;
			r.quality=qual2;
			r.obj=this;
			trimMatch(r);
		}
	}
	
//	/** Trim the left end of the read, from left to right */
//	private int trim(final boolean trimLeft, final boolean trimRight, final int minlen){
//		final int a, b;
//		if(optimalMode){
//			long packed=testOptimal(bases1, qual1, QualityTools.PROB_ERROR[trimq]);
//			a=trimLeft ? (int)((packed>>32)&0xFFFFFFFFL) : 0;
//			b=trimRight ? (int)((packed)&0xFFFFFFFFL) : 0;
//		}else{
//			a=(trimLeft ? testLeft(bases1, qual1, (byte)trimq) : 0);
//			b=(trimRight ? testRight(bases1, qual1, (byte)trimq) : 0);
//		}
//		return trim(a, b, minlen);
//	}
	
	/** Trim the left end of the read, from left to right */
	private int trim(int trimLeft, int trimRight, final int minlen){
		assert(trimLeft>=0 && trimRight>=0) : "trimLeft="+trimLeft+", trimRight="+trimRight+", minlen="+minlen+", len="+bases1.length;
		assert(trimLeft>0 || trimRight>0) : "trimLeft="+trimLeft+", trimRight="+trimRight+", minlen="+minlen+", len="+bases1.length;
		final int maxTrim=Tools.min(bases1.length, bases1.length-minlen);
		if(trimLeft+trimRight>maxTrim){
			int excess=trimLeft+trimRight-maxTrim;
			if(trimLeft>0 && excess>0){
				trimLeft=Tools.max(0, trimLeft-excess);
				excess=trimLeft+trimRight-maxTrim;
			}
			if(trimRight>0 && excess>0){
				trimRight=Tools.max(0, trimRight-excess);
				excess=trimLeft+trimRight-maxTrim;
			}
			
		}
		
		leftTrimmed=trimLeft;
		rightTrimmed=trimRight;
		final int sum=leftTrimmed+rightTrimmed;
		
		if(verbose){
			System.err.println("leftTrimmed="+leftTrimmed+", rightTrimmed="+rightTrimmed+", sum="+sum);
		}
		
		if(sum==0){
			bases2=bases1;
			qual2=qual1;
		}else{
			bases2=KillSwitch.copyOfRange(bases1, trimLeft, bases1.length-trimRight);
			qual2=((qual1==null || (trimLeft+trimRight>=qual1.length)) ? null : KillSwitch.copyOfRange(qual1, trimLeft, qual1.length-trimRight));
		}
		return sum;
	}
	
	/** Trim bases outside of leftLoc and rightLoc, excluding leftLoc and rightLoc */
	public static int trimToPosition(Read r, int leftLoc, int rightLoc, int minResultingLength){
		final int len=r.length();
		return trimByAmount(r, leftLoc, len-rightLoc-1, minResultingLength);
	}
	
	/** Remove non-genetic-code from reads */
	public static int trimBadSequence(Read r){
		final byte[] bases=r.bases, quals=r.quality;
		if(bases==null){return 0;}
		final int minGenetic=20;
		int lastNon=-1;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			if(!AminoAcid.isACGTN(b)){lastNon=i;}
			if(i-lastNon>minGenetic){break;}
		}
		if(lastNon>=0){
			r.bases=KillSwitch.copyOfRange(bases, lastNon+1, bases.length);
			if(quals!=null){
				r.quality=KillSwitch.copyOfRange(quals, lastNon+1, quals.length);
			}
		}
		return lastNon+1;
	}
	
	/** Trim this many bases from each end */
	public static int trimByAmount(Read r, int leftTrimAmount, int rightTrimAmount, int minResultingLength){

		leftTrimAmount=Tools.max(leftTrimAmount, 0);
		rightTrimAmount=Tools.max(rightTrimAmount, 0);
		
		//These assertions are unnecessary if the mapping information will never be used or output.
//		assert(r.match==null) : "TODO: Handle trimming of reads with match strings.";
		assert(r.sites==null) : "TODO: Handle trimming of reads with SiteScores.";
		
		if(r.match!=null){
			return trimReadWithMatch(r, (SamLine)r.obj, leftTrimAmount, rightTrimAmount, minResultingLength, Integer.MAX_VALUE);
		}
		
		final byte[] bases=r.bases, qual=r.quality;
		final int len=(bases==null ? 0 : bases.length), qlen=(qual==null ? 0 : qual.length);
		if(len<1){return 0;}
		minResultingLength=Tools.min(len, Tools.max(minResultingLength, 0));
		if(leftTrimAmount+rightTrimAmount+minResultingLength>len){
			rightTrimAmount=Tools.max(1, len-minResultingLength);
			leftTrimAmount=0;
		}
		
		final int total=leftTrimAmount+rightTrimAmount;
		if(total>0){
			r.bases=KillSwitch.copyOfRange(bases, leftTrimAmount, len-rightTrimAmount);
			r.quality=(leftTrimAmount+rightTrimAmount>=qlen ? null : KillSwitch.copyOfRange(qual, leftTrimAmount, qlen-rightTrimAmount));
			trimMatch(r);
			if(r.stop>r.start){ //TODO:  Fixing mapped coordinates needs more work.
				r.start+=leftTrimAmount;
				r.stop-=rightTrimAmount;
			}
		}
		
		if(verbose){
			System.err.println("leftTrimmed="+leftTrimAmount+", rightTrimmed="+rightTrimAmount+
					", sum="+total+", final length="+r.length());
		}
		
		return total;
	}
	
	/** Count number of bases that need trimming on each side, and pack into a long */
	private static long testOptimal(byte[] bases, byte[] qual, float avgErrorRate){
		if(optimalBias>=0){avgErrorRate=optimalBias;}//Override
		assert(avgErrorRate>0 && avgErrorRate<=1) : "Average error rate ("+avgErrorRate+") must be between 0 (exclusive) and 1 (inclusive)";
		if(bases==null || bases.length==0){return 0;}
		if(qual==null){return avgErrorRate>=1 ? 0 : ((((long)testLeftN(bases))<<32) | (((long)testRightN(bases))&0xFFFFFFFFL));}
		
		float maxScore=0;
		float score=0;
		int maxLoc=-1;
		int maxCount=-1;
		int count=0;
		
		final float nprob=Tools.max(Tools.min(avgErrorRate*1.1f, 1), NPROB);
		
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			byte q=qual[i];
//			float probError=(b=='N' ? nprob : ADJUST_QUALITY ? CalcTrueQuality.estimateErrorProbAvg(qual, bases, i) : QualityTools.PROB_ERROR[q]);
//			float probError=(b=='N' ? nprob : ADJUST_QUALITY ? CalcTrueQuality.estimateErrorProbGeoAvg(qual, bases, i) : QualityTools.PROB_ERROR[q]);
//			float probError=(b=='N' ? nprob : ADJUST_QUALITY ? CalcTrueQuality.estimateErrorProb2(qual, bases, i) : QualityTools.PROB_ERROR[q]);

//			float probError=(b=='N' ? nprob : q==1 ? PROB1 : QualityTools.PROB_ERROR[q]);
			float probError=(b=='N' ? nprob : QualityTools.PROB_ERROR[q]);
			
//			assert(q>0 || b=='N') : "index "+i+": q="+q+", b="+(char)b+"\n"+new String(bases)+"\n"+Arrays.toString(qual)+"\n";
			
			float delta=avgErrorRate-probError;
			score=score+delta;
			if(score>0){
				count++;
				if(score>maxScore || (score==maxScore && count>maxCount)){
					maxScore=score;
					maxCount=count;
					maxLoc=i;
				}
			}else{
				score=0;
				count=0;
			}
		}
		
		final int left, right;
		if(maxScore>0){
			assert(maxLoc>=0);
			assert(maxCount>0);
			left=maxLoc-maxCount+1;
			assert(left>=0 && left<=bases.length);
			right=bases.length-maxLoc-1;
		}else{
			left=0;
			right=bases.length;
		}
		final long packed=((((long)left)<<32) | (((long)right)&0xFFFFFFFFL));
		
		if(verbose){
			System.err.println(Arrays.toString(qual));
			System.err.println("After testLocal: maxScore="+maxScore+", maxLoc="+maxLoc+", maxCount="+maxCount+
					", left="+left+", right="+right+", returning "+Long.toHexString(packed));
		}
		return packed;
	}
	
	/** Count number of bases that need trimming on left side */
	private static int testLeft(byte[] bases, byte[] qual, final byte trimq){
		if(bases==null || bases.length==0){return 0;}
		if(qual==null){return trimq<0 ? 0 : testLeftN(bases);}
		int good=0;
		int lastBad=-1;
		int i=0;
		for(; i<bases.length && good<minGoodInterval; i++){
			final byte q=qual[i];
			final byte b=bases[i];
			assert(q>0 || b=='N') : "index "+i+": q="+q+", b="+(char)b+"\n"+new String(bases)+"\n"+Arrays.toString(qual)+"\n";
			if(q>trimq){good++;}
			else{good=0; lastBad=i;}
		}
		if(verbose){
//			System.err.println(Arrays.toString(qual));
			System.err.println("After testLeft: good="+good+", lastBad="+lastBad+", i="+i+", returning "+(lastBad+1));
//			assert(false);
		}
		return lastBad+1;
	}
	
	/** Count number of bases that need trimming on right side using a sliding window */
	private static int testRightWindow(byte[] bases, byte[] qual, final byte trimq, final int window){
		if(bases==null || bases.length==0){return 0;}
		if(qual==null || qual.length<window){return trimq>0 ? 0 : testRightN(bases);}
		final int thresh=Tools.max(window*trimq, 1);
		int sum=0;
		for(int i=0, j=-window; i<qual.length; i++, j++){
			final byte q=qual[i];
			sum+=q;
			if(j>=-1){
				if(j>=0){sum-=qual[j];}
				if(sum<thresh){
					return qual.length-j-1;
				}
			}
		}
		return 0;
	}
	
	/** Count number of bases that need trimming on right side */
	private static int testRight(byte[] bases, byte[] qual, final byte trimq){
		if(bases==null || bases.length==0){return 0;}
		if(qual==null){return trimq<0 ? 0 : testRightN(bases);}
		int good=0;
		int lastBad=bases.length;
		int i=bases.length-1;
		for(; i>=0 && good<minGoodInterval; i--){
			final byte q=qual[i];
			final byte b=bases[i];
			assert(q>0 || b=='N') : "index "+i+": q="+q+", b="+(char)b+"\n"+new String(bases)+"\n"+Arrays.toString(qual)+"\n";
			if(q>trimq){good++;}
			else{good=0; lastBad=i;}
		}
		if(verbose){
			System.err.println("After trimLeft: good="+good+", lastBad="+lastBad+", i="+i+", returning "+(bases.length-lastBad));
		}
		return bases.length-lastBad;
	}
	
	/** Count number of bases that need trimming on left side, considering only N as bad */
	public static int testLeftN(byte[] bases){
		if(bases==null || bases.length==0){return 0;}
		int good=0;
		int lastBad=-1;
		for(int i=0; i<bases.length && good<minGoodInterval; i++){
			final byte b=bases[i];
			//if(dna.AminoAcid.isFullyDefined(b)){good++;}
			if(b!=((byte)'N')){good++;}
			else{good=0; lastBad=i;}
		}
		return lastBad+1;
	}
	
	/** Count number of bases that need trimming on right side, considering only N as bad */
	public static int testRightN(byte[] bases){
		if(bases==null || bases.length==0){return 0;}
		int good=0;
		int lastBad=bases.length;
		for(int i=bases.length-1; i>=0 && good<minGoodInterval; i--){
			final byte b=bases[i];
			//if(dna.AminoAcid.isFullyDefined(b)){good++;}
			if(b!=((byte)'N')){good++;}
			else{good=0; lastBad=i;}
		}
		return bases.length-lastBad;
	}
	
	public boolean untrim(){
		if(leftTrimmed==0 && rightTrimmed==0){return false;}
		r.setPerfect(false);
		
		final int lt, rt;
		if(r.strand()==Gene.PLUS){
			lt=leftTrimmed;
			rt=rightTrimmed;
		}else{
			lt=rightTrimmed;
			rt=leftTrimmed;
		}
		
		boolean returnToShort=false;
		if(verbose){System.err.println("Untrimming");}
		if(r.match!=null){
			if(r.shortmatch()){
				r.toLongMatchString(false);
				r.setShortMatch(false);
				returnToShort=true;
			}
			byte[] match2=new byte[r.match.length+lt+rt];
			int i=0;
			for(; i<lt; i++){
				match2[i]='C';
			}
			for(int j=0; j<r.match.length; i++, j++){
				match2[i]=r.match[j];
			}
			for(; i<match2.length; i++){
				match2[i]='C';
			}
			r.match=match2;
		}
		r.bases=bases1;
		r.quality=qual1;
		r.start-=lt;
		r.stop+=rt;
		if(returnToShort){
			r.match=Read.toShortMatchString(r.match);
			r.setShortMatch(true);
		}
		
		if(r.sites!=null){
			for(SiteScore ss : r.sites){untrim(ss);}
		}
		
		return true;
	}
	
	private boolean untrim(SiteScore ss){
		if(ss==null){return false;}
		if(leftTrimmed==0 && rightTrimmed==0){return false;}
		ss.perfect=ss.semiperfect=false;
		
		final int lt, rt;
		if(ss.strand==Gene.PLUS){
			lt=leftTrimmed;
			rt=rightTrimmed;
		}else{
			lt=rightTrimmed;
			rt=leftTrimmed;
		}
		
		boolean returnToShort=false;
		if(verbose){System.err.println("Untrimming ss "+ss);}
		if(ss.match!=null){
			
			boolean shortmatch=false;
			for(byte b : ss.match){
				if(Character.isDigit(b)){shortmatch=true; break;}
			}
			
			if(shortmatch){
				ss.match=Read.toLongMatchString(ss.match);
				returnToShort=true;
			}
			byte[] match2=new byte[ss.match.length+lt+rt];
			int i=0;
			for(; i<lt; i++){
				match2[i]='C';
			}
			for(int j=0; j<ss.match.length; i++, j++){
				match2[i]=ss.match[j];
			}
			for(; i<match2.length; i++){
				match2[i]='C';
			}
			ss.match=match2;
		}
		ss.setLimits(ss.start-lt, ss.stop+rt);
		if(returnToShort){ss.match=Read.toShortMatchString(ss.match);}
		return true;
	}
	
	private static boolean trimMatch(Read r){
		if(r.match==null && r.sites==null){return false;}
		
		//Easy mode!
		r.match=null;
		if(r.sites!=null){
			for(SiteScore ss : r.sites){
				if(ss!=null){ss.match=null;}
			}
		}
		return true;
		
		//TODO - need to adjust read start and stop based on this information.  Also check strand!
//		byte[] match=r.match;
//		if(r.shortmatch()){match=Read.toLongMatchString(match);}
//		byte[] match2=new byte[match.length-leftTrimmed-rightTrimmed];
//		for(int mpos=0, bpos=0; bpos<leftTrimmed; mpos++){
//			byte m=match[mpos];
//			if(m=='D'){
//				//do nothing
//			}else{
//				bpos++;
//			}
//		}
	}
	

	
	public static int trimReadWithMatch(final Read r, final SamLine sl, final int leftTrimAmount, final int rightTrimAmount, int minFinalLength, int scafLen){
		if(leftTrimAmount<1 && rightTrimAmount<1){return 0;}
		if(leftTrimAmount+rightTrimAmount>r.length()){return 0;}
		if(r.match==null || r.bases==null){return 0;}

		r.toLongMatchString(false);
		byte[] match=r.match;
		
		assert(!r.shortmatch());

		
		ByteBuilder bb=new ByteBuilder(match.length);
		
		int leftTrim=leftTrimAmount>0 ? r.length() : 0, rightTrim=rightTrimAmount>0 ? r.length() : 0;
		{
			final int mlen=match.length;
			boolean keep=leftTrimAmount<1;
			int rpos=(sl==null ? 1 : sl.pos);
			for(int mpos=0, cpos=0; mpos<mlen; mpos++){
				final byte m=match[mpos];
				if(m=='m' || m=='C' || m=='S' || m=='V' || m=='N'){
					cpos++;
					rpos++;
				}else if(m=='D'){
					rpos++;
				}else if(m=='I'){
					cpos++;
				}else{
					assert(false) : "Unknown symbol "+(char)m;
				}

				if(keep){
					bb.append(m);
				}else if(cpos>=leftTrimAmount){
					byte next=(mpos<mlen-1 ? match[mpos+1] : (byte)'m');
					if(m=='m' || m=='S' ||  m=='V' ||  m=='N' || next=='m' || next=='S' || next=='V' || next=='N'){
						keep=true;
						if(sl!=null){sl.pos=rpos;}
						leftTrim=cpos;
					}
				}
			}
		}
		match=bb.toBytes();
		Tools.reverseInPlace(match);
		bb.clear();
		{
			final int mlen=match.length;
			boolean keep=rightTrimAmount<1;
			for(int mpos=0, cpos=0; mpos<mlen; mpos++){
				final byte m=match[mpos];
				if(m=='m' || m=='C' || m=='S' || m=='V' || m=='N'){
					cpos++;
				}else if(m=='D'){
				}else if(m=='I'){
					cpos++;
				}else{
					assert(false) : "Unknown symbol "+(char)m;
				}

				if(keep){
					bb.append(m);
				}else if(cpos>=rightTrimAmount){
					byte next=(mpos<mlen-1 ? match[mpos+1] : (byte)'m');
//					if(m=='m' || m=='S' ||  m=='V' ||  m=='N' || next=='m' || next=='S' || next=='V' || next=='N'){
					if(next!='I' && next!='D'){
						keep=true;
						rightTrim=cpos;
					}
//					}
				}
			}
		}
		match=bb.toBytes();
		Tools.reverseInPlace(match);
		
		if(leftTrim+rightTrim>=r.length()){return 0;}
		
		r.match=null;
		if(sl.strand()==Gene.MINUS){
			int temp=leftTrim;
			leftTrim=rightTrim;
			rightTrim=temp;
		}
		final int trimmed=trimByAmount(r, leftTrim, rightTrim, minFinalLength);
		r.match=match;
		
		if(sl!=null){
			int start=sl.pos-1;
			int stop=start+Read.calcMatchLength(match)-1;
			sl.cigar=SamLine.toCigar14(match, start, stop, scafLen, r.bases);
			sl.seq=r.bases;
			sl.qual=r.quality;
		}
		return trimmed;
	}
	
	
	public final Read r;
	
	/** untrimmed bases */
	public final byte[] bases1;
	/** untrimmed qualities */
	public final byte[] qual1;
	/** trimmed bases */
	public byte[] bases2;
	/** trimmed qualities */
	public byte[] qual2;
	
	
	public final byte trimq;
	public int leftTrimmed;
	public int rightTrimmed;
	
	/** Require this many consecutive good bases to stop trimming.  Minimum is 1. 
	 * This is for the old trimming mode and not really used anymore */
	public static int minGoodInterval=2;
	
	public static boolean verbose=false;
	public static boolean optimalMode=true;
	public static boolean windowMode=false;
	public static float optimalBias=-1f;
	
	public static int windowLength=4;
	
	private static final float NPROB=0.75f;
//	public static float PROB1=QualityTools.PROB_ERROR[1];
	
	
}

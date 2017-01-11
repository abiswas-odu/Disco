package stream;

import dna.AminoAcid;

/**
 * @author Brian Bushnell
 * @date May 5, 2016
 * 
 *
 */
public class MDWalker {
	
	MDWalker(String tag, byte[] longmatch_){
		mdTag=tag;
		longmatch=longmatch_;
		pos=(mdTag.startsWith("MD:Z:") ? 5 : 0);
		
		mpos=0;
		bpos=0;
		rpos=0;
		sym=0;
		current=0;
		mode=0;
		
		while(longmatch[mpos]=='C'){
			mpos++;
			bpos++;
		}
	}
	
	void fixMatch(byte[] bases){
		sym=0;
		while(pos<mdTag.length()){
			char c=mdTag.charAt(pos);
			pos++;
			
			if(Character.isDigit(c)){
				current=(current*10)+(c-'0');
				mode=NORMAL;
			}else{
				int mpos2=mpos;
				if(current>0){
					mpos2=mpos+current;
//					System.err.println(mpos+", "+current+", "+mpos2);
					assert(mode==NORMAL) : mode+", "+current;
					current=0;	
				}
				
				while(mpos<mpos2 || (mpos<longmatch.length && longmatch[mpos]=='I')){
					if(longmatch[mpos]=='I'){
//						System.err.println("I: mpos="+mpos+", bpos="+bpos);
						mpos2++;
					}else{
//						System.err.println("M: mpos="+mpos+", bpos="+bpos);
						rpos++;
					}
					mpos++;
					bpos++;
				}
				
//				while(mpos<longmatch.length && longmatch[mpos]=='I'){
//					System.err.println("I2: mpos="+mpos+", bpos="+bpos);
//					mpos++;
//					bpos++;
//				}
				
				if(c=='^'){
					mode=DEL;
//					System.err.println("c="+((char)c)+", mpos="+mpos+", rpos="+rpos+", bpos="+bpos+", mode="+mode+(mode==NORMAL ? "" : ", match="+(char)longmatch[mpos-1])+"\n"+new String(longmatch));
				}else if(mode==DEL){
//					System.err.println("c="+((char)c)+", mpos="+mpos+", rpos="+rpos+", bpos="+bpos+", mode="+mode+(mode==NORMAL ? "" : ", match="+(char)longmatch[mpos-1])+"\n"+new String(longmatch));
					rpos++;
					mpos++;
					sym=c;
				}
//				else if(longmatch[mpos]=='I'){
//					mode=INS;
//					bpos++;
//					mpos++;
//					sym=c;
//				}
				else if(mode==NORMAL || mode==SUB){
					assert(longmatch[mpos]!='D') : mpos+"\n"+new String(longmatch)+"\n"+mdTag;
					longmatch[mpos]=(byte)'S';
					if((bases!=null && !AminoAcid.isFullyDefined(bases[bpos])) || !AminoAcid.isFullyDefined(c)){longmatch[mpos]='N';}
					mode=SUB;
//					System.err.println("c="+((char)c)+", mpos="+mpos+", rpos="+rpos+", bpos="+bpos+", mode="+mode+(mode==NORMAL ? "" : ", match="+(char)longmatch[mpos-1])+"\n"+new String(longmatch));
					bpos++;
					rpos++;
					mpos++;
					sym=c;
				}else{
					assert(false);
				}
				
			}
			
		}
//		System.err.println();
//		assert((bases==null || Read.calcMatchLength(longmatch)==bases.length)) : 
//			bases.length+", "+Read.calcMatchLength(longmatch)+"\n"+new String(longmatch)+"\n"
//					+ new String(Read.toShortMatchString(longmatch))+"\n"+mdTag;
	}
	
	boolean nextSub(){
		sym=0;
		while(pos<mdTag.length()){
			char c=mdTag.charAt(pos);
			pos++;
			
			if(Character.isDigit(c)){
				current=(current*10)+(c-'0');
				mode=NORMAL;
			}else{
				if(current>0){
					bpos+=current;
					rpos+=current;
					mpos+=current; 
					assert(mode==NORMAL) : mode+", "+current;
					current=0;	
				}
				if(c=='^'){mode=DEL;}
				else if(mode==DEL){
					rpos++;
					mpos++;
					sym=c;
				}else if(longmatch[mpos]=='I'){
					mode=INS;
					bpos++;
					mpos++;
					sym=c;
				}else if(mode==NORMAL || mode==SUB || mode==INS){
					mode=SUB;
					bpos++;
					rpos++;
					mpos++;
					sym=c;
//					System.err.println("c="+((char)c)+", mpos="+mpos+", rpos="+rpos+", bpos="+bpos+", mode="+mode+(mode==NORMAL ? "" : ", match="+(char)longmatch[mpos-1])+"\n"+new String(longmatch));
					return true;
				}
			}
			
//			System.err.println("c="+((char)c)+", mpos="+mpos+", rpos="+rpos+", bpos="+bpos+", mode="+mode+(mode==NORMAL ? "" : ", match="+(char)longmatch[mpos-1])+"\n"+new String(longmatch));
		}
		return false;
	}
	
	public int matchPosition(){
		return mpos-1;
	}
	
	public int basePosition(){
		return bpos-1;
	}
	
	public int refPosition(){
		return rpos-1;
	}
	
	public char symbol(){
		assert(sym!=0);
		return sym;
	}

	/** Position in match string (excluding clipping and insertions) */
	private int mpos;
	/** Position in read bases (excluding clipping and insertions) */
	private int bpos;
	/** Position in reference bases (excluding clipping) */
	private int rpos;
	private char sym;
	
	private String mdTag;
	private byte[] longmatch;
	private int pos;
	private int current;
	private int mode;
	
//	private int dels=0, subs=0, normals=0;
	private static final int NORMAL=0, SUB=1, DEL=2, INS=3;
	
}

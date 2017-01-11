package var2;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import align2.QualityTools;
import dna.AminoAcid;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ByteBuilder;
import stream.Read;
import stream.SamLine;

/**
 * Tracks data for a variation.
 * Uses half-open coordinates.
 * @author Brian Bushnell
 * @date November 4, 2016
 *
 */
public class Var implements Comparable<Var> {
	
	public static void main(String[] args){
		Timer t=new Timer();
		VarMap map=VarMap.loadVars(args[0], null);
		t.stop("Loaded "+map.size()+" variants.\nTime: \t");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	public Var(int scafnum_, int start_, int stop_, int allele_){
		this(scafnum_, start_, stop_, AL_MAP[allele_]);
	}

	public Var(Var v) {
		this(v.scafnum, v.start, v.stop, v.allele);
	}
	
	public Var(int scafnum_, int start_, int stop_, byte[] allele_){
		scafnum=scafnum_;
		start=start_;
		stop=stop_;
		allele=allele_;
		hashcode=hash();
//		stamp=stamper.getAndIncrement();
		assert(allele.length>1 || allele==AL_0 || 
				allele==AL_A || allele==AL_C || allele==AL_G || allele==AL_T || allele==AL_N);
		assert(start<=stop) : "\n"+Var.toBasicHeader()+"\n"+this+"\n";
	}
	
	public Var(final byte[] line, final byte delimiter){
		int a=0, b=0;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 0: "+new String(line);
		scafnum=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 1: "+new String(line);
		start=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 2: "+new String(line);
		stop=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 3: "+new String(line);
//		type=
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>=a) : "Missing field 4: "+new String(line);
		if(b==a){allele=AL_0;}
		else if(b==a+1){allele=AL_MAP[line[a]];}
		else{allele=Arrays.copyOfRange(line, a, b);}
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 5: "+new String(line);
		r1plus=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 6: "+new String(line);
		r1minus=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 7: "+new String(line);
		r2plus=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 8: "+new String(line);
		r2minus=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 9: "+new String(line);
		properPairCount=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 10: "+new String(line);
		lengthSum=Tools.parseLong(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 11: "+new String(line);
		mapQSum=Tools.parseLong(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 12: "+new String(line);
		mapQMax=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 13: "+new String(line);
		baseQSum=Tools.parseLong(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 14: "+new String(line);
		baseQMax=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 15: "+new String(line);
		endDistSum=Tools.parseLong(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 16: "+new String(line);
		endDistMax=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 17: "+new String(line);
		idSum=Tools.parseLong(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 18: "+new String(line);
		idMax=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 19: "+new String(line);
		coverage=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 20: "+new String(line);
		minusCoverage=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		if(b>a){
//			phredScore=Tools.parseFloat(line, a, b);
			b++;
			a=b;
		}
		
		hashcode=hash();
		assert(allele.length>1 || allele==AL_0 || 
				allele==AL_A || allele==AL_C || allele==AL_G || allele==AL_T || allele==AL_N);
		assert(start<=stop) : this.toString();
	}

	//#CHROM POS    ID        REF  ALT     QUAL
	public static Var fromVCF(byte[] line, ScafMap scafMap) {
		int a=0, b=0;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(line);
		String scaf=new String(line, a, b-a);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 1: "+new String(line);
		int pos=Tools.parseInt(line, a, b);
		b++;
		a=b;

		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 2: "+new String(line);
		//String id=new String(line, a, b-a);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 3: "+new String(line);
		//byte[] ref=Arrays.copyOf(line, a, b);
		int reflen=line[a]=='.' ? 0 : b-a;
		b++;
		a=b;

		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 4: "+new String(line);
		byte[] alt;
		if(b<=a+1){alt=AL_MAP[line[a]];}
		else{alt=Arrays.copyOfRange(line, a, b);}
		b++;
		a=b;

//		while(b<line.length && line[b]!='\t'){b++;}
//		assert(b>a) : "Missing field 5: "+new String(line);
//		int qual=Tools.parseInt(line, a, b);
//		b++;
//		a=b;
		
		final int start;
		final int readlen;
		if(alt.length!=reflen && alt.length>0){
			alt=Arrays.copyOfRange(alt, 1, alt.length);
			start=pos;
		}else{
			start=pos-1;
		}
		readlen=alt.length;
		final int stop=start+reflen;
		assert(scaf!=null);
		assert(scafMap!=null);
		final int scafNum=scafMap.getNumber(scaf);
		assert(scafNum>=0) : scaf+"\n"+scafMap.keySet()+"\n"+scafMap.altKeySet()+"\n";
//		final Scaffold scaffold=scafMap.getScaffold(scafNum);
		
		return new Var(scafNum, start, stop, alt);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Mutators           ----------------*/
	/*--------------------------------------------------------------*/
	
	public void add(Var b){
		final int oldReads=count();
		
//		assert(oldReads>0) : this;
		assert(oldReads==0 || baseQSum/oldReads<=60) : this;
		
		assert(this.equals(b));
		r1plus+=b.r1plus;
		r1minus+=b.r1minus;
		r2plus+=b.r2plus;
		r2minus+=b.r2minus;
		properPairCount+=b.properPairCount;
		lengthSum+=b.lengthSum;
		
		mapQSum+=b.mapQSum;
		mapQMax=Tools.max(mapQMax, b.mapQMax);
		baseQSum+=b.baseQSum;
		baseQMax=Tools.max(baseQMax, b.baseQMax);

		endDistSum+=b.endDistSum;
		endDistMax=Tools.max(endDistMax, b.endDistMax);

		idSum+=b.idSum;
		idMax=Tools.max(idMax, b.idMax);

//		assert(count()>0 && count()>oldReads) : "\n"+this+"\n"+b;
		assert(count()>=oldReads) : "\n"+this+"\n"+b;
		assert(count()==oldReads+b.count()) : "\n"+this+"\n"+b;
		assert(count()==0 || baseQSum/count()<=60) : "\n"+this+"\n"+b;
	}
	
	public void add(Read r){
		SamLine sl=(SamLine)r.obj;
		final int bstart=calcBstart(r, sl);
		final int bstop=calcBstop(bstart, r);
		add(r, bstart, bstop);
	}
		
	public void add(Read r, final int bstart, final int bstop){

		final int oldReads=count();
		
		SamLine sl=(SamLine)r.obj;
		
		if(sl.strand()==0){
			if(sl.pairnum()==0){
				r1plus++;
			}else{
				r2plus++;
			}
		}else{
			if(sl.pairnum()==0){
				r1minus++;
			}else{
				r2minus++;
			}
		}
		
		lengthSum+=r.length();
		properPairCount+=(sl.properPair() ? 1 : 0);
		mapQSum+=sl.mapq;
		mapQMax=Tools.max(mapQMax, sl.mapq);
		
		int baseQ=calcBaseQ(bstart, bstop, r, sl);
		baseQSum+=baseQ;
		baseQMax=Tools.max(baseQMax, baseQ);
		
		int endDist=calcEndDist(bstart, bstop, r);
		endDistSum+=endDist;
		endDistMax=Tools.max(endDistMax, endDist);
		
		int id=(int)(1000*Read.identitySkewed(r.match, false, false, false, true));
		idSum+=id;
		idMax=Tools.max(idMax, id);
		
		assert(count()>0) : this;
		assert(count()==oldReads+1) : this;
		assert(baseQSum/count()<=60) : this;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static ArrayList<Var> toVars(Read r, SamLine sl, boolean callNs, ScafMap scafMap){
		if(!r.containsVariants()){return null;}
		final int scafnum=scafMap.getNumber(sl.rnameS());
		return toVars(r, sl, callNs, scafnum);
	}
	
	/**
	 * @TODO This crashes on indels in the last position in the match string.
	 * @param r
	 * @param sl
	 * @param callNs
	 * @param scafnum
	 * @return
	 */
	public static ArrayList<Var> toVars(Read r, SamLine sl, boolean callNs, final int scafnum){
		if(!r.containsVariants()){return null;}
		ArrayList<Var> list=new ArrayList<Var>();
		
		r.toLongMatchString(false);
		if(sl.strand()==1 && !r.swapped()){
			r.reverseComplement();
			r.setSwapped(true);
		}
		final byte[] match=r.match;
		final byte[] bases=r.bases;
		final int rpos0=sl.pos-1;

		int bstart=-1, bstop=-1;
		int rstart=-1, rstop=-1;
		int mode=-1;
		
		int mpos=0, bpos=0, rpos=rpos0;
		for(; mpos<match.length; mpos++){
			byte m=match[mpos];
			
			if(m!=mode){
				if(mode=='D'){
					bstop=bpos;
					rstop=rpos;
//					assert(false) : (char)m+", "+(char)mode+", "+rstart+", "+bstart;
					Var v=new Var(scafnum, rstart, rstop, 0);
					v.add(r, bstart, bstop);
					list.add(v);
					bstart=bstop=rstart=rstop=-1;
				}else if(mode=='I'){
					bstop=bpos;
					rstop=rpos;
					int blen=bstop-bstart;
					Var v;
					if(blen==1){
						v=new Var(scafnum, rstart, rstop, bases[bstart]);
					}else{
						v=new Var(scafnum, rstart, rstop, Arrays.copyOfRange(bases, bstart, bstop));
					}
					v.add(r, bstart, bstop);
					list.add(v);
					bstart=bstop=rstart=rstop=-1;
				}
			}
			
			if(m=='C'){
				bpos++;
			}else if(m=='m' || m=='S' || m=='N'){
				if(m=='S' || (m=='N' && callNs)){
					Var v=new Var(scafnum, rpos, rpos+1, bases[bpos]);
					v.add(r, bpos, bpos+1);
					list.add(v);
				}
				bpos++;
				rpos++;
			}else if(m=='D'){
				if(mode!=m){
					rstart=rpos;
					bstart=bpos;
				}
				rpos++;
			}else if(m=='I'){
				if(mode!=m){
					rstart=rpos;
					bstart=bpos;
				}
				bpos++;
			}else{
				assert(false) : "Unhandled symbol "+(char)m;
			}
			mode=m;
		}
		
		if(mode=='D'){
			bstop=bpos;
			rstop=rpos;
			Var v=new Var(scafnum, rstart, rstop, 0);
			v.add(r, bstart, bstop);
			list.add(v);
			bstart=bstop=rstart=rstop=-1;
		}else if(mode=='I'){
			bstop=bpos;
			rstop=rpos-1;
			int blen=bstop-bstart;
			Var v;
			assert(rstart<=rstop) : "\n"+rstart+", "+rstop+", "+rpos+
									"\n"+bstart+", "+bstop+", "+bpos+
									"\n"+r+"\n"+sl;
			if(blen==1){
				v=new Var(scafnum, rstart, rstop, bases[bstart]);
			}else{
				v=new Var(scafnum, rstart, rstop, Arrays.copyOfRange(bases, bstart, bstop));
			}
			v.add(r, bstart, bstop);
			list.add(v);
			bstart=bstop=rstart=rstop=-1;
		}
		
//		for(Var v : list){v.add(r);}
		
		return list;
	}

	public int calcBstart(Read r, SamLine sl){
		r.toLongMatchString(false);
		byte[] match=r.match;
		final int rstart=sl.pos-1;
		final int type=type();
		
		int bstart=-1;
		
		for(int mpos=0, rpos=rstart, bpos=0; mpos<match.length; mpos++){
			byte m=match[mpos];
			if(m=='C'){
				bpos++;
			}else if(m=='m' || m=='S' || m=='N'){
				if(rpos==rstart){
					assert(type==SUB || type==NOCALL) : type+", "+bpos+", "+rpos+"\n"+new String(match);
					bstart=bpos;
					break;
				}
				bpos++;
				rpos++;
			}else if(m=='D'){
				if(rpos==rstart){
					assert(type==DEL) : type+", "+rpos+"\n"+new String(match);
					bstart=bpos;
					break;
				}
				rpos++;
			}else if(m=='I'){
				if(rpos==rstart && type==INS){
					bstart=bpos;
					break;
				}
				bpos++;
			}else{
				assert(false) : "Unhandled symbol "+(char)m;
			}
		}
		assert(bstart>=0);
		return bstart;
	}
	
	public int calcBstop(int bstart, Read r){
		assert(bstart>=0);
		int bstop=bstart+readlen();
		assert(bstop<=r.length());
		return bstop;
	}
	
	public int calcEndDist(int bstart, int bstop, Read r){
		int dist=Tools.min(bstart, r.length()-bstop);
		assert(dist>=0 && dist<=r.length()/2) : dist+", "+r.length()+", "+bstart+", "+bstop+"\n"+this+"\n"+
			"\n"+new String(r.match)+"\n"+r.obj+"\n";
		assert(dist<=(r.length()-readlen())/2);
		return dist;
	}
	
	public int calcBaseQ(final int bstart0, final int bstop0, Read r, SamLine sl){
		final byte[] quals=r.quality;
		if(quals==null){return Shared.FAKE_QUAL;}
		final int type=type();
		final int bstart, bstop;
		final int len=r.length();
		
		if(sl.strand()==0 || (sl.strand()==1 && r.swapped())){
			bstart=bstart0;
			bstop=bstop0;
		}else{
			bstart=len-bstop0-1;
			bstop=len-bstart0-1;
			assert(bstop-bstart==bstop0-bstart0);
		}
		
		int sum=0, avg=0;
		if(type==DEL){
			if(bstart==0){
				sum=avg=quals[0];
			}else if(bstop>=len-1){
				sum=avg=quals[len-1];
			}else{
				assert(bstop==bstart) : bstart0+", "+bstop0+", "+bstart+", "+bstop+"\n"+
						r.length()+", "+r.swapped()+", "+type()+", "+readlen()+", "+reflen()+
						"\n"+this+"\n"+new String(r.match)+"\n"+r.obj+"\n";
				
//				-1, 73, -1, 73
//				151, true, 2, 0, 1
				
				sum=quals[bstart]+quals[bstop+1];
				avg=sum/2;
			}
		}else{
			for(int i=bstart; i<bstop; i++){
				sum+=quals[i];
			}
			avg=Math.round(sum/(float)(bstop-bstart));
		}
		return avg;
	}

	int reflen(){
		return stop-start;
	}
	
	int readlen(){
		return allele.length;
	}
	
	int type(){
		int reflen=reflen(), readlen=readlen();
		if(reflen==0){return INS;}
		if(readlen==0){return DEL;}
//		assert(start<=stop) : start+", "+stop;
//		assert(reflen==readlen) : reflen+", "+readlen+", "+new String(allele)+", "+start+", "+stop;
		for(byte b : allele){
			if(b!='N'){return SUB;}
		}
		return NOCALL;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Contract Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public boolean equals(Object b){
		return equals((Var)b);
	}
	
	public boolean equals(Var b){
		return hashcode==b.hashcode && compareTo(b)==0;
	}
	
	@Override
	public int hashCode(){
		return hashcode;
	}
	
	public long toKey() {
		long key=Long.rotateLeft(start, 30)^Long.rotateLeft(scafnum, 10)^hash(allele);
		return key&0x3FFFFFFFFFFFFFFFL;
	}
	
	@Override
	public int compareTo(Var v){
		if(scafnum!=v.scafnum){return scafnum-v.scafnum;}
		if(start!=v.start){return start-v.start;}
		if(stop!=v.stop){return v.stop-stop;}
		return compare(allele, v.allele);
	}
	
	public int compare(byte[] a, byte[] b){
		if(a==b){return 0;}
		if(a.length!=b.length){return b.length-a.length;}
		for(int i=0; i<a.length; i++){
			byte ca=a[i], cb=b[i];
			if(ca!=cb){return ca-cb;}
		}
		return 0;
	}
	
	public String toString(){
		return toText(new ByteBuilder(), 0.99f, 30, 30, 1, 2, null).toString();
	}
	
	public static String toHeader(float properPairRate, float totalQualityAvg, float mapqAvg, int ploidy){
		StringBuilder sb=new StringBuilder();
		sb.append("#ploidy\t"+ploidy+"\n");
		sb.append("#properPairRate\t"+properPairRate+"\n");
		sb.append("#totalQualityAvg\t"+totalQualityAvg+"\n");
		sb.append("#mapqAvg\t"+mapqAvg+"\n");
		sb.append("#scaf\tstart\tstop\ttype\tcall\tr1p\tr1m\tr2p\tr2m\tpaired\tlengthSum");
		sb.append("\tmapq\tmapqMax\tbaseq\tbaseqMax\tedist\tedistMax\tid\tidMax");
		sb.append("\tcov\tminusCov");
		sb.append("\tphredScore");
		if(extendedText){
			sb.append("\treadCount\talleleFraction\tstrandRatio\tbaseqAvg\tmapqAvg\tedistAvg\tidentityAvg");
			sb.append("\tedistScore\tidentityScore\tqualityScore\tpairedScore\tbiasScore\tcoverageScore\thomopolymerScore\tscore");
		}
		return sb.toString();
	}
	
	public static String toBasicHeader(){
		StringBuilder sb=new StringBuilder();
		sb.append("#scaf\tstart\tstop\ttype\tcall\tr1p\tr1m\tr2p\tr2m\tpaired\tlengthSum");
		sb.append("\tmapq\tmapqMax\tbaseq\tbaseqMax\tedist\tedistMax\tid\tidMax");
		sb.append("\tcov\tminusCov");
		sb.append("\tphredScore");
		return sb.toString();
	}
	
	public ByteBuilder toText(ByteBuilder bb, float properPairRate, float totalQualityAvg, float totalMapqAvg, float rarity, int ploidy, ScafMap map){
		useIdentity=true;
		bb.append(scafnum).append('\t');
		bb.append(start).append('\t');
		bb.append(stop).append('\t');
		bb.append(typeArray[type()]).append('\t');
		for(byte b : allele){bb.append(b);}
		bb.append('\t');
		
		bb.append(r1plus).append('\t');
		bb.append(r1minus).append('\t');
		bb.append(r2plus).append('\t');
		bb.append(r2minus).append('\t');
		bb.append(properPairCount).append('\t');
		bb.append(lengthSum).append('\t');

		bb.append(mapQSum).append('\t');
		bb.append(mapQMax).append('\t');
		bb.append(baseQSum).append('\t');
		bb.append(baseQMax).append('\t');
		bb.append(endDistSum).append('\t');
		bb.append(endDistMax).append('\t');
		bb.append(idSum).append('\t');
		bb.append(idMax).append('\t');

		bb.append(coverage).append('\t');
		bb.append(minusCoverage).append('\t');

//		bb.append(prevBase<0 ? 'N' : (char)prevBase).append('\t');
		
		final float score=score(properPairRate, totalQualityAvg, totalMapqAvg, rarity, ploidy, map);
		bb.append(String.format("%.2f", toPhredScore(score))).append('\t');
		
		if(extendedText){
			
			bb.append(count()).append('\t');
			bb.append(alleleFraction()).append('\t');
			bb.append(strandRatio()).append('\t');
			bb.append(baseQAvg()).append('\t');
			bb.append(mapQAvg()).append('\t');
			bb.append(edistAvg()).append('\t');
			bb.append(identityAvg()).append('\t');
			
			bb.append(edistScore()).append('\t');
			bb.append(identityScore()).append('\t');
			bb.append(qualityScore(totalQualityAvg)).append('\t');
			bb.append(pairedScore(properPairRate)).append('\t');
			bb.append(biasScore(properPairRate, map)).append('\t');
			bb.append(coverageScore(ploidy, rarity)).append('\t');
			bb.append(homopolymerSubScore(map)).append('\t');
			bb.append(score).append('\t');
		}
		
		bb.length--;
		
		return bb;
	}

	public static String toVcfHeader(float properPairRate, float totalQualityAvg, float mapqAvg, float rarity, float minAlleleFraction, int ploidy, 
			long reads, long pairs, long properPairs, ScafMap map, String sampleName, String ref, boolean trimWhitespace) {
		StringBuilder sb=new StringBuilder();
		
		sb.append("##fileformat=VCFv4.2\n");
		sb.append("##BBMapVersion="+Shared.BBMAP_VERSION_STRING+"\n");
		sb.append("##ploidy="+ploidy+"\n");
		sb.append(String.format("##rarity=%.5f\n", rarity));
		sb.append(String.format("##minallelefraction=%.5f\n", minAlleleFraction));
		sb.append("##reads="+reads+"\n");
		sb.append("##pairedReads="+pairs+"\n");
		sb.append("##properlyPairedReads="+properPairs+"\n");
		sb.append(String.format("##properPairRate=%.5f\n", properPairRate));
		sb.append(String.format("##totalQualityAvg=%.3f\n", totalQualityAvg));
		sb.append(String.format("##mapqAvg=%.3f\n", mapqAvg));
		if(ref!=null){sb.append("##reference="+ref+"\n");}
		
		for(Scaffold scaf : map.list){
			String name=scaf.name;
			if(trimWhitespace){name=Tools.trimWhitespace(name);}
			sb.append("##contig=<ID="+name+",length="+scaf.length+">\n");
		}
		
		{	
			sb.append("##FORMAT=<ID=PASS,Number=1,Type=String,Description=\"Pass\">\n");
			sb.append("##FORMAT=<ID=FAIL,Number=1,Type=String,Description=\"Fail\">\n");
			
			sb.append("##INFO=<ID=SN,Number=1,Type=Integer,Description=\"Scaffold Number\">\n");
			sb.append("##INFO=<ID=STA,Number=1,Type=Integer,Description=\"Start\">\n");
			sb.append("##INFO=<ID=STO,Number=1,Type=Integer,Description=\"Stop\">\n");
			sb.append("##INFO=<ID=TYP,Number=1,Type=Integer,Description=\"Type\">\n");
			
			sb.append("##INFO=<ID=R1P,Number=1,Type=Integer,Description=\"Read1 Plus Count\">\n");
			sb.append("##INFO=<ID=R1M,Number=1,Type=Integer,Description=\"Read1 Minus Count\">\n");
			sb.append("##INFO=<ID=R2P,Number=1,Type=Integer,Description=\"Read2 Plus Count\">\n");
			sb.append("##INFO=<ID=R2M,Number=1,Type=Integer,Description=\"Read2 Minus Count\">\n");
			sb.append("##INFO=<ID=PPC,Number=1,Type=Integer,Description=\"Paired Count\">\n");
			sb.append("##INFO=<ID=LS,Number=1,Type=Integer,Description=\"Length Sum\">\n");

			sb.append("##INFO=<ID=MQS,Number=1,Type=Integer,Description=\"MAPQ Sum\">\n");
			sb.append("##INFO=<ID=MQM,Number=1,Type=Integer,Description=\"MAPQ Max\">\n");
			sb.append("##INFO=<ID=BQS,Number=1,Type=Integer,Description=\"Base Quality Sum\">\n");
			sb.append("##INFO=<ID=BQM,Number=1,Type=Integer,Description=\"Base Quality Max\">\n");
			sb.append("##INFO=<ID=EDS,Number=1,Type=Integer,Description=\"End Distance Sum\">\n");
			sb.append("##INFO=<ID=EDM,Number=1,Type=Integer,Description=\"End Distance Max\">\n");
			sb.append("##INFO=<ID=IDS,Number=1,Type=Integer,Description=\"Identity Sum\">\n");
			sb.append("##INFO=<ID=IDM,Number=1,Type=Integer,Description=\"Identity Max\">\n");
			sb.append("##INFO=<ID=COV,Number=1,Type=Integer,Description=\"Coverage\">\n");
			sb.append("##INFO=<ID=MCOV,Number=1,Type=Integer,Description=\"Minus Coverage\">\n");
			sb.append("##INFO=<ID=HMP,Number=1,Type=Integer,Description=\"Homopolymer Count\">\n");

			sb.append("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
			sb.append("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Fraction\">\n");
			sb.append("##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Ref+, Ref-, Alt+, Alt-\">\n");
			
			sb.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
			sb.append("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
			sb.append("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">\n");
			sb.append("##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele Fraction\">\n");
			sb.append("##FORMAT=<ID=SC,Number=1,Type=Float,Description=\"Score\">\n");
			sb.append("##FORMAT=<ID=PF,Number=1,Type=String,Description=\"Pass Filter\">\n");
			
		}
		
		sb.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
		if(sampleName!=null){sb.append('\t').append(sampleName);}
		return sb.toString();
	}
	
	//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	1323_1066121
	//ChrIII_A_nidulans_FGSC_A4	133	.	T	C	204	PASS
	//	DP=27;VDB=1.614640e-01;RPB=9.947863e-01;AF1=0.5;AC1=1;DP4=5,8,6,8;MQ=49;FQ=180;PV4=1,1,0.055,1
	//	GT:PL:DP:SP:GQ	0/1:234,0,207:27:0:99
	public ByteBuilder toVCF(ByteBuilder bb, float properPairRate, float totalQualityAvg, float mapqAvg, int ploidy, ScafMap map, VarFilter filter, boolean trimWhitespace){
		
		final Scaffold scaf=map.getScaffold(scafnum);
		final byte[] bases=scaf.bases;
		final int reflen=reflen(), readlen=readlen(), type=type();
		final double score=phredScore(properPairRate, totalQualityAvg, mapqAvg, filter.rarity, ploidy, map);
		final boolean pass=(filter==null ? true : 
			filter.passesFilter(this, properPairRate, totalQualityAvg, mapqAvg, ploidy, map));
		
		bb.append(trimWhitespace ? Tools.trimWhitespace(scaf.name) : scaf.name).append('\t');
		boolean indel=(type==INS || type==DEL);
		boolean addPrevBase=true;
		bb.append(start+(indel && addPrevBase ? 0 : 1)).append('\t');
		bb.append('.').append('\t');
		
		final byte prevBase=(bases==null ? (byte)'N' : bases[Tools.mid(start-1, 0, bases.length-1)]);
		
		if(addPrevBase){
			if(reflen==0 || allele.length<1){bb.append(prevBase);}
			for(int i=0, rpos=start; i<reflen; i++, rpos++){
				bb.append(bases==null || rpos<0 || rpos>=bases.length ? (char)'N' : (char)bases[rpos]);
			}
			bb.append('\t');

			if(reflen==0 || allele.length<1){bb.append(prevBase);}
			bb.append(allele).append('\t');
		}else{
			if(reflen==0){
				bb.append('.');
			}else{
				for(int i=0, rpos=start; i<reflen; i++, rpos++){
					bb.append(bases==null || rpos<0 || rpos>=bases.length ? (char)'N' : (char)bases[rpos]);
				}
			}
			bb.append('\t');

			if(allele.length<1){
				bb.append('.').append('\t');
			}else{
				bb.append(allele).append('\t');
			}
		}

		bb.append(String.format("%.2f\t", score));
		bb.append(pass ? "PASS\t" : "FAIL\t");

		final int count=count();
		{
			assert(Scaffold.trackStrand==(minusCoverage>=0)) : Scaffold.trackStrand+", "+minusCoverage;
//			final int covMinus=(Scaffold.trackStrand ? scaf.minusCoverage(this) : coverage/2);
			final int covMinus=(Scaffold.trackStrand ? minusCoverage : coverage/2);
			final int covPlus=coverage-covMinus;
			final int refMinus=covMinus-minusCount();
			final int refPlus=covPlus-plusCount();
			bb.append("SN=").append(scafnum).append(';');
			bb.append("STA=").append(start).append(';');
			bb.append("STO=").append(stop).append(';');
			bb.append("TYP=").append(typeArray[type()]).append(';');
			
			bb.append("R1P=").append(r1plus).append(';');
			bb.append("R1M=").append(r1minus).append(';');
			bb.append("R2P=").append(r2plus).append(';');
			bb.append("R2M=").append(r2minus).append(';');
			bb.append("PPC=").append(properPairCount).append(';');
			bb.append("LS=").append(lengthSum).append(';');

			bb.append("MQS=").append(mapQSum).append(';');
			bb.append("MQM=").append(mapQMax).append(';');
			bb.append("BQS=").append(baseQSum).append(';');
			bb.append("BQM=").append(baseQMax).append(';');
			bb.append("EDS=").append(endDistSum).append(';');
			bb.append("EDM=").append(endDistMax).append(';');
			bb.append("IDS=").append(idSum).append(';');
			bb.append("IDM=").append(idMax).append(';');
			bb.append("COV=").append(coverage).append(';');
			bb.append("MCOV=").append(minusCoverage).append(';');
			bb.append("HMP=").append(homopolymerCount(map)).append(';');

			bb.append("DP=").append(Tools.max(coverage, count)).append(';');
			bb.append("AF=").append(String.format("%.4f",alleleFraction())).append(';');
			bb.append("DP4=").append(refPlus).append(',').append(refMinus).append(',').append(plusCount()).append(',').append(minusCount()).append(';');
			
			bb.length--;
		}
		{
			bb.append('\t');
			bb.append("GT:DP:AD:AF:SC:PF");
			bb.append('\t');

			bb.append(genotype(ploidy));
			bb.append(':');
			bb.append(Tools.max(coverage, count));
			bb.append(':');
			bb.append(count);
			bb.append(':');
			bb.append(String.format("%.4f",alleleFraction()));
			bb.append(':');
			bb.append(String.format("%.2f",score));
			bb.append(':');
			bb.append(pass ? "PASS" : "FAIL");
		}
		
		return bb;
	}
	
	//TODO: Actually, I should also track the ref coverage. 
	private String genotype(int ploidy) {
		if(ploidy==1){return "1";}
		final float af=alleleFraction();
		final int count=count();
		if(ploidy==2){
			if(af<0.2){return "0/0";}
			if(af<0.8){return "0/1";}
			return "1/1";
		}
		StringBuilder sb=new StringBuilder();
		
		int copies=Math.round(ploidy*af);
		if(af>=0.5){copies=Tools.max(copies, 1);}
		int refCopies=ploidy-copies;
		
		for(int i=0; i<refCopies; i++){
			sb.append(0).append('/');
		}
		for(int i=0; i<copies; i++){
			sb.append(1).append('/');
		}
		sb.setLength(sb.length()-1);
		return sb.toString();
	}

	private int hash(){
		return scafnum^Integer.rotateLeft(start, 9)^Integer.rotateRight(stop, 9)^hash(allele);
	}
	
	public static final int hash(byte[] a){
		int code=123456789;
		for(byte b : a){
			code=Integer.rotateLeft(code, 3)^codes[b];
		}
		return code&Integer.MAX_VALUE;
	}
	
	public int calcCoverage(ScafMap map){
		if(coverage>=0){return coverage;}
		
		Scaffold scaf=map.getScaffold(scafnum);
		coverage=scaf.calcCoverage(this);
		if(Scaffold.trackStrand){minusCoverage=scaf.minusCoverage(this);}
		return coverage;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Scoring Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	public static float toPhredScore(float score){
		if(score==0){return 0;}
		score=score*0.998f;
		return 2.5f*(float)QualityTools.probErrorToPhredDouble(1-score);
	}
	
	public float phredScore(float properPairRate, float totalQualityAvg, float totalMapqAvg, float rarity, int ploidy, ScafMap map){
		float score=score(properPairRate, totalQualityAvg, totalMapqAvg, rarity, ploidy, map);
		return toPhredScore(score);
	}
	
	public float score(float properPairRate, float totalQualityAvg, float totalMapqAvg, float rarity, int ploidy, ScafMap map){
		float cs=coverageScore(ploidy, rarity);
		if(cs==0){return 0;}
		float es=(useEdist ? edistScore() : 1);
		float qs=qualityScore(totalQualityAvg);
		float ps=(usePairing ? pairedScore(properPairRate) : 1);
		float bs=(useBias ? biasScore(properPairRate, map) : 1);
		float is=(useIdentity ? identityScore() : 1);
		float hs=(useHomopolymer ? homopolymerSubScore(map) : 1);
		return (float)Math.pow(es*qs*ps*bs*cs*is*hs, 0.2);
	}
	
	public float edistScore(){
		float lengthAvg=lengthAvg();
		float edistAvg=((edistAvg()*2+endDistMax))*0.333333333f;
		float constant=5+Tools.min(20, lengthAvg*0.1f)+lengthAvg*0.01f;
		float weighted=Tools.max(0.05f, edistAvg-Tools.min(constant, edistAvg*0.95f));
		weighted=weighted*weighted;
		return weighted/(weighted+4);
	}
	
	public float identityScore(){
		float lengthAvg=lengthAvg();
		float idAvg=0.001f*(((identityAvg()+idMax))*0.5f);
		float weighted=Tools.min(1, (idAvg*lengthAvg+(0.65f*Tools.max(1, readlen())))/lengthAvg); //Diminish impact of this var on itself
		weighted=0.75f+0.25f*weighted; //Compress the range
		return weighted;
	}
	
	public float qualityScore(float totalBaseqAvg){
		float bqAvg=baseQAvg();
		float mqAvg=0.5f*(mapQAvg()+mapQMax);
		
		final float delta=totalBaseqAvg-bqAvg;
		if(delta>1){
			bqAvg=Tools.max(bqAvg*0.5f, bqAvg-0.5f*delta);
		}
		
		float mult=0.25f;
		float thresh=12;
		if(bqAvg>thresh){
			bqAvg=bqAvg-thresh+(thresh*mult);
		}else{
			bqAvg=bqAvg*mult;
		}
		
		double baseProbAvg=1-Math.pow(10, 0-.1*bqAvg);
		double mapProbAvg=1-Math.pow(10, 0-.1*(mqAvg+2));
		double d=baseProbAvg*baseProbAvg*mapProbAvg;
		return (float)d;
	}
	
	public float pairedScore(float properPairRate){
		if(properPairRate<0.5){return 0.98f;}
		final float count=count();
		if(count==0){return 0;}
		float rate=properPairCount/count;
		rate=rate*(count/(0.1f+count));
		if(rate>=properPairRate){return Tools.max(rate, 1-0.001f*properPairRate);}
		float score=rate/properPairRate;
		return Tools.max(0.1f, score);
	}
	
	public float coverageScore(int ploidy, float rarity){
		int count=count();
		if(count==0){return 0;}
		float rawScore=count/(lowCoveragePenalty+count); //This may be too severe...
		
//		float ratio=alleleFraction();
		
		float ratio=0.98f;
		if(coverage>0){
			float dif=coverage-count;
			if(dif>0){
				dif=dif-coverage*.01f-Tools.min(0.5f, coverage*.1f);
				dif=Tools.max(0.1f, dif);
			}
			ratio=(coverage-dif)/coverage;
			if(rarity<1 && ratio>rarity){
				float minExpected=1f/ploidy;
				if(ratio<minExpected){
					ratio=minExpected-((minExpected-ratio)*0.1f);
				}
			}
		}
		
		float ratio2=Tools.min(1, ploidy*ratio);
		return rawScore*ratio2;
	}
	
	public float homopolymerSubScore(ScafMap map){
		if(type()!=SUB || map==null){return 1;}
		
		int count=homopolymerCount(map);
		if(count==0){return 1;}
		return 1f-(count*0.1f/9);
	}
	
	public int homopolymerCount(ScafMap map){
		if(type()!=SUB || map==null){return 0;}
		final byte[] bases=map.getScaffold(scafnum).bases;
		if(bases==null){return 1;}
		final byte base=allele[0];
		
		if(start>=bases.length || stop<0){return 0;}
		
		int count1=0;
		for(int i=start-1, lim=Tools.max(0, start-4); i>=lim; i--){
			if(bases[i]==base){count1++;}
			else{break;}
		}
		int count2=0;
		for(int i=stop, lim=Tools.min(bases.length, stop+4); i<lim; i++){
			if(bases[i]==base){count2++;}
			else{break;}
		}
		assert(count1+count2<=8) : count1+", "+count2;
		
		return count1+count2+(count1>0 && count2>0 ? 1 : 0);
	}
	
	public float biasScore(float properPairRate, ScafMap map){
		float strandBias=strandBiasScore(map);
		float readBias=readBiasScore(properPairRate);
		return (float)Math.sqrt(strandBias*readBias);
	}
	
	public float strandBiasScore(ScafMap map){
		float x=eventProb(plusCount(), minusCount());
		if(x>=0.99 || !doNscan){return x;}
		final int nScan=500;
		int scafEndDist=(map==null ? start : scafEndDist(map, nScan));
		if(scafEndDist>=nScan){return x;}
		
		float delta=1-x;
		delta=delta*(scafEndDist*scafEndDist)/(float)(nScan*nScan);
//		delta=delta*(scafEndDist)/(float)(nScan);
//		System.err.println("scafEndDist="+scafEndDist+"; score: "+x+" -> "+(1-delta));
		return 1-delta;
	}
	
	public float readBiasScore(float properPairRate){
		if(properPairRate<0.5){return 0.95f;}
		
		return eventProb(r1Count(), r2Count());
	}
	
	/** Adjusted probability of a binomial event being at least this lopsided. */
	public static float eventProb(int a, int b){
		
		double allowedBias=0.75;
		double slopMult=0.95;
		
		double n=a+b;
		double k=Tools.min(a, b);
		
		double slop=n*(allowedBias*0.5);
		double dif=n-k*2;
		dif=dif-(Tools.min(slop, dif)*slopMult);
		n=k*2+dif;
//		k=n*0.5-dif;
		assert(k<=n*0.5) : a+", "+b+", "+n+", "+k+", "+slop+", "+dif;
		
		if(n>PROBLEN){
			double mult=PROBLEN/(double)n;
			n=PROBLEN;
			k=(int)(k*mult);
		}

		int n2=(int)Math.round(n);
		int k2=Tools.min(n2/2, (int)(k+1));
		

//		if(a+b>3){
//			System.err.println(n+", "+k+", "+n2+", "+k2);
//		}
		
		float result=(float)prob[n2][k2];
		if(result<1 || a==b || a+1==b || a==b+1){return result;}
		
		float slope=Tools.min(a, b)/(float)Tools.max(a, b);
		return (float)(0.998+slope*0.002);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	public int plusCount(){return r1plus+r2plus;}
	public int minusCount(){return r1minus+r2minus;}
	public int r1Count(){return r1plus+r1minus;}
	public int r2Count(){return r2plus+r2minus;}
	public int count(){return r1plus+r1minus+r2plus+r2minus;}

	public float alleleFraction(){
		int count=count();
		int cov=Tools.max(count, coverage, 1);
		return count/(float)cov;
	}
	
	public float strandRatio(){
		int plus=plusCount();
		int minus=minusCount();
		if(plus==minus){return 1;}
		return (Tools.min(plus,  minus)+1)/(float)Tools.max(plus, minus);
	}
	public float baseQAvg(){return baseQSum/(float)count();}
	public float mapQAvg(){return mapQSum/(float)count();}
	public float edistAvg(){return endDistSum/(float)count();}
	public float identityAvg(){return idSum/(float)count();}
	public float lengthAvg(){return lengthSum/(float)count();}
	public float properPairRate(){return properPairCount/(float)count();}
	
	
	public void setCoverage(int coverage_, int minusCoverage_){
		coverage=coverage_;
		minusCoverage=minusCoverage_;
	}
	
	public int coverage(){
		assert(coverage>-1) : coverage;
		return coverage;
	}
	
	public int scafEndDist(ScafMap map, int nScan){
		Scaffold scaf=map.getScaffold(scafnum);
		int len=scaf.length;
		byte[] bases=scaf.bases;
		
		int scafEndDist=Tools.max(0, Tools.min(start, len-stop));
		if(bases==null || nScan<1){return scafEndDist;}
		int limit=Tools.min(nScan, scafEndDist);
		int contigEndDist=leftContigEndDist(bases, limit);
		limit=Tools.min(limit, contigEndDist);
		contigEndDist=rightContigEndDist(bases, limit);
		return Tools.min(scafEndDist, contigEndDist);
	}
	
	public int leftContigEndDist(byte[] bases, int maxDist){
		if(start>=bases.length){return Tools.min(bases.length, maxDist+1);}
		int ns=0;
		for(int i=start, lim=Tools.max(0, start-maxDist); i>=lim; i--){
			if(AminoAcid.isFullyDefined(bases[i])){
				ns=0;
			}else{
				ns++;
				if(ns>=10){
					int x=start-i-ns+1;
					assert(x>=0);
					return x;
				}
			}
		}
		return maxDist+1;
	}
	
	public int rightContigEndDist(byte[] bases, int maxDist){
		if(stop<0){return Tools.min(bases.length, maxDist+1);}
		int ns=0;
		for(int i=stop, lim=Tools.min(bases.length-1, stop+maxDist); i<=lim; i++){
			if(AminoAcid.isFullyDefined(bases[i])){
				ns=0;
			}else{
				ns++;
				if(ns>=10){
					int x=i-stop-ns+1;
					assert(x>=0);
					return x;
				}
			}
		}
		return maxDist+1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public final int scafnum;
	public final int start;
	public final int stop; //Half-open, so stop is always after start except for insertions
	public final byte[] allele;
	public final int hashcode;
	
	/*--------------------------------------------------------------*/
	/*----------------        Mutable Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	int r1plus;
	int r1minus;
	int r2plus;
	int r2minus;
	int properPairCount;

	private int coverage=-1;
	private int minusCoverage=-1;
	
	long mapQSum;
	public int mapQMax;
	
	long baseQSum;
	public int baseQMax;
	
	long endDistSum;
	public int endDistMax;
	
	long idSum;
	int idMax;
	
	long lengthSum;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static boolean extendedText=false;
	
	public static boolean useHomopolymer=true;
	public static boolean useIdentity=true;
	public static boolean usePairing=true;
	public static boolean useBias=true;
	public static boolean useEdist=true;
	public static boolean doNscan=true;
	public static float lowCoveragePenalty=0.8f;
	
	final static int PROBLEN=100;
	
//	static final AtomicLong stamper=new AtomicLong(0);
	
	static final String[] typeArray=new String[] {"NOCALL","SUB","DEL","INS"};
	static final int NOCALL=0, SUB=1, DEL=2, INS=3;
	
	static final byte[] AL_0=new byte[0];
	static final byte[] AL_A=new byte[] {(byte)'A'};
	static final byte[] AL_C=new byte[] {(byte)'C'};
	static final byte[] AL_G=new byte[] {(byte)'G'};
	static final byte[] AL_T=new byte[] {(byte)'T'};
	static final byte[] AL_N=new byte[] {(byte)'N'};
	static final byte[][] AL_MAP=makeMap();
	static final int[] codes=makeCodes();

//	public static ScafMap scafMap;
	
	static final int[] makeCodes(){
		Random randy=new Random(1);
		int[] array=new int[256];
		for(int i=0; i<array.length; i++){
			array[i]=randy.nextInt();
		}
		return array;
	}
	
	static final byte[][] makeMap(){
		byte[][] map=new byte[128][];
		map[0]=map['.']=map['\t']=AL_0;
		map['A']=map['a']=AL_A;
		map['C']=map['c']=AL_C;
		map['G']=map['g']=AL_G;
		map['T']=map['t']=AL_T;
		map['N']=map['n']=AL_N;
		return map;
	}
	
	/** factorial[n]=n! */
	private static final double[] factorial=makeFactorialArray(PROBLEN+1);
	/** binomial[n][k] = combinations in n pick k */
	private static final double[][] binomial=makeBinomialMatrix(PROBLEN+1);
	/** prob[n][k] = probability of an event this lopsided or worse. */
	private static final double[][] prob=makeProbMatrix(PROBLEN+1);

	private static double[] makeFactorialArray(int len) {
		double[] x=new double[len];
		x[0]=1;
		for(int i=1; i<len; i++){
			x[i]=x[i-1]*i;
		}
		return x;
	}

	private static double[][] makeBinomialMatrix(int len) {
		double[][] matrix=new double[len][];
		for(int n=0; n<len; n++){
			final int kmax=n/2;
			final double nf=factorial[n];
			matrix[n]=new double[kmax+1];
			for(int k=0; k<=kmax; k++){
				final double kf=factorial[k];
				final double nmkf=factorial[n-k];
				double combinations=nf/kf;
				combinations=combinations/nmkf;
				matrix[n][k]=combinations;
			}
		}
		return matrix;
	}

	private static double[][] makeProbMatrix(int len) {
		double[][] matrix=new double[len][];
		double mult=2;
		for(int n=0; n<len; n++){
			final int kmax=n/2;
			final double[] array=matrix[n]=new double[kmax+1];
			for(int k=0; k<=kmax; k++){
				final double combinations=binomial[n][k];
				array[k]=combinations*mult;
			}
//			if(n<=12){System.err.println(Arrays.toString(array));}
			for(int k=0; k<=kmax; k++){
				array[k]=Tools.min(1, (k==0 ? 0 : array[k-1])+array[k]);
			}
//			if(n<=12){System.err.println(Arrays.toString(array));}
//			assert(array[kmax]==1) : Arrays.toString(array);
			mult*=0.5;
//			if(n<=12){System.err.println();}
		}
		return matrix;
	}
	
}

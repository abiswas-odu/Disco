package shared;

import java.util.ArrayList;
import java.util.Arrays;

import align2.QualityTools;
import dna.AminoAcid;

import stream.Read;
import stream.SamLine;
import structures.LongList;
import fileIO.TextStreamWriter;


/**
 * @author Brian Bushnell
 * @date Mar 18, 2013
 *
 */
public class ReadStats {
	
	public ReadStats(){this(true);}
		
	public ReadStats(boolean addToList){
		if(addToList){
			synchronized(ReadStats.class){
				objectList.add(this);
			}
		}

		if(COLLECT_QUALITY_STATS){
			aqualArray=new long[2][127];
			qualLength=new long[2][MAXLEN];
			qualSum=new long[2][MAXLEN];
			qualSumDouble=new double[2][MAXLEN];
			bqualHistOverall=new long[127];
		}else{
			aqualArray=null;
			qualLength=null;
			qualSum=null;
			qualSumDouble=null;
			bqualHistOverall=null;
		}
		
		if(BQUAL_HIST_FILE!=null){
			bqualHist=new long[2][MAXLEN][127];
		}else{
			bqualHist=null;
		}
		
		if(QUAL_COUNT_HIST_FILE!=null){
			qcountHist=new long[2][127];
		}else{
			qcountHist=null;
		}

		if(COLLECT_MATCH_STATS){
			matchSum=new long[2][MAXLEN];
			delSum=new long[2][MAXLEN];
			insSum=new long[2][MAXLEN];
			subSum=new long[2][MAXLEN];
			nSum=new long[2][MAXLEN];
			clipSum=new long[2][MAXLEN];
			otherSum=new long[2][MAXLEN];
		}else{
			matchSum=null;
			delSum=null;
			insSum=null;
			subSum=null;
			nSum=null;
			clipSum=null;
			otherSum=null;
		}
		
		if(COLLECT_QUALITY_ACCURACY){
			qualMatch=new long[99];
			qualSub=new long[99];
			qualIns=new long[99];
			qualDel=new long[99];
		}else{
			qualMatch=null;
			qualSub=null;
			qualIns=null;
			qualDel=null;
		}

		if(COLLECT_INSERT_STATS){
			insertHist=new LongList(MAXLEN);
		}else{
			insertHist=null;
		}

		if(COLLECT_BASE_STATS){
			baseHist=new LongList[2][5];
			for(int i=0; i<baseHist.length; i++){
				for(int j=0; j<baseHist[i].length; j++){
					baseHist[i][j]=new LongList(400);
				}
			}
		}else{
			baseHist=null;
		}
		

		if(COLLECT_INDEL_STATS){
			insHist=new LongList(100);
			delHist=new LongList(100);
			delHist2=new LongList(100);
		}else{
			insHist=null;
			delHist=null;
			delHist2=null;
		}
		
		if(COLLECT_GC_STATS){
			gcHist=new long[GC_BINS+1];
		}else{
			gcHist=null;
		}
		
		if(COLLECT_ERROR_STATS){
			errorHist=new LongList(100);
		}else{
			errorHist=null;
		}
		
		if(COLLECT_LENGTH_STATS){
			lengthHist=new LongList(501);
		}else{
			lengthHist=null;
		}
		
		if(COLLECT_IDENTITY_STATS){
			idHist=new long[ID_BINS+1];
			idBaseHist=new long[ID_BINS+1];
		}else{
			idHist=null;
			idBaseHist=null;
		}
		
		if(COLLECT_TIME_STATS){
			timeHist=new LongList(1001);
		}else{
			timeHist=null;
		}
		
	}
	
	public static ReadStats mergeAll(){
		if(objectList==null || objectList.isEmpty()){return merged=null;}
		if(objectList.size()==1){return merged=objectList.get(0);}
		
		ReadStats x=new ReadStats(false);
		for(ReadStats rs : objectList){
			x.read2Count+=rs.read2Count;
			if(COLLECT_QUALITY_STATS){
				for(int i=0; i<MAXLEN; i++){
					x.qualLength[0][i]+=rs.qualLength[0][i];
					x.qualLength[1][i]+=rs.qualLength[1][i];
					x.qualSum[0][i]+=rs.qualSum[0][i];
					x.qualSum[1][i]+=rs.qualSum[1][i];
					x.qualSumDouble[0][i]+=rs.qualSumDouble[0][i];
					x.qualSumDouble[1][i]+=rs.qualSumDouble[1][i];
				}
				for(int i=0; i<x.aqualArray[0].length; i++){
					x.aqualArray[0][i]+=rs.aqualArray[0][i];
					x.aqualArray[1][i]+=rs.aqualArray[1][i];
				}
				for(int i=0; i<x.bqualHistOverall.length; i++){
					x.bqualHistOverall[i]+=rs.bqualHistOverall[i];
				}
				if(BQUAL_HIST_FILE!=null){
					for(int i=0; i<x.bqualHist.length; i++){
						for(int j=0; j<x.bqualHist[i].length; j++){
							for(int k=0; k<x.bqualHist[i][j].length; k++){
								x.bqualHist[i][j][k]+=rs.bqualHist[i][j][k];
							}
						}
					}
				}
				if(QUAL_COUNT_HIST_FILE!=null){
					for(int i=0; i<x.qcountHist.length; i++){
						for(int j=0; j<x.qcountHist[i].length; j++){
							x.qcountHist[i][j]+=rs.qcountHist[i][j];
						}
					}
				}
			}
			
			if(COLLECT_MATCH_STATS){
				for(int i=0; i<MAXLEN; i++){
					x.matchSum[0][i]+=rs.matchSum[0][i];
					x.matchSum[1][i]+=rs.matchSum[1][i];
					x.delSum[0][i]+=rs.delSum[0][i];
					x.delSum[1][i]+=rs.delSum[1][i];
					x.insSum[0][i]+=rs.insSum[0][i];
					x.insSum[1][i]+=rs.insSum[1][i];
					x.subSum[0][i]+=rs.subSum[0][i];
					x.subSum[1][i]+=rs.subSum[1][i];
					x.nSum[0][i]+=rs.nSum[0][i];
					x.nSum[1][i]+=rs.nSum[1][i];
					x.clipSum[0][i]+=rs.clipSum[0][i];
					x.clipSum[1][i]+=rs.clipSum[1][i];
					x.otherSum[0][i]+=rs.otherSum[0][i];
					x.otherSum[1][i]+=rs.otherSum[1][i];
				}
			}
			if(COLLECT_INSERT_STATS){
				x.insertHist.add(rs.insertHist);
				x.pairedCount+=rs.pairedCount;
				x.unpairedCount+=rs.unpairedCount;
			}
			if(COLLECT_BASE_STATS){
				for(int i=0; i<rs.baseHist.length; i++){
					for(int j=0; j<rs.baseHist[i].length; j++){
						x.baseHist[i][j].add(rs.baseHist[i][j]);
					}
				}
			}
			if(COLLECT_QUALITY_ACCURACY){
				for(int i=0; i<x.qualMatch.length; i++){
					x.qualMatch[i]+=rs.qualMatch[i];
					x.qualSub[i]+=rs.qualSub[i];
					x.qualIns[i]+=rs.qualIns[i];
					x.qualDel[i]+=rs.qualDel[i];
				}
			}
			

			if(COLLECT_INDEL_STATS){
				x.delHist.add(rs.delHist);
				x.delHist2.add(rs.delHist2);
				x.insHist.add(rs.insHist);
			}

			if(COLLECT_LENGTH_STATS){
				x.lengthHist.add(rs.lengthHist);
			}
			

			if(COLLECT_ERROR_STATS){
				x.errorHist.add(rs.errorHist);
			}
			
			if(COLLECT_GC_STATS){
				for(int i=0; i<rs.gcHist.length; i++){
					x.gcHist[i]+=rs.gcHist[i];
				}
			}
			
			if(COLLECT_IDENTITY_STATS){
				for(int i=0; i<rs.idHist.length; i++){
					x.idHist[i]+=rs.idHist[i];
					x.idBaseHist[i]+=rs.idBaseHist[i];
				}
			}
			
			if(COLLECT_TIME_STATS){
				x.timeHist.add(rs.timeHist);
			}

			x.gcMaxReadLen=Tools.max(x.gcMaxReadLen, rs.gcMaxReadLen);
			x.idMaxReadLen=Tools.max(x.idMaxReadLen, rs.idMaxReadLen);
		}
		
		merged=x;
		return x;
	}
	
	public void addToQualityHistogram(final Read r){
		if(r==null){return;}
		addToQualityHistogram2(r);
		if(r.mate!=null){addToQualityHistogram2(r.mate);}
	}
	
	private void addToQualityHistogram2(final Read r){
		int pairnum=r.pairnum();
		if(r==null || r.quality==null || r.quality.length<1){return;}
		byte[] quals=r.quality, bases=r.bases;
		final Object obj=r.obj;
		if(obj!=null){
			if(obj.getClass()==SamLine.class){
				pairnum=((SamLine)obj).pairnum();
			}else if(obj.getClass()==TrimRead.class){
				quals=(pairnum==0 ? ((TrimRead)obj).qual1 : ((TrimRead)obj).qual2);
				bases=(pairnum==0 ? ((TrimRead)obj).bases1 : ((TrimRead)obj).bases2);
			}
		}
		if(pairnum==1){read2Count++;}
		addToQualityHistogram(quals, pairnum);
		int x=Read.avgQualityByProbability(bases, quals, true, 0);
		aqualArray[pairnum][x]++;
		if(BQUAL_HIST_FILE!=null){
			addToBQualityHistogram(quals, pairnum);
		}
		if(QUAL_COUNT_HIST_FILE!=null){
			addToQCountHistogram(quals, pairnum);
		}
	}
	
	public void addToQualityHistogram(byte[] qual, int pairnum){
		if(qual==null || qual.length<1){return;}
		final int limit=Tools.min(qual.length, MAXLEN);
		final long[] ql=qualLength[pairnum], qs=qualSum[pairnum];
		final double[] qsd=qualSumDouble[pairnum];
		ql[limit-1]++;
		for(int i=0; i<limit; i++){
			qs[i]+=qual[i];
			qsd[i]+=QualityTools.PROB_ERROR[qual[i]];
		}
		for(byte q : qual){
			bqualHistOverall[q]++;
		}
	}
	
	private void addToBQualityHistogram(byte[] qual, int pairnum){
		if(qual==null || qual.length<1){return;}
		final int limit=Tools.min(qual.length, MAXLEN);
		final long[][] bqh=bqualHist[pairnum];
		for(int i=0; i<limit; i++){
			bqh[i][qual[i]]++;
		}
	}
	
	private void addToQCountHistogram(byte[] qual, int pairnum){
		if(qual==null || qual.length<1){return;}
		final long[] qch=qcountHist[pairnum];
		for(byte q : qual){
			qch[q]++;
		}
	}
	
	public void addToQualityAccuracy(final Read r){
		if(r==null){return;}
		addToQualityAccuracy(r, 0);
		if(r.mate!=null){addToQualityAccuracy(r.mate, 1);}
	}
	
	public void addToQualityAccuracy(final Read r, int pairnum){
		if(r==null || r.quality==null || r.quality.length<1 || !r.mapped() || r.match==null/* || r.discarded()*/){return;}
		final byte[] bases=r.bases;
		final byte[] qual=r.quality;
		byte[] match=r.match;
		
		if(match!=null && r.shortmatch()){match=Read.toLongMatchString(match);}

		final boolean plus=(r.strand()==0);
		int rpos=0;
		byte lastm='A';
		for(int mpos=0; mpos<match.length/* && rpos<limit*/; mpos++){
			byte b=bases[rpos];
			byte q=qual[rpos];
			byte m=match[plus ? mpos : match.length-mpos-1];
			
			{
				if(m=='m'){
					qualMatch[q]++;
				}else if(m=='S'){
					qualSub[q]++;
				}else if(m=='I'){
					if(AminoAcid.isFullyDefined(b)){qualIns[q]++;}
				}else if(m=='N'){
					//do nothing
				}else if(m=='C'){
					//do nothing
				}else if(m=='V'){
					//do nothing
				}else if(m=='D'){
					if(lastm!=m){
						int x=rpos, y=rpos-1;
						if(x<qual.length){
							if(AminoAcid.isFullyDefined(bases[x])){
								qualDel[qual[x]]++;
							}
						}
						if(y>=0){
							if(AminoAcid.isFullyDefined(bases[y])){
								qualDel[qual[y]]++;
							}
						}
					}
					rpos--;
				}else{
					assert(!Character.isDigit(m)) : ((char)m);
				}
			}

			rpos++;
			lastm=m;
		}
		
	}
	
	public void addToErrorHistogram(final Read r){
		if(r==null){return;}
		addToErrorHistogram(r, 0);
		if(r.mate!=null){addToErrorHistogram(r.mate, 1);}
	}
	
	private void addToErrorHistogram(final Read r, int pairnum){
		if(r==null || r.bases==null || r.length()<1 || !r.mapped() || r.match==null/* || r.discarded()*/){return;}
		r.toLongMatchString(false);
		int x=r.countSubs();
		errorHist.increment(x, 1);
	}
	
	public void addToLengthHistogram(final Read r){
		if(r==null){return;}
		addToLengthHistogram(r, 0);
		if(r.mate!=null){addToLengthHistogram(r.mate, 1);}
	}
	
	private void addToLengthHistogram(final Read r, int pairnum){
		if(r==null || r.bases==null){return;}
		int x=Tools.min(r.length(), MAXLENGTHLEN); 
		lengthHist.increment(x, 1);
	}
	
	public void addToGCHistogram(final Read r1){
		if(r1==null){return;}
		final Read r2=r1.mate;
		final int len1=r1.length(), len2=r1.mateLength();
		
		final float gc1=(len1>0 ? r1.gc() : -1);
		final float gc2=(len2>0 ? r2.gc() : -1);
		if(usePairGC){
			final float gc;
			if(r2==null){
				gc=gc1;
			}else{
				gc=(gc1*len1+gc2*len2)/(len1+len2);
			}
			addToGCHistogram(gc, len1+len2);
		}else{
			addToGCHistogram(gc1, len1);
			addToGCHistogram(gc2, len2);
		}
	}
	
	private void addToGCHistogram(final float gc, final int len){
		if(gc<0 || len<1){return;}
		gcHist[Tools.min(GC_BINS, (int)(gc*(GC_BINS+1)))]++;
		gcMaxReadLen=Tools.max(len, gcMaxReadLen);
	}
	
	public void addToIdentityHistogram(final Read r){
		if(r==null){return;}
		addToIdentityHistogram(r, 0);
		if(r.mate!=null){addToIdentityHistogram(r.mate, 1);}
	}
		
	private void addToIdentityHistogram(final Read r, int pairnum){
		if(r==null || r.bases==null || r.length()<1 || !r.mapped() || r.match==null/* || r.discarded()*/){return;}
		float id=r.identity();
		idHist[(int)(id*ID_BINS)]++;
		idBaseHist[(int)(id*ID_BINS)]+=r.length();
		idMaxReadLen=Tools.max(r.length(), idMaxReadLen);
	}
	
	public void addToTimeHistogram(final Read r){
		if(r==null){return;}
		addToTimeHistogram(r, 0);//Time for pairs is the same.
	}
	
	private void addToTimeHistogram(final Read r, int pairnum){
		if(r==null){return;}
		assert(r.obj!=null && r.obj.getClass()==Long.class);
		int x=(int)Tools.min(((Long)r.obj).longValue(), MAXTIMELEN); 
		timeHist.increment(x, 1);
	}
	
	public void addToIndelHistogram(final Read r){
		if(r==null){return;}
		addToIndelHistogram(r, 0);
		if(r.mate!=null){addToIndelHistogram(r.mate, 1);}
	}
		
	private void addToIndelHistogram(final Read r, int pairnum){
		if(r==null || r.bases==null || r.length()<1 || !r.mapped() || r.match==null/* || r.discarded()*/){return;}
		final byte[] bases=r.bases, match=r.match;
		final int limit=Tools.min(bases.length, MAXLEN);
		
		int rpos=0;
		int streak=0;
		byte lastm='A';
		for(int mpos=0; mpos<match.length && rpos<limit; mpos++){
			byte m=match[mpos];
			
			{	
				if(lastm!=m){
					if(lastm=='D'){
						streak=Tools.min(streak, MAXDELLEN2);
						if(streak<MAXDELLEN){delHist.increment(streak, 1);}
						delHist2.increment(streak/100, 1);
					}else if(lastm=='I'){
						streak=Tools.min(streak, MAXINSLEN);
						insHist.increment(streak, 1);
					}
					streak=0;
				}
			}
			streak++;
			rpos++;
			lastm=m;
		}
		if(lastm=='D'){
			streak=Tools.min(streak, MAXDELLEN2);
			if(streak<MAXDELLEN){delHist.increment(streak, 1);}
			delHist2.increment(streak/100, 1);
		}else if(lastm=='I'){
			streak=Tools.min(streak, MAXINSLEN);
			insHist.increment(streak, 1);
		}
	}
	
	public void addToMatchHistogram(final Read r){
		if(r==null){return;}
		addToMatchHistogram2(r);
		if(r.mate!=null){addToMatchHistogram2(r.mate);}
	}
	
	private void addToMatchHistogram2(final Read r){
		if(r==null || r.bases==null || r.length()<1 || !r.mapped() || r.match==null/* || r.discarded()*/){return;}
		int pairnum=r.pairnum();
		if(r.obj!=null){
			if(r.obj.getClass()==SamLine.class){
				pairnum=((SamLine)r.obj).pairnum();
			}
		}
		if(pairnum==1){read2Count++;}
		final byte[] bases=r.bases;
		final int limit=Tools.min(bases.length, MAXLEN);
		final long[] ms=matchSum[pairnum], ds=delSum[pairnum], is=insSum[pairnum],
				ss=subSum[pairnum], ns=nSum[pairnum], cs=clipSum[pairnum], os=otherSum[pairnum];
		
		byte[] match=r.match;
		if(r.shortmatch() && match!=null){match=Read.toLongMatchString(match);}
		
		if(match==null){
			for(int i=0; i<limit; i++){
				byte b=bases[i];
				if(b=='N'){ns[i]++;}
				else{os[i]++;}
			}
		}else{
			final boolean plus=(r.strand()==0);
			int rpos=0;
			byte lastm='A';
			for(int mpos=0; mpos<match.length && rpos<limit; mpos++){
				byte b=bases[rpos];//bases[plus ? rpos : bases.length-rpos-1];
				byte m=match[plus ? mpos : match.length-mpos-1];//match[mpos];
				if(b=='N'){
					if(m=='D'){
						if(lastm!=m){ds[rpos]++;}
						rpos--;
					}else{ns[rpos]++;}
				}else{
					if(m=='m'){
						ms[rpos]++;
					}else if(m=='S'){
						ss[rpos]++;
					}else if(m=='I'){
						is[rpos]++;
					}else if(m=='N' || m=='V'){
//						assert(false) : "\n"+r+"\n"+new String(Data.getChromosome(r.chrom).getBytes(r.start, r.stop))+"\nrpos="+rpos+", mpos="+mpos;
						os[rpos]++;
					}else if(m=='C'){
//						assert(false) : r;
						cs[rpos]++;
					}else if(m=='D'){
						if(lastm!=m){ds[rpos]++;}
						rpos--;
					}else{
						os[rpos]++;
						assert(false) : "For read "+r.numericID+", unknown symbol in match string: ASCII "+m+" = "+(char)m;
					}
				}
				rpos++;
				lastm=m;
			}
		}
	}
	
	public void addToInsertHistogram(final Read r, boolean ignoreMappingStrand){
		if(verbose){
			System.err.print(r.numericID);
			if(r==null || r.mate==null || !r.mapped() || !r.mate.mapped() || !r.paired()){
				System.err.println("\n");
			}else{
				System.err.println("\t"+r.strand()+"\t"+r.insertSizeMapped(ignoreMappingStrand)+"\t"+r.mate.insertSizeMapped(ignoreMappingStrand));
			}
		}
		if(r==null || r.mate==null || !r.mapped() || !r.mate.mapped() || !r.paired()){
			unpairedCount++;
			return;
		}
		int x=Tools.min(MAXINSERTLEN, r.insertSizeMapped(ignoreMappingStrand));
		if(x>0){
			insertHist.increment(x, 1);
			pairedCount++;
		}else{
			unpairedCount++;
		}
//		assert(x!=1) : "\n"+r+"\n\n"+r.mate+"\n";
//		System.out.println("Incrementing "+x);
	}
	
	public void addToInsertHistogram(final SamLine r1){
		int x=r1.tlen;
		x=Tools.min(MAXINSERTLEN, x);
		if(x>0){
			insertHist.increment(x, 1);
			pairedCount++;
		}else{
			unpairedCount++;
		}
	}
	
	public void addToInsertHistogram(final SamLine r1, final SamLine r2){
		if(r1==null){return;}
		int x=insertSizeMapped(r1, r2, REQUIRE_PROPER_PAIR);
		if(verbose){
			System.err.println(r1.qname+"\t"+x);
		}
		x=Tools.min(MAXINSERTLEN, x);
		if(x>0){
			insertHist.increment(x, 1);
			pairedCount++;
		}else{
			unpairedCount++;
		}
	}
	
	/** This is untested and only gives approximate answers when overlapping reads contain indels.
	 * It may give incorrect answers for same-strange pairs that are shorter than read length. 
	 * It might give negative answers but that would be a bug. */
	public static int insertSizeMapped(SamLine r1, SamLine r2, boolean requireProperPair){
		if(r2==null){return r1.length();}
		if(!r1.mapped() || !r2.mapped() || !r1.pairedOnSameChrom() || (requireProperPair && !r1.properPair())){
			return -1;
		}
		
		int a1=r1.start(true, false);
		int a2=r2.start(true, false);
		
		if(r1.strand()!=r2.strand()){
			if(r1.strand()==1){return insertSizeMapped(r2, r1, requireProperPair);}
		}else if(a1>a2){
			return insertSizeMapped(r2, r1, requireProperPair);
		}
		
		int b1=r1.stop(a1, true, false);
		int b2=r2.stop(a2, true, false);

		int clen1=r1.calcCigarLength(true, false);
		int clen2=r2.calcCigarLength(true, false);
		
		int mlen1=b1-a1+1;
		int mlen2=b2-a2+1;
		
		int dif1=mlen1-clen1;
		int dif2=mlen2-clen2;
		
		int mlen12=b2-a1+1;
		
		if(Tools.overlap(a1, b1, a2, b2)){//hard case
			return mlen12-Tools.max(dif1, dif2); //Approximate
		}else{//easy case
			return mlen12-dif1-dif2;
		}
	}
	
	public void addToBaseHistogram(final Read r){
		addToBaseHistogram2(r);
		if(r.mate!=null){addToBaseHistogram2(r.mate);}
	}
	
	public void addToBaseHistogram2(final Read r){
		if(r==null || r.bases==null){return;}
		int pairnum=r.pairnum();
		if(r.obj!=null){
			if(r.obj.getClass()==SamLine.class){
				pairnum=((SamLine)r.obj).pairnum();
			}
		}
		if(pairnum==1){read2Count++;}
		final byte[] bases=r.bases;
		final LongList[] lists=baseHist[pairnum];
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b]+1;
			lists[x].increment(i, 1);
		}
	}
	
	public static boolean testFiles(boolean allowDuplicates){
		return Tools.testOutputFiles(overwrite, append, allowDuplicates, 
				AVG_QUAL_HIST_FILE, QUAL_HIST_FILE, BQUAL_HIST_FILE, BQUAL_HIST_OVERALL_FILE, QUAL_COUNT_HIST_FILE, 
				MATCH_HIST_FILE, INSERT_HIST_FILE, BASE_HIST_FILE, QUAL_ACCURACY_FILE, INDEL_HIST_FILE, ERROR_HIST_FILE, LENGTH_HIST_FILE, 
				GC_HIST_FILE, IDENTITY_HIST_FILE, TIME_HIST_FILE);
	}
	
	public static boolean writeAll(){
		if(collectingStats()){
			ReadStats rs=mergeAll();
			boolean paired=rs.read2Count>0;
			
			if(AVG_QUAL_HIST_FILE!=null){rs.writeAverageQualityToFile(AVG_QUAL_HIST_FILE, paired);}
			if(QUAL_HIST_FILE!=null){rs.writeQualityToFile(QUAL_HIST_FILE, paired);}
			if(BQUAL_HIST_FILE!=null){rs.writeBQualityToFile(BQUAL_HIST_FILE, paired);}
			if(BQUAL_HIST_OVERALL_FILE!=null){rs.writeBQualityOverallToFile(BQUAL_HIST_OVERALL_FILE);}
			if(QUAL_COUNT_HIST_FILE!=null){rs.writeQCountToFile(QUAL_COUNT_HIST_FILE, paired);}
			if(MATCH_HIST_FILE!=null){rs.writeMatchToFile(MATCH_HIST_FILE, paired);}
			if(INSERT_HIST_FILE!=null){rs.writeInsertToFile(INSERT_HIST_FILE);}
			if(BASE_HIST_FILE!=null){rs.writeBaseContentToFile(BASE_HIST_FILE, paired);}
			if(QUAL_ACCURACY_FILE!=null){rs.writeQualityAccuracyToFile(QUAL_ACCURACY_FILE);}
			
			if(INDEL_HIST_FILE!=null){rs.writeIndelToFile(INDEL_HIST_FILE);}
			if(ERROR_HIST_FILE!=null){rs.writeErrorToFile(ERROR_HIST_FILE);}
			if(LENGTH_HIST_FILE!=null){rs.writeLengthToFile(LENGTH_HIST_FILE);}
			if(GC_HIST_FILE!=null){rs.writeGCToFile(GC_HIST_FILE, true);}
			if(IDENTITY_HIST_FILE!=null){rs.writeIdentityToFile(IDENTITY_HIST_FILE, true);}
			if(TIME_HIST_FILE!=null){rs.writeTimeToFile(TIME_HIST_FILE);}
			
			return rs.errorState;
		}
		return false;
	}
	
	public void writeAverageQualityToFile(String fname, boolean writePaired){
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, append, false);
		tsw.start();
		tsw.print("#Quality\tcount1\tfraction1"+(writePaired ? "\tcount2\tfraction2" : "")+"\n");

		long sum1=Tools.sum(aqualArray[0]);
		long sum2=Tools.sum(aqualArray[1]);
		double mult1=1.0/Tools.max(1, sum1);
		double mult2=1.0/Tools.max(1, sum2);
		
		long y=sum1+sum2;
		for(int i=0; i<aqualArray[0].length; i++){
			long x1=aqualArray[0][i];
			long x2=aqualArray[1][i];
			y-=x1;
			y-=x2;
			tsw.print(String.format("%d\t%d\t%.5f", i, x1, x1*mult1));
			if(writePaired){
				tsw.print(String.format("\t%d\t%.5f", x2, x2*mult2));
			}
			tsw.print("\n");
			if(y<=0){break;}
		}
		tsw.poison();
		tsw.waitForFinish();
		errorState|=tsw.errorState;
	}
	
	public void writeQCountToFile(String fname, boolean writePaired){
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, append, false);
		tsw.start();
		tsw.print("#Quality\tcount1\tfraction1"+(writePaired ? "\tcount2\tfraction2" : "")+"\n");

		long sum1=Tools.sum(qcountHist[0]);
		long sum2=Tools.sum(qcountHist[1]);
		double mult1=1.0/Tools.max(1, sum1);
		double mult2=1.0/Tools.max(1, sum2);
		
		long y=sum1+sum2;
		for(int i=0; i<qcountHist[0].length; i++){
			long x1=qcountHist[0][i];
			long x2=qcountHist[1][i];
			y-=x1;
			y-=x2;
			tsw.print(String.format("%d\t%d\t%.5f", i, x1, x1*mult1));
			if(writePaired){
				tsw.print(String.format("\t%d\t%.5f", x2, x2*mult2));
			}
			tsw.print("\n");
			if(y<=0){break;}
		}
		tsw.poison();
		tsw.waitForFinish();
		errorState|=tsw.errorState;
	}
	
	public void writeQualityToFile(String fname, boolean writePaired){
		StringBuilder sb=new StringBuilder();
		final boolean measure=(matchSum!=null);
		if(measure){
			if(writePaired){
				sb.append("#BaseNum\tRead1_linear\tRead1_log\tRead1_measured\tRead2_linear\tRead2_log\tRead2_measured\n");
			}else{
				sb.append("#BaseNum\tRead1_linear\tRead1_log\tRead1_measured\n");
			}
		}else{
			sb.append("#BaseNum\tRead1_linear\tRead1_log"+(writePaired ? "\tRead2_linear\tRead2_log" : "")+"\n");
		}
		
		final long[] qs1=qualSum[0], qs2=qualSum[1], ql1=qualLength[0], ql2=qualLength[1];
		final double[] qsd1=qualSumDouble[0], qsd2=qualSumDouble[1];
		
		for(int i=MAXLEN-2; i>=0; i--){
			ql1[i]+=ql1[i+1];
			ql2[i]+=ql2[i+1];
		}

		double div1sum=0;
		double div2sum=0;
		double deviation1sum=0;
		double deviation2sum=0;
		
		if(writePaired){
			for(int i=0; i<MAXLEN && (ql1[i]>0 || ql2[i]>0); i++){
				final int a=i+1;
				double blin, clin, blog, clog;
				final double div1=(double)Tools.max(1, ql1[i]);
				final double div2=(double)Tools.max(1, ql2[i]);
				
				blin=qs1[i]/div1;
				clin=qs2[i]/div2;
				blog=qsd1[i]/div1;
				clog=qsd2[i]/div2;
				blog=QualityTools.probErrorToPhredDouble(blog);
				clog=QualityTools.probErrorToPhredDouble(clog);
				if(measure){
					double bcalc=calcQualityAtPosition(i, 0);
					double ccalc=calcQualityAtPosition(i, 1);

					div1sum+=div1;
					div2sum+=div2;
					deviation1sum+=Math.abs(blog-bcalc)*div1;
					deviation2sum+=Math.abs(clog-ccalc)*div2;
					
					sb.append(String.format("%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", a, blin, blog, bcalc, clin, clog, ccalc));
				}else{
					sb.append(String.format("%d\t%.3f\t%.3f\t%.3f\t%.3f\n", a, blin, blog, clin, clog));
				}
			}
		}else{
			for(int i=0; i<MAXLEN && ql1[i]>0; i++){
				final int a=i+1;
				double blin, blog;
				final double div1=(double)Tools.max(1, ql1[i]);
				
				blin=qs1[i]/div1;
				blog=qsd1[i]/div1;
				blog=QualityTools.probErrorToPhredDouble(blog);
				if(measure){
					double bcalc=calcQualityAtPosition(i, 0);

					div1sum+=div1;
					deviation1sum+=Math.abs(blog-bcalc)*div1;
					
					sb.append(String.format("%d\t%.3f\t%.3f\t%.3f\n", a, blin, blog, bcalc));
				}else{
					sb.append(String.format("%d\t%.3f\t%.3f\n", a, blin, blog));
				}
			}
		}
		
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, append, false);
		tsw.start();
		if(measure){
			if(writePaired){
				tsw.print(String.format("#Deviation\t%.4f\t%.4f\n", deviation1sum/div1sum, deviation2sum/div2sum));
			}else{
				tsw.print(String.format("#Deviation\t%.4f\n", deviation1sum/div1sum));
			}
		}
		tsw.print(sb);
		tsw.poison();
		tsw.waitForFinish();
		errorState|=tsw.errorState;
	}
	
	private double calcQualityAtPosition(int pos, int pairnum){
		boolean includeNs=true;
		long m=matchSum[pairnum][pos];
		long d=delSum[pairnum][pos]; //left-adjacent deletion
		long d2=delSum[pairnum][Tools.min(pos, delSum[pairnum].length-1)]; //right-adjacent deletion
		long i=insSum[pairnum][pos];
		long s=subSum[pairnum][pos];
		long n=(includeNs ? nSum[pairnum][pos] : 0); //This only tracks no-calls, not no-refs.
		long good=Tools.max(0, m*2-d-d2+n/4);
		long total=Tools.max(0, m*2+i*2+s*2+n); //not d
		long bad=total-good;
		if(total<1){return 0;}
		double error=bad/(double)total;
		return QualityTools.probErrorToPhredDouble(error);
	}
	
	public void writeBQualityOverallToFile(String fname){
		final long[] cp30=Arrays.copyOf(bqualHistOverall, bqualHistOverall.length);
		for(int i=0; i<30; i++){cp30[i]=0;}
		
		final long sum=Tools.sum(bqualHistOverall);
		final long median=Tools.percentile(bqualHistOverall, 0.5);
		final double mean=Tools.averageHistogram(bqualHistOverall);
		final double stdev=Tools.standardDeviationHistogram(bqualHistOverall);
		final double mean30=Tools.averageHistogram(cp30);
		final double stdev30=Tools.standardDeviationHistogram(cp30);
		final double mult=1.0/Tools.max(1, sum);
		long y=sum;
		
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, append, false);
		tsw.start();
		tsw.print("#Median\t"+median+"\n");
		tsw.print("#Mean\t"+String.format("%.3f", mean)+"\n");
		tsw.print("#STDev\t"+String.format("%.3f", stdev)+"\n");
		tsw.print("#Mean_30\t"+String.format("%.3f", mean30)+"\n");
		tsw.print("#STDev_30\t"+String.format("%.3f", stdev30)+"\n");
		tsw.print("#Quality\tbases\tfraction\n");
		
		for(int i=0; i<bqualHistOverall.length; i++){
			long x=bqualHistOverall[i];
			y-=x;
			tsw.print(String.format("%d\t%d\t%.5f", i, x, x*mult));
			tsw.print("\n");
			if(y<=0){break;}
		}
		tsw.poison();
		tsw.waitForFinish();
		errorState|=tsw.errorState;
	}
	
	public void writeBQualityToFile(String fname, boolean writePaired){
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, append, false);
		tsw.start();
		tsw.print("#BaseNum\tcount_1\tmin_1\tmax_1\tmean_1\tQ1_1\tmed_1\tQ3_1\tLW_1\tRW_1");
		if(writePaired){tsw.print("\tcount_2\tmin_2\tmax_2\tmean_2\tQ1_2\tmed_2\tQ3_2\tLW_2\tRW_2");}
		tsw.print("\n");

		for(int i=0; i<MAXLEN; i++){
			final long[] a1=bqualHist[0][i], a2=bqualHist[1][i];
			final long sum1=Tools.sum(a1), sum2=Tools.sum(a2);
			if(sum1<1 && sum2<1){break;}
			
			{
				final long a[]=a1, sum=sum1;
				
				final long weightedSum=Tools.sumHistogram(a);
				final long med=Tools.median(a), min=Tools.minHistogram(a), max=Tools.maxHistogram(a);
				final long firstQuart=Tools.percentile(a, 0.25);
				final long thirdQuart=Tools.percentile(a, 0.75);
				final long leftWhisker=Tools.percentile(a, 0.02);
				final long rightWhisker=Tools.percentile(a, 0.98);
				final double mean=weightedSum*1.0/Tools.max(sum, 0);
				tsw.print(String.format("%d\t%d\t%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d", i, sum, min, max, mean, firstQuart, med, thirdQuart, leftWhisker, rightWhisker));
			}

			if(writePaired){
				final long a[]=a2, sum=sum2;
				
				final long weightedSum=Tools.sumHistogram(a);
				final long med=Tools.median(a), min=Tools.minHistogram(a), max=Tools.maxHistogram(a);
				final long firstQuart=Tools.percentile(a, 0.25);
				final long thirdQuart=Tools.percentile(a, 0.75);
				final long leftWhisker=Tools.percentile(a, 0.02);
				final long rightWhisker=Tools.percentile(a, 0.98);
				final double mean=weightedSum*1.0/Tools.max(sum, 0);
				tsw.print(String.format("\t%d\t%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d", sum, min, max, mean, firstQuart, med, thirdQuart, leftWhisker, rightWhisker));
			}
			tsw.print("\n");
		}
		tsw.poison();
		tsw.waitForFinish();
		errorState|=tsw.errorState;
	}
	
	public void writeQualityAccuracyToFile(String fname){
		
		int max=qualMatch.length;
		for(int i=max-1; i>=0; i--){
			if(qualMatch[i]+qualSub[i]+qualIns[i]+qualDel[i]>0){break;}
			max=i;
		}

		double devsum=0;
		double devsumSub=0;
		long observations=0;
		for(int i=0; i<max; i++){
			long qm=qualMatch[i]*2;
			long qs=qualSub[i]*2;
			long qi=qualIns[i]*2;
			long qd=qualDel[i];
			
			double phred=-1;
			double phredSub=-1;
			
			long sum=qm+qs+qi+qd;
			if(sum>0){
				double mult=1.0/sum;
				double subRate=(qs)*mult;
				double errorRate=(qs+qi+qd)*mult;
				
				phredSub=QualityTools.probErrorToPhredDouble(subRate);
				phred=QualityTools.probErrorToPhredDouble(errorRate);
				double deviation=phred-i;
				double deviationSub=phredSub-i;
				if(i==Read.MIN_CALLED_QUALITY && deviation<0){deviation=0;}
				else if(i==Read.MAX_CALLED_QUALITY && max==Read.MAX_CALLED_QUALITY+1 && deviation>0){deviation=0;}
				if(i==Read.MIN_CALLED_QUALITY && deviationSub<0){deviationSub=0;}
				else if(i==Read.MAX_CALLED_QUALITY && max==Read.MAX_CALLED_QUALITY+1 && deviationSub>0){deviationSub=0;}
				devsum+=(Math.abs(deviation)*sum);
				devsumSub+=(Math.abs(deviationSub)*sum);
				observations+=sum;
			}
		}

		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, append, false);
		tsw.start();
		tsw.print(String.format("#Deviation\t%.3f\n", devsum/observations));
		tsw.print(String.format("#DeviationSub\t%.3f\n", devsumSub/observations));
		tsw.print("#Quality\tMatch\tSub\tIns\tDel\tTrueQuality\tTrueQualitySub\n");
		for(int i=0; i<max; i++){
			long qm=qualMatch[i]*2;
			long qs=qualSub[i]*2;
			long qi=qualIns[i]*2;
			long qd=qualDel[i];
			
			double phred=-1;
			double phredSub=-1;
			
			long sum=qm+qs+qi+qd;
			if(sum>0){
				double mult=1.0/sum;
				double subRate=(qs)*mult;
				double errorRate=(qs+qi+qd)*mult;
				
				phredSub=QualityTools.probErrorToPhredDouble(subRate);
				phred=QualityTools.probErrorToPhredDouble(errorRate);
				
//				System.err.println("sub: "+qs+"/"+sum+" -> "+subRate+" -> "+phredSub);
			}
			
			tsw.print(i+"\t"+qm+"\t"+qs+"\t"+qi+"\t"+qd);
			tsw.print(phred>=0 ? String.format("\t%.2f", phred) : "\t");
			tsw.print(phredSub>=0 ? String.format("\t%.2f\n", phredSub) : "\t\n");
			
//			System.err.println(qm+"\t"+qs+"\t"+qi+"\t"+qd);
		}
		
		tsw.poison();
		tsw.waitForFinish();
		errorState|=tsw.errorState;
	}
	
	public void writeMatchToFile(String fname, boolean writePaired){
		if(!writePaired){
			writeMatchToFileUnpaired(fname);
			return;
		}
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, false);
		tsw.start();
		tsw.print("#BaseNum\tMatch1\tSub1\tDel1\tIns1\tN1\tOther1\tMatch2\tSub2\tDel2\tIns2\tN2\tOther2\n");
		
		final long[] ms1=matchSum[0], ds1=delSum[0], is1=insSum[0],
				ss1=subSum[0], ns1=nSum[0], cs1=clipSum[0], os1=otherSum[0];
		final long[] ms2=matchSum[1], ds2=delSum[1], is2=insSum[1],
				ss2=subSum[1], ns2=nSum[1], cs2=clipSum[1], os2=otherSum[1];
		
		for(int i=0; i<MAXLEN; i++){
			int a=i+1;
			long sum1=ms1[i]+is1[i]+ss1[i]+ns1[i]+cs1[i]+os1[i]; //no deletions
			long sum2=ms2[i]+is2[i]+ss2[i]+ns2[i]+cs2[i]+os2[i]; //no deletions
			if(sum1==0 && sum2==0){break;}
			double inv1=1.0/(double)Tools.max(1, sum1);
			double inv2=1.0/(double)Tools.max(1, sum2);

			tsw.print(String.format("%d", a));
			tsw.print(String.format("\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f", 
					ms1[i]*inv1, ss1[i]*inv1, ds1[i]*inv1, is1[i]*inv1, ns1[i]*inv1, (os1[i]+cs1[i])*inv1));
			tsw.print(String.format("\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f", 
					ms2[i]*inv2, ss2[i]*inv2, ds2[i]*inv2, is2[i]*inv2, ns2[i]*inv2, (os2[i]+cs2[i])*inv2)
//					+", "+ms2[i]+", "+is2[i]+", "+ss2[i]+", "+ns2[i]+", "+cs2[i]+", "+os2[i]
					);
			tsw.print("\n");
		}
		tsw.poison();
		tsw.waitForFinish();
		errorState|=tsw.errorState;
	}
	
	public void writeMatchToFileUnpaired(String fname){
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, false);
		tsw.start();
		tsw.print("#BaseNum\tMatch1\tSub1\tDel1\tIns1\tN1\tOther1\n");
		
		final long[] ms1=matchSum[0], ds1=delSum[0], is1=insSum[0],
				ss1=subSum[0], ns1=nSum[0], cs1=clipSum[0], os1=otherSum[0];
		
		for(int i=0; i<MAXLEN; i++){
			int a=i+1;
			long sum1=ms1[i]+is1[i]+ss1[i]+ns1[i]+cs1[i]+os1[i]; //no deletions
			if(sum1==0){break;}
			double inv1=1.0/(double)Tools.max(1, sum1);

			tsw.print(String.format("%d", a));
			tsw.print(String.format("\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f", 
					ms1[i]*inv1, ss1[i]*inv1, ds1[i]*inv1, is1[i]*inv1, ns1[i]*inv1, (os1[i]+cs1[i])*inv1)
//					+", "+ms1[i]+", "+is1[i]+", "+ss1[i]+", "+ns1[i]+", "+cs1[i]+", "+os1[i]
					);
			tsw.print("\n");
		}
		tsw.poison();
		tsw.waitForFinish();
		errorState|=tsw.errorState;
	}
	
	public void writeInsertToFile(String fname){
		StringBuilder sb=new StringBuilder();
		sb.append("#Mean\t"+String.format("%.3f", Tools.averageHistogram(insertHist.array))+"\n");
		sb.append("#Median\t"+Tools.percentile(insertHist.array, 0.5)+"\n");
		sb.append("#Mode\t"+Tools.calcMode(insertHist.array)+"\n");
		sb.append("#STDev\t"+String.format("%.3f", Tools.standardDeviationHistogram(insertHist.array))+"\n");
		double percent=pairedCount*100.0/(pairedCount+unpairedCount);
		sb.append("#PercentOfPairs\t"+String.format("%.3f", percent)+"\n");
//		sb.append("#PercentOfPairs\t"+String.format("%.3f", matedPercent)+"\n");
		sb.append("#InsertSize\tCount\n");
		writeHistogramToFile(fname, sb.toString(), insertHist, !skipZeroInsertCount);
	}
	
	public void writeBaseContentToFile(String fname, boolean paired){
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, false);
		tsw.start();
//		if(paired){
			tsw.print("#Pos\tA\tC\tG\tT\tN\n");
//		}
		
		LongList[] lists;
		
		int max=writeBaseContentToFile2(tsw, baseHist[0], 0);
		if(paired){
			writeBaseContentToFile2(tsw, baseHist[1], max);
		}
		
		tsw.poison();
		tsw.waitForFinish();
		errorState|=tsw.errorState;
	}
	
	private static int writeBaseContentToFile2(TextStreamWriter tsw, LongList[] lists, int offset){
		int max=0;
		StringBuilder sb=new StringBuilder(100);
		for(LongList ll : lists){max=Tools.max(max, ll.size);}
		for(int i=0; i<max; i++){
			long a=lists[1].get(i);
			long c=lists[2].get(i);
			long g=lists[3].get(i);
			long t=lists[4].get(i);
			long n=lists[0].get(i);
			double mult=1.0/(a+c+g+t+n);

			sb.setLength(0);
			sb.append(i+offset).append('\t');
			sb.append(String.format("%.5f\t", a*mult));
			sb.append(String.format("%.5f\t", c*mult));
			sb.append(String.format("%.5f\t", g*mult));
			sb.append(String.format("%.5f\t", t*mult));
			sb.append(String.format("%.5f\n", n*mult));
			
			tsw.print(sb.toString());
		}
		return max;
	}
	
	public void writeIndelToFile(String fname){
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, false);
		tsw.start();
		tsw.print("#Length\tDeletions\tInsertions\n");
		
		int max=Tools.max(insHist.size, delHist.size);

		for(int i=0; i<max; i++){
			long x=delHist.get(i);
			long y=insHist.get(i);
			if(x>0 || y>0 || !skipZeroIndel){
				tsw.print(i+"\t"+x+"\t"+y+"\n");
			}
		}
		
		//TODO: Disabled because it was irritating when graphing.  Should write to a different file.
//		tsw.print("#Length_bin\tDeletions\n");
//		max=delHist2.size;
//		for(int i=0; i<max; i++){
//			long x=delHist2.get(i);
//			if(x>0 || !skipZeroIndel){
//				tsw.print((i*DEL_BIN)+"\t"+x+"\n");
//			}
//		}
		
		tsw.poison();
		tsw.waitForFinish();
		errorState|=tsw.errorState;
	}
	
	public void writeErrorToFile(String fname){
		writeHistogramToFile(fname, "#Errors\tCount\n", errorHist, false);
	}
	
	public void writeLengthToFile(String fname){
		writeHistogramToFile(fname, "#Length\tCount\n", lengthHist, false);
	}
	
	public void writeTimeToFile(String fname){
		writeHistogramToFile(fname, "#Time\tCount\n", timeHist, false);
	}
	
	public void writeHistogramToFile(String fname, String header, LongList hist, boolean printZeros){
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, false);
		tsw.start();
		tsw.print(header);
		
		int max=hist.size;

		for(int i=0; i<max; i++){
			long x=hist.get(i);
			if(x>0 || printZeros){
				tsw.print(i+"\t"+x+"\n");
			}
		}
		tsw.poison();
		tsw.waitForFinish();
		errorState|=tsw.errorState;
	}
	
	public void writeGCToFile(String fname, boolean printZeros){
		final long[] hist;
		if(GC_BINS_AUTO && gcMaxReadLen+1<gcHist.length){
			hist=Tools.downsample(gcHist, gcMaxReadLen+1);
		}else{
			hist=gcHist;
		}
		final int bins=hist.length;
		final double gcMult=100.0/Tools.max(1, bins-1);
		final long total=Tools.sum(hist);
		final long max=Tools.max(hist);
		final double countsPerX=Tools.max(1, ((max*1000.0)/40));
		final double fractionMult=1.0/Tools.max(1, total);
		long sum=0;
		
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, false);
		tsw.start();
		tsw.print("#Mean\t"+String.format("%.3f", Tools.averageHistogram(hist)*gcMult)+"\n");
		tsw.print("#Median\t"+String.format("%.3f", Tools.percentile(hist, 0.5)*gcMult)+"\n");
		tsw.print("#Mode\t"+String.format("%.3f", Tools.calcMode(hist)*gcMult)+"\n");
		tsw.print("#STDev\t"+String.format("%.3f", Tools.standardDeviationHistogram(hist)*gcMult)+"\n");
		if(GC_PLOT_X){
			tsw.print("#GC\tCount\tCumulative\tPlot\n");
		}else{
			tsw.print("#GC\tCount\n");
		}
		
		for(int i=0; i<bins; i++){
			long x=hist[i];
			sum+=x;
			if(x>0 || printZeros){
				//This assumes GC_BINS==100
//				tsw.print(i+"\t"+x+"\n");
				if(GC_PLOT_X){
					tsw.print(String.format("%.1f\t%d\t", i*gcMult, x));
					
					tsw.print(String.format("%.3f\t", sum*fractionMult));
					
					int len=(int)((x*1000)/countsPerX);
					for(int j=0; j<len; j++){tsw.print("X");}
					if(len<1 && x>0){
						if((x*1000f)/countsPerX>0.1f){tsw.print("x");}
						else{tsw.print(".");}
					} 
					
					tsw.print("\n");
				}else{
					tsw.print(String.format("%.1f\t%d\n", i*gcMult, x));
				}
			}
		}
		tsw.poison();
		tsw.waitForFinish();
		errorState|=tsw.errorState;
	}
	
	public void writeIdentityToFile(String fname, boolean printZeros){
		final long[] hist, histb;
		if(ID_BINS_AUTO && idMaxReadLen+1<idHist.length){
			hist=Tools.downsample(idHist, idMaxReadLen+1);
			histb=Tools.downsample(idBaseHist, idMaxReadLen+1);
		}else{
			hist=idHist;
			histb=idBaseHist;
		}
		final int max=hist.length;
		final double mult=100.0/(max-1);
		
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, false);
		tsw.start();
		
		
		tsw.print("#Mean_reads\t"+String.format("%.3f", (Tools.averageHistogram(hist)*mult))+"\n");
		tsw.print("#Mean_bases\t"+(String.format("%.3f", Tools.averageHistogram(histb)*mult))+"\n");
		tsw.print("#Median_reads\t"+(int)Math.round(Tools.percentile(hist, 0.5)*mult)+"\n");
		tsw.print("#Median_bases\t"+(int)Math.round(Tools.percentile(histb, 0.5)*mult)+"\n");
		tsw.print("#Mode_reads\t"+(int)Math.round(Tools.calcMode(hist)*mult)+"\n");
		tsw.print("#Mode_bases\t"+(int)Math.round(Tools.calcMode(histb)*mult)+"\n");
		tsw.print("#STDev_reads\t"+String.format("%.3f", (Tools.standardDeviationHistogram(hist)*mult))+"\n");
		tsw.print("#STDev_bases\t"+String.format("%.3f", (Tools.standardDeviationHistogram(histb)*mult))+"\n");
		tsw.print("#Identity\tReads\tBases\n");
		
		for(int i=0; i<max; i++){
			long x=hist[i], x2=histb[i];
			if(x>0 || printZeros){
				tsw.print(String.format("%.1f", i*mult)+"\t"+x+"\t"+x2+"\n");
			}
		}
		tsw.poison();
		tsw.waitForFinish();
		errorState|=tsw.errorState;
	}
	
	public static boolean collectingStats(){
		return COLLECT_QUALITY_STATS || COLLECT_QUALITY_ACCURACY || COLLECT_MATCH_STATS || COLLECT_INSERT_STATS || COLLECT_BASE_STATS
				|| COLLECT_INDEL_STATS || COLLECT_GC_STATS || COLLECT_ERROR_STATS || COLLECT_LENGTH_STATS || COLLECT_IDENTITY_STATS || COLLECT_TIME_STATS;
	}
	
	//Tracks to see if read2s have been encountered, for displaying stats.
	private long read2Count=0;

	public long pairedCount=0;
	public long unpairedCount=0;
	
	public final long[][] aqualArray;
	public final long[][] qualLength;
	public final long[][] qualSum;
	
	public final long[][][] bqualHist;
	public final long[] bqualHistOverall;
	
	public final long[][] qcountHist;
	
	public final double[][] qualSumDouble;
	
	public final long[][] matchSum;
	public final long[][] delSum;
	public final long[][] insSum;
	public final long[][] subSum;
	public final long[][] nSum;
	public final long[][] clipSum;
	public final long[][] otherSum;

	public final long[] qualMatch;
	public final long[] qualSub;
	public final long[] qualIns;
	public final long[] qualDel;
	
	public final long[] gcHist;
	public final long[] idHist;
	public final long[] idBaseHist;
	private int gcMaxReadLen=1;
	private int idMaxReadLen=1;

	public final LongList[][] baseHist;
	
	/** Insert size */
	public final LongList insertHist;
	/** Read length */
	public final LongList lengthHist;
	/** Number errors per read */
	public final LongList errorHist;
	/** Insertion length */
	public final LongList insHist;
	/** Deletion length */
	public final LongList delHist;
	/** Deletion length, binned  */
	public final LongList delHist2;
	/** Time */
	public final LongList timeHist;
	
	public static boolean REQUIRE_PROPER_PAIR=true;
	public static final int MAXLEN=6000;
	public static final int MAXINSERTLEN=40000;
	public static final int MAXLENGTHLEN=80000;
	public static final int MAXTIMELEN=80000;
	public static final int MAXINSLEN=1000;
	public static final int MAXDELLEN=1000;
	public static final int MAXDELLEN2=1000000;
	public static final int DEL_BIN=100;
	public static int GC_BINS=100;
	public static int ID_BINS=100;
	public static boolean ID_BINS_AUTO=false;
	public static boolean GC_BINS_AUTO=false;
	public static boolean GC_PLOT_X=false;
	
	public boolean errorState=false;
	
	public static ReadStats merged=null;
	
//	public static double matedPercent=0;
	
	public static ArrayList<ReadStats> objectList=new ArrayList<ReadStats>();
	public static boolean COLLECT_QUALITY_STATS=false;
	public static boolean COLLECT_QUALITY_ACCURACY=false;
	public static boolean COLLECT_MATCH_STATS=false;
	public static boolean COLLECT_INSERT_STATS=false;
	public static boolean COLLECT_BASE_STATS=false;
	public static boolean COLLECT_INDEL_STATS=false;
	public static boolean COLLECT_GC_STATS=false;
	public static boolean COLLECT_ERROR_STATS=false;
	public static boolean COLLECT_LENGTH_STATS=false;
	public static boolean COLLECT_IDENTITY_STATS=false;
	public static boolean COLLECT_TIME_STATS=false;
	
	public static boolean usePairGC=true;
	
	public static String AVG_QUAL_HIST_FILE=null;
	public static String QUAL_HIST_FILE=null;
	public static String BQUAL_HIST_FILE=null;
	public static String QUAL_COUNT_HIST_FILE=null;
	public static String BQUAL_HIST_OVERALL_FILE=null;
	public static String QUAL_ACCURACY_FILE=null;
	public static String MATCH_HIST_FILE=null;
	public static String INSERT_HIST_FILE=null;
	public static String BASE_HIST_FILE=null;
	public static String INDEL_HIST_FILE=null;
	public static String ERROR_HIST_FILE=null;
	public static String LENGTH_HIST_FILE=null;
	public static String GC_HIST_FILE=null;
	public static String IDENTITY_HIST_FILE=null;
	public static String TIME_HIST_FILE=null;
	
	public static boolean overwrite=false;
	public static boolean append=false;
	public static final boolean verbose=false;

	public static boolean skipZeroInsertCount=true;
	public static boolean skipZeroIndel=true;
	
}

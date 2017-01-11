package jgi;

import java.util.ArrayList;
import java.util.Arrays;

import stream.ConcurrentLegacyReadInputStream;
import stream.RTextInputStream;
import stream.Read;
import stream.SiteScore;
import structures.CoverageArray;
import structures.CoverageArray2;
import structures.ListNum;
import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import dna.Gene;
import dna.Parser;
import fileIO.ReadWrite;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Jul 16, 2012
 *
 */
public class MakeCoverageHistogram {
	
	public static void main(String[] args){
		System.err.println("Executing "+(new Object() { }.getClass().getEnclosingClass().getName())+" "+Arrays.toString(args)+"\n");
		Timer t=new Timer();
		
		Data.GENOME_BUILD=-1;
		for(int i=2; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.isJavaFlag(arg)){
				//jvm argument; do nothing
			}else if(a.equals("genome") || a.equals("build")){
				Data.setGenome(Integer.parseInt(b));
			}else if(a.equals("maxdepth")){
				MAX_DEPTH=Integer.parseInt(b);
			}
		}
//		assert(false) : "MAX_DEPTH="+MAX_DEPTH;
		assert(Data.GENOME_BUILD>-1);
		
		calc(args[0], args[1]);
		t.stop();
		System.out.println("Time: \t"+t);
	}
	
	public static void calc(String fname1, String fname2){
		RTextInputStream rtis=new RTextInputStream(fname1, (fname2==null || fname2.equals("null") ? null : fname2), -1);
		ConcurrentLegacyReadInputStream cris=new ConcurrentLegacyReadInputStream(rtis, -1);
		
		cris.start();
		System.err.println("Started cris");
		boolean paired=cris.paired();
		System.err.println("Paired: "+paired);

		ArrayList<CoverageArray> pcov=new ArrayList<CoverageArray>(8);
		pcov.add(new CoverageArray2(0,1000));
		ArrayList<CoverageArray> cov=new ArrayList<CoverageArray>(8);
		cov.add(new CoverageArray2(0,1000));
		
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(paired==(r.mate!=null));
			}
			
			while(reads!=null && reads.size()>0){
				//System.err.println("reads.size()="+reads.size());
				for(Read r : reads){
					readsProcessed++;
					
//					System.out.println("Processing read "+r.numericID);
					
					if(r!=null){
						if(r.sites!=null){
							
							for(int x=0; x<r.sites.size(); x++){
								SiteScore ss=r.sites.get(x);
								
								if(PROCESS_ALL_SITES || x==0 || ss.semiperfect){
									sitesProcessed++;

									boolean b=false;
									if(ss.perfect || ss.semiperfect){
										b=true;
									}else{//Check for no-refs
										int len=ss.stop-ss.start+1;
										if(len==r.length()){
											b=checkPerfection(ss.start, ss.stop, r.bases, Data.getChromosome(ss.chrom), ss.strand==Gene.MINUS, 0.5f);
										}
									}
									if(b){
										while(pcov.size()<=ss.chrom){
											pcov.add(new CoverageArray2(pcov.size(), 500));
										}
										CoverageArray ca=pcov.get(ss.chrom);
										for(int i=ss.start; i<=ss.stop; i++){
											ca.increment(i);
										}
									}
									
									while(cov.size()<=ss.chrom){
										cov.add(new CoverageArray2(cov.size(), 500));
									}
									CoverageArray ca=cov.get(ss.chrom);
									for(int i=ss.start; i<=ss.stop; i++){
										ca.increment(i);
									}
								}
							}
//							System.out.println(sitesProcessed);
						}
					}
					
					if(r.mate!=null){
						Read r2=r.mate;
						if(r2.sites!=null){
							
							for(int x=0; x<r2.sites.size(); x++){
								SiteScore ss=r2.sites.get(x);

								if(PROCESS_ALL_SITES || x==0 || ss.semiperfect){
									sitesProcessed++;

									boolean b=false;
									if(ss.perfect || ss.semiperfect){
										b=true;
									}else{//Check for no-refs
										int len=ss.stop-ss.start+1;
										if(len==r2.length()){
											b=checkPerfection(ss.start, ss.stop, r2.bases, Data.getChromosome(ss.chrom), ss.strand==Gene.MINUS, 0.5f);
										}
									}
									if(b){
										while(pcov.size()<=ss.chrom){
											pcov.add(new CoverageArray2(pcov.size(), 500));
										}
										CoverageArray ca=pcov.get(ss.chrom);
										for(int i=ss.start; i<=ss.stop; i++){
											ca.increment(i);
										}
									}
									
									while(cov.size()<=ss.chrom){
										cov.add(new CoverageArray2(cov.size(), 500));
									}
									CoverageArray ca=cov.get(ss.chrom);
									for(int i=ss.start; i<=ss.stop; i++){
										ca.increment(i);
									}
								}
							}
						}
					}
					
//					System.out.println(r.toString());
//					assert(r.list!=null);
//					assert(r.list.size()>0);
					
				}
				//System.err.println("returning list");
				cris.returnList(ln.id, ln.list.isEmpty());
				//System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			System.err.println("Finished reading");
			cris.returnList(ln.id, ln.list.isEmpty());
			System.err.println("Returned list");
			ReadWrite.closeStream(cris);
			System.err.println("Closed stream");
			System.err.println("Processed "+readsProcessed+" reads.");
			System.err.println("Processed "+sitesProcessed+" sites.");
		}
		
		int max=MAX_DEPTH;
		long[] hist=new long[max+1];
		long[] phist=new long[max+1];
		double[] histF=new double[max+1];
		double[] phistF=new double[max+1];
		long[] histC=new long[max+1];
		long[] phistC=new long[max+1];
		double[] histCF=new double[max+1];
		double[] phistCF=new double[max+1];
		
		for(int chrom=1; chrom<=Data.numChroms; chrom++){
			ChromosomeArray cha=Data.getChromosome(chrom);
			if(pcov.size()>chrom){
				CoverageArray ca=pcov.get(chrom);
				for(int i=0; i<=cha.maxIndex; i++){
					int x=ca.get(i);
					byte b=cha.get(i);
					if(b!='N'){
						phist[Tools.min(max, x)]++;
					}
				}
			}
		}
		
		for(int chrom=1; chrom<=Data.numChroms; chrom++){
			ChromosomeArray cha=Data.getChromosome(chrom);
			if(cov.size()>chrom){
				CoverageArray ca=cov.get(chrom);
				for(int i=0; i<=cha.maxIndex; i++){
					int x=ca.get(i);
					byte b=cha.get(i);
					if(b!='N'){
						hist[Tools.min(max, x)]++;
					}
				}
			}
		}
		
		phistC[max]=phist[max];
		histC[max]=hist[max];
		for(int i=max; i>0; i--){
			phistC[i-1]=phistC[i]+phist[i-1];
			histC[i-1]=histC[i]+hist[i-1];
		}
		for(int i=0; i<=max; i++){
			phistCF[i]=phistC[i]*100d/phistC[0];
			phistF[i]=phist[i]*100d/phistC[0];
			histCF[i]=histC[i]*100d/histC[0];
			histF[i]=hist[i]*100d/histC[0];
		}
		
		System.out.println("\nTotal Coverage:");
		for(int i=0; i<=max; i++){
			System.out.println(i+"\t"+hist[i]+String.format("\t%.3f%%", histF[i])+"\t"+histC[i]+String.format("\t%.3f%%", histCF[i]));
		}

		
		System.out.println("\nPerfect Coverage:");
		for(int i=0; i<=max; i++){
			System.out.println(i+"\t"+phist[i]+String.format("\t%.3f%%", phistF[i])+"\t"+phistC[i]+String.format("\t%.3f%%", phistCF[i]));
		}
		
	}
	
	private static boolean checkPerfection(int start, int stop, byte[] bases, ChromosomeArray cha, boolean rcomp, float f) {
		
		int noref=0;
		if(rcomp){
			for(int i=0; i<bases.length; i++){
				byte a=AminoAcid.baseToComplementExtended[bases[bases.length-i-1]];
				byte b=cha.get(start+i);
				if(b=='N'){noref++;}
				else if(a!=b){return false;}
			}
		}else{
			for(int i=0; i<bases.length; i++){
				byte a=bases[i];
				byte b=cha.get(start+i);
				if(b=='N'){noref++;}
				else if(a!=b){return false;}
			}
		}
		return bases.length-noref>=f*bases.length;
	}
	
	public static long readsProcessed=0;
	public static long sitesProcessed=0;
	public static boolean PROCESS_ALL_SITES=false;
	public static int MAX_DEPTH=100;
	
}

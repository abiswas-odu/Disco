package jgi;

import java.util.ArrayList;
import java.util.Arrays;

import align2.MSA;
import stream.ByteBuilder;
import stream.ConcurrentReadInputStream;
import stream.Read;
import stream.SamLine;
import stream.SiteScore;
import structures.ListNum;
import dna.AminoAcid;
import dna.Gene;
import dna.Parser;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Oct 6, 2014
 *
 */
public class FindPrimers {

	public static void main(String[] args){
		Timer t=new Timer();
		FindPrimers as=new FindPrimers(args);
		as.process(t);
	}
	
	public FindPrimers(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		String query_=null;
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("primer") || a.equals("query") || a.equals("literal")){
				query_=b;
			}else if(a.equals("msa")){
				msaType=b;
			}else if(a.equals("columns")){
				columns=Integer.parseInt(b);
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			in1=parser.in1;
			out1=parser.out1;
		}
		
		
		if(query_==null){
			querys=rquerys=null;
			maxqlen=0;
		}else{
			int max=0;
			String[] s2=query_.split(",");
			querys=new byte[s2.length][];
			rquerys=new byte[s2.length][];
			for(int i=0; i<s2.length; i++){
				max=Tools.max(max, s2[i].length());
				querys[i]=s2[i].getBytes();
				rquerys[i]=AminoAcid.reverseComplementBases(querys[i]);
			}
			maxqlen=max;
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, true, false, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
	}
	
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			//Distributed version
//			ConcurrentReadInputStream cris0=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
//			cris=new ConcurrentReadInputStreamD(cris0, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
		
		final ByteStreamWriter bsw;
		if(out1!=null){

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			bsw=new ByteStreamWriter(ffout1);
			bsw.start();
		}else{bsw=null;}
		
		
		MSA msa=MSA.makeMSA(maxqlen+3, columns, msaType);
		
		final Read queryRead=new Read(querys[0], null, 1);
		queryRead.sites=new ArrayList<SiteScore>();
		
		long readsProcessed=0;
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(reads!=null && reads.size()>0){
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r=reads.get(idx);
					
					if(r.length()+2>msa.maxColumns){
						msa=MSA.makeMSA(maxqlen+3, r.length()+2+r.length()/2, "MultiStateAligner11ts");
					}
					
					SiteScore ss=null;
					int score1=-999999, score2=-999999;
					
					final int a=0, b=r.length()-1;
					int[] max;
					
					SiteScore ssf=null;
					for(int qnum=0; qnum<querys.length; qnum++){
						final byte[] query=querys[qnum];
						max=msa.fillLimited(query, r.bases, a, b, -9999, null);
						if(max!=null){
							int[] score=msa.score(query, r.bases, a, b, max[0], max[1], max[2], false);
							ss=new SiteScore(1, (byte)0, score[1], score[2], 1, score[0]);
							if(ssf==null || ss.quickScore>ssf.quickScore){
								ss.setSlowScore(ss.quickScore);
								score1=ss.score=ss.quickScore;
								ss.match=msa.traceback(query, r.bases, a, b, score[3], score[4], score[5], false);
								ss.hits=qnum;
								ssf=ss;
							}
						}
					}
					
					SiteScore ssr=null;
					for(int qnum=0; qnum<rquerys.length; qnum++){
						final byte[] rquery=rquerys[qnum];
						max=msa.fillLimited(rquery, r.bases,a, b, -9999, null);
						if(max!=null){
							int[] score=msa.score(rquery, r.bases, a, b, max[0], max[1], max[2], false);
							ss=new SiteScore(1, (byte)1, score[1], score[2], 1, score[0]);
							if(ssr==null || ss.quickScore>ssr.quickScore){
								ss.setSlowScore(ss.quickScore);
								score1=ss.score=ss.quickScore;
								ss.match=msa.traceback(rquery, r.bases, a, b, score[3], score[4], score[5], false);
								ss.hits=qnum;
								ssr=ss;
							}
						}
					}
					
					if(ssf==null && ssr==null){}
					else{
						if(score1>=score2 && ssf!=null){
							bsw.println(toBytes(null, r, ssf));
						}
						if(score2>score1 && ssr!=null){
							bsw.println(toBytes(null, r, ssr));
						}
					}
					
					readsProcessed++;
				}

				cris.returnList(ln.id, ln.list.isEmpty());
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		ArrayList<Read> rosList=new ArrayList<Read>();
		rosList.add(queryRead);
		if(!queryRead.sites.isEmpty()){
			queryRead.setFromTopSite();
			queryRead.setMapped(true);
		}
		
		bsw.poisonAndWait();
		ReadWrite.closeStreams(cris);
		if(verbose){outstream.println("Finished.");}
		
		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+String.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
	}
	

	/*--------------------------------------------------------------*/
	
	private ByteBuilder toBytes(ByteBuilder bb, Read r, SiteScore ss){
		
		final byte[] query=querys[ss.hits], rquery=rquerys[ss.hits];
		
		if(bb==null){bb=new ByteBuilder(80);}
		bb.append("query").append('\t');
		bb.append(makeFlag(ss)).append('\t');
		bb.append(r.id.replace('\t', '_')).append('\t');
		bb.append(ss.start+1).append('\t');
		bb.append(Tools.max(ss.score/query.length, 4)).append('\t');
		String cigar=SamLine.toCigar14(ss.match, ss.start, ss.stop, r.length(), query);
		if(cigar==null){bb.append('*').append('\t');}else{bb.append(cigar).append('\t');}
		bb.append('0').append('\t');
		bb.append('*').append('\t');
		bb.append('0').append('\t');
		
		bb.append(ss.strand()==Gene.MINUS ? rquery : query).append('\t');
		bb.append('*').append('\t');
		
		float f=Read.identity(ss.match);
		bb.append(String.format("YI:f:%.2f", (100*f)));
		
		return bb;
	}
	
	public static int makeFlag(SiteScore ss){
		int flag=0;
		if(ss.strand()==Gene.MINUS){flag|=0x10;}
//		if(r.secondary()){flag|=0x100;}
//		if(r.discarded()){flag|=0x200;}
		return flag;
	}
	
	
	
	/*--------------------------------------------------------------*/
	
	private void printOptions(){
		throw new RuntimeException("printOptions: TODO");
	}
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	private final byte[][] querys, rquerys;
	private final int maxqlen;
	private int columns=2000;
	private String msaType="MultiStateAligner11ts";
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}

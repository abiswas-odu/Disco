package var2;

import java.util.ArrayList;
import java.util.Arrays;

import shared.Tools;
import stream.ByteBuilder;

public class VCFLine {
	
	public VCFLine(byte[] line) {
		int a=0, b=0;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(line);
		scaf=new String(line, a, b-a);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 1: "+new String(line);
		pos=Tools.parseInt(line, a, b);
		b++;
		a=b;

		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 2: "+new String(line);
		id=line[a]=='.' ? DOT : Arrays.copyOfRange(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 3: "+new String(line);
		if(b<=a+1){ref=Var.AL_MAP[line[a]];}
		else{ref=Arrays.copyOfRange(line, a, b);}
		reflen=line[a]=='.' ? 0 : b-a;
		b++;
		a=b;

		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 4: "+new String(line);
		if(b<=a+1){alt=Var.AL_MAP[line[a]];}
		else{alt=Arrays.copyOfRange(line, a, b);}
		b++;
		a=b;

		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 5: "+new String(line);
		qual=Tools.parseFloat(line, a, b);
		b++;
		a=b;

		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 6: "+new String(line);
		filter=Arrays.copyOfRange(line, a, b);;
		b++;
		a=b;

		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 7: "+new String(line);
		info=Arrays.copyOfRange(line, a, b);;
		b++;
		a=b;

		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 8: "+new String(line);
		format=Arrays.copyOfRange(line, a, b);;
		b++;
		a=b;
		
		while(b<line.length){
			while(b<line.length && line[b]!='\t'){b++;}
			if(b<=a){
				break;
			}
			byte[] sample=Arrays.copyOfRange(line, a, b);
			samples.add(sample);
			b++;
			a=b;
		}
	}
	
	public Var toVar(){
		return makeVar(info, alt);
	}
	
	public static Var makeVar(byte[] info, byte[] alt){
		int a=0, b=0;
		
		//SN=0;STA=547693;STO=547694;TYP=SUB;
		assert(Tools.startsWith(info, "SN", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 0: "+new String(info);
		int scaf=Tools.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 1: "+new String(info);
		int start=Tools.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 2: "+new String(info);
		int stop=Tools.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 3: "+new String(info);
		int type=Var.SUB;
		if(Tools.contains(info, SUB, a)){type=Var.SUB;}
		else if(Tools.contains(info, DEL, a)){type=Var.DEL;}
		else if(Tools.contains(info, INS, a)){type=Var.INS;}
		else if(Tools.contains(info, NOCALL, a)){type=Var.NOCALL;}
		else{assert(false) : new String(info);}
		b++;
		a=b;
		
		//R1P=20;R1M=29;R2P=25;R2M=19;
		assert(Tools.startsWith(info, "R1P", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 4: "+new String(info);
		int r1p=Tools.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 5: "+new String(info);
		int r1m=Tools.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 6: "+new String(info);
		int r2p=Tools.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 7: "+new String(info);
		int r2m=Tools.parseInt(info, a, b);
		b++;
		a=b;
		
		//PPC=93;LS=13113;MQS=3975;MQM=44;
		assert(Tools.startsWith(info, "PPC", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 8: "+new String(info);
		int pc=Tools.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 9: "+new String(info);
		long ls=Tools.parseLong(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 10: "+new String(info);
		long mqs=Tools.parseLong(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 11: "+new String(info);
		int mqm=Tools.parseInt(info, a, b);
		b++;
		a=b;
		
		//BQS=3042;BQM=37;EDS=3496;EDM=71;
		assert(Tools.startsWith(info, "BQS", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 12: "+new String(info);
		long bqs=Tools.parseLong(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 13: "+new String(info);
		int bqm=Tools.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 14: "+new String(info);
		long eds=Tools.parseLong(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 15: "+new String(info);
		int edm=Tools.parseInt(info, a, b);
		b++;
		a=b;
		
		//IDS=91828;IDM=992;COV=94;MCOV=48;
		assert(Tools.startsWith(info, "IDS", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 16: "+new String(info);
		long ids=Tools.parseLong(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 17: "+new String(info);
		int idm=Tools.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 18: "+new String(info);
		int cov=Tools.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 19: "+new String(info);
		int mcov=Tools.parseInt(info, a, b);
		b++;
		a=b;
		
		//HMP=2;DP=94;AF=0.989;DP4=1,0,45,48
		assert(Tools.startsWith(info, "HMP", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 20: "+new String(info);
		int hmp=Tools.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 21: "+new String(info);
//		int dp=Tools.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 22: "+new String(info);
//		float af=Tools.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 23: "+new String(info);
//		String dp4=new String(info, a, b-a);
		b++;
		a=b;
		
		if(type==Var.DEL || type==Var.INS){
			if(alt.length<=1){alt=Var.AL_0;}
			else if(alt.length==2){alt=Var.AL_MAP[alt[1]];}
			else{alt=Arrays.copyOfRange(alt, 1, alt.length);}
		}
		
		//HMP=2;DP=94;AF=0.989;DP4=1,0,45,48	GT:DP:AD:AF	1:94:93:0.989

		Var v=new Var(scaf, start, stop, alt);
		v.r1plus=r1p;
		v.r1minus=r1m;
		v.r2plus=r2p;
		v.r2minus=r2m;
		v.properPairCount=pc;
		v.lengthSum=ls;
		v.mapQSum=mqs;
		v.mapQMax=mqm;
		v.baseQSum=bqs;
		v.baseQMax=bqm;
		v.endDistSum=eds;
		v.endDistMax=edm;
		v.idSum=ids;
		v.idMax=idm;
		v.setCoverage(cov, mcov);
//		v.homopolymerCount=hmp; //derived
		
		return v;
	}
	
	public String toString(){
		ByteBuilder bb=new ByteBuilder();
		return toText(bb).toString();
	}
	
	public ByteBuilder toText(ByteBuilder bb){
		bb.append(scaf).append('\t');
		bb.append(pos).append('\t');
		bb.append(id).append('\t');
		bb.append(ref).append('\t');
		bb.append(alt).append('\t');
		bb.append(String.format("%.3f",qual)).append('\t');
		bb.append(filter).append('\t');
		bb.append(info).append('\t');
		bb.append(format);
		for(byte[] sample : samples){
			bb.append('\t').append(sample);
		}
		return bb;
	}
	
	public String scaf;
	public int pos;
	public byte[] id;
	public byte[] ref;
	public int reflen;
	public byte[] alt;
	public float qual;
	public byte[] filter;
	public byte[] info;
	public byte[] format;
	public ArrayList<byte[]> samples=new ArrayList<byte[]>();

	private static final byte[] NOCALL="NOCALL".getBytes();
	private static final byte[] SUB="SUB".getBytes();
	private static final byte[] DEL="DEL".getBytes();
	private static final byte[] INS="INS".getBytes();
	private static final byte[] DOT=".".getBytes();

}

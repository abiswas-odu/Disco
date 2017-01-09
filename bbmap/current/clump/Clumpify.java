package clump;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import dna.Parser;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import hiseq.FlowcellCoordinate;
import jgi.BBMerge;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sort.SortByName;
import stream.FASTQ;
import stream.Read;
import structures.Quantizer;

/**
 * @author Brian Bushnell
 * @date Nov 6, 2015
 *
 */
public class Clumpify {

	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 6);
		BBMerge.changeQuality=Read.CHANGE_QUALITY=false;
		Clumpify cl=new Clumpify(args);
		cl.process(t);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public Clumpify(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		System.err.println("\nClumpify version "+Shared.BBMAP_VERSION_STRING);
		
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		args2=new ArrayList<String>();
		args2.add("in");
		args2.add("out");
		args2.add("groups");
		args2.add("ecco=f");
		args2.add("rename=f");
		args2.add("shortname=f");
		args2.add("unpair=f");
		args2.add("repair=f");
		args2.add("namesort=f");
		args2.add("overwrite=t");
		
		String gString="auto";
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens

			if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("out") || a.equals("out1")){
				out1=b;
			}else if(a.equals("groups") || a.equals("g") || a.equals("sets") || a.equals("ways")){
				gString=b;
			}else if(a.equals("delete")){
				delete=Tools.parseBoolean(b);
			}else if(a.equals("usetmpdir")){
				useTmpdir=Tools.parseBoolean(b);
			}else if(a.equals("ecco")){
				ecco=Tools.parseBoolean(b);
			}else if(a.equals("compresstemp") || a.equals("ct")){
				if(b!=null & b.equalsIgnoreCase("auto")){forceCompressTemp=forceRawTemp=false;}
				else{
					forceCompressTemp=Tools.parseBoolean(b);
					forceRawTemp=!forceCompressTemp;
				}
			}else if(a.equals("tmpdir")){
				Shared.setTmpdir(b);
			}else if(a.equals("rename") || a.equals("addname")){
				addName=Tools.parseBoolean(b);
			}else if(a.equals("shortname") || a.equals("shortnames")){
				shortName=b;
			}else if(a.equals("seed")){
				KmerComparator.defaultSeed=Long.parseLong(b);
			}else if(a.equals("hashes")){
				KmerComparator.setHashes(Integer.parseInt(b));
			}else if(a.equals("passes")){
				passes=Integer.parseInt(b);
				args2.add(arg);
//			}else if(a.equals("k")){
//				k=Integer.parseInt(b);
//				args2.add(arg);
			}else if(a.equals("border")){
				KmerComparator.defaultBorder=Integer.parseInt(b);
			}

			else if(a.equals("unpair")){
				unpair=Tools.parseBoolean(b);
			}else if(a.equals("repair")){
				repair=Tools.parseBoolean(b);
			}else if(a.equals("namesort") || a.equals("sort")){
				namesort=Tools.parseBoolean(b);
			}else if(a.equals("overwrite")){
				overwrite=Tools.parseBoolean(b);
			}else if(a.equals("v1") || a.equals("kmersort1")){
				boolean x=Tools.parseBoolean(b);
				if(x){V2=V3=false;}
			}else if(a.equals("v2") || a.equals("kmersort2")){
				V2=Tools.parseBoolean(b);
				if(V2){V3=false;}
			}else if(a.equals("v3") || a.equals("kmersort3")){
				V3=Tools.parseBoolean(b);
				if(V3){V2=false;}
			}
			
			else if(a.equals("comparesequence")){
				KmerComparator.compareSequence=Tools.parseBoolean(b);
			}else if(a.equals("allowadjacenttiles") || a.equals("spantiles")){
				FlowcellCoordinate.spanTiles=Tools.parseBoolean(b);
			}
			
//			else if(a.equals("repair")){
//				repair=Tools.parseBoolean(b);
//			}else if(a.equals("namesort") || a.equals("sort")){
//				namesort=Tools.parseBoolean(b);
//			}
			
			else if(a.equals("interleaved") || a.equals("int")){
				if("auto".equalsIgnoreCase(b)){FASTQ.FORCE_INTERLEAVED=!(FASTQ.TEST_INTERLEAVED=true);}
				else{
					FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=Tools.parseBoolean(b);
					System.err.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}else if(a.equals("cq") || a.equals("changequality")){
				BBMerge.changeQuality=Read.CHANGE_QUALITY=Tools.parseBoolean(b);
			}else if(a.equals("quantize") || a.equals("quantizesticky")){
				quantizeQuality=Quantizer.parse(arg, a, b);
			}
			
			else if(Clump.parseStatic(arg, a, b)){
				//Do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//Do nothing
			}
			
			else{
				args2.add(arg);
			}
		}
		
		KmerSplit.quantizeQuality=KmerSort.quantizeQuality=quantizeQuality;
		
		Parser.processQuality();
		
		assert(!unpair || !KmerComparator.mergeFirst) : "Unpair and mergefirst may not be used together.";
		
		if(in1==null){
			throw new RuntimeException("\nOne input file is required.\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read to some input files.\n");
		}
		
		autoSetGroups(gString);
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	public void process(Timer t){
		String[] args=args2.toArray(new String[0]);
		args[2]="groups="+groups;
		
		useSharedHeader=(FileFormat.hasSamOrBamExtension(in1) && out1!=null
				&& FileFormat.hasSamOrBamExtension(out1));
		
		if(groups==1){
			args[0]="in="+in1;
			args[1]="out="+out1;
			args[3]="ecco="+ecco;
			args[4]="rename="+addName;
			args[5]="shortname="+shortName;
			args[6]="unpair="+unpair;
			args[7]="repair="+repair;
			args[8]="namesort="+namesort;
			args[9]="ow="+overwrite;
			KmerSort.main(args);
		}else{
			String in=in1, out;
			final int conservativePasses=Clump.conservativeFlag ? passes : Tools.max(1, passes/2);
			if(passes>1){Clump.setConservative(true);}
			long fileMem=-1;
			for(int pass=1; pass<=passes; pass++){
				if(/*passes>1 &&*/ (V2 || V3)){
//					System.err.println("Running pass with fileMem="+fileMem);
					out=(pass==passes ? out1 : getTempFname("clumpify_p"+(pass+1)+"_temp%_"));
					fileMem=runOnePass_v2(args, pass, in, out, fileMem);
//					System.err.println("New fileMem="+fileMem);
				}else{
					out=(pass==passes ? out1 : getTempFname("clumpify_temp_pass"+pass+"_"));
					runOnePass(args, pass, in, out);
				}
				in=out;
				KmerComparator.defaultBorder=Tools.max(0, KmerComparator.defaultBorder-1);
				KmerComparator.defaultSeed++;
				if(pass>=conservativePasses){Clump.setConservative(false);}
			}
		}
		t.stop();
		System.err.println("Total time: \t"+t);
	}
	
	private void runOnePass(String[] args, int pass, String in, String out){
		assert(groups>1);
		if(pass>1){
			ecco=false;
			shortName="f";
			addName=false;
		}

		String temp=getTempFname("clumpify_p"+pass+"_temp%_");
		
		String temp2=temp.replace("%", "FINAL");
		final boolean externalSort=(pass==passes && (repair || namesort));
		
		args[0]="in="+in;
		args[1]="out="+temp;
		args[3]="ecco="+ecco;
		args[4]="addname=f";
		args[5]="shortname="+shortName;
		args[6]="unpair="+unpair;
		args[7]="repair=f";
		args[8]="namesort=f";
		args[9]="ow="+overwrite;
		KmerSplit.maxZipLevel=2;
		KmerSplit.main(args);
		
		FASTQ.DETECT_QUALITY=FASTQ.DETECT_QUALITY_OUT=false;
		FASTQ.ASCII_OFFSET=FASTQ.ASCII_OFFSET_OUT;
		
		args[0]="in="+temp;
		args[1]="out="+(externalSort ? temp2 : out);
		args[3]="ecco=f";
		args[4]="addname="+addName;
		args[5]="shortname=f";
		args[6]="unpair=f";
		args[7]="repair="+(repair && externalSort);
		args[8]="namesort="+(namesort && externalSort);
		args[9]="ow="+overwrite;
		if(unpair){
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		KmerSort.main(args);
		
		if(delete){
			for(int i=0; i<groups; i++){
				new File(temp.replaceFirst("%", ""+i)).delete();
			}
			if(pass>1){new File(in).delete();}
		}
		
		if(externalSort){
			outstream.println();
			SortByName.main(new String[] {"in="+temp2, "out="+out, "ow="+overwrite});
			if(delete){new File(temp2).delete();}
		}
	}
	
	private long runOnePass_v2(String[] args, int pass, String in, String out, long fileMem){
		assert(groups>1);
		if(pass>1){
			ecco=false;
			shortName="f";
			addName=false;
		}
		
		String temp=getTempFname("clumpify_p"+pass+"_temp%_");
		
		String temp2=temp.replace("%", "FINAL");
		String namesorted=temp.replace("%", "namesorted_%");
		final boolean externalSort=(pass==passes && (repair || namesort));
		
		if(pass==1){
			args[0]="in="+in;
			args[1]="out="+temp;
			args[3]="ecco="+ecco;
			args[4]="addname=f";
			args[5]="shortname="+shortName;
			args[6]="unpair="+unpair;
			args[7]="repair=f";
			args[8]="namesort=f";
			args[9]="ow="+overwrite;
			KmerSplit.maxZipLevel=2;
			KmerSplit.main(args);
			fileMem=KmerSplit.lastMemProcessed;
			
			FASTQ.DETECT_QUALITY=FASTQ.DETECT_QUALITY_OUT=false;
			FASTQ.ASCII_OFFSET=FASTQ.ASCII_OFFSET_OUT;
		}
		
		args[0]="in="+(pass==1 ? temp : in);
//		args[1]="out="+(externalSort ? temp2 : out);
		args[1]="out="+(externalSort ? namesorted : out);
		args[3]="ecco=f";
		args[4]="addname="+addName;
		args[5]="shortname=f";
		args[6]="unpair=f";
		args[7]="repair="+(repair && externalSort);
		args[8]="namesort="+(namesort && externalSort);
		args[9]="ow="+overwrite;
		if(unpair){
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		if(externalSort){
			KmerSort2.doHashAndSplit=KmerSort3.doHashAndSplit=false;
		}
		if(V3){
			KmerSort3.main(fileMem, args);
			if(fileMem<1){fileMem=KmerSort3.lastMemProcessed;}
		}else{KmerSort2.main(args);}
		
		if(delete){
			for(int i=0; i<groups; i++){
				new File((pass==1 ? temp : in).replaceFirst("%", ""+i)).delete();
			}
		}
		
		if(externalSort){
			outstream.println();
			
			ArrayList<String> names=new ArrayList<String>();
			for(int i=0; i<groups; i++){
				names.add(namesorted.replaceFirst("%", ""+i));
			}
			ReadWrite.MAX_ZIP_THREADS=Shared.threads();
			
			ReadWrite.USE_PIGZ=true;
			ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 6);
			FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
			FileFormat dest=FileFormat.testOutput(out, FileFormat.FASTQ, null, true, overwrite, false, false);
			SortByName.mergeAndDump(names, dest, null, delete, useSharedHeader);
		}
		
//		if(externalSort){
//			outstream.println();
//			SortByName.main(new String[] {"in="+temp2, "out="+out, "ow="+overwrite});
//			if(delete){new File(temp2).delete();}
//		}
		return fileMem;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** This is called if the program runs with no parameters */
	private void printOptions(){
		throw new RuntimeException("TODO");
	}
	
	private void autoSetGroups(String s) {
		if(s==null || s.equalsIgnoreCase("null")){return;}
		if(Character.isDigit(s.charAt(0))){
			groups=Integer.parseInt(s);
			return;
		}
		assert(s.equalsIgnoreCase("auto")) : "Unknown groups setting: "+s;
		
		final long maxMem=Shared.memAvailable(1);
		FileFormat ff=FileFormat.testInput(in1, FileFormat.FASTQ, null, false, false);
		if(ff==null || ff.stdio()){return;}
		
//		outstream.println("in1="+in1+", overhead="+(0.5*(ReadKey.overhead+Clump.overhead)));
		
		double[] estimates=Tools.estimateFileMemory(in1, 200, 0.5*(ReadKey.overhead+Clump.overhead), true);
		
//		outstream.println(Arrays.toString(estimates));
		
		double memEstimate=estimates==null ? 0 : estimates[0];
		double diskEstimate=estimates==null ? 0 : estimates[1];
		double worstCase=memEstimate*1.5;

//		outstream.println("Raw Disk Size Estimate: "+(long)(diskEstimate/(1024*1024))+" MB");
		outstream.println("Memory Estimate:        "+(long)(memEstimate/(1024*1024))+" MB");
		outstream.println("Memory Available:       "+(maxMem/(1024*1024))+" MB");
		
		if(maxMem>worstCase){
			groups=1;
		}else{
			groups=Tools.max(11, 3+(int)(4*worstCase/maxMem))|1;
		}
		outstream.println("Set groups to "+groups);
	}
	
	private String getTempFname(String core){
//		outstream.println(core);
		String temp;
		String path="", extension=".fq";
		if(out1!=null){
			core=ReadWrite.stripToCore(out1)+"_"+core;
			path=ReadWrite.getPath(out1);
			extension=ReadWrite.getExtension(out1);
		}
		
		if(useTmpdir && Shared.tmpdir()!=null){
			temp=Shared.tmpdir()+core+Long.toHexString((randy.nextLong()&Long.MAX_VALUE))+extension;
		}else{
			temp=path+core+Long.toHexString((randy.nextLong()&Long.MAX_VALUE))+extension;
		}
		
		final String comp=ReadWrite.compressionType(temp);
		if(forceCompressTemp && comp==null){
			temp+=".gz";
		}else if(comp!=null && forceRawTemp){
			temp=temp.substring(0, temp.lastIndexOf('.'));
		}

//		outstream.println(temp);
		return temp;
	}
	
	public static void shrinkName(Read r) {
		if(r==null){return;}
		String s=r.id;
		if(s.contains("HISEQ")){s=s.replace("HISEQ", "H");}
		if(s.contains("MISEQ")){
			s=s.replace("MISEQ", "M");		
		}
		if(s.contains(":000000000-")){
			s=s.replace(":000000000-", ":");
		}
		r.id=s;
	}
	
	public static void shortName(Read r) {
		StringBuilder sb=new StringBuilder(14);
		long x=r.numericID|1;
		
		while(x<1000000000L){
			x*=10;
			sb.append('0');
		}
		sb.append(r.numericID);
		
//		while(x<0x10000000L){
//			x*=16;
//			sb.append('0');
//		}
//		sb.append(Long.toHexString(r.numericID));
		
		sb.append(r.pairnum()==0 ? " 1:" : " 2:");
		r.id=sb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	
	private boolean quantizeQuality=false;
	private Random randy=new Random();
	private int groups=31;
	private int passes=1;
//	private int k=31;
	private boolean ecco=false;
	private boolean addName=false;
	private String shortName="f";
	private boolean useTmpdir=false;
	private boolean delete=true;
	private boolean useSharedHeader=false;
	private boolean forceCompressTemp=false;
	private boolean forceRawTemp=false;
	private boolean overwrite=true;

	private boolean unpair=false;
	private boolean repair=false;
	private boolean namesort=false;
	private boolean V2=false;
	private boolean V3=true;
	
	private String in1=null;
	private String out1=null;
	
	ArrayList<String> args2=new ArrayList<String>();
	private PrintStream outstream=System.err;
	
}

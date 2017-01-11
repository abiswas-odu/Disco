package jgi;

import java.io.File;
import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.TimeZone;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ListNum;
import dna.Data;
import dna.Parser;
import driver.RenameAndMux;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ReadWrite;
import fileIO.FileFormat;
import fileIO.TextFile;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import align2.BBMap;
import assemble.Tadpole;

/**
 * @author Brian Bushnell
 * @date October 9, 2014
 *
 */
public class DecontaminateByNormalization {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		Timer t=new Timer();
		DecontaminateByNormalization dbn=new DecontaminateByNormalization(args);
		dbn.process(t);
	}
	
	public DecontaminateByNormalization(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		for(String s : args){if(s.startsWith("out=standardout") || s.startsWith("out=stdout")){outstream=System.err;}}
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		boolean setInterleaved=false; //Whether it was explicitly set.

		
		
		Shared.READ_BUFFER_LENGTH=Tools.min(200, Shared.READ_BUFFER_LENGTH);
//		Shared.READ_BUFFER_NUM_BUFFERS=Shared.READ_BUFFER_NUM_BUFFERS;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		CoveragePileup.USE_WINDOW=true;

		final ArrayList<String> readNameFiles=new ArrayList<String>();
		final ArrayList<String> refNameFiles=new ArrayList<String>();
		int onlyProcessFirstN=-1;
		
		Parser parser=new Parser();
		parser.overwrite=true;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens
			
			if(Parser.isJavaFlag(arg)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(parser.parseMapping(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
//				align2.FastaReadInputStream2.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("minc")){
				minc=Float.parseFloat(b);
			}else if(a.equals("minp")){
				minp=Float.parseFloat(b);
			}else if(a.equals("minr")){
				minr=Integer.parseInt(b);
			}else if(a.equals("minl")){
				minl=Integer.parseInt(b);
			}else if(a.equals("mind") || a.equals("mindepth")){
				minDepth=Integer.parseInt(b);
			}else if(a.equals("dp") || a.equals("depthpercentile")){
				depthPercentile=Float.parseFloat(b);
				assert(depthPercentile>=0 && depthPercentile<=1) : "depthPercentile must be between 0 and 1.";
			}else if(a.equals("minprob")){
				minprob=Float.parseFloat(b);
				assert(minprob>=0 && minprob<=1) : "minprob must be between 0 and 1.";
			}else if(a.equals("minratio") || a.equals("ratio")){
				minRatio=Double.parseDouble(b);
			}else if(a.equals("basesundermin")){
				basesUnderMin=Integer.parseInt(b);
			}else if(a.equals("onlyprocessfirstn") || a.equals("opfn")){
				onlyProcessFirstN=Integer.parseInt(b);
			}else if(a.equals("window")){
				CoveragePileup.LOW_COV_WINDOW=Integer.parseInt(b);
			}else if(a.equals("windowcov")){
				CoveragePileup.LOW_COV_DEPTH=Double.parseDouble(b);
			}else if(a.equals("mapraw")){
				mapRawReads=Tools.parseBoolean(b);
			}
			
			/* Tadpole parameters */
			else if(a.equals("tadpole") || a.equals("ecctadpole") || a.equals("ecct") || a.equals("ecc")){
				ecct=Tools.parseBoolean(b);
			}else if(a.equals("tadpoleaggressive") || a.equals("aggressive") || a.equals("aecc")){
				tadpoleAggressive=Tools.parseBoolean(b);
				if(tadpoleAggressive){tadpoleConservative=false;}
			}else if(a.equals("tadpoleconservative") || a.equals("conservative") || a.equals("cecc")){
				tadpoleConservative=Tools.parseBoolean(b);
				if(tadpoleConservative){tadpoleAggressive=false;}
			}else if(a.equals("kt") || a.equals("ktadpole") || a.equals("tadpolek")){
				tadpoleK=Integer.parseInt(b);
			}else if(a.equals("tadpoleprefilter") || a.equals("tadpre")){
				tadpolePrefilter=Integer.parseInt(b);
			}
			
			else if(a.equals("k")){
				normK=Integer.parseInt(b);
			}else if(a.equals("target")){
				normTarget=Integer.parseInt(b);
			}else if(a.equals("hashes")){
				normHashes=Integer.parseInt(b);
			}else if(a.equals("passes")){
				normPasses=Integer.parseInt(b);
			}else if(a.equals("kfilter")){
				kfilter=Integer.parseInt(b);
			}else if(a.equals("ambig") || a.equals("ambiguous")){
				ambigMode=b;
			}
			
			//Deprecated
//			else if(a.equals("ecc")){
//				ecc=Tools.parseBoolean(b);
//			}else if(a.equals("cecc")){
//				cecc=Tools.parseBoolean(b);
//				if(cecc){ecc=true;aecc=false;}
//			}else if(a.equals("aecc")){
//				aecc=Tools.parseBoolean(b);
//				if(aecc){ecc=true;cecc=false;}
//			}
			
			else if(a.equals("prefilter")){
				prefilter=Tools.parseBoolean(b);
			}else if(a.equals("filterbits") || a.equals("fbits")){
				filterBits=Integer.parseInt(b);
			}else if(a.equals("prefilterbits") || a.equals("prebits") || a.equals("pbits")){
				prefilterBits=Integer.parseInt(b);
			}else if(a.equals("logname") || a.equals("log")){
				logName=b;
			}else if(a.equals("resultsname") || a.equals("results") || a.equals("summary")){
				resultsName=b;
			}else if(a.equals("tempdir") || a.equals("tmpdir")){
				tempdir=b;
			}else if(a.equals("outdir") || a.equals("out")){
				outdir=b;
			}else if(a.equals("ref") || a.equals("refs")){
				String[] split2=b.split(",");
				for(String name : split2){
					refNames.add(name);
				}
			}else if(a.equals("read") || a.equals("reads") || a.equals("data")){
				String[] split2=b.split(",");
				for(String name : split2){
					readNames.add(name);
				}
			}else if(a.equals("refnamefile") || a.equals("refnamelist")){
				String[] split2=b.split(",");
				for(String name : split2){
					refNameFiles.add(name);
				}
			}else if(a.equals("readnamefile") || a.equals("readnamelist")){
				String[] split2=b.split(",");
				for(String name : split2){
					readNameFiles.add(name);
				}
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;

			setInterleaved=parser.setInterleaved;
		}

		parseStringsFromFiles(readNameFiles);
		parseStringsFromFiles(refNameFiles);
		
		readNames.addAll(readNameFiles);
		refNames.addAll(refNameFiles);
		
		assert(readNames.size()==refNames.size()) : "Must have same number of read and assembly files. "+readNames.size()+"!="+refNames.size();
		
		while(onlyProcessFirstN>0 && onlyProcessFirstN<readNames.size()){ //For testing.
			readNames.remove(readNames.size()-1);
			refNames.remove(refNames.size()-1);
		}
//		System.err.println("\n************ 5\n"+readNames+"\n\n"+refNames+"\n\n"+readNameFiles+"\n\n"+refNameFiles);
		
		if(outdir!=null && outdir.length()>0 && !outdir.endsWith("/")){outdir=outdir+"/";}
		if(tempdir!=null && tempdir.length()>0 && !tempdir.endsWith("/")){tempdir=tempdir+"/";}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	void process(Timer t){
		log("Decontaminate start", false);
		String mergePath=tempdir+"_merged.fq.gz";
		String tadpolePath=tempdir+"_corrected.fq.gz";
		String normPath=tempdir+"_normalized.fq.gz";
		
		if(mapRawReads){
			map(readNames, refNames, 0);
		}
		renameAndMux_MT(readNames, mergePath);
		if(ecct){
			eccTadpole(mergePath, tadpolePath);
			if(deleteFiles){delete(mergePath);}
			mergePath=tadpolePath;
		}
		normalize(mergePath, normPath, normK, minDepth, normTarget, normHashes, normPasses, ecc, prefilter, normalizeByLowerDepth);
		if(deleteFiles){delete(mergePath);}
		demux(normPath, readNames);
		if(deleteFiles){delete(normPath);}
		map(readNames, refNames, 1);
		filter(readNames, refNames);
		
		t.stop();
		
		outstream.println("Time: \t"+t);
		log("Decontaminate finish", true);
	}
	
	public static void parseStringsFromFiles(ArrayList<String> list){
		String[] x=list.toArray(new String[list.size()]);
		list.clear();
		for(String s : x){
			File f=new File(s);
			if(f.exists() && f.isFile()){
				TextFile tf=new TextFile(s);
				String[] lines=tf.toStringLines();
				for(String s2 : lines){
					list.add(s2);
				}
			}else{
				list.add(s);
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private void renameAndMux_ST(ArrayList<String> readPaths, String fnameOut){
		log("renameAndMux start", true);
		System.err.println("\nRename/Merge Phase Start");

		FileFormat ffout=FileFormat.testOutput(fnameOut, FileFormat.FASTQ, null, true, true, false, false);
		
		long readsProcessed=0;
		long basesProcessed=0;
		
		for(String in : readPaths){
			
			final FileFormat ffin=FileFormat.testInput(in, FileFormat.FASTQ, null, true, true);
			final ConcurrentReadInputStream cris;
			final String core=ReadWrite.stripToCore(in);
			{
				cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin, null);
				if(verbose){outstream.println("Started cris");}
				cris.start(); //4567
			}
			boolean paired=cris.paired();
			if(!ffin.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
			final ConcurrentReadOutputStream ros;
			
			final int buff=4;

			if(cris.paired()){
				outstream.println("Writing interleaved.");
			}			

			assert(!in.equalsIgnoreCase(fnameOut)) : "Input file and output file have same name.";

			ros=ConcurrentReadOutputStream.getStream(ffout, null, buff, null, false);
			ros.start();
			
			{
				
				ListNum<Read> ln=cris.nextList();
				ArrayList<Read> reads=(ln!=null ? ln.list : null);
				
				if(reads!=null && !reads.isEmpty()){
					Read r=reads.get(0);
					assert((r.mate!=null)==cris.paired());
				}

				while(reads!=null && reads.size()>0){
					
					for(int idx=0; idx<reads.size(); idx++){
						final Read r1=reads.get(idx);
						final Read r2=r1.mate;
						
						final int initialLength1=r1.length();
						final int initialLength2=(r1.mateLength());
						
						{
							readsProcessed++;
							basesProcessed+=initialLength1;
							r1.id=core+"_"+r1.numericID+" 1:";
						}
						if(r2!=null){
							readsProcessed++;
							basesProcessed+=initialLength2;
							r2.id=core+"_"+r1.numericID+" 2:";
						}
						
						
					}
					
					final ArrayList<Read> listOut=reads;
					
					if(ros!=null){ros.add(listOut, ln.id);}

					cris.returnList(ln.id, ln.list.isEmpty());
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
				if(ln!=null){
					cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
				}
			}
			
			errorState|=ReadWrite.closeStreams(cris, ros);
			
			ffout=FileFormat.testOutput(fnameOut, FileFormat.FASTQ, null, true, false, true, false);
		}
		log("renameAndMux finish", true);
	}
	
	private void renameAndMux_MT(ArrayList<String> readPaths, String fnameOut){
		log("renameAndMux start", true);
		System.err.println("\nRename/Merge Phase Start");
		
//		String dir=ReadWrite.parseRoot(fnameIn);
		
		ArrayList<String> argList=new ArrayList<String>();
		
		{
			argList.add("out="+fnameOut);
			
			StringBuilder sb=new StringBuilder("in=");
			String comma="";
			for(String s : readPaths){
				sb.append(comma);
				sb.append(s);
				comma=",";
			}
			argList.add(sb.toString());
		}
		
		String[] args=argList.toArray(new String[0]);
		
		{//Run Mux
			int oldZT=ReadWrite.MAX_ZIP_THREADS;
			boolean oldUP=ReadWrite.USE_PIGZ;
			boolean oldUU=ReadWrite.USE_UNPIGZ;
			try {
				RenameAndMux.main(args);
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
			ReadWrite.MAX_ZIP_THREADS=oldZT;
			ReadWrite.USE_PIGZ=oldUP;
			ReadWrite.USE_UNPIGZ=oldUU;
		}

		log("renameAndMux finish", true);
	}
	
	private void eccTadpole(String fnameIn, String fnameOut){
		log("tadpole start", true);
		System.err.println("\nTadpole Phase Start");
		
//		String dir=ReadWrite.parseRoot(fnameIn);
		
		ArrayList<String> argList=new ArrayList<String>();
		
		{
			argList.add("in="+fnameIn);
			argList.add("out="+fnameOut);
			argList.add("k="+tadpoleK);
			argList.add("mode=correct");
			if(tadpoleAggressive){
				argList.add("aggressive");
			}else if(tadpoleConservative){
				argList.add("conservative");
			}
			if(tadpolePrefilter>=0){
				argList.add("prefilter="+tadpolePrefilter);
			}
			argList.add("prealloc="+tadpolePrealloc);
		}
		
		String[] args=argList.toArray(new String[0]);
		
		{//Run Mux
			int oldZT=ReadWrite.MAX_ZIP_THREADS;
			boolean oldUP=ReadWrite.USE_PIGZ;
			boolean oldUU=ReadWrite.USE_UNPIGZ;
			try {
				Tadpole.main(args);
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
			ReadWrite.MAX_ZIP_THREADS=oldZT;
			ReadWrite.USE_PIGZ=oldUP;
			ReadWrite.USE_UNPIGZ=oldUU;
		}

		log("tadpole finish", true);
	}
	
	private void normalize(String fnameIn, String fnameOut, int k, int min, int target, int hashes, int passes, boolean ecc, boolean prefilter, boolean uselowerdepth){
		log("normalization start", true);
		System.err.println("\nNormalization/Error Correction Phase Start");
		
		ArrayList<String> argList=new ArrayList<String>();
		
		{
			argList.add("ecc="+ecc);
			if(ecc && aecc){argList.add("aecc="+aecc);}
			if(ecc && cecc){argList.add("cecc="+cecc);}
			argList.add("prefilter="+prefilter);
			argList.add("bits="+filterBits);
			argList.add("prebits="+prefilterBits);
			argList.add("hashes="+hashes);
			argList.add("passes="+passes);
			argList.add("target="+target);
			argList.add("k="+k);
			argList.add("mindepth="+min);
			argList.add("maxdepth="+target);
			argList.add("minprob="+minprob);
			argList.add("dp="+depthPercentile);
			argList.add("in="+fnameIn);
			argList.add("out="+fnameOut);
			argList.add("uld="+uselowerdepth);
		}
		
		String[] normargs=argList.toArray(new String[0]);
		
		{//Run BBNorm
			try {
				KmerNormalize.main(normargs);
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
		}
		
		log("normalization finish", true);
	}
	
	private void demux(String fnameIn, ArrayList<String> readPaths){
		log("demux start", true);
		System.err.println("\nDemux Phase Start");
		
//		String dir=ReadWrite.parseRoot(fnameIn);
		
		ArrayList<String> argList=new ArrayList<String>();
		
		{

			argList.add("in="+fnameIn);
			argList.add("out="+tempdir+"%_demuxed.fq.gz");
			
			StringBuilder sb=new StringBuilder("names=");
			String comma="";
			for(String s : readPaths){
				sb.append(comma);
				sb.append(ReadWrite.stripToCore(s));
				comma=",";
			}
			argList.add(sb.toString());
		}
		
		String[] args=argList.toArray(new String[0]);
		
		{//Run Demux
			int oldZT=ReadWrite.MAX_ZIP_THREADS;
			boolean oldUP=ReadWrite.USE_PIGZ;
			boolean oldUU=ReadWrite.USE_UNPIGZ;
			ReadWrite.MAX_ZIP_THREADS=Tools.mid(ReadWrite.MAX_ZIP_THREADS, 2, 4);
			try {
				DemuxByName.main(args);
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
			ReadWrite.MAX_ZIP_THREADS=oldZT;
			ReadWrite.USE_PIGZ=oldUP;
			ReadWrite.USE_UNPIGZ=oldUU;
		}
		
		log("demux finish", true);
	}
	
	private void map(ArrayList<String> readPaths, ArrayList<String> refnames, int pass){
		log("map start", true);
		System.err.println("\nMapping Phase Start");
		
		for(int i=0; i<readPaths.size(); i++){
			final String path=readPaths.get(i);
			final String ref=refnames.get(i);
			final String core=ReadWrite.stripToCore(path);
			final String demuxed=tempdir+core+"_demuxed.fq.gz";
			final String dir=(outdir==null ? "" : outdir);
			
			final String infile=(pass==0 ? path : demuxed);
			
			ArrayList<String> argList=new ArrayList<String>();

			{
				argList.add("in="+infile);
				argList.add("ref="+ref);
				argList.add("covstats="+dir+core+"_covstats"+pass+".txt");
				argList.add("arrays=t");
				argList.add("nodisk");
				argList.add("ambig="+ambigMode);
				if(kfilter>1){argList.add("kfilter="+kfilter);}
				argList.add("fast");
				argList.add("ow="+overwrite);
				argList.add("minscaf=0");
			}

			String[] args=argList.toArray(new String[0]);

			{//Run BBMap
				try {
					BBMap.main(args);
				} catch (Exception e) {
					e.printStackTrace();
					log("failed", true);
					System.exit(1);
				}
			}
			Data.unloadAll();
		}
		
		log("map finish", true);
	}
	
	private void filter(ArrayList<String> readPaths, ArrayList<String> refnames){
		//filterbycoverage.sh -Xmx2g in=PUCA.fasta cov=PUCA_stats.txt out=PUCA_clean.fa minc=5 minp=40
		log("filter start", true);
		System.err.println("\nFiltering Phase Start");
		
		for(int i=0; i<readPaths.size(); i++){
			String path=readPaths.get(i);
			String ref=refnames.get(i);
			String core=ReadWrite.stripToCore(path);
			String dir=(outdir==null ? /*ReadWrite.parseRoot(ref)*/ "" : outdir);
			String stats0=mapRawReads ? (outdir==null ? "" : outdir)+core+"_covstats0.txt" : null;
			String stats1=(outdir==null ? "" : outdir)+core+"_covstats1.txt";
			String results=(resultsName==null ? "" : (resultsName.contains("/") || resultsName.contains("\\") ? resultsName : (outdir==null ? "" : outdir)+resultsName));
			
			ArrayList<String> argList=new ArrayList<String>();

			{
				argList.add("log="+results);
				argList.add("appendlog="+(i>0));
				argList.add("logheader="+(i==0));
				if(stats0!=null){argList.add("cov0="+stats0);}
				argList.add("cov1="+stats1);
				argList.add("in="+ref);
				argList.add("out="+dir+core+"_clean.fasta");
				argList.add("outd="+dir+core+"_dirty.fasta");
				argList.add("minc="+minc);
				argList.add("minp="+minp);
				argList.add("minr="+minr);
				argList.add("minl="+minl);
				argList.add("basesundermin="+basesUnderMin);
				if(stats0!=null){argList.add("minratio="+minRatio);}
				argList.add("ow="+overwrite);
			}

			String[] args=argList.toArray(new String[0]);

			{//Run filtering
				try {
					FilterByCoverage.main(args);
				} catch (Exception e) {
					e.printStackTrace();
					log("failed", true);
					System.exit(1);
				}
			}
		}
		
		log("filter finish", true);
	}
	
	/**
	 * Log a message in the log file
	 * @param message Message to log
	 * @param append True to append, false to overwrite
	 */
	private void log(String message, boolean append){
		if(logName!=null){
			ReadWrite.writeString(message+", "+timeString()+"\n", logName, append);
		}
	}
	
	/**
	 * TODO:  Some machines are set to UTC rather than PST
	 * @return Timestamp in RQC's format
	 */
	public static String timeString(){
		SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
//		sdf.setTimeZone(TimeZone.getTimeZone("PST"));
		sdf.setTimeZone(TimeZone.getDefault());
		return sdf.format(new Date());
	}
	
	/**
	 * Delete all non-null filenames.
	 * @param prefix Append this prefix to filenames before attempting to delete them
	 * @param names Filenames to delete
	 */
	private void delete(String path){
		if(path==null){return;}
		if(verbose){System.err.println("Trying to delete "+path);}
		File f=new File(path);
		if(f.exists()){
			try {
				f.delete();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	
	private void printOptions(){
		assert(false) : "printOptions: TODO";
//		outstream.println("Syntax:\n");
//		outstream.println("java -ea -Xmx512m -cp <path> jgi.ReformatReads in=<infile> in2=<infile2> out=<outfile> out2=<outfile2>");
//		outstream.println("\nin2 and out2 are optional.  \nIf input is paired and there is only one output file, it will be written interleaved.\n");
//		outstream.println("Other parameters and their defaults:\n");
//		outstream.println("overwrite=false  \tOverwrites files that already exist");
//		outstream.println("ziplevel=4       \tSet compression level, 1 (low) to 9 (max)");
//		outstream.println("interleaved=false\tDetermines whether input file is considered interleaved");
//		outstream.println("fastawrap=70     \tLength of lines in fasta output");
//		outstream.println("qin=auto         \tASCII offset for input quality.  May be set to 33 (Sanger), 64 (Illumina), or auto");
//		outstream.println("qout=auto        \tASCII offset for output quality.  May be set to 33 (Sanger), 64 (Illumina), or auto (meaning same as input)");
//		outstream.println("outsingle=<file> \t(outs) Write singleton reads here, when conditionally discarding reads from pairs.");
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private final ArrayList<String> readNames=new ArrayList<String>();
	private final ArrayList<String> refNames=new ArrayList<String>();
	
	private String logName=null;
	private String resultsName="results.txt";
	private String tempdir=(Shared.tmpdir() == null ? "" : Shared.tmpdir());
	private String outdir=null;
	
	/*--------------------------------------------------------------*/


	private int kfilter=55;
	private String ambigMode="random";
	
	private long maxReads=-1;

	private float minc=3.5f;
	private float minp=20;
	private int minr=20;
	private int minl=500;
	
	/** Delete temp files after use */
	private boolean deleteFiles=true;

	/** Error correct with Tadpole */
	private boolean ecct=false;
	private boolean tadpoleAggressive=false;
	private boolean tadpoleConservative=false;
	private boolean tadpolePrealloc=true;
	private int tadpoleK=42;
	private int tadpolePrefilter=1;
	
	/** Scaffolds will be discarded if there are at least this many bases in windows below a coverage cutoff. */
	private int basesUnderMin=-1;
	
	private float depthPercentile=0.75f;
	private float minprob=0.5f;
	private int minDepth=2;
	private int normK=31;
	private int normTarget=20;
	private int normHashes=4;
	private int normPasses=1;
	private int filterBits=32;
	private int prefilterBits=2;
	private boolean ecc=false;
	private boolean cecc=false;
	private boolean aecc=false;
	private boolean prefilter=true;
	private boolean normalizeByLowerDepth=false;

	private double minRatio=1.2f;
	private boolean mapRawReads=true;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	
}

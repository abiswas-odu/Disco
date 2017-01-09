package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import dna.Parser;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ByteBuilder;
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;
import structures.StringPair;
import var2.CallVariants2.Sample;

/**
 * @author Brian Bushnell
 * @date December 18, 2016
 *
 */
public class MergeSamples {
	
	public static void main(String[] args){
		Timer t=new Timer();
		MergeSamples ms=new MergeSamples(args);
		//ms.process(t);
	}
	
	public MergeSamples(){}
	
	public MergeSamples(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
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
			}else if(a.equals("invalid")){
				outInvalid=b;
			}else if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			overwrite=parser.overwrite;
			append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	public void mergeSamples(ArrayList<Sample> list, ScafMap scafMap, String outVcf){
		map=scafMap;
		ArrayList<StringPair> vcfList=new ArrayList<StringPair>(list.size());
		for(Sample s : list){vcfList.add(new StringPair(s.name, s.vcfName));}
		mergeFiles(vcfList, scafMap, outVcf);
	}
	
	public void mergeFiles(ArrayList<StringPair> list, ScafMap scafMap, String outVcf){
		System.err.println("Merging "+list);
		final int ways=list.size();
		ByteFile[] bfa=new ByteFile[ways];
		for(int i=0; i<ways; i++){
			StringPair pair=list.get(i);
			FileFormat ff=FileFormat.testInput(pair.b, FileFormat.TEXT, null, false, false);
			bfa[i]=ByteFile.makeByteFile(ff, false);
		}
//		System.err.println("Made byte files.");
		
		ByteStreamWriter bswVcf=null;
		if(outVcf!=null){
			bswVcf=new ByteStreamWriter(outVcf, true, false, true, FileFormat.TEXT);
			bswVcf.start();
		}
//		System.err.println("Started writer.");
		
		ByteBuilder bb=new ByteBuilder(34000);
		VCFLine[] row=processRow(bfa, bb);
		while(row!=null){
//			System.err.println("Processed a line.");
			if(row[0]!=null){
				VCFLine merged=merge(row);
				merged.toText(bb);
				bb.append('\n');
				if(bb.length>32000){
					bswVcf.print(bb);
					bb=new ByteBuilder(34000);
				}
			}
			row=processRow(bfa, bb);
		}
//		System.err.println("Finished processing.");
		
		if(bswVcf!=null){
			if(bb.length>0){bswVcf.print(bb);}
			bswVcf.poisonAndWait();
		}

//		System.err.println("Closed stream.");
	}
	
	VCFLine[] processRow(ByteFile[] bfa, ByteBuilder bb){
		byte[][] lines=new byte[bfa.length][];
		for(int i=0; i<bfa.length; i++){
			byte[] line=bfa[i].nextLine();
			if(line==null){return null;}
			lines[i]=line;
		}
		
		VCFLine[] row=new VCFLine[bfa.length];
		if(lines[0][0]=='#'){
			processHeader(lines, bb);
			return row;
		}
		for(int i=0; i<lines.length; i++){
			byte[] line=lines[i];
			row[i]=new VCFLine(line);
			if(i>0){assert(row[i].pos==row[0].pos) : "\n"+row[0]+"\n"+row[i];}
		}
		return row;
	}
	
	void processHeader(byte[][] lines, ByteBuilder bb){
		String[][] matrix=new String[lines.length][];
		for(int i=0; i<lines.length; i++){
			matrix[i]=new String(lines[i]).split("=");
		}
		
		if(matrix[0][0].equals("##ploidy")){
			ploidy=Integer.parseInt(matrix[0][1]);
			bb.append("##ploidy="+ploidy+"\n");
		}else if(matrix[0][0].equals("##reads")){
			for(String[] split : matrix){
				reads+=Long.parseLong(split[1]);
			}
			bb.append("##reads="+reads+"\n");
		}else if(matrix[0][0].equals("##pairedReads")){
			for(String[] split : matrix){
				pairedReads+=Long.parseLong(split[1]);
			}
			bb.append("##pairedReads="+pairedReads+"\n");
		}else if(matrix[0][0].equals("##properlyPairedReads")){
			for(String[] split : matrix){
				properlyPairedReads+=Long.parseLong(split[1]);
			}
			properPairRate=(float)(properlyPairedReads*1.0/(Tools.max(1, reads)));
			bb.append("##properlyPairedReads="+properlyPairedReads+"\n");
			bb.append("##properPairRate="+String.format("%.5f\n", properPairRate));
		}else if(matrix[0][0].equals("##properPairRate")){
			//do nothing
		}else if(matrix[0][0].equals("##totalQualityAvg")){
			for(String[] split : matrix){
				totalQualityAvg+=Float.parseFloat(split[1]);
			}
			totalQualityAvg/=lines.length;
			bb.append("##totalQualityAvg="+String.format("%.3f\n", totalQualityAvg));
		}else if(matrix[0][0].equals("##mapqAvg")){
			for(String[] split : matrix){
				mapqAvg+=Float.parseFloat(split[1]);
			}
			mapqAvg/=lines.length;
			bb.append("##mapqAvg="+String.format("%.3f\n", mapqAvg));
		}else if(matrix[0][0].startsWith("#CHROM\tPOS\t")){
			bb.append(lines[0]);
			for(int i=1; i<lines.length; i++){
				String[] split=new String(lines[i]).split("\t");
				bb.append('\t').append(split[split.length-1]);
			}
			bb.append('\n');
		}else{
			bb.append(lines[0]);
			bb.append('\n');
		}
	}
	
	VCFLine merge(VCFLine[] row){

//		System.err.println(row.length);
//		System.err.println(row[0]);
		
		Var sum=null;
		VCFLine best=null;
		for(VCFLine line : row){
			if(best==null || line.qual>best.qual){best=line;}
			Var v=line.toVar();
			assert(v!=null);
			if(sum==null){sum=v;}
			else{sum.add(v);}
		}
		
		assert(sum!=null) : row.length+", "+row[0];
		
		//ByteBuilder bb, float properPairRate, float totalQualityAvg, float mapqAvg, int ploidy, ScafMap map, VarFilter filter, boolean trimWhitespace
		ByteBuilder vcfBuilder=sum.toVCF(new ByteBuilder(), properPairRate, totalQualityAvg, mapqAvg, ploidy, map, filter, trimWhitespace);
		VCFLine merged=new VCFLine(vcfBuilder.toBytes());
		merged.samples.clear();
		for(VCFLine line : row){
			merged.samples.addAll(line.samples);
		}
		if(merged.qual<best.qual){
			merged.qual=best.qual;
			merged.filter=best.filter;
		}
		return merged;
	}
	
	/*--------------------------------------------------------------*/
	
	private void printOptions(){
		assert(false) : "printOptions: TODO";
	}
	
	
	/*--------------------------------------------------------------*/

	long readsSum;
	long pairsSum;
	int ploidy=1;
	
	float properPairRate;
	float totalQualityAvg;
	float mapqAvg;
	
	long reads;
	long pairedReads;
	long properlyPairedReads;
	
	VarFilter filter;
	ScafMap map;
	boolean trimWhitespace=true;
	
	private String in1=null;
	private String out1=null;
	private String outInvalid=null;
	
	/*--------------------------------------------------------------*/
	
	private long linesProcessed=0;
	private long linesValid=0;
	private long bytesProcessed=0;
	
	private long maxLines=Long.MAX_VALUE;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	
}

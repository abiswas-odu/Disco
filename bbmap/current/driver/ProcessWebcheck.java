package driver;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

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
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;
import structures.IntList;

/**
 * @author Brian Bushnell
 * @date May 9, 2016
 *
 */
public class ProcessWebcheck {
	
	public static void main(String[] args){
		Timer t=new Timer();
		ProcessWebcheck sample=new ProcessWebcheck(args);
		sample.process(t);
	}
	
	public ProcessWebcheck(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
//		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
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
			}else if(a.equals("extendedstats") || a.equals("es")){
				extendedStats=Tools.parseBoolean(b);
			}else if(a.equals("invalid") || a.equals("outinvalid")){
				outInvalid=b;
			}else if(a.equals("fail") || a.equals("outfail")){
				outFail=b;
			}else if(a.equals("ms") || a.equals("millis")){
				boolean x=Tools.parseBoolean(b);
				ms=(x ? "ms" : "");
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
			}else if(!arg.contains("=")){
				for(String s : arg.split(",")){
					in1.add(s);
				}
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			overwrite=parser.overwrite;
			append=parser.append;
			
			if(parser.in1!=null){
				for(String s : parser.in1.split(",")){
					in1.add(s);
				}
			}
			
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

		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);
		ffoutInvalid=FileFormat.testOutput(outInvalid, FileFormat.TXT, null, true, overwrite, append, false);
		ffoutFail=FileFormat.testOutput(outFail, FileFormat.TXT, null, true, overwrite, append, false);
		
		ffin1=new ArrayList<FileFormat>(in1.size());
		for(String s : in1){
			FileFormat ff=FileFormat.testInput(s, FileFormat.TXT, null, true, true);
			assert(ff!=null) : "Cannot read file "+s;
			ffin1.add(ff);
		}
		assert(ffin1.size()>0) : "No input files.";
	}
	
	void process(Timer t){

		ByteStreamWriter bsw=null;
		if(ffout1!=null){
			bsw=new ByteStreamWriter(ffout1);
			bsw.start();
		}

		ByteStreamWriter bswInvalid=null;
		if(ffoutInvalid!=null){
			bswInvalid=new ByteStreamWriter(ffoutInvalid);
			bswInvalid.start();
		}

		ByteStreamWriter bswFail=null;
		if(ffoutFail!=null){
			bswFail=new ByteStreamWriter(ffoutFail);
			bswFail.start();
		}
		
		for(FileFormat ff : ffin1){
			ByteFile bf=ByteFile.makeByteFile(ff, false);
			process2(bf, bswFail, bswInvalid);
		}

		//		double rpnano=linesProcessed/(double)(t.elapsed);
		//		double bpnano=bytesProcessed/(double)(t.elapsed);
		//
		//		String rpstring=(linesProcessed<100000 ? ""+linesProcessed : linesProcessed<100000000 ? (linesProcessed/1000)+"k" : (linesProcessed/1000000)+"m");
		//		String bpstring=(bytesProcessed<100000 ? ""+bytesProcessed : bytesProcessed<100000000 ? (bytesProcessed/1000)+"k" : (bytesProcessed/1000000)+"m");
		//
		//		while(rpstring.length()<8){rpstring=" "+rpstring;}
		//		while(bpstring.length()<8){bpstring=" "+bpstring;}
		//		
		//		outstream.println("Time:                         \t"+t);
		//		outstream.println("Lines Processed:    "+rpstring+" \t"+String.format("%.2fk lines/sec", rpnano*1000000));
		//		outstream.println("Bytes Processed:    "+bpstring+" \t"+String.format("%.2fm bytes/sec", bpnano*1000));

		passLatency.shrink();
		failLatency.shrink();
		failCode.sort();
		failCode.shrinkToUnique();

		StringBuilder sb=new StringBuilder();

		ArrayList<String> list=new ArrayList<String>();
		list.addAll(map.keySet());
		Shared.sort(list);
		for(String s : list){
			sb.append(s+"\t"+map.get(s)[0]+"\n");
		}

		if(extendedStats){
			sb.append('\n');
			sb.append("Lines_Processed:\t"
					+ linesProcessed).append('\n');
			sb.append("Invalid_Lines:\t"
					+ (linesProcessed-linesValid)).append('\n');
			sb.append("Passing:\t"
					+ passLatency.size).append('\n');
			sb.append("Failing:\t"
					+ failLatency.size).append('\n');
			sb.append("Avg_Pass_Latency:\t"
					+ (passLatency.size>0 ? Tools.average(passLatency.array) : 0)+ms).append('\n');
			sb.append("Max_Pass_Latency:\t"
					+ (passLatency.size>0 ? Tools.max(passLatency.array) : 0)+ms).append('\n');
			sb.append("Avg_Fail_Latency:\t"
					+ (failLatency.size>0 ? Tools.average(failLatency.array) : 0)+ms).append('\n');
			sb.append("Max_Fail_Latency:\t"
					+ (failLatency.size>0 ? Tools.max(failLatency.array) : 0)+ms).append('\n');
			sb.append("Observed_Fail_Codes:");
			for(int i=0; i<failCode.size; i++){
				sb.append('\t').append(failCode.get(i));
			}
			sb.append('\n');
		}

		outstream.print(sb);

		if(bsw!=null){
			bsw.print(sb);
			errorState|=bsw.poisonAndWait();
		}
		if(bswInvalid!=null){errorState|=bswInvalid.poisonAndWait();}
		if(bswFail!=null){errorState|=bswFail.poisonAndWait();}

		t.stop();

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
		
	private void process2(ByteFile bf, ByteStreamWriter bswFail, ByteStreamWriter bswInvalid){
		
		byte[] line=bf.nextLine();
		
		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=line.length;
				
				boolean valid=(line[0]!='#' && Character.isDigit(line[line.length-1]));
				String[] split=null;
				if(valid){
					split=new String(line).split("\\|");
					valid=(split.length==4);
				}
//				assert(false) : (char)line[0]+", "+(char)line[line.length-1]+", "+new String(line).split("\\|").length;
				
				int code=-1;
				float latency=-1;
				if(valid){
					try {
						code=Integer.parseInt(split[2].substring(0, split[2].indexOf(' ')));
						latency=Float.parseFloat(split[3]);
					} catch (Exception e) {
						valid=false;
					}
				}
				
				if(valid){
					linesValid++;
					
					long[] cnt=map.get(split[2]);
					if(cnt==null){
						cnt=new long[1];
						map.put(split[2], cnt);
					}
					cnt[0]++;
					
					int latency2=(int)(latency*1000);
					
					if(code==200){
						passLatency.add(latency2);
					}else{
						failLatency.add(latency2);
						failCode.add(code);
						if(bswFail!=null){
							bswFail.println(line);
						}
					}
				}else{
					if(bswInvalid!=null){
						bswInvalid.println(line);
					}
				}
			}
			line=bf.nextLine();
		}
		
		errorState|=bf.close();
	}
	
	/*--------------------------------------------------------------*/
	
	private void printOptions(){
		assert(false) : "printOptions: TODO";
	}
	
	
	/*--------------------------------------------------------------*/
	
	private ArrayList<String> in1=new ArrayList<String>();
	private String out1=null;
	private String outInvalid=null;
	private String outFail=null;
	
	private String ms="ms";
	
	/*--------------------------------------------------------------*/
	
	private long linesProcessed=0;
	private long linesValid=0;
	private long bytesProcessed=0;

	private IntList passLatency=new IntList();
	private IntList failLatency=new IntList();
	private IntList failCode=new IntList();
	
	private HashMap<String, long[]> map=new HashMap<String, long[]>();
	
	private long maxLines=Long.MAX_VALUE;
	private boolean extendedStats=false;
	
	/*--------------------------------------------------------------*/
	
	private final ArrayList<FileFormat> ffin1;
	private final FileFormat ffout1;
	private final FileFormat ffoutInvalid;
	private final FileFormat ffoutFail;
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.out;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	
}

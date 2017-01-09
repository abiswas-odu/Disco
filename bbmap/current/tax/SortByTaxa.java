package tax;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;

import stream.ByteBuilder;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import dna.Parser;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Mar 11, 2015
 *
 */
public class SortByTaxa {
	
	public static void main(String[] args){
		Timer t=new Timer();
		SortByTaxa mb=new SortByTaxa(args);
		mb.process(t);
	}
	
	public SortByTaxa(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		
		
		Shared.READ_BUFFER_LENGTH=Tools.min(200, Shared.READ_BUFFER_LENGTH);
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		
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
			}else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				tableFile=b;
				if("auto".equalsIgnoreCase(b)){tableFile=TaxTree.defaultTableFile();}
			}else if(a.equals("tree") || a.equals("taxtree")){
				treeFile=b;
				if("auto".equalsIgnoreCase(b)){treeFile=TaxTree.defaultTreeFile();}
			}else if(a.equals("fuse")){
				fuse=Tools.parseBoolean(b);
			}else if(a.equals("dummyreads") || a.equals("adddummies") || a.equals("dummy")){
				addDummyReads=Tools.parseBoolean(b);
			}else if(a.equals("dummylevel")){
				dummyLevel=Integer.parseInt(b);
			}else if(a.equals("promote")){
				promote=Integer.parseInt(b);
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
//				align2.FastaReadInputStream2.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(parser.in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				parser.in1=arg;
			}else if(parser.out1==null && i==1 && !arg.contains("=")){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
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
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FA, null, true, overwrite, append, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FA, null, true, true);
		
		if(tableFile!=null){
			outstream.println("Loading gi table.");
			GiToNcbi.initialize(tableFile);
		}
		if(treeFile!=null){
			outstream.println("Loading tree.");
			tree=ReadWrite.read(TaxTree.class, treeFile, true);
		}else{
			tree=null;
		}
	}
	
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null, null, null);
			cris.start();
			if(verbose){outstream.println("Started cris");}
		}
		
		long readsProcessed=0;
		long basesProcessed=0;
		
		long dummiesAdded=0;
		
		ArrayList<Read> all=new ArrayList<Read>();
		
		{
			outstream.println("Loading sequences.");
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(reads!=null && reads.size()>0){
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					
					final int initialLength1=r1.length();
					
					if(tree!=null){
						TaxNode tn=tree.getNode(r1.id);
						if(tn!=null){tn.incrementRaw(1);}
					}
					
					readsProcessed++;
					basesProcessed+=initialLength1;
				}
				
				all.addAll(reads);

				cris.returnList(ln.id, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadWrite.closeStream(cris);
		
		if(addDummyReads){
			outstream.println("Adding dummies.");
			for(TaxNode n : tree.nodes){
				if(n!=null && n.level>=dummyLevel && n.countRaw<1){
					Read r=new Read(dummyBases, -1, -1, -1, "ncbi|"+n.id, null, all.size(), 0);
					all.add(r);
					dummiesAdded++;
				}
			}
		}
		
		{
			outstream.println("Sorting.");
			Shared.sort(all, taxaComparator);
		}
		
		if(fuse){
			outstream.println("Fusing.");
			ArrayList<Read> fused=new ArrayList<Read>();
			ArrayList<Read> current=new ArrayList<Read>();
			final ByteBuilder bb=new ByteBuilder(1000000);
			
			int taxid=-2;
			int segment=0;
			long currentLength=0;
			for(int i=0; i<all.size(); i++){
				Read r=all.remove(i);
				int tax=GiToNcbi.getID(r.id);
				if(promote>-1){
					TaxNode n=tree.getNode(tax);
					assert(n!=null) : "Can't find node for "+r.id+", "+r.numericID+", "+r.length();
					while(n.level<promote){
						n=tree.getNode(n.pid);
					}
					tax=n.id;
				}
				if(tax!=taxid || r.length()+currentLength>MAX_FUSE_LENGTH){
					Read x=fuse(current, taxid, segment, bb);
					current.clear();
					currentLength=0;
					if(tax==taxid){segment++;}
					else{segment=0;}
					if(x!=null){
						fused.add(x);
					}
				}
				current.add(r);
				currentLength+=(r.length()+padding);
				taxid=tax;
			}
			{
				Read x=fuse(current, taxid, segment, bb);
				current.clear();
				if(x!=null){
					fused.add(x);
				}
			}
			all=fused;
		}
		
		if(out1!=null){
			outstream.println("Writing output.");
			final ConcurrentReadOutputStream ros;
			final int buff=4;		

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";

			ros=ConcurrentReadOutputStream.getStream(ffout1, null, null, null, buff, null, false);
			ros.start();
			
			ArrayList<Read> list=new ArrayList<Read>();
			long num=0;
			for(Read r : all){
				list.add(r);
				if(list.size()>=200){
					ros.add(list, num);
					num++;
					list=new ArrayList<Read>();
				}
			}
			if(list.size()>0){
				ros.add(list, num);
				num++;
			}
			
			errorState|=ReadWrite.closeStream(ros);
		}
		
		t.stop();
		
		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);

		String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
		String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		outstream.println();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		
		outstream.println();
		if(dummiesAdded>0){
			outstream.println("Dummies Added:      "+dummiesAdded);
		}
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * @param current
	 * @param taxid
	 * @return
	 */
	private Read fuse(ArrayList<Read> current, int taxid, int segment, ByteBuilder bb) {
//		System.err.println("Calling fuse("+current+")");
		
		if(current.isEmpty()){return null;}
		if(current.size()==1){
			Read x=current.get(0);
			x.id="ncbi|"+taxid+"|"+segment;
//			System.err.println("a) Returning "+x);
			return x;
		}
		Read x=current.get(0);
		bb.setLength(0);

		long lensum=0;
		
		Read last=null;
		try {
			for(int i=0; i<current.size(); i++){
				Read r=current.remove(i);
				last=r;
				lensum+=r.length();
				if(bb.length()>0 && r.length()>0){
					for(int j=0; j<padding; j++){
						bb.append('N');
						lensum++;
					}
				}
				bb.append(r.bases);
			}
		} catch (Throwable e) {
			System.err.println(lensum+", "+last.length()+", "+last.id+"\n"+current.size()+", "+taxid+", "+tree.getNode(taxid));
			System.err.println(e);
		}
		x.bases=bb.toBytes();
		x.quality=null;
		x.id="ncbi|"+taxid;
//		System.err.println("b) Returning "+x);
		return x;
	}
	
	@SuppressWarnings("unused")
	private void validate(ArrayList<Read> all){
		for(Read a : all){
			for(Read b : all){
				assert(taxaComparator.compare(a, b)==-taxaComparator.compare(b, a)) : (verbose=true)+"\n"+a.id+", "+b.id+"\n"+
						taxaComparator.compare(a, b)+"\n"+taxaComparator.compare(b, a);
				for(Read c : all){
					int ab=taxaComparator.compare(a, b);
					int bc=taxaComparator.compare(b, c);
					int ca=taxaComparator.compare(c, a);
					int zero=(ab==0 ? 1 : 0)+(bc==0 ? 1 : 0)+(ca==0 ? 1 : 0);
					int more=(ab>0 ? 1 : 0)+(bc>0 ? 1 : 0)+(ca>0 ? 1 : 0);
					int less=(ab<0 ? 1 : 0)+(bc<0 ? 1 : 0)+(ca<0 ? 1 : 0);
					assert(zero+more+less==3) : a.id+", "+b.id+", "+c.id+"; "+ab+", "+bc+", "+ca;
					assert(zero==0 || zero==1 || zero==3) : a.id+", "+b.id+", "+c.id+"; "+ab+", "+bc+", "+ca;
					if(ab==0 && bc==0){assert(ca==0) : a.id+", "+b.id+", "+c.id+"; "+ab+", "+bc+", "+ca;}
					if(zero==0){
						assert(less>0 && more>0) : a.id+", "+b.id+", "+c.id+"; "+ab+", "+bc+", "+ca;
					}else if(zero==1){
//						assert(less==2 || more==2) : a.id+", "+b.id+", "+c.id+"; "+ab+", "+bc+", "+ca;
					}
				}
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
	
	public final class TaxaComparator implements Comparator<Read>{
		
		@Override
		public int compare(Read a, Read b) {
			//if(verbose){System.err.println("Comparing "+a.id+", "+b.id);}
			int atax=GiToNcbi.getID(a.id);
			int btax=GiToNcbi.getID(b.id);
			
			if(tree!=null){
				int x=compareWithTree(atax, btax);
				if(x!=0){
					//if(verbose){System.err.println("\na) returning "+x);}
					return x;
				}
			}
			
//			System.err.println("Comparing "+atax+" to "+btax+" for reads "+a.numericID+", "+b.numericID);
			if(atax!=btax){
				//if(verbose){System.err.println("b) returning "+(atax-btax));}
				return atax-btax;
			}
			if(a.length()!=b.length()){
				//if(verbose){System.err.println("c) returning "+(b.length()-a.length()));}
				return b.length()-a.length();
			}
			//if(verbose){System.err.println("d) returning "+(a.id.compareTo(b.id)));}
			return a.id.compareTo(b.id);
		}
		
//		public int compareWithTree(int a, int b){
//			if(a==b){return 0;}
//			if(a==-1){return 1;}
//			if(b==-1){return -1;}
//			TaxNode na=tree.getNode(a);
//			TaxNode nb=tree.getNode(b);
//			if(na==nb){return 0;}
//			if(na==null){return 1;}
//			if(nb==null){return -1;}
//			while(na.level<promote){na=tree.getNode(na.pid);}
//			while(nb.level<promote){nb=tree.getNode(nb.pid);}
//			while(na.pid!=nb.pid){
//				if(na.pid==nb.id){
//					return 1;
//				}else if(nb.pid==na.id){
//					return -1;
//				}
//				if(na.level<=nb.level){
//					na=tree.getNode(na.pid);
//				}else{
//					nb=tree.getNode(nb.pid);
//				}
//			}
//			assert(na.id>-1 && nb.id>-1);
//			return na.id-nb.id;
//		}
		
//		public int compareWithTree(int a, int b){
//			if(verbose){System.err.print("e");}
//			if(a==b){return 0;}
//			if(verbose){System.err.print("f");}
//			if(a==-1){return 1;}
//			if(verbose){System.err.print("g");}
//			if(b==-1){return -1;}
//			if(verbose){System.err.print("h");}
//			TaxNode na=tree.getNode(a);
//			TaxNode nb=tree.getNode(b);
//			if(na==nb){return 0;}
//			if(verbose){System.err.print("i");}
//			if(na==null){return 1;}
//			if(verbose){System.err.print("j");}
//			if(nb==null){return -1;}
//			if(verbose){System.err.print("k");}
//			while(na.level<promote){na=tree.getNode(na.pid);}
//			while(nb.level<promote){nb=tree.getNode(nb.pid);}
//			while(na.pid!=nb.pid && na.pid!=nb.id && nb.pid!=na.id){
//				TaxNode pa=tree.getNode(na.pid);
//				TaxNode pb=tree.getNode(nb.pid);
//				if(pa.level<=pb.level){
//					if(verbose){System.err.println(na.id+", lv "+na.level+" promoted to "+na.pid+", level "+pa.level);}
//					na=pa;
//				}else{
//					if(verbose){System.err.println(nb.id+", lv "+nb.level+" promoted to "+nb.pid+", level "+pb.level);}
//					nb=pb;
//				}
//			}
//			if(na==nb){return 0;}
//			if(na.pid==nb.id){
//				if(verbose){System.err.println("\na -> b");
//				System.err.println(a+" -> "+b);
//				System.err.println(na+" -> "+nb);}
//				return 1;
//			}else if(nb.pid==na.id){
//				if(verbose){System.err.println("\nb -> a");
//				System.err.println(b+" -> "+a);
//				System.err.println(nb+" -> "+na);}
//				return -1;
//			}
//				if(verbose){System.err.print("n");}
//			assert(na.id>-1 && nb.id>-1);
//			return na.id-nb.id;
//		}
		
		public int compareWithTree(int a, int b){
			if(a==b){return 0;}
			if(a==-1){return 1;}
			if(b==-1){return -1;}
			TaxNode na=tree.getNode(a);
			TaxNode nb=tree.getNode(b);
			if(na==nb){return 0;}
			if(na==null){return 1;}
			if(nb==null){return -1;}
			while(na.level<promote){na=tree.getNode(na.pid);}
			while(nb.level<promote){nb=tree.getNode(nb.pid);}
			while(na.pid!=nb.pid && na.pid!=nb.id && nb.pid!=na.id){
				TaxNode pa=tree.getNode(na.pid);
				TaxNode pb=tree.getNode(nb.pid);
				if(pa.level<=pb.level){
					na=pa;
				}else{
					nb=pb;
				}
			}
			if(na==nb){return 0;}
			if(na.pid==nb.id){
				return 1;
			}else if(nb.pid==na.id){
				return -1;
			}
			assert(na.id>-1 && nb.id>-1);
			return na.id-nb.id;
		}
		
	}
	
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;

	private String tableFile=null;
	private String treeFile=null;
	
	private boolean fuse=false;
	private int promote=0;
	private int padding=3;
	
	private boolean addDummyReads=true;
	private int dummyLevel=TaxTree.stringToLevel("species");
	
	private final TaxTree tree;
	private final TaxaComparator taxaComparator=new TaxaComparator();
	private final byte[] dummyBases=new byte[] {'N'};
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;

	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	
	private static int MAX_FUSE_LENGTH=500000000;
	
}

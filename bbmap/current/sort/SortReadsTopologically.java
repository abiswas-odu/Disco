package sort;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.zip.ZipOutputStream;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentLegacyReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.FastqReadInputStream;
import stream.RTextInputStream;
import stream.Read;
import stream.ReadStreamStringWriter;
import stream.ReadStreamWriter;
import structures.ListNum;
import dna.AminoAcid;
import dna.Data;
import dna.Parser;
import fileIO.ReadWrite;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import fileIO.FileFormat;

public class SortReadsTopologically {
	
	
	public static void main(String[] args){
		
		Parser parser=new Parser();
		String in1=null;
		String in2=null;
		String out="raw_tsorted#.txt.gz";
		int prefix=4;
		
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if("null".equalsIgnoreCase(b)){b=null;}
//			System.err.println("Processing "+args[i]);
			
			if(Parser.isJavaFlag(arg)){
				//jvm argument; do nothing
			}else if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(a.equals("i") || a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
				in1=b;
				if(b.indexOf('#')>=0){
					in1=b.replaceFirst("#", "1");
					in2=b.replaceFirst("#", "2");
				}
			}else if(a.equals("in2") || a.equals("input2")){
				in2=b;
			}else if(a.equals("o") || a.equals("out") || a.equals("output")){
				out=b;
			}else if(a.endsWith("merge")){
				MERGE_DUPLICATES=Tools.parseBoolean(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Tools.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
				Data.sysout.println("Set overwrite to "+overwrite);
			}else if(a.equals("prefix")){
				prefix=Integer.parseInt(b);
			}else if(a.endsWith("blocksize")){
				MAX_BLOCKSIZE_TO_SORT=Integer.parseInt(b);
			}else{
				throw new RuntimeException("Unknown parameter: "+args[i]);
			}
		}
		
		if(in1==null){throw new RuntimeException("Please specify input file.");}
		if(out==null){throw new RuntimeException("Please specify output file.");}
		if(in1.equalsIgnoreCase(in2) || in1.equalsIgnoreCase(out) || (in2!=null && in2.equalsIgnoreCase(out))){
			throw new RuntimeException("Duplicate filenames.");
		}
		
		FileFormat ff=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, false);
		boolean fastq=ff.fastq();
		boolean fasta=ff.fasta();
		boolean bread=ff.bread();

		if(out!=null && !out.contains("#")){
			throw new RuntimeException("Output filename must contain '#' symbol.");
		}
		
		SortReadsTopologically srt;
		if(fasta){
			FastaReadInputStream fris1=new FastaReadInputStream(in1, (FASTQ.FORCE_INTERLEAVED && in2==null), false, true, in2==null ? Shared.READ_BUFFER_MAX_DATA : -1);
			FastaReadInputStream fris2=(in2==null ? null : new FastaReadInputStream(in2, false, false, true, -1));
			ConcurrentGenericReadInputStream cris=new ConcurrentGenericReadInputStream(fris1, fris2, -1);
			srt=new SortReadsTopologically(cris, out, prefix);
		}else if(fastq){
			FastqReadInputStream fris1=new FastqReadInputStream(in1, true);
			FastqReadInputStream fris2=(in2==null ? null : new FastqReadInputStream(in2, true));
			ConcurrentGenericReadInputStream cris=new ConcurrentGenericReadInputStream(fris1, fris2, -1);
			srt=new SortReadsTopologically(cris, out, prefix);
		}else{
			srt=new SortReadsTopologically(in1, in2, out, prefix);
		}
		
		srt.processMT();
		if(MERGE_DUPLICATES){
			Data.sysout.println("Merged "+srt.merged+" duplicates of "+srt.processed+" total.");
			if(srt.correctMerged>0 || srt.incorrectMerged>0){
				Data.sysout.println("Merged "+srt.correctMerged+" reads from same origin (correct).");
				Data.sysout.println("Merged "+srt.incorrectMerged+" reads from different origin (incorrect).");
			}
		}
	}
	
	public SortReadsTopologically(String fname1, String fname2, String outname_, int prefix_){
		assert(fname2==null || !fname1.equals(fname2)) : "Error - input files have same name.";
		
		RTextInputStream rtis=new RTextInputStream(fname1, fname2, -1);
		outname=outname_;
		paired=rtis.paired();
		cris=new ConcurrentLegacyReadInputStream(rtis, -1);
		prefix=prefix_;
		assert(prefix<=5);

		blockwriter1=(fname1==null ? null : new ReadStreamStringWriter(null, true, 4, false));
		blockwriter2=(fname2==null ? null : new ReadStreamStringWriter(null, false, 4, false));
	}
	
	public SortReadsTopologically(ConcurrentReadInputStream cris_, String outname_, int prefix_){
		cris=cris_;
		outname=outname_;
		paired=cris.paired();
		prefix=prefix_;
		assert(prefix<=5);

		blockwriter1=(new ReadStreamStringWriter(null, true, 4, false));
		blockwriter2=(!paired ? null : new ReadStreamStringWriter(null, false, 4, false));
	}

	public void processMT(){

		final String fname1=outname.replaceFirst("#", "1");
		final String fname2=(!paired ? null : outname.replaceFirst("#", "2"));
		if(fname1!=null && new File(fname1).exists()){
			if(overwrite){new File(fname1).delete();}
			else{throw new RuntimeException("Destination file "+fname1+" already exists.");}
		}
		if(fname2!=null && new File(fname2).exists()){
			if(overwrite){new File(fname1).delete();}
			else{throw new RuntimeException("Destination file "+fname2+" already exists.");}
		}
		
		Timer t=new Timer();
		Timer total=new Timer();
		t.start();
		total.start();
		
//		assert(false) : fname1+", "+fname2+", "+outname+", "+prefix;

		cris.start();
		System.err.println("Started cris");
		
		Thread bwt1=null, bwt2=null;
		if(fname1!=null){
			bwt1=new Thread(blockwriter1);
			bwt1.start();
		}
		if(fname2!=null){
			bwt2=new Thread(blockwriter2);
			bwt2.start();
		}
		System.err.println("Started blockwriters");
		
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(paired==(r.mate!=null));
				
				assert(prefix<=5);
			}
			
			while(reads!=null && reads.size()>0){
				//System.err.println("reads.size()="+reads.size());
				for(Read r : reads){addRead(r);}
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
		}

		
		synchronized(this){this.notifyAll();}
		System.err.println("Notified all");
		
		finishWritingBlocks();
		System.err.println("Wrote blocks");
		
		if(bwt1!=null){blockwriter1.poison();}
		if(bwt2!=null){blockwriter2.poison();}
//		if(bwt1!=null){blockwriter1.addList(null);}
//		if(bwt2!=null){blockwriter2.addList(null);}
		
		if(bwt1!=null){
			while(bwt1.isAlive()){
				try {
					bwt1.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		if(bwt2!=null){
			while(bwt2.isAlive()){
				try {
					bwt2.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		
		t.stop();
		Data.sysout.println("Temp Write Time: "+t);
		t.start();

		if(ReadWrite.ZIPLEVEL<2){ReadWrite.ZIPLEVEL=2;}
		ReadStreamWriter wt1=(fname1==null ? null : new ReadStreamStringWriter(fname1, true, 4, false));
		ReadStreamWriter wt2=(fname2==null ? null : new ReadStreamStringWriter(fname2, false, 4, false));

		Thread wtt1=(wt1==null ? null : new Thread(wt1));
		Thread wtt2=(wt2==null ? null : new Thread(wt2));

		if(wtt1!=null){wtt1.start();}
		if(wtt2!=null){wtt2.start();}
		
		ArrayList<String> keys=new ArrayList<String>(table.size());
		keys.addAll(table.keySet());
		Shared.sort(keys);
		
		ReadComparatorTopological tcomp=new ReadComparatorTopological();
		
		for(String key : keys){
			Block b=table.get(key);
			table.remove(key);
			processed+=b.added;
			
			if(b.added>MAX_BLOCKSIZE_TO_SORT){
				System.err.println("Skipping sorting for key "+key+" of size "+b.added);
				RTextInputStream temp=new RTextInputStream(b.fname1, b.fname2, -1);
				ArrayList<Read> reads=temp.nextList();
				while(reads!=null && reads.size()>0){
					if(reads!=null && reads.size()>0){
						if(wt1!=null){wt1.addList(reads);}
						if(wt2!=null){wt2.addList(reads);}
					}
					b.numRead+=reads.size();
					reads=temp.nextList();
				}
				temp.close();
				temp=null;
				
//				Data.sysout.println(key+"\t"+b.added);
				b.delete();
			}else{
				ArrayList<Read> list=b.readBlock();
				if(PRINT_BLOCKS){Data.sysout.println(key+"\t"+list.size());}
				b.delete();

				Shared.sort(list, tcomp);
				if(MERGE_DUPLICATES){
					int count;
					count=mergeDuplicates(list, 0, 0, (byte)-99);
					if(count>0){
						Tools.condense(list);
						Shared.sort(list, tcomp);
					}
					count=mergeDuplicates(list, 1, 0, (byte)-99);
//					if(count>0){
//						Tools.condense(list);
//						Shared.sort(list, tcomp);
//					}
//					count=mergeDuplicates(list, 0, 1, (byte)2);
					
					Tools.condense(list);
					Shared.sort(list, tcomp);
				}
				if(list!=null && list.size()>0){
					if(wt1!=null){wt1.addList(list);}
					if(wt2!=null){wt2.addList(list);}
				}
			}
		}
		
		//Add poison
//		if(wt1!=null){wt1.addList(null);}
//		if(wt2!=null){wt2.addList(null);}
		if(wt1!=null){wt1.poison();}
		if(wt2!=null){wt2.poison();}
		
		if(wtt1!=null){
			while(wtt1.isAlive()){
				try {
					wtt1.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		
		if(wtt2!=null){
			while(wtt2.isAlive()){
				try {
					wtt2.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		
		t.stop();
		total.stop();
		Data.sysout.println("Final Sort + Write Time: "+t);
		Data.sysout.println("Total Time: "+total);
		
	}

	
	private int mergeDuplicates(ArrayList<Read> list, int nmax, int mmax, byte qmax){
		if(list==null || list.size()<2){return 0;}
		Read current=list.get(0);
		
		int correct=0;
		int incorrect=0;
		
		int count=0;
		for(int i=1; i<list.size(); i++){
			Read r=list.get(i);
			boolean keep=false;
			if(current.length()==r.length() && ((current.mate==null && r.mate==null) || (current.mateLength()==r.mateLength())) && 
					current.isDuplicateByBases(r, nmax, mmax, qmax, true, false)){
				Read r2=r.mate;
				Read c2=current.mate;
				if(c2==null || c2.isDuplicateByBases(r2, nmax, mmax, qmax, true, false)){
//					assert(r.synthetic()) : r;
					if(r.synthetic()){
						if(r.originalSite!=null && current.originalSite!=null){
							if(r.originalSite.equals(current.originalSite)){
								correct++;
							}else{
								incorrect++;
							}
						}else if(r.chrom>0){
							if(r.chrom==current.chrom && r.start==current.start && r.stop==current.stop && r.strand()==current.strand()){
								correct++;
							}else{
								incorrect++;
							}
						}
						if(r2!=null && c2!=null && r2.originalSite!=null && c2.originalSite!=null){
							if(r2.originalSite.equals(c2.originalSite)){
								correct++;
							}else{
								incorrect++;
							}
						}else if(r2!=null && c2!=null && r2.chrom>0){
							if(r2.chrom==c2.chrom && r2.start==c2.start && r2.stop==c2.stop && r2.strand()==c2.strand()){
								correct++;
							}else{
								incorrect++;
							}
						}
					}
					current.merge(r, true, true);
					list.set(i, null);
					count++;
					keep=true;
				}
			}
			if(!keep){current=r;}
		}
		merged+=count;
		correctMerged+=correct;
		incorrectMerged+=incorrect;
		return count;
	}
	
	
	private void addRead(Read r){
		StringBuilder sb=new StringBuilder(prefix);
		boolean bad=false;
		for(int i=0; i<prefix && i<r.length(); i++){
			byte b=r.bases[i];
			
			if(b>=0 && b<=3){
				sb.append((int)b);
			}else{

				if(AminoAcid.isFullyDefined(b)){
					sb.append((char)b);
				}else{
					bad=true;
					sb.append('N');
				}
			}
			
		}
		
		String key=bad ? "ZN" : sb.toString();
//		String key=sb.toString();
		Block b=table.get(key);
		if(b==null){
			//System.err.println("Created block "+key);
			b=new Block(key, outname);
			table.put(key, b);
		}
		b.add(r);
	}
	
	
	public void finishWritingBlocks(){
		System.err.println("Called finishWritingBlocks()");
		for(String key : table.keySet()){
			Block b=table.get(key);
			b.finishWritingBuffer();
		}
	}
	

	
	private static final void finishWriting(PrintWriter writer, OutputStream outStream){
		writer.flush();
		if(outStream.getClass()==ZipOutputStream.class){
			ZipOutputStream zos=(ZipOutputStream)outStream;
			try {
				zos.closeEntry();
				zos.finish();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		writer.close();
		try {
			outStream.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private class Block{
		
		public Block(String name_, String fname_){
			
			if(DONT_COMPRESS_TEMP_FILES){
				while(fname_.endsWith(".gz") || fname_.endsWith(".zip") || fname_.endsWith(".bz2")){
					fname_=fname_.substring(0, fname_.lastIndexOf('.'));
				}
			}
			
			name=name_;
			fname1=fname_.replaceFirst("#", "_tsort_tempBlock_"+name+"_1");
			fname2=(!paired ? null : fname_.replaceFirst("#", "_tsort_tempBlock_"+name+"_2"));
			
			if(fname1==null){
				assert(false);
				outStream1=null;
				writer1=null;
			}else{
				outStream1=ReadWrite.getOutputStream(fname1, append, true, false);
				writer1=new PrintWriter(outStream1);
			}
			
			if(fname2==null){
				outStream2=null;
				writer2=null;
			}else{
				outStream2=ReadWrite.getOutputStream(fname2, append, true, false);
				writer2=new PrintWriter(outStream2);
			}
		}

		public void add(Read r){
			buffer.add(r);
			added++;
			if(buffer.size()>=WRITE_BUFFER){
				writeBuffer(false);
			}
		}
		
		public void writeBuffer(boolean close){
			
			written+=buffer.size();
			ArrayList<Read> temp=buffer;
			buffer=(close ? null : new ArrayList<Read>(WRITE_BUFFER));
			
			if(close){
//				System.err.println("Closing "+name+": "+ fname1+", "+fname2);
				if(blockwriter1!=null){blockwriter1.addList(temp, writer1, outStream1, close);}
				if(blockwriter2!=null){blockwriter2.addList(temp, writer2, outStream2, close);}
			}else{
				if(blockwriter1!=null && temp!=null && !temp.isEmpty()){blockwriter1.addList(temp, writer1, outStream1, close);}
				if(blockwriter2!=null && temp!=null && !temp.isEmpty()){blockwriter2.addList(temp, writer2, outStream2, close);}
			}
			
			assert(added==written);
//			buffer.clear();
		}
		
		public void finishWritingBuffer(){
			//System.err.println("Writing block "+name);
			writeBuffer(true);

//			finishWriting(writer1, outStream1);
//			if(fname2!=null){
//				finishWriting(writer2, outStream2);
//			}
			
		}
		
		public synchronized ArrayList<Read> readBlock(){
			RTextInputStream temp=new RTextInputStream(fname1, fname2, -1);
			ArrayList<Read> out=new ArrayList<Read>((int)written);
			ArrayList<Read> reads=temp.nextList();
			while(reads!=null && reads.size()>0){
				out.addAll(reads);
				numRead+=reads.size();
				reads=temp.nextList();
			}
			temp.close();
			temp=null;
			assert(numRead==written);
			
			return out;
		}
		
		public synchronized void delete() {
			if(fname1!=null){new File(fname1).delete();}
			if(fname2!=null){new File(fname2).delete();}
		}
		
		public final String name;
		public final String fname1, fname2;
		
		public final OutputStream outStream1, outStream2;
		public final PrintWriter writer1, writer2;
		private ArrayList<Read> buffer=new ArrayList<Read>(WRITE_BUFFER);
		
		public long added=0, written=0, numRead=0;
	}
	
	public final String outname;
	private final ConcurrentReadInputStream cris;
	public long merged=0;
	public long processed=0;

	public long correctMerged=0;
	public long incorrectMerged=0;
	
	private final HashMap<String, Block> table=new HashMap<String, Block>(4096);
	
	public final boolean paired;
	public final int prefix;
	
	public static boolean USE_CRIS=true; //Similar speed either way.  "true" may be better with many threads.
	
	public static final int WRITE_BUFFER=1000; //Bigger number uses more memory, for less frequent writes.
	public static int MAX_BLOCKSIZE_TO_SORT=16000000;

	public static final boolean DONT_COMPRESS_TEMP_FILES=false;
	public static boolean MERGE_DUPLICATES=false;
	public static boolean overwrite=false;
	public static boolean append=false;
	public static boolean PRINT_BLOCKS=false;
	

	private final ReadStreamWriter blockwriter1;
	private final ReadStreamWriter blockwriter2;
	
	
}

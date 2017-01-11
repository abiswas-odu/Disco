package shared;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Random;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.concurrent.atomic.AtomicLongArray;
import java.util.regex.Pattern;

import stream.KillSwitch;
import stream.Read;
import stream.ReadInputStream;
import stream.SamLine;
import stream.SiteScore;
import structures.CoverageArray;
import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;

public final class Tools {
	
	public static long estimateFileSize(String fname){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, false, false);
		if(ff==null || ff.stdio()){return -1;}

		double mult=1;
		if(ff.compressed()){
			if(ff.bam()){mult=6;}
			if(ff.fasta()){mult=5;}
			else{mult=4;}
		}

		File f=new File(fname);
		long size=f.length();
		double diskEstimate=size*mult;
		return (long)diskEstimate;
	}
	
	/**
	 * @param fname
	 * @param readsToFetch
	 * @param extraOverheadPerRead
	 * @param earlyExit
	 * @return {memEstimate, diskEstimate, memRatio, diskRatio};
	 */
	public static double[] estimateFileMemory(String fname, int readsToFetch, double extraOverheadPerRead, boolean earlyExit){
		
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, false, false);
		if(ff==null || ff.stdio()){return null;}

		File f=new File(fname);
		final long size=f.length();
		
		if(earlyExit){
			
			long available=Shared.memFree();
			
			
			double memRatio, diskRatio;
			if(ff.compressed()){
				memRatio=40;
				diskRatio=5;
			}else{
				memRatio=8;
				diskRatio=1;
			}
			
			long memEstimate=(long)(memRatio*size);
			long diskEstimate=(long)(diskRatio*size);
			
//			System.err.println(memEstimate+", "+available+", "+(memEstimate*1.5));
			
			if(memEstimate*2.1<available){
				return new double[] {memEstimate, diskEstimate, memRatio, diskRatio};
			}
		}
		
		readsToFetch=Tools.max(readsToFetch, 200);
		ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, ff.compressed() && !ff.gzip(), true);
		ArrayList<Read> reads=ReadInputStream.toReads(ff, readsToFetch);
		long minBytes=Integer.MAX_VALUE;
		long maxBytes=1;
		long sumBytes=0;
		long minMem=Integer.MAX_VALUE;
		long maxMem=1;
		long sumMem=0;
		long minLen=Integer.MAX_VALUE;
		long maxLen=1;
		long sumLen=0;
		long minQLen=Integer.MAX_VALUE;
		long maxQLen=1;
		long sumQLen=0;
		long minHdr=Integer.MAX_VALUE;
		long maxHdr=1;
		long sumHdr=0;
		long readCount=0;
		
		if(reads==null || reads.size()<1){
			minBytes=maxBytes=minLen=maxLen=minMem=maxMem=minHdr=maxHdr=1;
		}else{
			for(Read r1 : reads){
				long x;
				readCount++;
				
				x=r1.length();
				minLen=min(minLen, x);
				maxLen=max(maxLen, x);
				sumLen+=x;
				
				x=r1.qlength();
				minQLen=min(minQLen, x);
				maxQLen=max(maxQLen, x);
				sumQLen+=x;
				
				x=r1.countBytes();
				minMem=min(minMem, x);
				maxMem=max(maxMem, x);
				sumMem+=x;
				
				x=r1.id.length();
				minHdr=min(minHdr, x);
				maxHdr=max(maxHdr, x);
				sumHdr+=x;
				
				x=r1.countFastqBytes();
				minBytes=min(minBytes, x);
				maxBytes=max(maxBytes, x);
				sumBytes+=x;
				
				Read r2=r1.mate;
				if(r2!=null){
					readCount++;
					
					x=r2.length();
					minLen=min(minLen, x);
					maxLen=max(maxLen, x);
					sumLen+=x;
					
					x=r2.qlength();
					minQLen=min(minQLen, x);
					maxQLen=max(maxQLen, x);
					sumQLen+=x;
					
					x=r2.countBytes();
					minMem=min(minMem, x);
					maxMem=max(maxMem, x);
					sumMem+=x;
					
					x=r2.id.length();
					minHdr=min(minHdr, x);
					maxHdr=max(maxHdr, x);
					sumHdr+=x;
					
					x=r2.countFastqBytes();
					minBytes=min(minBytes, x);
					maxBytes=max(maxBytes, x);
					sumBytes+=x;
				}
			}
		}
		
		double compressedSize;
		if(ff.compressed()){
			compressedSize=0.125*( //bytes per bit
					0.2*sumHdr //Repetitive header characters
					+1.5*sumLen
					+2.5*sumQLen
					+(4*6*readCount) //@, +, and newline
					+(8*2*readCount) //Actual information content of a typical header
					);
			if(ff.bz2() || ff.fqz()){
				compressedSize*=0.83;
			}
		}else{
			compressedSize=sumBytes;
		}
		double memRatio=(sumMem+readCount*extraOverheadPerRead)/compressedSize;
		double diskRatio=sumBytes/compressedSize;
		
		long memEstimate=(long)(memRatio*size);
		long diskEstimate=(long)(diskRatio*size);

//		System.err.println("compressedSize="+compressedSize);
//		System.err.println("memRatio="+memRatio);
//		System.err.println("diskRatio="+diskRatio);
		
//		assert(false) : diskEstimate+", "+minBytes+", "+
//		double worstCase=estimate*1.75;
		return new double[] {memEstimate, diskEstimate, memRatio, diskRatio};
	}

	public static int secondHighestPosition(int[] array) {
		int maxP, maxP2;
		if(array[0]>=array[1]){
			maxP=0;
			maxP2=1;
		}else{
			maxP=1;
			maxP2=0;
		}
		for(int i=2; i<array.length; i++){
			int x=array[i];
			if(x>array[maxP2]){
				if(x>=array[maxP]){
					maxP2=maxP;
					maxP=i;
				}else{
					maxP2=i;
				}
			}
		}
		return maxP2;
	}
	
	public static final boolean nextBoolean(Random randy){
		return randy.nextBoolean();
//		int r=randy.nextInt()&0x7FFFFFFF;
//		return r%294439>=147219;
	}
	
	/**
	 * @param probError
	 * @return
	 */
	public static float[] inverse(float[] array) {
		float[] out=new float[array.length];
		for(int i=0; i<array.length; i++){
			out[i]=1/max(array[i], 1000000000f);
		}
		return out;
	}
	
	public static boolean checkHeader(String s){
		if(s==null){return true;}
		boolean ok=true;
		for(int i=0; i<s.length() && ok; i++){
			char c=s.charAt(i);
			ok=(c>=32 && c<=126);
		}
		return ok;
	}
	
	public static String fixHeader(String s){
//		assert(false) : new String(specialChars);
		if(checkHeader(s)){return s;}
		StringBuilder sb=new StringBuilder(s.length());
		for(int i=0; i<s.length(); i++){
			final char c=s.charAt(i), d;
			
			if(c>=0 && c<=255){
				d=specialChars[c];
			}else{
				d='X';
			}
//			System.err.println(c+"="+(int)c);
			sb.append(d);
		}
		return sb.toString();
	}
	
	
	/**
	 * Returns this file name if it is a file, or all the files in the directory if it is a directory.
	 * @param b
	 * @param fasta
	 * @param fastq
	 * @param sam
	 * @param any
	 * @return
	 */
	public static ArrayList<String> getFileOrFiles(String b, ArrayList<String> list, boolean fasta, boolean fastq, boolean sam, boolean any){
		if(list==null){list=new ArrayList<String>();}
		String[] split=b.split(",");
		for(String s : split){
			File f=new File(s);
			if(f.isDirectory()){
				for(File f2 : f.listFiles()){
					if(f2.isFile()){
						String name=f2.getName().toLowerCase();
						String ext=ReadWrite.rawExtension(name);
						
						boolean pass=any || (fasta && FileFormat.isFasta(ext)) || (fastq && FileFormat.isFastq(ext)) || (sam && FileFormat.isSamOrBam(ext));
						
						if(pass){
							String s2=f2.getAbsolutePath();
							list.add(s2);
						}
					}
				}
			}else{
				list.add(s);
			}
		}
		return list;
	}
	
	public static ArrayList<byte[]> toAdapterList(String name, int maxLength){
		if(maxLength<1){maxLength=Integer.MAX_VALUE;}
		if(name==null){return null;}
		String[] split;
		ArrayList<byte[]> list=new ArrayList<byte[]>();
		if(new File(name).exists()){
			split=new String[] {name};
		}else{
			split=name.split(",");
		}
		for(String s : split){
			if(new File(s).exists()){
				ArrayList<Read> reads=ReadInputStream.toReads(s, FileFormat.FASTA, -1);
				for(Read r : reads){
					if(r!=null && r.length()>0){
						byte[] array=checkAdapter(r.bases, maxLength);
						if(array.length>0){list.add(array);}
					}
				}
			}else{
				byte[] array=checkAdapter(s.getBytes(), maxLength);
				if(array.length>0){list.add(array);}
			}
		}
		return list==null || list.isEmpty() ? null : list;
	}
	
	private static byte[] checkAdapter(byte[] array, int maxLength){
		if(array.length>maxLength){array=Arrays.copyOf(array, maxLength);}
		
		for(int i=0; i<array.length; i++){
			byte b=array[i];
			int x=AminoAcid.baseToNumberExtended[b];
			if(x<0 || !Character.isLetter(b) || !Character.isUpperCase(b)){
				throw new RuntimeException("Invalid nucleotide "+(char)b+" in literal sequence "+new String(array)+"\n"
					+ "If this was supposed to be a filename, the file was not found.");
			}
			if(AminoAcid.baseToNumber[b]<0){array[i]='N';}//Degenerate symbols become N
		}
		
		int trailingNs=0;
		for(int i=array.length-1; i>=0; i--){
			if(array[i]=='N'){trailingNs++;}
		}
		if(trailingNs>0){
			array=Arrays.copyOf(array, array.length-trailingNs);
		}
		return array;
	}
	
	public static byte[][] toAdapters(String name, final int maxLength){
		ArrayList<byte[]> list=toAdapterList(name, maxLength);
		return list==null ? null : list.toArray(new byte[list.size()][]);
	}
	
	/** Add names to a collection.
	 * This can be a literal name, or a text file with one name per line, 
	 * or a fastq, fasta, or sam file, in which case the read names will be added.
	 * @param s
	 * @param names
	 * @return
	 */
	public static final int addNames(String s, Collection<String> names, boolean allowSubprocess){
		int added=0;
		if(new File(s).exists()){

			int[] vector=FileFormat.testFormat(s, false, false);
			final int type=vector[0];
			ByteFile bf=ByteFile.makeByteFile(s, false, allowSubprocess);
			
			if(type==FileFormat.FASTQ){
				int num=0;
				for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine(), num++){
					if((num&3)==0 && line.length>0){
						names.add(new String(line, 1, line.length-1));
					}
				}
			}else if(type==FileFormat.FASTA){
				for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
					if(line.length>0 && line[0]=='>'){
						names.add(new String(line, 1, line.length-1));
					}
				}
			}else if(type==FileFormat.SAM){
				for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
					if(line.length>0 && line[0]!='@'){
						String name=SamLine.parseNameOnly(line);
						if(name!=null && name.length()>0){names.add(name);}
					}
				}
			}else{
				for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
					if(line.length>0){
						names.add(new String(line));
					}
				}
			}
			bf.close();
		}else{
			added++;
			names.add(s);
		}
		return added;
	}
	
	/**
	 * Make copies of any read with ambiguous bases to represent all possible non-ambiguous representations.
	 * @param reads A list of reads
	 * @param k minimum length of reads to replicate
	 * @return A list of reads with no ambiguity codes.
	 */
	public static ArrayList<Read> replicateAmbiguous(ArrayList<Read> reads, int minlen) {
		ArrayList<Read> out=new ArrayList<Read>();
		for(Read r1 : reads){
			final Read r2=r1.mate;
			r1.mate=null;
			
			if(r1.containsUndefined() && r1.length()>=minlen){
				ArrayList<Read> temp=makeReplicates(r1);
				out.addAll(temp);
			}else{
				out.add(r1);
			}
			if(r2!=null){
				r2.mate=null;
				if(r2.containsUndefined() && r2.length()>=minlen){
					ArrayList<Read> temp=makeReplicates(r2);
					out.addAll(temp);
				}else{
					out.add(r2);
				}
			}
		}
		return out;
	}
	
	/**
	 * Make copies of this read to represent all possible non-ambiguous representations.
	 * Return a list of all fully-defined versions.
	 * @param r A read to replicate
	 * @return A list of reads with no ambiguity codes.
	 */
	public static ArrayList<Read> makeReplicates(final Read r) {
//		System.err.println("\n***Called makeReplicates("+new String(r.bases)+")");
		ArrayList<Read> temp=null;
		if(!r.containsUndefined()){
			temp=new ArrayList<Read>();
			temp.add(r);
			return temp;
		}
		final byte[] bases=r.bases;
		for(int i=0; i<r.bases.length; i++){
			byte b=bases[i];
			if(!AminoAcid.isFullyDefined(b)){
				temp=replicateAtPosition(r, i);
				break;
			}
		}
		assert(temp!=null);
		final ArrayList<Read> out;
		if(temp.get(0).containsUndefined()){
			out=new ArrayList<Read>();
			for(Read rr : temp){
				out.addAll(makeReplicates(rr));
			}
		}else{
			out=temp;
		}
		return out;
	}
	
	/**
	 * @param r A read
	 * @param pos The position of an ambiguous base
	 * @param out A list of replicates
	 */
	private static ArrayList<Read> replicateAtPosition(final Read r, final int pos) {
//		System.err.println("Called replicateAtPosition("+new String(r.bases)+", "+pos+")");
		if(r.quality!=null){
			r.quality[pos]=Shared.FAKE_QUAL;
		}
		final byte[] bases=r.bases;
		final byte b=bases[pos];
		final int num=AminoAcid.baseToNumberExtended[b]&0xF;
		assert(num>0 && Integer.bitCount(num)>1 && Integer.bitCount(num)<=4) : b+", "+num;
		ArrayList<Read> out=new ArrayList<Read>(4);
		for(int i=0; i<4; i++){
			int mask=(1<<i);
			if((num&mask)==mask){
				Read rr=r.clone();
				rr.bases=rr.bases.clone();
				rr.bases[pos]=AminoAcid.numberToBase[i];
//				System.err.println("Added clone ("+new String(rr.bases)+")");
				out.add(rr);
			}
		}
		return out;
	}
	
	/** Checks for permission to read files, and input name collisions. */
	public static boolean testOutputFiles(boolean overwrite, boolean append, boolean allowDuplicates, ArrayList<String>...args){
		if(args==null || args.length==0){return true;}
		ArrayList<String> list=new ArrayList<String>();
		for(ArrayList<String> als : args){
			if(als!=null){
				list.addAll(als);
			}
		}
		return testOutputFiles(overwrite, append, allowDuplicates, list.toArray(new String[list.size()]));
	}
	
	/** Checks for permission to overwrite files, and output name collisions. */
	public static boolean testOutputFiles(boolean overwrite, boolean append, boolean allowDuplicates, String...args){
		if(args==null || args.length==0){return true;}
		HashSet<String> set=new HashSet<String>(args.length*2);
		int terms=0;
		for(String s : args){
			if(s!=null){
				if(isOutputFileName(s)){
					terms++;
					
					if(!overwrite && !append && new File(s).exists()){
						assert(overwrite) : "File "+s+" exists and overwrite=false";
						return false;
					}
					
					if(!allowDuplicates && set.contains(s)){
						assert(false) : "Duplicate file "+s+" was specified for multiple output streams.";
						return false;
					}
					
					set.add(s);
				}		
			}
		}
		return true;
	}
	
	/** Checks for permission to read files, and input name collisions. */
	public static boolean testInputFiles(boolean allowDuplicates, boolean throwException, ArrayList<String>...args){
		if(args==null || args.length==0){return true;}
		ArrayList<String> list=new ArrayList<String>();
		for(ArrayList<String> als : args){
			if(als!=null){
				list.addAll(als);
			}
		}
		return testInputFiles(allowDuplicates, throwException, list.toArray(new String[list.size()]));
	}
	
	/** Checks for permission to read files, and input name collisions. */
	public static boolean testInputFiles(boolean allowDuplicates, boolean throwException, String[]...args){
		if(args==null || args.length==0){return true;}
		for(String[] s : args){
			if(!testInputFiles(allowDuplicates, throwException, s)){return false;}
		}
		return true;
	}
	
	/** Checks for permission to read files, and input name collisions. */
	public static boolean testInputFiles(boolean allowDuplicates, boolean throwException, String...args){
		if(args==null || args.length==0){return true;}
		HashSet<String> set=new HashSet<String>(args.length*2);
		int terms=0;
		for(String s : args){
			if(s!=null){
				String s2=s.toLowerCase();
				if(canRead(s)){
					terms++;
				}else{
					if(throwException){throw new RuntimeException("Can't read file '"+s+"'");}
					return false;
				}

				if(!allowDuplicates && set.contains(s2)){
					if(throwException){throw new RuntimeException("Duplicate file "+s+" was specified for multiple input streams.");}
					return false;
				}

				set.add(s2);
			}
		}
		return true;
	}
	
	/** Checks for permission to overwrite files, and output name collisions. 
	 * @return True if no problems are detected */
	public static boolean testForDuplicateFiles(boolean throwException, String...args){
		if(args==null || args.length==0){return true;}
		HashSet<String> set=new HashSet<String>(args.length*2);
		int terms=0;
		for(String s0 : args){
			if(s0!=null){
				String s=s0.toLowerCase();
				terms++;
				if(set.contains(s) && !s.equals("stdout") && !s.startsWith("stdout.")){
					if(throwException){throw new RuntimeException("File '"+s0+"' was specified multiple times.");}
					return false;
				}
				set.add(s);
			}
		}
		return true;
	}
	
	public static final boolean canWrite(String s, boolean overwrite){
		if(isNullFileName(s) || isSpecialOutputName(s)){return true;}
		File f=new File(s);
		if(f.exists()){return overwrite && f.canWrite();}
		return true;
	}
	
//	public static final boolean outputDestinationExists(String s){
//		if(isNullFileName(s)){return false;}
//		if(isSpecialOutputName(s)){return false;}
//		File f=new File(s);
//		return f.exists();
//	}
	
	public static final boolean isOutputFileName(String s){
		return !(isNullFileName(s) || isSpecialOutputName(s)); 
	}
	
	public static final boolean isNullFileName(String s){
		if(s==null || s.equalsIgnoreCase("null") || s.equalsIgnoreCase("none")){return true;}
		for(int i=0; i<s.length(); i++){
			if(!Character.isWhitespace(s.charAt(i))){return false;}
		}
		return true;
	}
	
	public static final boolean isSpecialOutputName(String s){
		if(s==null){return false;}
		s=s.toLowerCase();
		return s.equals("stdout") || s.equals("stderr") || s.equals("standardout") || s.equals("standarderr")
				|| s.equals("/dev/null") || s.startsWith("stdout.") || s.startsWith("stderr.");
	}
	
	public static final boolean isSpecialInputName(String s){
		if(s==null){return false;}
		s=s.toLowerCase();
		return s.equals("stdin") || s.equals("standardin") || s.startsWith("stdin.");
	}
	
	public static final boolean canRead(String s){
		if(s==null){return false;}
		if(isSpecialInputName(s)){return true;}
		File f=new File(s);
		return f.canRead();
	}
	
	/** Returns index of first matching location */
	public static final int locationOf(final byte[] big, final byte[] small, final int maxMismatches){
		int x=containsForward(big, small, maxMismatches);
		return x>=0 ? x : containsReverse(big, small, maxMismatches);
	}
	
	/** Returns index of first matching location */
	public static final int containsForward(final byte[] big, final byte[] small, final int maxMismatches){
		final int ilimit=big.length-small.length;
//		System.err.println("Entering: ilimit="+ilimit+", maxMismatches="+maxMismatches+", small.length="+small.length);
		for(int i=0; i<=ilimit; i++){
			int mismatches=0;
			for(int j=0; j<small.length && mismatches<=maxMismatches; j++){
				final byte b=big[i+j];
				final byte s=small[j];
				if(b!=s){mismatches++;}
			}
			if(mismatches<=maxMismatches){
//				System.err.println("Returning "+i+", mismatches="+mismatches);
				return i;
			}
		}
		return -1;
	}
	
	/** Returns index of first matching location */
	public static final int containsReverse(final byte[] big, final byte[] small, final int maxMismatches){
		final int ilimit=big.length-small.length;
		for(int i=0; i<=ilimit; i++){
			int mismatches=0;
			for(int j=0, k=small.length-1; j<small.length && mismatches<=maxMismatches; j++, k--){
				final byte b=big[i+j];
				final byte s=AminoAcid.baseToComplementExtended[small[k]];
				if(b!=s){mismatches++;}
			}
			if(mismatches<=maxMismatches){return i;}
		}
		return -1;
	}
	
	/** Removes null elements by shrinking the list.  May change list order. */
	public static final <X> int condense(ArrayList<X> list){
		if(list==null || list.size()==0){return 0;}
		int removed=0;
		
		for(int i=list.size()-1; i>0; i--){
			if(list.get(i)==null){
				removed++;
				X last=list.get(list.size()-1);
				list.set(i, last);
				list.remove(list.size()-1);
			}
		}
		return removed;
	}
	
	/** Removes null elements by shrinking the list.  Will not change list order. */
	public static final <X> int condenseStrict(ArrayList<X> list){
		if(list==null || list.size()==0){return 0;}
		int removed=0;
		
		int insertPos=0;
		for(int i=0; i<list.size(); i++){
			X x=list.get(i);
			if(x!=null){
				if(insertPos!=i){
					assert(insertPos<i);
					while(list.get(insertPos)!=null){insertPos++;}
					assert(insertPos<i && list.get(insertPos)==null) : insertPos+", "+i; //slow, temporary
					list.set(i, null);
					list.set(insertPos, x);
				}
				insertPos++;
			}else{
				removed++;
			}
		}
		for(int i=0; i<removed; i++){
			X x=list.remove(list.size()-1);
			assert(x==null);
		}
		return removed;
	}
	
	/** Removes null elements by shrinking the array.  Will not change array order. */
	public static final <X> X[] condenseStrict(X[] array){
		if(array==null){return array;}
		int nulls=0;
		for(X x : array){if(x==null){nulls++;}}
		if(nulls==0){return array;}
		X[] array2=Arrays.copyOf(array, array.length-nulls);
		
		int j=0;
		for(X x : array){
			if(x!=null){
				array2[j]=x;
				j++;
			}
		}
		return array2;
	}
	
	/** Creates a new list without null elements. */
	public static final <X> ArrayList<X> condenseNew(ArrayList<X> list){
		ArrayList<X> temp=new ArrayList<X>(list.size());
		for(X x : list){
			if(x!=null){temp.add(x);}
		}
		return temp;
	}
	
	//This should also be correct.  I'm not sure which is faster.
//	/** Removes null elements by shrinking the list.  Will not change list order. */
//	public static final <X> int condenseStrict(ArrayList<X> list){
//		if(list==null || list.size()==0){return 0;}
//		int removed=0;
//		int last=0;
//		
//		for(int i=0; i<list.size(); i++){
//			X x=list.get(i);
//			if(x==null){
//				removed++;
//			}else{
//				while(last<i && list.get(last)!=null){last++;}
//				assert(last==i || list.get(last)==null);
//				if(last!=i){
//					assert(last<i);
//					list.set(last, x);
//					list.set(i, null);
//				}
//			}
//		}
//		for(int i=0; i<removed; i++){
//			X x=list.remove(list.size()-1);
//			assert(x==null);
//		}
//		return removed;
//	}
	

//	public static final int trimSiteList(ArrayList<SiteScore> ssl, float fractionOfMax, boolean retainPaired){
////		assert(false);
//		if(ssl==null || ssl.size()==0){return -999999;}
//		if(ssl.size()==1){return ssl.get(0).score;}
//		int maxScore=-999999;
//		for(SiteScore ss : ssl){
//			maxScore=Tools.max(maxScore, ss.score);
//		}
//		
//		int cutoff=(int) (maxScore*fractionOfMax);
//		trimSitesBelowCutoff(ssl, cutoff, retainPaired);
////		trimSitesBelowCutoffInplace(ssl, cutoff);
//		return maxScore;
//	}
	
	/** minSitesToRetain should be set to 1 if the list is not sorted by score (for efficiency of removal).  Otherwise, it can be higher. */
	public static final int trimSiteList(ArrayList<SiteScore> ssl, float fractionOfMax, boolean retainPaired, boolean retainSemiperfect, 
			int minSitesToRetain, int maxSitesToRetain){
//		assert(false);
		if(ssl==null || ssl.size()==0){return -999999;}
		if(ssl.size()==1){return ssl.get(0).score;}
		int maxScore=-999999;
		
		if(minSitesToRetain>1 && minSitesToRetain<ssl.size()){
			assert(inOrder(ssl));
			maxScore=ssl.get(0).score;
		}else{
			for(SiteScore ss : ssl){
				maxScore=Tools.max(maxScore, ss.score);
			}
		}
		
		int cutoff=(int) (maxScore*fractionOfMax);
		trimSitesBelowCutoff(ssl, cutoff, retainPaired, retainSemiperfect, minSitesToRetain, maxSitesToRetain);
		return maxScore;
	}
	
	/** minSitesToRetain should be set to 1 if the list is not sorted by score.  Otherwise, it can be higher. */
	public static final void trimSiteListByMax(ArrayList<SiteScore> ssl, int cutoff, boolean retainPaired, boolean retainSemiperfect, 
			int minSitesToRetain, int maxSitesToRetain){
//		assert(false);
		if(ssl==null || ssl.size()==0){return;}
		if(ssl.size()==1){return;}
		
		trimSitesBelowCutoff(ssl, cutoff, retainPaired, retainSemiperfect, minSitesToRetain, maxSitesToRetain);
	}
	
	public static final <X extends Comparable<? super X>> boolean inOrder(ArrayList<X> list){
		if(list==null || list.size()<2){return true;}
		for(int i=1; i<list.size(); i++){
			X xa=list.get(i-1);
			X xb=list.get(i);
			if(xa.compareTo(xb)>0){return false;}
		}
		return true;
	}
	

	
	public static final int mergeDuplicateSites(ArrayList<SiteScore> list, boolean doAssertions, boolean mergeDifferentGaps){
		if(list==null || list.size()<2){return 0;}
		Shared.sort(list, SiteScore.PCOMP);
		
		int removed=0;
		
		SiteScore a=list.get(0);
		for(int i=1; i<list.size(); i++){
			SiteScore b=list.get(i);
			if(a.positionalMatch(b, true)){
				
				if(doAssertions){
					if(!(a.perfect==b.perfect || 
							(a.perfect && (a.score>b.score || a.slowScore>b.slowScore)))){
						throw new RuntimeException("\n"+SiteScore.header()+"\n"+a.toText()+"\n"+b.toText()+"\n");
					}

					assert(a.perfect==b.perfect || 
							(a.perfect && (a.score>b.score || a.slowScore>b.slowScore))) : 
								"\n"+SiteScore.header()+"\n"+a.toText()+"\n"+b.toText()+"\n";
				}
				
				a.setSlowScore(max(a.slowScore, b.slowScore));
//				a.setPairedScore(a.pairedScore<=0 && b.pairedScore<=0 ? 0 : max(a.slowScore+1, a.pairedScore, b.pairedScore));
				a.setPairedScore((a.pairedScore<=a.slowScore && b.pairedScore<=a.slowScore) ? 0 : max(0, a.pairedScore, b.pairedScore));
				a.setScore(max(a.score, b.score));
				a.perfect=(a.perfect || b.perfect);
				a.semiperfect=(a.semiperfect || b.semiperfect);
				
				removed++;
				list.set(i, null);
			}else if(mergeDifferentGaps && a.positionalMatch(b, false)){ //Same outermost boundaries, different gaps
				
				SiteScore better=null;
				if(a.score!=b.score){
					better=(a.score>b.score ? a : b);
				}else if(a.slowScore!=b.slowScore){
					better=(a.slowScore>b.slowScore ? a : b);
				}else if(a.pairedScore!=b.pairedScore){
					better=(a.pairedScore>b.pairedScore ? a : b);
				}else{
					better=a;
				}
				
				a.setSlowScore(max(a.slowScore, b.slowScore));
				a.setPairedScore((a.pairedScore<=a.slowScore && b.pairedScore<=a.slowScore) ? 0 : max(0, a.pairedScore, b.pairedScore));
				a.setScore(max(a.score, b.score));
				a.perfect=(a.perfect || b.perfect);//TODO: This is not correct.  And perfect sites should not have gaps anyway.
				a.semiperfect=(a.semiperfect || b.semiperfect);
				a.gaps=better.gaps;
				
				removed++;
				list.set(i, null);
			}
			else{
				a=b;
			}
		}

//		if(removed>0){condense(list);}
		if(removed>0){condenseStrict(list);}
		return removed;
	}
	

	
	public static final int subsumeOverlappingSites(ArrayList<SiteScore> list, boolean subsumeIfOnlyStartMatches, boolean subsumeInexact){
		if(list==null || list.size()<2){return 0;}
		Shared.sort(list, SiteScore.PCOMP);
		
		int removed=0;
		
		
		for(int i=0; i<list.size(); i++){
			SiteScore a=list.get(i);
			
			assert(a==null || !a.perfect || a.semiperfect);
			
			boolean overlappingA=true;
			if(a!=null){
				for(int j=i+1; overlappingA && j<list.size(); j++){
					SiteScore b=list.get(j);
					assert(b==null || !b.perfect || b.semiperfect);
					if(b!=null){
						overlappingA=(a.chrom==b.chrom && b.start<a.stop && b.stop>a.start);
						if(overlappingA && a.strand==b.strand){
							
							SiteScore better=null;
							if(a.perfect!=b.perfect){
								better=a.perfect ? a : b;
							}else if(a.semiperfect!=b.semiperfect){
								better=a.semiperfect ? a : b;
							}else if(a.score!=b.score){
								better=(a.score>b.score ? a : b);
							}else if(a.slowScore!=b.slowScore){
								better=(a.slowScore>b.slowScore ? a : b);
							}else if(a.pairedScore!=b.pairedScore){
								better=(a.pairedScore>b.pairedScore ? a : b);
							}else if(a.pairedScore!=b.pairedScore){
								better=(a.quickScore>b.quickScore ? a : b);
							}else{
								better=a;
							}
							
//							if((a.perfect && b.perfect) || (a.semiperfect && b.semiperfect)){
							if(a.semiperfect && b.semiperfect){
								if(a.start==b.start || a.stop==b.stop){
									list.set(i, better);
									list.set(j, null);
									removed++;
									a=better;
								}else{
									//retain both of them
								}
							}else if(a.perfect || b.perfect){
								list.set(i, better);
								list.set(j, null);
								removed++;
								a=better;
							}else if(a.semiperfect || b.semiperfect){
								if(a.start==b.start && a.stop==b.stop){
									list.set(i, better);
									list.set(j, null);
									removed++;
									a=better;
								}else{
									//retain both of them
								}
							}else if(subsumeInexact || (a.start==b.start && (subsumeIfOnlyStartMatches || a.stop==b.stop))){
								assert(!a.semiperfect && !a.perfect && !b.semiperfect && !b.perfect);
								a.setLimits(min(a.start, b.start), max(a.stop, b.stop));
								a.setSlowScore(max(a.slowScore, b.slowScore));
								a.setPairedScore(a.pairedScore<=0 && b.pairedScore<=0 ? 0 : max(a.slowScore+1, a.pairedScore, b.pairedScore));
								a.quickScore=max(a.quickScore, b.quickScore);
								a.setScore(max(a.score, b.score, a.pairedScore));
								a.gaps=better.gaps;//Warning!  Merging gaps would be better; this could cause out-of-bounds.
								//TODO: Test for a subsumption length limit.
								list.set(j, null);
								removed++;
							}
						}
					}
				}
			}
		}
		
//		if(removed>0){condense(list);}
		if(removed>0){condenseStrict(list);}
		return removed;
	}
	

	
	public static final int removeOverlappingSites(ArrayList<SiteScore> list, boolean requireAMatchingEnd){
		if(list==null || list.size()<2){return 0;}
		Shared.sort(list, SiteScore.PCOMP);
		
		int removed=0;
		
		
		for(int i=0; i<list.size(); i++){
			SiteScore a=list.get(i);
			boolean overlappingA=true;
			if(a!=null){
				for(int j=i+1; overlappingA && j<list.size(); j++){
					SiteScore b=list.get(j);
					if(b!=null){
						overlappingA=(a.chrom==b.chrom && b.start<a.stop && b.stop>a.start);
						if(overlappingA && a.strand==b.strand){
							
							SiteScore better=null;
							if(a.perfect!=b.perfect){
								better=a.perfect ? a : b;
							}else if(a.score!=b.score){
								better=(a.score>b.score ? a : b);
							}else if(a.slowScore!=b.slowScore){
								better=(a.slowScore>b.slowScore ? a : b);
							}else if(a.pairedScore!=b.pairedScore){
								better=(a.pairedScore>b.pairedScore ? a : b);
							}else if(a.pairedScore!=b.pairedScore){
								better=(a.quickScore>b.quickScore ? a : b);
							}else{
								better=a;
							}
							
							if(a.start==b.start && a.stop==b.stop){
								list.set(i, better);
								list.set(j, null);
								a=better;
								removed++;
							}else if(a.start==b.start || a.stop==b.stop){ //In this case they cannot both be perfect
								list.set(i, better);
								list.set(j, null);
								a=better;
								removed++;
							}else if(!requireAMatchingEnd && a.score!=b.score){
								list.set(i, better);
								list.set(j, null);
								a=better;
								removed++;
							}
						}
					}
				}
			}
		}
		
//		if(removed>0){condense(list);}
		if(removed>0){condenseStrict(list);}
		return removed;
	}
	

	
	/** Returns the number of sitescores in the list within "thresh" of the top score.  Assumes list is sorted descending. 
	 * This is used to determine whether a mapping is ambiguous. */
	public static final int countTopScores(ArrayList<SiteScore> list, int thresh){
		assert(thresh>=0) : thresh;
		if(list==null || list.isEmpty()){return 0;}
		int count=1;
		final SiteScore ss=list.get(0);
		final int limit=ss.score-thresh;
		
		for(int i=1; i<list.size(); i++){
			SiteScore ss2=list.get(i);
			if(ss2.score<limit){break;}
			if(ss.start!=ss2.start && ss.stop!=ss2.stop){ //Don't count mappings to the same location
				count++;
			}
		}
		return count;
	}
	

	
	/** Assumes list is sorted by NON-PAIRED score. 
	 * Returns number removed. */
	public static final int removeLowQualitySitesPaired(ArrayList<SiteScore> list, int maxSwScore, float multSingle, float multPaired){
		if(list==null || list.size()==0){return 0;}
		
		assert(multSingle>=multPaired);
		
		int initialSize=list.size();
		final int swScoreThresh=(int)(maxSwScore*multSingle); //Change low-quality alignments to no-hits.
		final int swScoreThreshPaired=(int)(maxSwScore*multPaired);
		if(list.get(0).score<swScoreThreshPaired){list.clear(); return initialSize;}
		
		for(int i=list.size()-1; i>=0; i--){
			SiteScore ss=list.get(i);
			assert(ss.score==ss.slowScore) : ss.quickScore+", "+ss.slowScore+", "+ss.pairedScore+", "+ss.score+"\n"+ss;
			assert(i==0 || ss.slowScore<=list.get(i-1).slowScore) : "List is not sorted by singleton score!";
			if(ss.pairedScore>0){
				assert(ss.pairedScore>ss.quickScore || ss.pairedScore>ss.slowScore) : ss;
				if(ss.slowScore<swScoreThreshPaired){list.remove(i);}
			}else{
				assert(ss.pairedScore<=0) : ss.toText();
				if(ss.slowScore<swScoreThresh){list.remove(i);}
			}
		}
		
		return initialSize-list.size();
	}
	

	
//	/** Assumes list is sorted by NON-PAIRED score. 
//	 * Returns number removed. */
//	public static final int removeLowQualitySitesUnpaired(ArrayList<SiteScore> list, int maxSwScore, float multSingle){
//		if(list==null || list.size()==0){return 0;}
//		
//		int initialSize=list.size();
//		final int swScoreThresh=(int)(maxSwScore*multSingle); //Change low-quality alignments to no-hits.
//		if(list.get(0).score<swScoreThresh){list.clear(); return initialSize;}
//		
////		for(int i=list.size()-1; i>=0; i--){
//		for(int i=list.size()-1; i>1; i--){
//			SiteScore ss=list.get(i);
//			assert(ss.score==ss.slowScore);
//			assert(i==0 || ss.slowScore<=list.get(i-1).slowScore) : "List is not sorted by singleton score!";
//			assert(ss.pairedScore==0) : ss.toText();
//			if(ss.slowScore<swScoreThresh){list.remove(i);}
//		}
//		
//		return initialSize-list.size();
//	}

	
	/** Assumes list is sorted by NON-PAIRED score. 
	 * Returns number removed. */
	public static final int removeLowQualitySitesUnpaired(ArrayList<SiteScore> list, int thresh){
		if(list==null || list.size()==0){return 0;}
		
		int initialSize=list.size();
		if(list.get(0).score<thresh){list.clear(); return initialSize;}
		
//		for(int i=list.size()-1; i>=0; i--){
		for(int i=list.size()-1; i>1; i--){
			SiteScore ss=list.get(i);
			assert(ss.score==ss.slowScore || (ss.score<=0 && ss.slowScore<=0)) : ss;
			assert(i==0 || ss.slowScore<=list.get(i-1).slowScore) : "List is not sorted by singleton score!";
			assert(ss.pairedScore<=0) : ss.toText();
			if(ss.slowScore<thresh){list.remove(i);}
		}
		
		return initialSize-list.size();
	}
	

	
	/** Assumes list is sorted by NON-PAIRED score. 
	 * Returns number removed. */
	public static final int removeLowQualitySitesPaired2(ArrayList<SiteScore> list, int maxSwScore, float multSingle, float multPaired, int expectedSites){
		if(list==null || list.size()==0){return 0;}
		
		assert(multSingle>=multPaired);
		
		int initialSize=list.size();
		final int swScoreThresh=(int)(maxSwScore*multSingle); //Change low-quality alignments to no-hits.
		final int swScoreThreshPaired=(int)(maxSwScore*multPaired);
		final int swScoreThresh2=(int)(maxSwScore*multSingle*1.2f);
		final int swScoreThreshPaired2=(int)(maxSwScore*multPaired*1.1f);
		if(list.get(0).score<swScoreThreshPaired){list.clear(); return initialSize;}
		final int nthBest=list.get(Tools.min(list.size(), expectedSites)-1).score-maxSwScore/64;
		
		for(int i=list.size()-1, min=expectedSites*2; i>min; i--){
			if(list.get(i).slowScore>=nthBest){break;}
			list.remove(i);
		}
		
		for(int i=list.size()-1; i>=0; i--){
			SiteScore ss=list.get(i);
			assert(ss.score==ss.slowScore);
			assert(i==0 || ss.slowScore<=list.get(i-1).slowScore) : "List is not sorted by singleton score!";
			if(ss.pairedScore>0){
				int thresh=(i>=expectedSites ? swScoreThreshPaired2 : swScoreThreshPaired);
				assert(ss.pairedScore>ss.quickScore || ss.pairedScore>ss.slowScore) : ss;
				if(ss.slowScore<thresh){list.remove(i);}
			}else{
				int thresh=(i>=expectedSites ? swScoreThresh2 : swScoreThresh);
//				assert(ss.pairedScore==0) : ss.toText(); //Disable in case of negative values
				if(ss.slowScore<thresh){list.remove(i);}
			}
		}
		
		return initialSize-list.size();
	}
	

	
	/** Assumes list is sorted by NON-PAIRED score. 
	 * Returns number removed. 
	 * This has a couple of changes (like potentially removing the second-best site) that make it applicable to SKIMMER not MAPPER.
	 * */
	public static final int removeLowQualitySitesUnpaired2(ArrayList<SiteScore> list, int maxSwScore, float multSingle, int expectedSites){
		if(list==null || list.size()==0){return 0;}
		
		for(int i=expectedSites/2; i<list.size(); i++){
			if(list.get(i).perfect){expectedSites++;}
		}
		
		int initialSize=list.size();
		final int swScoreThresh=(int)(maxSwScore*multSingle); //Change low-quality alignments to no-hits.
		final int swScoreThresh2=(int)(maxSwScore*multSingle*1.2f); //Change low-quality alignments to no-hits.
		if(list.get(0).score<swScoreThresh){list.clear(); return initialSize;}
		final int nthBest=list.get(Tools.min(list.size(), expectedSites)-1).score-maxSwScore/64;
		
		for(int i=list.size()-1, min=expectedSites*2; i>min; i--){
			if(list.get(i).slowScore>=nthBest){break;}
			list.remove(i);
		}
		
//		for(int i=list.size()-1; i>=0; i--){
		for(int i=list.size()-1; i>=1; i--){
			SiteScore ss=list.get(i);
			assert(ss.score==ss.slowScore);
			assert(i==0 || ss.slowScore<=list.get(i-1).slowScore) : "List is not sorted by singleton score!";
			assert(ss.pairedScore<=0) : ss.toText(); //This was "==0" but that makes the assertion fire for negative values. 
			int thresh=(i>=expectedSites ? swScoreThresh2 : swScoreThresh);
			if(ss.slowScore<thresh){list.remove(i);}
		}
		
		return initialSize-list.size();
	}
	
	
//	public static final void trimSitesBelowCutoff(ArrayList<SiteScore> ssl, int cutoff, boolean retainPaired){
//		trimSitesBelowCutoff(ssl, cutoff, retainPaired, 1);
//	}
	
	
//	public static final void trimSitesBelowCutoff(ArrayList<SiteScore> ssl, int cutoff, boolean retainPaired, int minSitesToRetain){
////		assert(false);
//		assert(minSitesToRetain>=1);
//		if(ssl==null || ssl.size()<minSitesToRetain){return;}
//		
//		ArrayList<SiteScore> ssl2=new ArrayList<SiteScore>(ssl.size());
//		for(SiteScore ss : ssl){
//			if(ss.score>=cutoff || (retainPaired && ss.pairedScore>0)){
//				ssl2.add(ss);
//			}
//		}
//		
////		Shared.sort(ssl2);
////		System.err.println("Cutoff: "+cutoff);
////		for(SiteScore ss : ssl2){
////			System.err.print("("+ss.chrom+", "+ss.score+"), ");
////		}
////		System.err.println();
//		
//		if(ssl2.size()==ssl.size()){return;}
////		System.err.println("cutoff: "+cutoff+",\tsize: "+ssl.size()+" -> "+ssl2.size());
//		ssl.clear();
//		ssl.addAll(ssl2);
//	}
	
	
	public static final void trimSitesBelowCutoff(ArrayList<SiteScore> ssl, int cutoff, boolean retainPaired, boolean retainSemiperfect,
			int minSitesToRetain, int maxSitesToRetain){
//		assert(false);
		assert(minSitesToRetain>=1);
		assert(maxSitesToRetain>minSitesToRetain) : maxSitesToRetain+", "+minSitesToRetain+"\nError - maxsites2 must be greater than "+minSitesToRetain+"!";
		if(ssl==null || ssl.size()<=minSitesToRetain){return;}
		while(ssl.size()>maxSitesToRetain){ssl.remove(ssl.size()-1);}
		
		int removed=0;
		final int maxToRemove=ssl.size()-minSitesToRetain;
		
		assert(minSitesToRetain==1 || inOrder(ssl));
		
		if(retainPaired){
			for(int i=ssl.size()-1; i>=0; i--){
				SiteScore ss=ssl.get(i);
				if(!retainSemiperfect || !ss.semiperfect){
					if(ss.score<cutoff && ss.pairedScore<=0){
						ssl.set(i, null);
						removed++;
						if(removed>=maxToRemove){
							assert(removed==maxToRemove);
							break;
						}
					}
				}
			}
		}else{
			for(int i=ssl.size()-1; i>=0; i--){
				SiteScore ss=ssl.get(i);
				if(!retainSemiperfect || !ss.semiperfect){
					if(ss.score<cutoff){
						ssl.set(i, null);
						removed++;
						if(removed>=maxToRemove){
							assert(removed==maxToRemove);
							break;
						}
					}
				}
			}
		}
		
		if(removed>0){
			condenseStrict(ssl);
		}
		assert(ssl.size()>=minSitesToRetain);
	}
	
	//Messes up order
//	public static final void trimSitesBelowCutoffInplace(ArrayList<SiteScore> ssl, int cutoff, boolean retainPaired){
////		assert(false);
//		if(ssl==null || ssl.size()<2){return;}
//		
//		for(int i=0; i<ssl.size(); i++){
//			SiteScore ss=ssl.get(i);
//			if(ss.score<cutoff && (!retainPaired || ss.pairedScore==0)){
//				SiteScore temp=ssl.remove(ssl.size()-1);
//				if(i<ssl.size()){
//					ssl.set(i, temp);
//					i--;
//				}
//			}
//		}
//	}
	
	public static CharSequence toStringSafe(byte[] array){
		if(array==null){return "null";}
		StringBuilder sb=new StringBuilder();
		sb.append(Arrays.toString(array));
		if(array.length<1 || array[0]<32 || array[0]>126){return sb;}
		sb.append('\n');
		for(int i=0; i<array.length; i++){
			byte b=array[i];
			if(b<32 || b>126){break;}
			sb.append((char)b);
		}
		return sb;
	}
	
	public static boolean equals(int[] a, int[] b){
		if(a==b){return true;}
		if(a==null || b==null){return false;}
		if(a.length!=b.length){return false;}
		for(int i=0; i<a.length; i++){
			if(a[i]!=b[i]){return false;}
		}
		return true;
	}
	
	public static boolean equals(byte[] a, byte[] b){
		if(a==b){return true;}
		if(a==null || b==null){return false;}
		if(a.length!=b.length){return false;}
		for(int i=0; i<a.length; i++){
			if(a[i]!=b[i]){return false;}
		}
		return true;
	}
	
	public static boolean equals(String a, byte[] b){
		if(a==null || b==null){
			return (a==null && b==null);
		}
		if(a.length()!=b.length){return false;}
		for(int i=0; i<b.length; i++){
			if(a.charAt(i)!=b[i]){return false;}
		}
		return true;
	}
	
	public static boolean contains(byte[] a, byte[] b, int start){
		if(a==null || b==null){
			return (a==null && b==null);
		}
		if(a.length<b.length+start){return false;}
		for(int i=start, j=0; j<b.length; i++, j++){
			if(a[i]!=b[j]){return false;}
		}
		return true;
	}
	
	/**
	 * @param array
	 * @param s
	 * @return True if the array starts with the String.
	 */
	public static boolean startsWith(byte[] array, String s) {
		return startsWith(array, s, 0);
	}
	
	/**
	 * @param array
	 * @param s
	 * @return True if the array starts with the String.
	 */
	public static boolean startsWith(byte[] array, String s, int initialPos) {
		if(array==null || s==null || array.length+initialPos<s.length()){return false;}
		for(int i=initialPos; i<s.length(); i++){
			if(array[i]!=s.charAt(i)){return false;}
		}
		return true;
	}

	public static int compare(byte[] a, byte[] b){
		if(a==b){return 0;}
		if(a==null){return -1;}
		if(b==null){return 1;}
		int lim=min(a.length, b.length);
		for(int i=0; i<lim; i++){
			if(a[i]!=b[i]){return a[i]-b[i];}
		}
		return a.length-b.length;
	}

	public static int sumInt(byte[] array){
		long x=0;
		for(byte y : array){x+=y;}
		assert(x<=Integer.MAX_VALUE && x>=Integer.MIN_VALUE) : x;
		return (int)x;
	}

	public static void add(long[] array, long[] incr) {
		for(int i=0; i<array.length; i++){
			array[i]+=incr[i];
		}
	}

	public static long sum(byte[] array){
		if(array==null){return 0;}
		long x=0;
		for(byte y : array){x+=y;}
		return x;
	}

	public static long sum(char[] array){
		long x=0;
		for(char y : array){x+=y;}
		return x;
	}
	
	public static long sum(short[] array){
		long x=0;
		for(short y : array){x+=y;}
		return x;
	}
	
	public static int cardinality(short[] array){
		int x=0;
		for(int y : array){if(y!=0){x++;}}
		return x;
	}
	
	public static long sum(int[] array){
		long x=0;
		for(int y : array){x+=y;}
		return x;
	}

	public static double sum(double[] array){
		double x=0;
		for(double y : array){x+=y;}
		return x;
	}
	
	public static double mean(int[] array){
		return sum(array)/(double)array.length;
	}
	
	public static double harmonicMean(int[] array){
		double sum=0;
		for(int x : array){
			if(x>0){sum+=1.0/x;}
		}
		return array.length/sum;
	}
	
	public static int cardinality(int[] array){
		int x=0;
		for(int y : array){if(y!=0){x++;}}
		return x;
	}
	
	public static double weightedAverage(long[] array){
		if(array.length<2){
			return array.length==1 ? array[0] : 0;
		}
		double wsum=0;
		long div=0;
		final int mid=array.length/2;
		for(int i=0; i<mid; i++){
			wsum+=(i+1)*(array[i]+array[array.length-1-1]);
			div+=(i+1)*2;
		}
		if((array.length&1)==1){
			wsum+=(mid+1)*array[mid];
			div+=(mid+1);
		}
		return wsum/div;
	}
	
	public static long sum(int[] array, int from, int to){
		long x=0;
		for(int i=from; i<=to; i++){x+=array[i];}
		return x;
	}
	
	public static long sum(long[] array){
		long x=0;
		for(long y : array){x+=y;}
		return x;
	}
	
	public static long sum(long[] array, int from, int to){
		long x=0;
		for(int i=from; i<=to; i++){x+=array[i];}
		return x;
	}
	
	public static long sumHistogram(long[] array){
		long x=0;
		for(int i=1; i<array.length; i++){
			x+=(i*array[i]);
		}
		return x;
	}
	
	public static long minHistogram(long[] array){
		for(int i=0; i<array.length; i++){
			if(array[i]>0){return i;}
		}
		return 0;
	}
	
	public static long maxHistogram(long[] array){
		for(int i=array.length-1; i>=0; i--){
			if(array[i]>0){return i;}
		}
		return 0;
	}
	
	public static long sum(AtomicIntegerArray array){
		long x=0;
		for(int i=0; i<array.length(); i++){x+=array.get(i);}
		return x;
	}
	
	public static long sum(AtomicLongArray array){
		long x=0;
		for(int i=0; i<array.length(); i++){x+=array.get(i);}
		return x;
	}
	
	public static long[] toArray(AtomicLongArray array){
		long[] x=new long[array.length()];
		for(int i=0; i<array.length(); i++){x[i]=array.get(i);}
		return x;
	}
	
	public static long[] toArray(CoverageArray array){
		long[] x=new long[array.maxIndex+1];
		for(int i=0; i<=array.maxIndex; i++){x[i]=array.get(i);}
		return x;
	}
	
	public static int min(int[] array){
		int min=Integer.MAX_VALUE;
		for(int y : array){if(y<min){min=y;}}
		return min;
	}
	
	public static byte min(byte[] array){
		byte min=Byte.MAX_VALUE;
		for(byte y : array){if(y<min){min=y;}}
		return min;
	}
	
	public static int intSum(int[] array){
		int x=0;
		for(int y : array){x+=y;}
		return x;
	}
	
	public static void reverseInPlace(final byte[] array){
		if(array==null){return;}
		final int max=array.length/2, last=array.length-1;
		for(int i=0; i<max; i++){
			byte temp=array[i];
			array[i]=array[last-i];
			array[last-i]=temp;
		}
	}
	
	public static void reverseInPlace(final char[] array){
		if(array==null){return;}
		final int max=array.length/2, last=array.length-1;
		for(int i=0; i<max; i++){
			char temp=array[i];
			array[i]=array[last-i];
			array[last-i]=temp;
		}
	}
	
	public static void reverseInPlace(final int[] array){
		if(array==null){return;}
		reverseInPlace(array, 0, array.length);
	}
	
	public static void reverseInPlace(final long[] array){
		if(array==null){return;}
		reverseInPlace(array, 0, array.length);
	}
	
	public static void reverseInPlace(final byte[] array, final int from, final int to){
		if(array==null){return;}
		final int max=to/2, last=to-1;
		for(int i=0; i<max; i++){
			byte temp=array[i];
			array[i]=array[last-i];
			array[last-i]=temp;
		}
	}
	
	public static void reverseInPlace(final int[] array, final int from, final int to){
		if(array==null){return;}
		final int max=to/2, last=to-1;
		for(int i=0; i<max; i++){
			int temp=array[i];
			array[i]=array[last-i];
			array[last-i]=temp;
		}
	}
	
	public static void reverseInPlace(final long[] array, final int from, final int to){
		if(array==null){return;}
		final int max=to/2, last=to-1;
		for(int i=0; i<max; i++){
			long temp=array[i];
			array[i]=array[last-i];
			array[last-i]=temp;
		}
	}
	
	public static byte[] reverseAndCopy(final byte[] array){
//		if(array==null){return null;}
//		byte[] copy=Arrays.copyOf(array, array.length);
//		reverseInPlace(copy);
//		return copy;
		return reverseAndCopy(array, null);
	}
	
	public static int[] reverseAndCopy(final int[] array){
//		if(array==null){return null;}
//		int[] copy=Arrays.copyOf(array, array.length);
//		reverseInPlace(copy);
//		return copy;
		return reverseAndCopy(array, null);
	}
	
	public static byte[] reverseAndCopy(final byte[] array, byte[] out){
		if(array==null){
			assert(out==null);
			return null;
		}
		if(out==null){out=new byte[array.length];}
		assert(array.length==out.length && array!=out);
		for(int i=0, last=array.length-1; i<array.length; i++){out[i]=array[last-i];}
		return out;
	}
	
	public static int[] reverseAndCopy(final int[] array, int[] out){
		if(array==null){
			assert(out==null);
			return null;
		}
		if(out==null){out=new int[array.length];}
		assert(array.length==out.length && array!=out);
		for(int i=0, last=array.length-1; i<array.length; i++){out[i]=array[last-i];}
		return out;
	}
	
	public static void cullHighFreqEntries(int[][] data, float fractionToExclude){
		if(fractionToExclude<=0){return;}
		int[] count=new int[data.length];
		
		long numBases=0;
		
		for(int i=0; i<data.length; i++){
			count[i]=(data[i]==null ? 0 : data[i].length);
			numBases+=count[i];
		}
		
		int numIndicesToRemove=((int)(numBases*fractionToExclude));
		
		Arrays.sort(count);
		
		for(int i=1; i<count.length; i++){
			assert(count[i]>=count[i-1]) : "\n\ncount["+i+"]="+count[i]+"\ncount["+(i-1)+"]="+count[i-1]+"\n";
		}
		
		int pos=count.length-1;
		for(int sum=0; pos>1 && sum<numIndicesToRemove; pos--){
			sum+=count[pos];
		}
		int maxLengthToKeep2=count[pos];
		
		for(int i=0; i<data.length; i++){
			if(data[i]!=null && data[i].length>maxLengthToKeep2){data[i]=null;}
		}
	}
	
	public static int findLimitForHighFreqEntries(int[][] data, float fractionToExclude){
		if(fractionToExclude<=0){return Integer.MAX_VALUE;}
		int[] count=new int[data.length];
		
		long numBases=0;
		
		for(int i=0; i<data.length; i++){
			count[i]=(data[i]==null ? 0 : data[i].length);
			numBases+=count[i];
		}
		
		int numIndicesToRemove=((int)(numBases*fractionToExclude));
		
		Arrays.sort(count);
		
		for(int i=1; i<count.length; i++){
			assert(count[i]>=count[i-1]) : "\n\ncount["+i+"]="+count[i]+"\ncount["+(i-1)+"]="+count[i-1]+"\n";
		}
		
		int pos=count.length-1;
		for(int sum=0; pos>1 && sum<numIndicesToRemove; pos--){
			sum+=count[pos];
		}
		int maxLengthToKeep2=count[pos];
		
		return maxLengthToKeep2;
	}
	
	public static void cullClumpyEntries(final int[][] data, final int maxDist, final int minLength, final float fraction){
		
		long total=0;
		long removedSites=0;
		long removedKeys=0;
		
		if(maxDist<=0){return;}
		for(int i=0; i<data.length; i++){
			int[] array=data[i];
			total+=(array==null ? 0 : array.length);
			if(array!=null && array.length>=minLength){
				if(isClumpy(array, maxDist, fraction)){
					removedSites+=array.length;
					removedKeys++;
					data[i]=null;
				}
			}
		}

//		System.err.println("Removed\t"+removedSites+"\t/ "+total+"\tsites," +
//				" or "+String.format("%.4f", (removedSites*100f/total))+"%");
//		System.err.println("Removed\t"+removedKeys+"\t/ "+data.length+"\tkeys," +
//				" or  "+String.format("%.4f", (removedKeys*100f/data.length))+"%");
		
	}
	
	public static HashSet<Integer> banClumpyEntries(final int[][] data, final int maxDist, final int minLength, final float fraction){
		
		HashSet<Integer> set=new HashSet<Integer>(128);
		
		long total=0;
		long removedSites=0;
		long removedKeys=0;
		
		if(maxDist<=0){return set;}
		
		for(int i=0; i<data.length; i++){
			int[] array=data[i];
			total+=(array==null ? 0 : array.length);
			if(array!=null && array.length>=minLength){
				if(isClumpy(array, maxDist, fraction)){
					removedSites+=array.length;
					removedKeys++;
					set.add(i);
				}
			}
		}

//		System.err.println("Banned\t"+removedSites+"\t/ "+total+"\tsites," +
//				" or "+String.format("%.4f", (removedSites*100f/total))+"%");
//		System.err.println("Banned\t"+removedKeys+"\t/ "+data.length+"\tkeys," +
//				" or  "+String.format("%.4f", (removedKeys*100f/data.length))+"%");
		
		return set;
		
	}
	
	public static final boolean isClumpy(final int[] array, final int maxDist, final float fraction){
		if(array==null){return false;}
		int count=0;
		for(int i=1; i<array.length; i++){
			int dif=array[i]-array[i-1];
			if(dif<=maxDist){count++;}
		}
		return count>=(array.length*fraction);
	}

	public static int[] makeLengthHistogram(int[][] x, int buckets) {
		int[] lengths=new int[x.length];
		long total=0;
		for(int i=0; i<x.length; i++){
			int[] list=x[i];
			if(list!=null){
				lengths[i]=list.length;
				total+=list.length;
			}
		}
		Arrays.sort(lengths);
		
		int[] hist=new int[buckets+1];
		
		long sum=0;
		int ptr=0;
		for(int i=0; i<buckets; i++){
			long nextLimit=((total*i)+buckets/2)/buckets;
			while(ptr<lengths.length && sum<nextLimit){
				sum+=lengths[ptr];
				ptr++;
			}
			
			hist[i]=lengths[Tools.max(0, ptr-1)];
		}
		hist[hist.length-1]=lengths[lengths.length-1];
		
//		System.out.println(Arrays.toString(hist));
//		assert(false);
		return hist;
	}
	
	public static String toKMG(long x){
		double div=1;
		String ext="";
		if(x>10000000000000L){
			div=1000000000000L;
			ext="T";
		}else if(x>10000000000L){
			div=1000000000L;
			ext="B";
		}else if(x>10000000){
			div=1000000;
			ext="M";
		}else if(x>100000){
			div=1000;
			ext="K";
		}
		return String.format("%.2f", x/div)+ext;
	}
	
	public static byte[] parseRemap(String b){
		final byte[] remap;
		if(b==null || ("f".equalsIgnoreCase(b) || "false".equalsIgnoreCase(b))){
			remap=null;
		}else{
			assert((b.length()&1)==0) : "Length of remap argument must be even.  No whitespace is allowed.";
			
			remap=new byte[128];
			for(int j=0; j<remap.length; j++){remap[j]=(byte)j;}
			for(int j=0; j<b.length(); j+=2){
				char x=b.charAt(j), y=b.charAt(j+1);
				remap[x]=(byte)y;
			}
		}
		return remap;
	}
	
	public static long parseKMG(String b){
		if(b==null){return 0;}
		char c=Character.toLowerCase(b.charAt(b.length()-1));
		boolean dot=b.indexOf('.')>=0;
		if(!Character.isLetter(c) && !dot){return Long.parseLong(b);}
		
		long mult=1;
		if(Character.isLetter(c)){
			if(c=='k'){mult=1000;}
			else if(c=='m'){mult=1000000;}
			else if(c=='g' || c=='b'){mult=1000000000;}
			else if(c=='t'){mult=1000000000000L;}
			else{throw new RuntimeException(b);}
			b=b.substring(0, b.length()-1);
		}
		
		if(!dot){return Long.parseLong(b)*mult;}
		
		return (long)(Double.parseDouble(b)*mult);
	}
	
	public static boolean isNumber(String s){
		if(s==null || s.length()==0){return false;}
		char c=s.charAt(0);
		return Character.isDigit(c) || c=='.' || c=='-';
	}
	
	public static boolean parseBoolean(String s){
		if(s==null || s.length()<1){return true;}
		if(s.length()==1){
			char c=Character.toLowerCase(s.charAt(0));
			return c=='t' || c=='1';
		}
		if(s.equalsIgnoreCase("null") || s.equalsIgnoreCase("none")){return false;}
		return Boolean.parseBoolean(s);
	}
	
	public static boolean parseYesNo(String s){
		if(s==null || s.length()<1){return true;}
		if(s.length()==1){
			char c=Character.toLowerCase(s.charAt(0));
			if(c=='y'){return true;}
			if(c=='n'){return false;}
			throw new RuntimeException(s);
		}

		if(s.equalsIgnoreCase("yes")){return true;}
		if(s.equalsIgnoreCase("no")){return false;}
		if(s.equalsIgnoreCase("unknown")){return false;} //Special case for IMG database
		
		throw new RuntimeException(s);
	}
	
	public static int[] parseIntArray(String s, String regex){
		if(s==null){return null;}
		String[] split=s.split(regex);
		int[] array=new int[split.length];
		for(int i=0; i<split.length; i++){
			array[i]=Integer.parseInt(split[i]);
		}
		return array;
	}
	
	public static byte[] parseByteArray(String s, String regex){
		if(s==null){return null;}
		String[] split=s.split(regex);
		byte[] array=new byte[split.length];
		for(int i=0; i<split.length; i++){
			array[i]=Byte.parseByte(split[i]);
		}
		return array;
	}
	
	public static int parseIntHexDecOctBin(final String s){
		if(s==null || s.length()<1){return 0;}
		int radix=10;
		if(s.length()>1 && s.charAt(1)=='0'){
			final char c=s.charAt(1);
			if(c=='x' || c=='X'){radix=16;}
			else if(c=='b' || c=='B'){radix=2;}
			else if(c=='o' || c=='O'){radix=8;}
		}
		return Integer.parseInt(s, radix);
	}
	
	public static float parseFloat(byte[] array, int a, int b){
		String s=new String(array, a, b-a);
		return Float.parseFloat(s);
	}
	
	public static int parseInt(byte[] array, int a, int b){
		assert(b>a);
		int r=0;
		final byte z='0';
		boolean negative=false;
		if(array[a]=='-'){negative=true; a++;}
		for(; a<b; a++){
			int x=(array[a]-z);
			assert(x<10 && x>=0) : x+" = "+(char)array[a]+"\narray="+new String(array)+", start="+a+", stop="+b;
			r=(r*10)+x;
		}
		if(negative){r*=-1;}
		return r;
	}
	
	public static long parseLong(byte[] array){return parseLong(array, 0, array.length);}
	
	public static long parseLong(byte[] array, int a, int b){
		assert(b>a);
		long r=0;
		final byte z='0';
		boolean negative=false;
		if(array[a]=='-'){negative=true; a++;}
		for(; a<b; a++){
			int x=(array[a]-z);
			assert(x<10 && x>=0) : x+" = "+(char)array[a]+"\narray="+new String(array)+", start="+a+", stop="+b;
			r=(r*10)+x;
		}
		if(negative){r*=-1;}
		return r;
	}
	
	/** TODO:  This (temporarily) uses a lot of memory.  Could be reduced by making an array of length max(x) and counting occurrences. */
	public static int[] makeLengthHistogram2(int[] x, int buckets, boolean verbose) {
		int[] lengths=KillSwitch.copyOf(x, x.length);
		long total=sum(x);
		Shared.sort(lengths);
		
		if(verbose){
			System.out.println("Length array size:\t"+x.length);
			System.out.println("Min value:        \t"+lengths[0]);
			System.out.println("Med value:        \t"+lengths[lengths.length/2]);
			System.out.println("Max value:        \t"+lengths[lengths.length-1]);
			System.out.println("Total:            \t"+total);
		}
		
		int[] hist=new int[buckets+1];
		
		long sum=0;
		int ptr=0;
		for(int i=0; i<buckets; i++){
			long nextLimit=((total*i)+buckets/2)/buckets;
			while(ptr<lengths.length && sum<nextLimit){
				sum+=lengths[ptr];
				ptr++;
			}
			
			hist[i]=lengths[Tools.max(0, ptr-1)];
		}
		hist[hist.length-1]=lengths[lengths.length-1];
		
//		System.out.println(Arrays.toString(hist));
//		assert(false);
		return hist;
	}
	
	public static int[] makeLengthHistogram3(int[] x, int buckets, boolean verbose) {
		int max=max(x);
		if(max>x.length){
			Data.sysout.println("Reverted to old histogram mode.");
			return makeLengthHistogram2(x, buckets, verbose);
		}
		
		int[] counts=new int[max+1];
		long total=0;
		for(int i=0; i<x.length; i++){
			int a=x[i];
			if(a>=0){
				counts[a]++;
				total+=a;
			}
		}
		
		return makeLengthHistogram4(counts, buckets, total, verbose);
	}
	
	/** Uses counts of occurrences of lengths rather than raw lengths */
	public static int[] makeLengthHistogram4(int[] counts, int buckets, long total, boolean verbose) {
		if(total<=0){
			total=0;
			for(int i=1; i<counts.length; i++){
				total+=(i*counts[i]);
			}
		}
		
		if(verbose){
//			System.out.println("Length array size:\t"+x.length);
//			System.out.println("Min value:        \t"+lengths[0]);
//			System.out.println("Med value:        \t"+lengths[lengths.length/2]);
//			System.out.println("Max value:        \t"+lengths[lengths.length-1]);
			System.err.println("Total:            \t"+total);
		}
		
		int[] hist=new int[buckets+1];
		
		long sum=0;
		int ptr=0;
		for(int i=0; i<buckets; i++){
			long nextLimit=((total*i)+buckets/2)/buckets;
			while(ptr<counts.length && sum<nextLimit){
				sum+=counts[ptr]*ptr;
				ptr++;
			}
			
			hist[i]=Tools.max(0, ptr-1);
		}
		hist[hist.length-1]=counts.length-1;
		
//		System.out.println(Arrays.toString(hist));
//		assert(false);
		return hist;
	}
	
	/**
	 * @param cov
	 * @return
	 */
	public static int average(short[] array) {
		return (int)(array==null || array.length==0 ? 0 : sum(array)/array.length);
	}
	
	/**
	 * @param cov
	 * @return
	 */
	public static int average(int[] array) {
		return (int)(array==null || array.length==0 ? 0 : sum(array)/array.length);
	}

	public static double averageDouble(int[] array) {
		return (int)(array==null || array.length==0 ? 0 : sum(array)*1.0/array.length);
	}
	
	public static int median(int[] array){return percentile(array, .5);}
	
	public static int median(long[] array){return percentile(array, .5);}
	
	public static int percentile(int[] array, double fraction){
		if(array==null || array.length<1){return 0;}
		long target=(long)(sum(array)*fraction);
		long sum=0;
		for(int i=0; i<array.length; i++){
			sum+=array[i];
			if(sum>=target){
				return i;
			}
		}
		return array.length-1;
	}
	
	public static int percentile(long[] array, double fraction){
		if(array==null || array.length<1){return 0;}
		long target=(long)(sum(array)*fraction);
		long sum=0;
		for(int i=0; i<array.length; i++){
			sum+=array[i];
			if(sum>=target){
				return i;
			}
		}
		return array.length-1;
	}
	
	public static int calcMode(long array[]){
		if(array==null || array.length<1){return 0;}
		int median=percentile(array, 0.5);
		int mode=0;
		long modeCount=array[mode];
		for(int i=1; i<array.length; i++){
			long count=array[i];
			if(count>modeCount || (count==modeCount && absdif(i, median)<absdif(mode, median))){
				mode=i;
				modeCount=count;
			}
		}
		return mode;
	}

	public static int absdif(int a, int b) {
		return a>b ? a-b : b-a;
	}

	public static float absdif(float a, float b) {
		return a>b ? a-b : b-a;
	}

	public static double absdif(double a, double b) {
		return a>b ? a-b : b-a;
	}
	
	/** Uses unsigned math */
	public static final int absdifUnsigned(int a, int b){
		return (a<0 == b<0) ? a>b ? a-b : b-a : Integer.MAX_VALUE;
	}
	
	public static final boolean overlap(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1 && b2>=a1; 
	}
	
	public static final int overlapLength(int a1, int b1, int a2, int b2){
		if(!overlap(a1,b1,a2,b2)){return 0;}
		if(a1<=a2){
			return b1>=b2 ? b2-a2+1 : b1-a2+1;
		}else{
			return b2>=b1 ? b1-a1+1 : b2-a1+1;
		}
	}
	
	/** Is (a1, b1) within (a2, b2) ? */
	public static final boolean isWithin(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a1>=a2 && b1<=b2; 
	}
	
	public static final int constrict(int point, int a, int b){
		assert(a<=b);
		return(point<a ? a : point>b ? b : point);
	}
	
	public static final int indexOf(byte[] array, byte b){
		int i=0;
		while(i<array.length && array[i]!=b){i++;}
		return (i==array.length ? -1 : i);
	}
	
	public static final int indexOf(char[] array, char b){
		int i=0;
		while(i<array.length && array[i]!=b){i++;}
		return (i==array.length ? -1 : i);
	}
	
	public static final int lastIndexOf(byte[] array, byte b){
		int i=array.length-1;
		while(i>=0 && array[i]!=b){i--;}
		return i;
	}
	
	public static final int stringLength(long x){
		if(x<0){
			if(x==Integer.MIN_VALUE){return 11;}
			return lengthOf(-x)+1;
		}
		return lengthOf(x);
	}
	
	public static final int stringLength(int x){
		if(x<0){
			if(x==Long.MIN_VALUE){return 20;}
			return lengthOf(-x)+1;
		}
		return lengthOf(x);
	}
	
	public static final int lengthOf(int x){
		assert(x>=0);
		int i=1;
		while(x>ilens[i]){i++;}
		return i;
	}
	
	public static final int lengthOf(long x){
		assert(x>=0);
		int i=1;
		while(x>llens[i]){i++;}
		return i;
	}

	public static final int max(int[] array){return array[maxIndex(array)];}
	
	public static final int maxIndex(int[] array){
		int max=array[0], maxIndex=0;
		for(int i=1; i<array.length; i++){
			if(array[i]>max){max=array[i];maxIndex=i;}
		}
		return maxIndex;
	}

	public static final long max(long[] array){return array[maxIndex(array)];}
	
	public static final int maxIndex(long[] array){
		long max=array[0];
		int maxIndex=0;
		for(int i=1; i<array.length; i++){
			if(array[i]>max){max=array[i];maxIndex=i;}
		}
		return maxIndex;
	}
	
	public static final double max(double[] array){return array[maxIndex(array)];}
	
	public static final int maxIndex(double[] array){
		double max=array[0];
		int maxIndex=0;
		for(int i=1; i<array.length; i++){
			if(array[i]>max){max=array[i];maxIndex=i;}
		}
		return maxIndex;
	}
	
	public static final double standardDeviation(long[] numbers){
		if(numbers==null || numbers.length<2){return 0;}
		long sum=sum(numbers);
		double avg=sum/(double)numbers.length;
		double sumdev2=0;
		for(int i=0; i<numbers.length; i++){
			long x=numbers[i];
			double dev=avg-x;
			sumdev2+=(dev*dev);
		}
		return Math.sqrt(sumdev2/numbers.length);
	}
	
	public static final double standardDeviation(int[] numbers){
		if(numbers==null || numbers.length<2){return 0;}
		long sum=sum(numbers);
		double avg=sum/(double)numbers.length;
		double sumdev2=0;
		for(int i=0; i<numbers.length; i++){
			long x=numbers[i];
			double dev=avg-x;
			sumdev2+=(dev*dev);
		}
		return Math.sqrt(sumdev2/numbers.length);
	}
	
	public static final double standardDeviation(char[] numbers){
		if(numbers==null || numbers.length<2){return 0;}
		long sum=sum(numbers);
		double avg=sum/(double)numbers.length;
		double sumdev2=0;
		for(int i=0; i<numbers.length; i++){
			long x=numbers[i];
			double dev=avg-x;
			sumdev2+=(dev*dev);
		}
		return Math.sqrt(sumdev2/numbers.length);
	}
	
	public static final double standardDeviation(short[] numbers){
		if(numbers==null || numbers.length<2){return 0;}
		long sum=sum(numbers);
		double avg=sum/(double)numbers.length;
		double sumdev2=0;
		for(int i=0; i<numbers.length; i++){
			long x=numbers[i];
			double dev=avg-x;
			sumdev2+=(dev*dev);
		}
		return Math.sqrt(sumdev2/numbers.length);
	}
	
	public static final double averageHistogram(long[] histogram){
		long sum=max(1, sum(histogram));
		long sum2=0;
		for(int i=0; i<histogram.length; i++){
			sum2+=(histogram[i]*i);
		}
		double avg=sum2/(double)sum;
		return avg;
	}
	
	public static final double standardDeviationHistogram(char[] histogram){
		long sum=max(1, sum(histogram));
		long sum2=0;
		for(int i=0; i<histogram.length; i++){
			sum2+=(histogram[i]*i);
		}
		double avg=sum2/(double)sum;
		double sumdev2=0;
		for(int i=0; i<histogram.length; i++){
			double dev=avg-i;
			double dev2=dev*dev;
			sumdev2+=(histogram[i]*dev2);
		}
		return Math.sqrt(sumdev2/sum);
	}
	
	public static final double standardDeviationHistogram(int[] histogram){
		long sum=max(1, sum(histogram));
		long sum2=0;
		for(int i=0; i<histogram.length; i++){
			sum2+=(histogram[i]*i);
		}
		double avg=sum2/(double)sum;
		double sumdev2=0;
		for(int i=0; i<histogram.length; i++){
			double dev=avg-i;
			double dev2=dev*dev;
			sumdev2+=(histogram[i]*dev2);
		}
		return Math.sqrt(sumdev2/sum);
	}
	
	public static final double standardDeviationHistogram(long[] histogram){
		long sum=max(1, sum(histogram));
		long sum2=0;
		for(int i=0; i<histogram.length; i++){
			sum2+=(histogram[i]*i);
		}
		double avg=sum2/(double)sum;
		double sumdev2=0;
		for(int i=0; i<histogram.length; i++){
			double dev=avg-i;
			double dev2=dev*dev;
			sumdev2+=(histogram[i]*dev2);
		}
		return Math.sqrt(sumdev2/sum);
	}
	
	/** Special version that calculates standard deviation based on unique kmers rather than overall events */
	public static final double standardDeviationHistogramKmer(long[] histogram){
		final long sum=sum(histogram);
		double sumU=0;
		for(int i=0; i<histogram.length; i++){
			long x=histogram[i];
			sumU+=(x/(double)max(i, 1));
		}
		double avg=sum/max(sumU, 1);
		double sumdev2=0;
		for(int i=1; i<histogram.length; i++){
			double dev=avg-i;
			double dev2=dev*dev;
			long x=histogram[i];
			sumdev2+=((x/(double)max(i, 1))*dev2);
		}
		return Math.sqrt(sumdev2/sumU);
	}
	
	public static final double standardDeviationHistogram(AtomicLongArray histogram){
		long sum=max(1, sum(histogram));
		long sum2=0;
		for(int i=0; i<histogram.length(); i++){
			sum2+=(histogram.get(i)*i);
		}
		double avg=sum2/(double)sum;
		double sumdev2=0;
		for(int i=0; i<histogram.length(); i++){
			double dev=avg-i;
			double dev2=dev*dev;
			sumdev2+=(histogram.get(i)*dev2);
		}
		return Math.sqrt(sumdev2/sum);
	}
	
	/** Special version that calculates standard deviation based on unique kmers rather than overall events */
	public static final double standardDeviationHistogramKmer(AtomicLongArray histogram){
		final long sum=sum(histogram);
		double sumU=0;
		for(int i=0; i<histogram.length(); i++){
			long x=histogram.get(i);
			sumU+=(x/(double)max(i, 1));
		}
		double avg=sum/max(sumU, 1);
		double sumdev2=0;
		for(int i=1; i<histogram.length(); i++){
			double dev=avg-i;
			double dev2=dev*dev;
			long x=histogram.get(i);
			sumdev2+=((x/(double)max(i, 1))*dev2);
		}
		return Math.sqrt(sumdev2/sumU);
	}
	
	public static final long[] downsample(long[] array, int bins){
		if(array==null || array.length==bins){return array;}
		assert(bins<=array.length);
		assert(bins>=0);
		long[] r=new long[bins];
		if(bins==0){return r;}
		double mult=bins/(double)array.length;
		for(int i=0; i<array.length; i++){
			int j=(int)(mult*i);
			r[j]+=array[i];
//			if(array[i]>0){System.err.println(i+"->"+j+": "+array[i]);}
		}
		return r;
	}

	
	public static final void pause(int millis){
		try {
			Thread.sleep(millis);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static final int min(int x, int y){return x<y ? x : y;}
	public static final int max(int x, int y){return x>y ? x : y;}
	public static final int min(int x, int y, int z){return x<y ? (x<z ? x : z) : (y<z ? y : z);}
	public static final int max(int x, int y, int z){return x>y ? (x>z ? x : z) : (y>z ? y : z);}
	public static final int min(int x, int y, int z, int z2){return min(min(x,y), min(z,z2));}
	public static final int max(int x, int y, int z, int z2){return max(max(x,y), max(z,z2));}
	
	//Median of 3
	public static final int mid(int x, int y, int z){return x<y ? (x<z ? min(y, z) : x) : (y<z ? min(x, z) : y);}

	public static final char min(char x, char y){return x<y ? x : y;}
	public static final char max(char x, char y){return x>y ? x : y;}

	public static final byte min(byte x, byte y){return x<y ? x : y;}
	public static final byte max(byte x, byte y){return x>y ? x : y;}
	public static final byte min(byte x, byte y, byte z){return x<y ? min(x, z) : min(y, z);}
	public static final byte max(byte x, byte y, byte z){return x>y ? max(x, z) : max(y, z);}
	public static final byte min(byte x, byte y, byte z, byte a){return min(min(x, y), min(z, a));}
	public static final byte max(byte x, byte y, byte z, byte a){return max(max(x, y), max(z, a));}

	public static final byte mid(byte x, byte y, byte z){return x<y ? (x<z ? min(y, z) : x) : (y<z ? min(x, z) : y);}
	
	public static final long min(long x, long y){return x<y ? x : y;}
	public static final long max(long x, long y){return x>y ? x : y;}
	public static final long min(long x, long y, long z){return x<y ? (x<z ? x : z) : (y<z ? y : z);}
	public static final long max(long x, long y, long z){return x>y ? (x>z ? x : z) : (y>z ? y : z);}
	public static final long min(long x, long y, long z, long z2){return min(min(x,y), min(z,z2));}
	public static final long max(long x, long y, long z, long z2){return max(max(x,y), max(z,z2));}
	public static final long mid(long x, long y, long z){return x<y ? (x<z ? min(y, z) : x) : (y<z ? min(x, z) : y);}
	public static final int longToInt(long x){return x<Integer.MIN_VALUE ? Integer.MIN_VALUE : x>Integer.MAX_VALUE ? Integer.MAX_VALUE : (int)x;}
	
	public static final double min(double x, double y){return x<y ? x : y;}
	public static final double max(double x, double y){return x>y ? x : y;}
	public static final double mid(double x, double y, double z){return x<y ? (x<z ? min(y, z) : x) : (y<z ? min(x, z) : y);}
	
	public static final float min(float x, float y){return x<y ? x : y;}
	public static final float max(float x, float y){return x>y ? x : y;}
	public static final float min(float x, float y, float z){return x<y ? (x<z ? x : z) : (y<z ? y : z);}
	public static final float max(float x, float y, float z){return x>y ? (x>z ? x : z) : (y>z ? y : z);}
	public static final float mid(float x, float y, float z){return x<y ? (x<z ? min(y, z) : x) : (y<z ? min(x, z) : y);}
	
	public static final int min(int[] array, int fromIndex, int toIndex){
		int min=array[fromIndex];
		for(int i=fromIndex+1; i<=toIndex; i++){
			min=min(min, array[i]);
		}
		return min;
	}
	
	public static final int max(int[] array, int fromIndex, int toIndex){
		int max=array[fromIndex];
		for(int i=fromIndex+1; i<=toIndex; i++){
			max=max(max, array[i]);
		}
		return max;
	}

	public static int minIndex(int[] array) {
		if(array==null || array.length<1){return -1;}
		float min=array[0];
		int index=0;
		for(int i=1; i<array.length; i++){
			if(array[i]<min){
				min=array[i];
				index=i;
			}
		}
		return index;
	}
	
	public static String trimWhitespace(String s){
		for(int i=0; i<s.length(); i++){
			if(Character.isWhitespace(s.charAt(i))){
				String s2=s.substring(0, i);
				return s2;
			}
		}
		return s;
	}
	
	public static double exponential(Random randy, double lamda){
//		for(int i=0; i<20; i++){
//			double p=randy.nextDouble();
//			double r=-Math.log(1-p)/lamda;
//			System.err.println(p+", "+lamda+"->"+"\n"+r);
//		}
//		assert(false);
		double p=randy.nextDouble();
		return -Math.log(1-p)/lamda;
	}
	
	public static double log2(double d){
		return Math.log(d)*invlog2;
	}
	
	public static double logRoot2(double d){
		return Math.log(d)*invlogRoot2;
	}
	
	public static double log1point2(double d){
		return Math.log(d)*invlog1point2;
	}

	private static final double log2=Math.log(2);
	private static final double invlog2=1/log2;
	private static final double logRoot2=Math.log(Math.sqrt(2));
	private static final double invlogRoot2=1/logRoot2;
	private static final double log1point2=Math.log(1.2);
	private static final double invlog1point2=1/log1point2;
	
	public static final char[] specialChars;
	
	public static final int[] ilens;
	public static final long[] llens;
	
	/** A single whitespace */
	public static final Pattern whitespace = Pattern.compile("\\s");
	/** Multiple whitespace */
	public static final Pattern whitespacePlus = Pattern.compile("\\s+");
	
	static{
		ilens=new int[Integer.toString(Integer.MAX_VALUE).length()+1];
		llens=new long[Long.toString(Long.MAX_VALUE).length()+1];
		for(int i=1, x=9; i<ilens.length; i++){
			ilens[i]=x;
			x=(x*10)+9;
		}
		ilens[ilens.length-1]=Integer.MAX_VALUE;
		for(long i=1, x=9; i<llens.length; i++){
			llens[(int)i]=x;
			x=(x*10)+9;
		}
		llens[llens.length-1]=Long.MAX_VALUE;
		
		specialChars=new char[256];
		Arrays.fill(specialChars, 'X');
		for(int i=0; i<32; i++){
			specialChars[i]=' ';
		}
		for(int i=32; i<127; i++){
			specialChars[i]=(char)i;
		}
		specialChars[127]=' ';
		specialChars[128]='C';
		specialChars[129]='u';
		specialChars[130]='e';
		specialChars[131]='a';
		specialChars[132]='a';
		specialChars[133]='a';
		specialChars[134]='a';
		specialChars[135]='c';
		specialChars[136]='e';
		specialChars[137]='e';
		specialChars[138]='e';
		specialChars[139]='i';
		specialChars[140]='i';
		specialChars[141]='i';
		specialChars[142]='S';
		specialChars[143]='S';
		specialChars[144]='E';
		specialChars[145]='a';
		specialChars[146]='a';
		specialChars[147]='o';
		specialChars[148]='o';
		specialChars[149]='o';
		specialChars[150]='u';
		specialChars[151]='u';
		specialChars[152]='y';
		specialChars[153]='O';
		specialChars[154]='U';
		specialChars[155]='c';
		specialChars[156]='L';
		specialChars[157]='Y';
		specialChars[158]='P';
		specialChars[159]='f';
		specialChars[160]='a';
		specialChars[161]='i';
		specialChars[162]='o';
		specialChars[163]='u';
		specialChars[164]='n';
		specialChars[165]='N';
		specialChars[166]='a';
		specialChars[167]='o';
		specialChars[168]='?';
		specialChars[224]='a';
		specialChars[224]='B';
		specialChars[230]='u';
		specialChars[252]='n';
		specialChars[253]='2';
	}
	
}

package sort;

import java.io.File;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.ArrayBlockingQueue;

import align2.TranslateColorspaceRead;
import stream.ConcurrentLegacyReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.RTextInputStream;
import stream.Read;
import stream.ReadStreamStringWriter;
import stream.ReadStreamWriter;
import structures.ListNum;
import dna.AminoAcid;
import dna.Data;
import dna.Gene;
import fileIO.ReadWrite;
import shared.Shared;
import shared.Timer;
import shared.Tools;

public class SortReadsByMapping {
	
	
	public static void main(String[] args){
		
		for(String s_ : args){
			String s=s_.toLowerCase();
			String split[]=(s.contains("=") ? s.split("=") : null);
			if(s.equalsIgnoreCase("merge")){MERGE_DUPLICATES=true;}
			else if(s.equalsIgnoreCase("regen")){REGENERATE_MATCH_STRING=true;}
			else if(s.equalsIgnoreCase("trim")){TRIM_LOW_QUALITY_TAILS=true;}
			else if(s.equalsIgnoreCase("fixshort")){FIX_SHORT_PAIRED_READS=true;}
			else if(s.equalsIgnoreCase("removesingletonduplicates")){REMOVE_SINGLETON_DUPLICATES_OF_PAIRS=true;}
			else if(s.equalsIgnoreCase("swaptoplus") || s.equalsIgnoreCase("swap")){SWAP_READ1_TO_PLUS=true;}
			else if(s.equalsIgnoreCase("mergeoppositestrand")){MERGE_OPPOSITE_STRAND_DUPLICATES=true;}
			else if(s.startsWith("merge=")){
				MERGE_DUPLICATES=(split[1].startsWith("t") || split[1].equals("1") ? true : false);
			}else if(s.startsWith("regen=")){
				REGENERATE_MATCH_STRING=(split[1].startsWith("t") || split[1].equals("1") ? true : false);
			}else if(s.startsWith("trim=")){
				TRIM_LOW_QUALITY_TAILS=(split[1].startsWith("t") || split[1].equals("1") ? true : false);
			}else if(s.startsWith("fixshort=")){
				FIX_SHORT_PAIRED_READS=(split[1].startsWith("t") || split[1].equals("1") ? true : false);
			}else if(s.startsWith("removesingletonduplicates=") || s.startsWith("removesingletonduplicatesofpairs=")){
				REMOVE_SINGLETON_DUPLICATES_OF_PAIRS=(split[1].startsWith("t") || split[1].equals("1") ? true : false);
			}else if(s.startsWith("minq=") || s.startsWith("minquality=") || s.startsWith("trimquality=")){
				TRIM_QUALITY=Byte.parseByte(split[1]);
			}else if(s.startsWith("window=") || s.startsWith("trimwindow=")){
				TRIM_WINDOW=Byte.parseByte(split[1]);
			}else if(s.startsWith("swaptoplus=") || s.startsWith("swap=")){
				SWAP_READ1_TO_PLUS=(split[1].startsWith("t") || split[1].equals("1") ? true : false);
			}else if(s.startsWith("mergeoppositestrand=")){
				MERGE_OPPOSITE_STRAND_DUPLICATES=(split[1].startsWith("t") || split[1].equals("1") ? true : false);;
			}else if(s.startsWith("readlimit=")){
				READ_LIMIT=Long.parseLong(split[1]);
				Data.sysout.println("Set READ_LIMIT to "+READ_LIMIT);
			}else if(s.startsWith("threads=")){
				REGEN_THREADS=Integer.parseInt(split[1]);
			}else if(s.startsWith("overwrite=")){
				overwrite=Tools.parseBoolean(split[1]);
			}
		}
		
		Read.DECOMPRESS_MATCH_ON_LOAD=true;
		
		SortReadsByMapping srt;
		
		String reads1=args[0];
		String reads2=args[1].equalsIgnoreCase("null") ?  null : args[1];
		String outname=args[2].equalsIgnoreCase("null") ?  ReadWrite.parseRoot(reads1)+"mapped_sorted#.txt.gz" : args[2];
		assert(outname.contains("#"));
		int blocksize=Integer.parseInt(args[3]);

		srt=new SortReadsByMapping(reads1, reads2, outname, blocksize);
		
		
		srt.process();

		double rmult=100d/(srt.processed);
		double bmult=100d/srt.basesInitiallyMapped;
		
		float pmult=(srt.paired ? 2 : 1);
		
		long remaining=srt.processed-srt.merged-srt.merged2-srt.removedSingletonDupe-srt.removedLQ-srt.removedShort;
		Data.sysout.println("Processed "+srt.processed+" reads; "+remaining+" remaining"+String.format(" (%.2f%%)", remaining*rmult));
		if(MERGE_DUPLICATES){
			Data.sysout.println("Merged "+srt.merged2+" strict duplicates"+String.format(" (%.2f%%)", srt.merged2*rmult));
			Data.sysout.println("Merged "+srt.merged+" duplicates"+String.format(" (%.2f%%)", srt.merged*rmult));
			if(srt.paired && REMOVE_SINGLETON_DUPLICATES_OF_PAIRS){
				Data.sysout.println("Removed "+srt.removedSingletonDupe+" singleton duplicates of pairs"+
						String.format(" (%.2f%%)", srt.removedSingletonDupe*rmult));
			}
		}
		if(FIX_SHORT_PAIRED_READS){
			Data.sysout.println("Removed "+srt.removedShort+" short reads"+String.format(" (%.2f%%)", srt.removedShort*rmult));
			Data.sysout.println("Trimmed "+srt.basesOverlapping+" overlapping bases of "+srt.basesInitiallyMapped+" initially mapped"+
					String.format(" (%.2f%%)", srt.basesOverlapping*bmult));
		}
		if(TRIM_LOW_QUALITY_TAILS){
			Data.sysout.println("Removed "+srt.removedLQ+" low-quality reads"+String.format(" (%.2f%%)", srt.removedLQ*rmult));
			Data.sysout.println("Trimmed "+srt.basesRemoved+" low-quality bases of "+srt.basesMapped+" mapped"+
					String.format(" (%.2f%%)", srt.basesRemoved*bmult));
		}
	}
	
	public SortReadsByMapping(String fname1, String fname2, String outname_, int blocksize_){
		assert(fname2==null || !fname1.equals(fname2)) : "Error - input files have same name.";
		
		long limit=READ_LIMIT;
		RTextInputStream rtis=new RTextInputStream(fname1, fname2, limit);
		outname=outname_;
		paired=rtis.paired();
		cris=new ConcurrentLegacyReadInputStream(rtis, limit);
		blocksize=blocksize_;
		assert(blocksize>200000);

		blockwriter1=(fname1==null ? null : new ReadStreamStringWriter(null, true, 4, false));
		blockwriter2=(fname2==null ? null : new ReadStreamStringWriter(null, false, 4, false));
	}
	
	public void process(){

		final String fname1=outname.replaceFirst("#", "1");
		final String fname2=(!paired ? null : outname.replaceFirst("#", "2"));
		if(!overwrite){
			if(fname1!=null && new File(fname1).exists()){throw new RuntimeException("Destination file "+fname1+" already exists.");}
			if(fname2!=null && new File(fname2).exists()){throw new RuntimeException("Destination file "+fname2+" already exists.");}
		}
		
		Timer t=new Timer();
		Timer total=new Timer();
		t.start();
		total.start();
		

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
				if(paired){
					asymmetricReads=(r.length()!=r.mateLength());
				}
			}
			
			while(reads!=null && reads.size()>0){
				//System.err.println("reads.size()="+reads.size());
				
				
				if(KILL_BAD_PAIRS && paired){
					for(Read r : reads){
						
						if(r.isBadPair(REQUIRE_CORRECT_STRANDS_PAIRS, SAME_STRAND_PAIRS, 20000)){
							int x=r.mapScore/r.length();
							int y=r.mate.mapScore/r.mateLength();
							if(x>=y){
								r.mate.clearAnswers(false);
							}else{
								r.clearAnswers(false);
							}
						}
						
						addRead(r);
					}
				}else{
					for(Read r : reads){addRead(r);}
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
		}

		
		synchronized(this){this.notifyAll();}
		System.err.println("Notified all");
		
		finishWritingBlocks();
		System.err.println("Wrote blocks");
		

		if(bwt1!=null){blockwriter1.poison();}
		if(bwt2!=null){blockwriter2.poison();}
		
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
		
		final ReadComparatorMapping mcomp=new ReadComparatorMapping();
		
		int lastChrom=-1;
		for(String key : keys){
			Block b=table.get(key);
			table.remove(key);
			processed+=b.added;
			
			if(UNLOAD_CHROMS_WHEN_DONE && lastChrom>-1 && b.chrom!=lastChrom){
				Data.unload(lastChrom, false); //Saves memory when regenerating match strings
			}
			lastChrom=b.chrom;
			
			if(b.added>MAX_BLOCKSIZE_TO_SORT){
				if(true){throw new RuntimeException("Skipping sorting for key "+key+" of size "+b.added);}
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
				
				Data.sysout.println(key+"\t"+b.added);
				b.delete();
			}else{
				ArrayList<Read> list=b.readBlock();
				Data.sysout.println(key+"\t"+list.size());
				b.delete();

				//Shared.sort(list, mcomp);
				if(MERGE_DUPLICATES){
					if(!paired){
						Shared.sort(list, mcomp);
						if(USE_STRICT_MERGE){findAndMergeDuplicatesStrict(list);}
						findAndMergeDuplicates(list, false);
					}else{
						
						//Possibly, doing two passes (unswap, merge, reswap, merge) is unnecessary...
						
						if(SWAP_READ1_TO_PLUS && paired){
							//Unswap
							for(int i=0; i<list.size(); i++){
								Read r=list.get(i);
								if(r!=null && r.mate!=null && r.swapped()){list.set(i, r.mate);}
							}
							//Merge
							doPairedSplitAndMergeSeries(list, mcomp, false);
							//Reswap
							for(int i=0; i<list.size(); i++){
								Read r=list.get(i);
								if(r!=null && r.mate!=null && r.swapped()){list.set(i, r.mate);}
							}
							
							if(MERGE_OPPOSITE_STRAND_DUPLICATES){
								//Merge
								doPairedSplitAndMergeSeries(list, mcomp, true);
							}else{
								Shared.sort(list, mcomp);
							}
						}else{
							doPairedSplitAndMergeSeries(list, mcomp, false);
						}
					}
				}
				
				//Unswap
				if(SWAP_READ1_TO_PLUS && paired){
					for(int i=0; i<list.size(); i++){
						Read r=list.get(i);
						if(r!=null && r.mate!=null && r.swapped()){
							list.set(i, r.mate);
						}
					}
				}
				
				if(FIX_SHORT_PAIRED_READS && paired){
					int[] rvector=new int[4];
					int removedTemp=fixShort(list, rvector);
					removedShort+=removedTemp;
					basesOverlapping+=rvector[1];
					basesInitiallyMapped+=rvector[2];
					int needRegen=rvector[3];
					
					if(REGENERATE_MATCH_STRING && needRegen>0){
						regenMatchStrings(list);
					}
				}else{
					for(Read r : list){
						if(r!=null){
							if(r.mapped()){basesInitiallyMapped+=r.length();}
							if(r.mate!=null && r.mate.mapped()){basesInitiallyMapped+=r.mate.length();}
						}
					}
				}
				
				if(TRIM_LOW_QUALITY_TAILS){
					int[] rvector=new int[4];
					int removedTemp=trimTails(list, TRIM_WINDOW, TRIM_QUALITY, rvector);
					removedLQ+=removedTemp;
					basesRemoved+=rvector[1];
					basesMapped+=rvector[2];
					int needRegen=rvector[3];
					
					if(REGENERATE_MATCH_STRING && needRegen>0){
						regenMatchStrings(list);
					}
				}else{
					for(Read r : list){
						if(r!=null){
							if(r.mapped() && !r.invalid()){basesMapped+=r.length();}
							if(r.mate!=null  && !r.mate.invalid() && r.mate.mapped()){basesMapped+=r.mateLength();}
						}
					}
				}
				
				//Reswap
				if(SWAP_READ1_TO_PLUS && paired){
					for(int i=0; i<list.size(); i++){
						Read r=list.get(i);
						if(r!=null && r.mate!=null && r.swapped()){
							list.set(i, r.mate);
						}
					}
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
		
		readsWritten=0;
		basesWritten=0;
		
		if(wtt1!=null){
			while(wtt1.isAlive()){
				try {
					wtt1.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			readsWritten+=wt1.readsWritten();
			basesWritten+=wt1.basesWritten();
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
			readsWritten+=wt2.readsWritten();
			basesWritten+=wt2.basesWritten();
		}
		
		t.stop();
		total.stop();
		Data.sysout.println("Final Sort + Write Time: "+t);
		Data.sysout.println("Total Time: "+total);
		
	}
	
	
	private void doPairedSplitAndMergeSeries(ArrayList<Read> list, final ReadComparatorMapping mcomp, boolean mergeDifferentLength){

		//This special section is probably not necessary.
		//Theoretically, keeping everything in a single list should work fine.
		
		int p=0, e1=0, e2=0, e12=0;
		for(Read r : list){
			if(r!=null){
				if(r.paired()){
					p++;
				}else if(r.mapped() && r.mate.mapped()){
					e12++;
				}else if(r.mapped()){
					e1++;
				}else if(r.mate.mapped()){
					e2++;
				}
			}
		}

		ArrayList<Read> listP=new ArrayList<Read>(p);
		ArrayList<Read> list1=new ArrayList<Read>(e1);
		ArrayList<Read> list2=new ArrayList<Read>(e2);
		ArrayList<Read> list12=new ArrayList<Read>(e12);
		
		for(Read r : list){
			if(r!=null){
				if(r.paired()){
					listP.add(r);
				}else if(r.mapped() && r.mate.mapped()){
					list12.add(r);
				}else if(r.mapped()){
					list1.add(r);
				}else if(r.mate.mapped()){
					list2.add(r);
				}
			}
		}
		list.clear();
		
		Shared.sort(listP, mcomp);
		if(USE_STRICT_MERGE){findAndMergeDuplicatesStrict(listP);}
		findAndMergeDuplicates(listP, mergeDifferentLength);
		list.addAll(listP);
		listP=null;
		
		Shared.sort(list1, mcomp);
		findAndMergeDuplicates(list1, mergeDifferentLength);
		list.addAll(list1);
		list1=null;
		
		Shared.sort(list2, mcomp);
		findAndMergeDuplicates(list2, mergeDifferentLength);
		list.addAll(list2);
		list2=null;
		
		Shared.sort(list12, mcomp);
		if(USE_STRICT_MERGE){findAndMergeDuplicatesStrict(list12);}
		findAndMergeDuplicates(list12, mergeDifferentLength);
		list.addAll(list12);
		list12=null;
		
		Tools.condense(list);
		Shared.sort(list, mcomp);
		if(REMOVE_SINGLETON_DUPLICATES_OF_PAIRS){
			findAndRemoveSingletonDuplicatesOfPairs(list);
			Tools.condense(list);
			Shared.sort(list, mcomp);
		}
	}
	
	
	private void doPairedSplitAndMergeSeries_old(ArrayList<Read> list, final ReadComparatorMapping mcomp){

		//This special section is probably not necessary.
		//Theoretically, keeping everything in a single list should work fine.
		
		int p=0, e1=0, e2=0, e12=0;
		for(Read r : list){
			if(r!=null){
				if(r.paired()){
					p++;
				}else if(r.mapped() && r.mate.mapped()){
					e12++;
				}else if(r.mapped()){
					e1++;
				}else if(r.mate.mapped()){
					e2++;
				}
			}
		}

		ArrayList<Read> listP=new ArrayList<Read>(p);
		ArrayList<Read> list1=new ArrayList<Read>(e1);
		ArrayList<Read> list2=new ArrayList<Read>(e2);
		ArrayList<Read> list12=new ArrayList<Read>(e12);
		
		for(Read r : list){
			if(r!=null){
				if(r.paired()){
					listP.add(r);
				}else if(r.mapped() && r.mate.mapped()){
					list12.add(r);
				}else if(r.mapped()){
					list1.add(r);
				}else if(r.mate.mapped()){
					list2.add(r);
				}
			}
		}
		list.clear();
		
		Shared.sort(listP, mcomp);
		if(USE_STRICT_MERGE){findAndMergeDuplicatesStrict(listP);}
		findAndMergeDuplicates(listP, false);
		if(asymmetricReads && SWAP_READ1_TO_PLUS){
			Tools.condense(listP);
			Shared.sort(listP, mcomp);
			findAndMergeDuplicates(listP, true);
		}
		list.addAll(listP);
		listP=null;
		
		Shared.sort(list1, mcomp);
		findAndMergeDuplicates(list1, false);
		if(asymmetricReads && SWAP_READ1_TO_PLUS){
			Tools.condense(list1);
			Shared.sort(list1, mcomp);
			findAndMergeDuplicates(list1, true);
		}
		list.addAll(list1);
		list1=null;
		
		Shared.sort(list2, mcomp);
		findAndMergeDuplicates(list2, false);
		if(asymmetricReads && SWAP_READ1_TO_PLUS){
			Tools.condense(list2);
			Shared.sort(list2, mcomp);
			findAndMergeDuplicates(list2, true);
		}
		list.addAll(list2);
		list2=null;
		
		Shared.sort(list12, mcomp);
		if(USE_STRICT_MERGE){findAndMergeDuplicatesStrict(list12);}
		findAndMergeDuplicates(list12, false);
		if(asymmetricReads && SWAP_READ1_TO_PLUS){
			Tools.condense(list12);
			Shared.sort(list12, mcomp);
			findAndMergeDuplicates(list12, true);
		}
		list.addAll(list12);
		list12=null;
		
		Tools.condense(list);
		Shared.sort(list, mcomp);
		if(REMOVE_SINGLETON_DUPLICATES_OF_PAIRS){
			findAndRemoveSingletonDuplicatesOfPairs(list);
			Tools.condense(list);
			Shared.sort(list, mcomp);
		}
	}
	
	
	private static int trimTails(ArrayList<Read> list, int thresh, byte minq, int[] rvector){
		
		int removed=0;
		int basesRemoved=0;
		int basesMapped=0;
		int needRegen=0;
		for(int i=0; i<list.size(); i++){
			final Read r=list.get(i);
			if(r!=null){
				final Read r2=r.mate;
				if(r.mapped() && !r.invalid()){
					basesMapped+=r.length();
					basesRemoved+=trimTail(r, thresh, minq);
					if(r.match==null && r.valid()){needRegen++;}
					else{
						assert(r.invalid() || TranslateColorspaceRead.verifyMatchString2(r, true));
					}
				}
				if(r2!=null && r2.mapped() && !r2.invalid()){
					basesMapped+=r2.length();
					basesRemoved+=trimTail(r2, thresh, minq);
					if(r2.match==null && r2.valid()){needRegen++;}
					else{
						assert(r2.invalid() || TranslateColorspaceRead.verifyMatchString2(r2, true));
					}
				}
				
				if((r.invalid() || !r.mapped()) && (r2==null || r2.invalid() || !r2.mapped())){
					removed++;
					list.set(i, null);
				}
			}
		}
		
		if(rvector!=null){
			rvector[0]=removed;
			rvector[1]=basesRemoved;
			rvector[2]=basesMapped;
			rvector[3]=needRegen;
		}
		
		return removed;
	}
	
	
	private static int fixShort(ArrayList<Read> list, int[] rvector){
		
		int removed=0;
		int basesRemoved=0;
		int basesMapped=0;
		int needRegen=0;
		for(int i=0; i<list.size(); i++){
			final Read r=list.get(i);
			final Read r2=(r==null ? null : r.mate);
			
			if(r!=null){
				if(r.mapped()){basesMapped+=r.length();}
				if(r2!=null && r2.mapped()){basesMapped+=r2.length();}
			}
			
			if(r!=null && r2!=null && r.mapped() && r.paired() && r2.mapped()){

				final int expectedMinLengthOuter=r.length()+r2.length();

				final int refLengthOuter;
				final int refLengthInner;
				
				if(SAME_STRAND_PAIRS){
					throw new RuntimeException("TODO");
				}else{
					if(r.strand()==Gene.PLUS && r2.strand()==Gene.MINUS){
						refLengthOuter=r2.stop-r.start+1;
						refLengthInner=r2.start-r.stop-1;
					}else if(r.strand()==Gene.MINUS && r2.strand()==Gene.PLUS){
						refLengthOuter=r.stop-r2.start+1;
						refLengthInner=r.start-r2.stop-1;
					}else{
						//Wrong strands - don't do anything
						refLengthOuter=expectedMinLengthOuter;
						refLengthInner=1;
					}
				}
				
				//TODO: Merge the qualities and calls of the trimmed portions
				if(false && refLengthOuter<expectedMinLengthOuter){
					//Short read type 1: outer distance is too short
					
					if(refLengthOuter<Tools.max(r.length()/2, 35)){
						r.setInvalid(true);
						r2.setInvalid(true);
						basesRemoved+=(r.length()+r2.length());
					}else{
						int dif=expectedMinLengthOuter-refLengthOuter;
						int rem1=dif/2;
						int rem2=dif-rem1;
						if(r2.length()-rem2<10){
							r2.setInvalid(true);
							r.setPaired(false);
							r2.setPaired(false);
							basesRemoved+=r2.length();
							rem1=dif;
							rem2=0;
						}else if(r.length()-rem1<10){
							r.setInvalid(true);
							r.setPaired(false);
							r2.setPaired(false);
							basesRemoved+=r.length();
							rem1=0;
							rem2=dif;
						}
						if(rem1>0){
							basesRemoved+=trimTailByXBases(r, rem1);
							if(r.match==null && r.valid()){needRegen++;}
							else{
								assert(r.invalid() || TranslateColorspaceRead.verifyMatchString2(r, true));
							}
						}
						if(rem2>0){
							basesRemoved+=trimTailByXBases(r2, rem2);
							if(r2.match==null && r2.valid()){needRegen++;}
							else{
								assert(r2.invalid() || TranslateColorspaceRead.verifyMatchString2(r2, true));
							}
						}
					}
				}else if(refLengthInner<=0 || refLengthOuter<expectedMinLengthOuter){
//					System.err.println(r.toText(false));
					//Short read type 2: inner distance is negative

					int cOverlap1=countCalledBasesOnOrAfterRefLoc(r, r.strand()==Gene.PLUS ? r2.start : r2.stop);
					int cOverlap2=countCalledBasesOnOrAfterRefLoc(r2, r2.strand()==Gene.PLUS ? r.start : r.stop);
					
					int overlap=Tools.max(cOverlap1, cOverlap2);
					
					if(overlap>0){
						int toRemain=r.length()+r2.length()-overlap;

						if(toRemain<Tools.max(r.length()/2, 34)){
							r.setInvalid(true);
							r2.setInvalid(true);
							basesRemoved+=(r.length()+r2.length());
//							System.err.println("Removed read. refLengthOuter="+refLengthOuter+", refLengthInner="+
//									refLengthInner+", cOverlap1="+cOverlap1+", cOverlap2="+cOverlap2+", toRemain="+toRemain+
//									"\n"+new String(r.match)+"\n"+new String(r2.match)+"\n"+
//									r.start+"~"+r.stop+"("+Gene.strandCodes[r.strand()]+")\n"+
//									r2.start+"~"+r2.stop+"("+Gene.strandCodes[r2.strand()]+")\n");
							
						}else{
//							System.err.println((cOverlap1==cOverlap2 ? "" : "****\t")+cOverlap1+", "+cOverlap2);
							int rem1=overlap/2;
							int rem2=overlap-rem1;
							if(r2.length()-rem2<10){
								r2.setInvalid(true);
								r.setPaired(false);
								r2.setPaired(false);
								basesRemoved+=r2.length();
								rem1=overlap;
								rem2=0;
							}else if(r.length()-rem1<10){
								r.setInvalid(true);
								r.setPaired(false);
								r2.setPaired(false);
								basesRemoved+=r.length();
								rem1=0;
								rem2=overlap;
							}
							if(rem1>0){
								basesRemoved+=trimTailByXBases(r, rem1);
								if(r.match==null && r.valid()){needRegen++;}
								else{
									assert(r.invalid() || TranslateColorspaceRead.verifyMatchString2(r, true));
								}
							}
							if(rem2>0){
								basesRemoved+=trimTailByXBases(r2, rem2);
								if(r2.match==null && r2.valid()){needRegen++;}
								else{
									assert(r2.invalid() || TranslateColorspaceRead.verifyMatchString2(r2, true));
								}
							}
						}
					}
				}
				
				if((r.invalid() || !r.mapped()) && (r2==null || r2.invalid() || !r2.mapped())){
					removed++;
					list.set(i, null);
				}
			}
		}
		
		if(rvector!=null){
			rvector[0]=removed;
			rvector[1]=basesRemoved;
			rvector[2]=basesMapped;
			rvector[3]=needRegen;
		}
		
		return removed;
	}

	
	//TODO: Add support for deletions
	/** thresh: Must see this many consecutive 'm' to stop. */ 
	private static int trimTail(Read r, int thresh, byte minq){
		byte[] bases=r.bases;
		byte[] match=r.match;
		byte[] quality=r.quality;
		
		assert(match!=null);
		if(r.strand()==Gene.MINUS){ //Remember to un-reverse later
			Tools.reverseInPlace(match);
		}
		

		int lastBadLoc=quality.length;
		int lastBadMLoc=match.length;
		int qloc=quality.length-1;
		int mloc=match.length-1;
		
		for(; mloc>=0 && qloc>=0; mloc--){

			assert(qloc<lastBadLoc) : "\n"+qloc+", "+lastBadLoc+", "+mloc+", "+lastBadMLoc+"\n"+r.toText(false)+"\n";
			assert(mloc<lastBadMLoc) : "\n"+qloc+", "+lastBadLoc+", "+mloc+", "+lastBadMLoc+"\n"+r.toText(false)+"\n";
			if(lastBadLoc-qloc>thresh){break;}
			
			byte m=match[mloc];
			byte q=quality[qloc];
			
			if(m=='D'){
				//do nothing
				lastBadLoc=qloc+1;
				lastBadMLoc=mloc;
			}else{
				if(q<minq || m!='m'){
					lastBadLoc=qloc;
					lastBadMLoc=mloc;
				}
				qloc--;
			}
		}
		
		if(lastBadLoc==quality.length){
			if(r.strand()==Gene.MINUS){Tools.reverseInPlace(match);}
			return 0;
		}
		
		if(lastBadLoc<6){
			if(r.strand()==Gene.MINUS){Tools.reverseInPlace(match);}
			r.setInvalid(true);
			return quality.length;
		}
		
//		{
//			if(r.strand()==Gene.MINUS){ //Remember to un-reverse later
//				Tools.reverseInPlace(match);
//			}
//			System.err.println("\nBefore:\n"+r.toText(false));
//			if(r.strand()==Gene.MINUS){ //Remember to un-reverse later
//				Tools.reverseInPlace(match);
//			}
//		}
		
		assert(lastBadLoc<quality.length);
		assert(lastBadMLoc<match.length);
		
		bases=Arrays.copyOf(bases, lastBadLoc);
		quality=Arrays.copyOf(quality, lastBadLoc);
		match=Arrays.copyOf(match, lastBadMLoc);
		
		if(r.strand()==Gene.MINUS){Tools.reverseInPlace(match);}
		
		boolean realign=false;
		int lengthOfMatchString=0;
		for(byte m : match){
			if(m=='m' || m=='N' || m=='s' || m=='S' || m=='D'){
				lengthOfMatchString++;
			}else if(m=='X' || m=='Y'){
				realign=true;
			}
		}
		
//		assert(!realign) : r.toText(false);
		
		if(realign){
			System.err.println("Killed match string while trimming this read:\n"+r.toText(false));
			r.match=null;
			match=null;
		}else{
			if(r.strand()==Gene.PLUS){
				r.stop=r.start+lengthOfMatchString-1;
			}else{
				r.start=r.stop-lengthOfMatchString+1;
			}
		}
		
		int trimmed=r.quality.length-quality.length;
		r.quality=quality;
		r.match=match;
		r.bases=bases;
		
		assert(trimmed>0);
		assert(r.quality.length==r.length());
		assert(r.match==null || r.match.length>=r.quality.length);
		
//		System.err.println("After:\n"+r.toText(false));
		
		return trimmed;
	}
	
	
	private static int trimTailByXBases(Read r, final int x){
		byte[] bases=r.bases;
		byte[] match=r.match;
		byte[] quality=r.quality;
		
		final int newLen=bases.length-x;
		
		if(newLen<6){
			r.setInvalid(true);
			return quality.length;
		}
		
		assert(match!=null);
		if(r.strand()==Gene.MINUS){ //Remember to un-reverse later
			Tools.reverseInPlace(match);
		}
		
		int qloc=quality.length-1;
		int mloc=match.length-1;
		
		for(; mloc>=0 && qloc>=newLen; mloc--){
			
			byte m=match[mloc];
//			byte q=quality[qloc];
			
			if(m=='D'){
				//do nothing
			}else{
				qloc--;
			}
		}
		
		while(mloc>=0 && match[mloc]=='D'){mloc--;}
		assert(qloc==newLen-1);
		
		bases=Arrays.copyOf(bases, newLen);
		quality=Arrays.copyOf(quality, newLen);
		match=Arrays.copyOf(match, mloc+1);
		
		if(r.strand()==Gene.MINUS){Tools.reverseInPlace(match);}
		
		boolean realign=false;
		int lengthOfMatchString=0;
		for(byte m : match){
			if(m=='m' || m=='N' || m=='s' || m=='S' || m=='D'){
				lengthOfMatchString++;
			}else if(m=='X' || m=='Y'){
				realign=true;
			}
		}
		
//		assert(!realign) : r.toText(false);
		
		if(realign){
			System.err.println("Killed match string while trimming this read:\n"+r.toText(false));
			r.match=null;
			match=null;
		}else{
			if(r.strand()==Gene.PLUS){
				r.stop=r.start+lengthOfMatchString-1;
			}else{
				r.start=r.stop-lengthOfMatchString+1;
			}
		}
		
		int trimmed=r.quality.length-quality.length;
		r.quality=quality;
		r.match=match;
		r.bases=bases;
		
		assert(trimmed>0);
		assert(r.quality.length==r.length());
		assert(r.match==null || r.match.length>=r.quality.length);
		
//		System.err.println("After:\n"+r.toText(false));
		
		return trimmed;
	}
	
	
	private static int countCalledBasesOnOrAfterRefLoc(Read r, final int rlimit){
		final int clen=r.length();
		byte[] match=r.match;
		
		if(r.strand()==Gene.PLUS){
			
//			final int rlimit=rlimit_0-1;
			
			int cloc=0;
			int mloc=0;
			int rloc=r.start;
			for(; mloc<match.length && rloc<rlimit; mloc++){
				byte m=match[mloc];

				if(m=='D'){
					rloc++;
				}else if(m=='X' || m=='Y' || m=='I'){
					cloc++;
				}else{
					cloc++;
					rloc++;
				}
			}
			
			if(rloc>=rlimit){
				
				if(rloc>rlimit){
					return clen+(rloc-rlimit);
				}
				
				int ret=clen-cloc;
				assert(rloc==rlimit) : "ret="+ret+", clen="+clen+", cloc="+cloc+",\n"+
					"rloc="+rloc+", rlimit="+rlimit+", mloc="+mloc+", mlen="+match.length+",\n"+
					"r.start="+r.start+", r.stop="+r.stop+", r2.start="+r.mate.start+", r2.stop="+r.mate.stop+"\n\n"+r.toText(false)+"\n\n";
				assert(ret>=0 && ret<=clen) : "ret="+ret+", clen="+clen+", cloc="+cloc+",\n"+
					"rloc="+rloc+", rlimit="+rlimit+", mloc="+mloc+", mlen="+match.length+",\n"+
					"r.start="+r.start+", r.stop="+r.stop+", r2.start="+r.mate.start+", r2.stop="+r.mate.stop;
				return ret;
			}else{
				assert(cloc==clen) : clen+", "+cloc+"\n"+r.toText(false)+"\n"; //Maybe cloc==clen
				return 0;
			}
		}else{

//			final int rlimit=rlimit_0+1;
			
			int cloc=clen-1;
			int mloc=match.length-1;
			int rloc=r.stop;
			for(; mloc>=0 && rloc>rlimit; mloc--){
				byte m=match[mloc];

				if(m=='D'){
					rloc--;
				}else if(m=='X' || m=='Y' || m=='I'){
					cloc--;
				}else{
					cloc--;
					rloc--;
				}
			}
			
			if(rloc<=rlimit){
				if(rloc<rlimit){
					return clen+(rlimit-rloc);
				}
				
				assert(rloc==rlimit);
				int ret=cloc+1;
				assert(ret>=0 && ret<=clen) : "ret="+ret+", clen="+clen+", cloc="+cloc+",\n"+
					"rloc="+rloc+", rlimit="+rlimit+", mloc="+mloc+", mlen="+match.length+",\n"+
					"r.start="+r.start+", r.stop="+r.stop+", r2.start="+r.mate.start+", r2.stop="+r.mate.stop;
				return ret;
			}else{
				assert(cloc==-1) : clen+", "+cloc; //Maybe cloc==-1
				return 0;
			}
		}
	}
	
	
	private void findAndMergeDuplicates(ArrayList<Read> list, boolean mergeDifferentLengthReads){
		if(list==null || list.size()<2){return;}
		Read current=list.get(0);
		
		for(int i=1; i<list.size(); i++){
			final Read r=list.get(i);
			final Read r2=r.mate;
			boolean merge=false;
			
			assert(paired==(current.mate!=null));
			
			final boolean dupeLoose=current.isDuplicateByMapping(r, false, false);
			final boolean mdupeLoose=(current.mate==null ? false : current.mate.isDuplicateByMapping(r.mate, false, false));
			final boolean lengthOK=(mergeDifferentLengthReads || 
					(r.length()==current.length() && (!paired || r2.length()==current.mateLength())));
			
			if(lengthOK && (dupeLoose || mdupeLoose)){
				if(paired){

					if(r.length()==current.length()){
						//Normal case

						if(dupeLoose && mdupeLoose){

							boolean dupeStrict=current.isDuplicateByMapping(r, true, true);
							boolean mdupeStrict=current.mate.isDuplicateByMapping(r2, true, true);
							if(dupeStrict && mdupeStrict){
								if(current.perfect() && current.mate.perfect() && r.perfect() && r2.perfect()){
									//									assert(!contains(current.match, 'N'));
									//									assert(!contains(r.match, 'N'));
									current.merge(r, true, true);
									//									assert(!contains(current.match, 'N'));
								}else{

									boolean N1=contains(current.bases, 'N');
									boolean N2=contains(current.mate.bases, 'N');
									boolean N3=contains(r.bases, 'N');
									boolean N4=contains(r2.bases, 'N');

									current.merge(r, false, true);

									if(!N1 || !N3){
										assert(!contains(current.bases, 'N')) : N1+", "+N2+", "+N3+", "+N4+"\n"+
										current.toText(false)+"\n"+r.toText(false)+"\n";
									}
									if(!N2 || !N4){
										assert(!contains(current.mate.bases, 'N')) : N1+", "+N2+", "+N3+", "+N4+"\n"+
										current.mate.toText(false)+"\n"+r2.toText(false)+"\n";
									}

								}
							}else{
								current.merge(r, false, false);
							}
							merge=true;
						}else if(current.paired() && r.paired()){
							//do nothing - not duplicates
						}else if(current.paired() && !r.paired()){
							//							if(dupe){
							//								r2.clearAnswers();
							//								current.merge(r, false, false);
							//								merge=true;
							//							}else{
							//								assert(mdupe);
							//								r.clearAnswers();
							//								current.merge(r, false, false);
							//								merge=true;
							//							}
						}else if(r.paired() && !current.paired()){ //This should not happen...
							//							if(dupe){
							//								current.mate.clearAnswers();
							//								r.merge(current, false, false);
							//								merge=false;
							//								merged++;
							//								list.set(cindex, null);
							//								current=r;
							//								cindex=i;
							//							}else{
							//								assert(mdupe);
							//								current.clearAnswers();
							//								r.merge(current, false, false);
							//								merge=false;
							//								merged++;
							//								list.set(cindex, null);
							//								current=r;
							//								cindex=i;
							//							}
						}else{//Neither is paired
							if(dupeLoose && !current.mate.mapped() && !r2.mapped()){
								boolean dupeStrict=current.isDuplicateByMapping(r, true, true);
								current.merge(r, false, dupeStrict);
								merge=true;
							}else if(mdupeLoose && !current.mapped() && !r.mapped()){
								boolean mdupeStrict=current.mate.isDuplicateByMapping(r2, true, true);
								current.merge(r, false, mdupeStrict);
								merge=true;
							}
						}
						
						assert(r2==null || r.numericID==r2.numericID);
						assert(current.mate==null || current.numericID==current.mate.numericID);
						assert(r.mate==null || r.mate.mate==r);
						assert(current.mate==null || current.mate.mate==current);
						
					}else{
						//Asymmetric paired reads, such as 50-35 from Solid, where one read got swapped
						
						if(dupeLoose && mdupeLoose){

							boolean dupeStrict=current.isDuplicateByMapping(r, false, true);
							boolean mdupeStrict=current.mate.isDuplicateByMapping(r2, false, true);
							if(dupeStrict && mdupeStrict){
								if(current.perfect() && current.mate.perfect() && r.perfect() && r2.perfect()){
									//									assert(!contains(current.match, 'N'));
									//									assert(!contains(r.match, 'N'));
									current.merge(r, true, true);
									//									assert(!contains(current.match, 'N'));
								}else{

//									boolean N1=contains(current.bases, 'N');
//									boolean N2=contains(current.mate.bases, 'N');
//									boolean N3=contains(r.bases, 'N');
//									boolean N4=contains(r2.bases, 'N');

									current.merge(r, false, true);

//									if(!N1 || !N3){
//										assert(!contains(current.bases, 'N')) : N1+", "+N2+", "+N3+", "+N4+"\n"+
//										current.toText(false)+"\n"+r.toText(false)+"\n";
//									}
//									if(!N2 || !N4){
//										assert(!contains(current.mate.bases, 'N')) : N1+", "+N2+", "+N3+", "+N4+"\n"+
//										current.mate.toText(false)+"\n"+r2.toText(false)+"\n";
//									}

								}
							}else{
								//Do nothing, but delete second copy
//								current.merge(r, false, false);
							}
							merge=true;
						}else if(current.paired() && r.paired()){
							//do nothing - not duplicates
						}else if(current.paired() && !r.paired()){
							//							if(dupe){
							//								r2.clearAnswers();
							//								current.merge(r, false, false);
							//								merge=true;
							//							}else{
							//								assert(mdupe);
							//								r.clearAnswers();
							//								current.merge(r, false, false);
							//								merge=true;
							//							}
						}else if(r.paired() && !current.paired()){ //This should not happen...
							//							if(dupe){
							//								current.mate.clearAnswers();
							//								r.merge(current, false, false);
							//								merge=false;
							//								merged++;
							//								list.set(cindex, null);
							//								current=r;
							//								cindex=i;
							//							}else{
							//								assert(mdupe);
							//								current.clearAnswers();
							//								r.merge(current, false, false);
							//								merge=false;
							//								merged++;
							//								list.set(cindex, null);
							//								current=r;
							//								cindex=i;
							//							}
						}else{//Neither is paired
							
							//In this case just remove the lower-quality read (generally the second one).
							//Should be very, very rare.
							
							merge=true; //Pretend to merge but really just delete the second copy.
							
//							if(dupeLoose && !current.mate.mapped() && !r2.mapped()){
//								boolean dupeStrict=current.isDuplicateByMapping(r, true, true);
//								current.merge(r, false, dupeStrict);
//								merge=true;
//							}else if(mdupeLoose && !current.mapped() && !r.mapped()){
//								boolean mdupeStrict=current.mate.isDuplicateByMapping(r2, true, true);
//								current.merge(r, false, mdupeStrict);
//								merge=true;
//							}
						}
					}

					assert(r2==null || r.numericID==r2.numericID);
					assert(current.mate==null || current.numericID==current.mate.numericID);
					assert(r.mate==null || r.mate.mate==r);
					assert(current.mate==null || current.mate.mate==current);
					
				}else{
					//Single-ended
					
					assert(dupeLoose);
					boolean dupeStrict=current.isDuplicateByMapping(r, true, true);
					if(current.perfect() && r.perfect()){
						current.merge(r, true, true);
					}else{
						current.merge(r, false, dupeStrict);
					}
					merge=true;
				}
			}
			
			if(merge){
				merged++;
				list.set(i, null);
			}else{
				current=r;
			}
		}
		
		if(REGENERATE_MATCH_STRING){regenMatchStrings(list);}
	}
	
	
	private void findAndMergeDuplicatesStrict(ArrayList<Read> list){
		if(list==null || list.size()<2){return;}
		
		int addIndex=0;
		
		ArrayList<Read> toMerge=new ArrayList<Read>();
		ArrayList<Read> toMerge2=new ArrayList<Read>();

		for(int i=0; i<list.size(); i++){
			Read r=list.set(i, null);

			if(!toMerge.isEmpty()){
				Read current=toMerge.get(0);
				assert(paired==(current.mate!=null));

				final boolean dupeStrict=current.isDuplicateByMapping(r, true, true);

				if(!paired){
					if(!dupeStrict){
						Read x=toMerge.get(0);
						if(toMerge.size()>1){
							merged2+=toMerge.size()-1;
							x=mergeReads(toMerge, true);
						}
						assert(list.get(addIndex)==null);
						list.set(addIndex, x);
						addIndex++;
						toMerge.clear();
					}
				}else{
					assert(toMerge.size()==toMerge2.size());
					final boolean mdupeStrict=current.mate.isDuplicateByMapping(r.mate, true, true);
					if(!dupeStrict || !mdupeStrict){
						Read x=toMerge.get(0);
						if(toMerge.size()>1){
							merged2+=toMerge.size()-1;
							x=mergeReads(toMerge, true);
							Read y=mergeReads(toMerge2, true);
							assert(x.mate==y);
							assert(y.mate==x);
							assert(x!=y);
						}
						assert(list.get(addIndex)==null);
						list.set(addIndex, x);
						addIndex++;
						toMerge.clear();
						toMerge2.clear();
					}
				}
			}
			
			toMerge.add(r);
			if(paired){toMerge2.add(r.mate);}
		}

		if(!toMerge.isEmpty()){
			Read x=toMerge.get(0);
			if(toMerge.size()>1){
				merged2+=toMerge.size()-1;
				x=mergeReads(toMerge, true);
				if(paired){
					Read y=mergeReads(toMerge2, true);
					assert(x.mate==y);
					assert(y.mate==x);
					assert(x!=y);
				}
			}
			assert(list.get(addIndex)==null);
			list.set(addIndex, x);
			addIndex++;
		}

		for(int i=list.size()-1; i>=0 && list.get(i)==null; i--){list.remove(i);}

		if(REGENERATE_MATCH_STRING){regenMatchStrings(list);}
	}
	
	
	private void findAndRemoveSingletonDuplicatesOfPairs(ArrayList<Read> list){
		if(list==null || list.size()<2){return;}
		assert(paired);
		
		Read current=null;
		for(int i=0; i<list.size(); i++){
			Read r=list.get(i);
			if(r.paired()){current=r;}
			else if(current!=null && !r.paired()){
				assert(r.length()==current.length()) : "Merging different-length reads is supported but seems to be not useful.";
				if(r.mapped() && !r.mate.mapped()){
					final boolean sameLength=current.length()==r.length();
					if(current.isDuplicateByMapping(r, sameLength, true)){
						if(current.mate.isDuplicateByBases(r.mate, r.mateLength()/5, r.mateLength()/5, (byte)14, false, true)){
							list.set(i, null);
							removedSingletonDupe++;
						}
					}
				}else{
					final boolean sameLength=current.mateLength()==r.mateLength();
					if(current.mate.isDuplicateByMapping(r.mate, sameLength, true)){
						if(!sameLength || current.isDuplicateByBases(r, r.length()/5, r.length()/5, (byte)14, false, true)){
							list.set(i, null);
							removedSingletonDupe++;
						}
					}
				}
			}
		}

		if(REGENERATE_MATCH_STRING){regenMatchStrings(list);}
	}

	
	/** Assumes all reads in the list are duplicates. 
	 * Only merges reads in list, not their mates. 
	 * Mutates the first input read. */
	public static Read mergeReads(ArrayList<Read> list, boolean retainPerfect){
		if(list==null || list.isEmpty()){return null;}
		if(list.size()==1){return list.get(0);}
		
		//This block prevents the destruction of perfect reads.
		{
			Read a=list.get(0);
			for(int i=1; i<list.size(); i++){
				Read b=list.get(i);
				a.copies+=b.copies;
				if(a.perfect() && !b.perfect()){
					list.set(i, null);
				}
			}
			Tools.condense(list);
			if(list.size()==1){return list.get(0);}
		}
		
		int len=0;
		for(Read r : list){len=Tools.max(len, r.length());}
		
		assert(len==list.get(0).length());
		
		int[][] count=new int[4][len];
		int[][] qual=new int[4][len];
		byte[][] maxQual=new byte[4][len];
		
		for(Read r : list){
			for(int i=0; i<r.length(); i++){
				byte b=r.bases[i];
				assert((b!='A' && b!='C' && b!='G' && b!='T'));
				byte q=r.quality[i];
				if(b>=0 && b<=3){
					count[b][i]+=r.copies;
					qual[b][i]+=q;
					maxQual[b][i]=Tools.max(q, maxQual[b][i]);
				}
			}
		}
		
		int[] carray=new int[4];
		int[] qarray=new int[4];
		byte[] marray=new byte[4];

		
		byte[] bases=new byte[len];
		byte[] quality=new byte[len];
		
		for(int i=0; i<len; i++){
			for(int j=0; j<4; j++){
				carray[j]=count[j][i];
				qarray[j]=qual[j][i];
				marray[j]=maxQual[j][i];
				
				qarray[j]=(qarray[j]+marray[j])/2;
			}
			byte best=findBest(carray, qarray, marray);
			if(best<0){
				bases[i]='N';
				quality[i]=0;
			}else{
				bases[i]=AminoAcid.numberToBase[best];
				int q=2*qarray[best];
				for(int j=0; j<4; j++){q-=qarray[j];}
				q=Tools.min(48, Tools.max(q, 1));
				quality[i]=(byte)q;
			}
		}
		
		final Read r=list.get(0);
		
//		if(r.match!=null){
//			for(int i=0; i<len; i++){
//				if(bases[i]!=r.bases[i]){
//					r.match=null;
//					break;
//				}
//			}
//		}
		
		//Uses the primary read as a template and merging the new data in.
		//That way, if the new data differs from the primary (best) read, but the new quality score is very low,
		//the old data can be retained (with quality reduced).
		//This may or may not be a good idea.
		
		boolean killMatch=false;
		final boolean retain=retainPerfect && r.perfect();
		for(int i=0; i<bases.length; i++){
			byte b1=r.bases[i];
			byte b2=bases[i];
			byte q1=r.quality[i];
			byte q2=quality[i];
			
			assert(q1>=0);
			assert(q2>=0);
			
			if(b1==b2){
				r.quality[i]=q2;
				if(b1=='N'){r.quality[i]=0;}
			}else{
				if(b2=='N'){
					r.quality[i]=Tools.min(r.quality[i], (byte)2);
				}else if(b1=='N'){
					r.bases[i]=b2;
					r.quality[i]=q2;
					killMatch=true;
				}else{
					if(retain){
						r.quality[i]=Tools.max((byte)2, (byte)(q1-q2));
					}else if(q2-q1>10){
						r.bases[i]=b2;
						r.quality[i]=q2;
						killMatch=true;
					}else if(q1<15 && q2>20){
						r.bases[i]=b2;
						r.quality[i]=q2;
						killMatch=true;
					}else{
						r.quality[i]=Tools.max((byte)2, (byte)(q1-q2));
					}
				}
			}
		}
		
		if(killMatch){r.match=null;}
		
		return r;
	}
	
	private static final byte findBest(int[] count, int[] qual, byte[] maxqual){
		byte best=-1;
		int bestScore=-1;
		for(byte i=0; i<count.length; i++){
			int score=count[i]*8+qual[i]*2+maxqual[i];
			if(score>0){
				if(score>bestScore){
					best=i;
					bestScore=score;
				}else if(score==bestScore){
					if(qual[i]>qual[best]){
						best=i;
						bestScore=score;
					}else if(qual[i]==qual[best] && count[i]>count[best]){
						best=i;
						bestScore=score;
					}
				}
			}
		}
		return best;
	}
	
	
	private void regenMatchStrings(ArrayList<Read> list){
		if(list==null || list.isEmpty()){return;}
		
		int needed=0;
		for(Read r : list){
			if(r!=null){
				if(r.mapped() && r.match==null){
					needed++;
				}else if(r.mate!=null && r.mate.mapped() && r.mate.match==null){
					needed++;
				}
			}
		}
		if(needed<1){return;}
		
		
		final int lim=100;
		
//		System.err.println("Starting RMTs");
		RegenMatchThread[] rmt=new RegenMatchThread[Tools.max(1, Tools.min(REGEN_THREADS, needed/lim))];
		for(int i=0; i<rmt.length; i++){
			rmt[i]=new RegenMatchThread();
			rmt[i].start();
		}
		ArrayList<Read> list2=new ArrayList<Read>(lim);
		for(Read r : list){
			if(r!=null){
				boolean flag=false;
				if(r.mapped() && r.match==null){
					flag=true;
				}else if(r.mate!=null && r.mate.mapped() && r.mate.match==null){
					flag=true;
				}
				if(flag){
					list2.add(r);
					if(list2.size()>=lim){
						while(list2!=null){
							try {
								REGEN_PIPE.put(list2);
								list2=null;
							} catch (InterruptedException e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
						}
						list2=new ArrayList<Read>(lim); 
					}
				}
			}
		}

		if(list2!=null && list2.size()>0){
			while(list2!=null){
				try {
					REGEN_PIPE.put(list2);
					list2=null;
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}

//		System.err.println("Poisoning RMTs");
		//Poison
		for(int i=0; i<rmt.length; i++){
			list2=new ArrayList<Read>(0);
			while(list2!=null){
				try {
					REGEN_PIPE.put(list2);
					list2=null;
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}

//		System.err.println("Joining RMTs");
		for(int i=0; i<rmt.length; i++){
			while(rmt[i].isAlive()){
				try {
					rmt[i].join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}


	}
	
	
	private static boolean contains(byte[] array, char c){
		if(array==null){return false;}
		for(byte b : array){if(b==c){return true;}}
		return false;
	}
	
	private String makeKey(Read r){
		String key;

		if(!r.mapped() && (r.mate==null || !r.mate.mapped())){
			key="Z";
		}else{
			if(r.mapped()){key=makeKey2(r, r.strand()).insert(0, 'C').toString();}
			else{
				byte strand=r.mate.strand();
				if(!SAME_STRAND_PAIRS){
					strand=(byte)(1^strand);
				}
				key=makeKey2(r.mate, strand).insert(0, MOVE_SINGLETONS_TO_END ? 'D' : 'C').toString();
			}
		}
		
		return key;
	}
	
	
	private StringBuilder makeKey2(Read r, byte strand){
		assert(r.mapped());

		StringBuilder sb=new StringBuilder(32);
		
		int num=Tools.max(r.start/blocksize, 0);
		if(r.chrom<=9){sb.append('0');}
		sb.append(r.chrom);
		sb.append("_").append(strand).append("_");
		for(int i=6-(""+num).length(); i>0; i--){
			sb.append('0');
		}
		sb.append(num);
		return sb;
	}
	
	
	private void addRead(Read r){
		Read r2=r.mate;
		assert(r2==null || r.numericID==r2.numericID);
		boolean swap=false;
		if(SWAP_READ1_TO_PLUS && r2!=null){
			
			if(r.paired() && r.mapped() && r.valid() && r2.mapped() && r2.valid()){ //Ideal pair
				if(r.strand()==Gene.MINUS && r2.strand()==Gene.PLUS){swap=true;}
			}else if(r.mapped() && r.valid() && r2.mapped() && r2.valid()){
				if(r.strand()==Gene.MINUS && r2.strand()==Gene.PLUS){swap=true;}
			}else if(r.mapped() && r.valid()){
				if(r.strand()==Gene.MINUS){swap=true;}
			}else if(r2.mapped() && r2.valid()){
				if(r2.strand()==Gene.PLUS){swap=true;}
			}
		}
		
		if(swap){
			r.setSwapped(true);
			r2.setSwapped(true);
			Read temp=r;
			r=r2;
			r2=temp;
		}
		assert(r2==null || (r.numericID==r2.numericID && r!=r2));
		
		String key=makeKey(r);
		
//		String key=sb.toString();
		Block b=table.get(key);
		if(b==null){
			//System.err.println("Created block "+key);
			b=new Block(key, outname, r.chrom);
			table.put(key, b);
		}
		b.add(r);
	}
	
	
	public void finishWritingBlocks(){
		System.err.println("Called finishWritingBlocks()");
		int numWritten=0;
		for(String key : table.keySet()){
			Block b=table.get(key);
			b.finishWritingBuffer();
			numWritten++;
		}
		assert(numWritten==table.size()) : "Only wrote "+numWritten+" of "+table.size();
	}
	
	
	private class Block{
		
		public Block(String name_, String fname_, int chrom_){
			
			if(DONT_COMPRESS_TEMP_FILES){
				if(fname_.endsWith(".gz") || fname_.endsWith(".zip") || fname_.endsWith(".bz2")){
					fname_=fname_.substring(0, fname_.lastIndexOf('.'));
				}
			}
			
			name=name_;
			fname1=fname_.replaceFirst("#", "_msort_tempBlock_"+name+"_1");
			fname2=(!paired ? null : fname_.replaceFirst("#", "_msort_tempBlock_"+name+"_2"));
			chrom=chrom_;
//			Data.sysout.println(fname1);
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
		public final int chrom; //Necessary for unloading data
		
		public final OutputStream outStream1, outStream2;
		public final PrintWriter writer1, writer2;
		private ArrayList<Read> buffer=new ArrayList<Read>(WRITE_BUFFER);
		
		public long added=0, written=0, numRead=0;
	}
	
	private class RegenMatchThread extends Thread{
		
		@Override
		public void run(){
			for(ArrayList<Read> list=take(); !list.isEmpty(); list=take()){
				for(Read r : list){
					if(r!=null){
						final Read r2=r.mate;
						if(r.mapped() && r.match==null && r.valid()){regenMatchString(r);}
						if(r2!=null && r2.mapped() && r2.match==null && r.valid()){regenMatchString(r2);}
					}
				}
			}
		}
		
		private void regenMatchString(Read r){
			assert(r.match==null);
			
			TranslateColorspaceRead.realign_new(r.topSite(), r.bases, tcr.msaBS, 4, 1, 0, false, true, r.numericID);
			assert(false) : "TODO: move ss locs back to read.";
			
			r.setPerfectFlag(Integer.MAX_VALUE);
			assert(!r.perfect() || r.stop-r.start==(r.length()-1)) : 
				"\n"+r.toText(false)+"\n"+new String(r.bases)+"\n"+new String(AminoAcid.reverseComplementBases(r.bases))+
				"\n"+Data.getChromosome(r.chrom).getString(r.topSite().start, r.topSite().stop)+"\n";
			
			if(r.match!=null){
//				boolean xy=TranslateColorspaceRead.containsXY(r.match);
				assert(TranslateColorspaceRead.verifyMatchString2(r, true)) : r.toText(false);
			}
		}
		
		private ArrayList<Read> take(){
			ArrayList<Read> list=null;
			while(list==null){
				try {
					list=REGEN_PIPE.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			return list;
		}
		
		private final TranslateColorspaceRead tcr=null; /*new TranslateColorspaceRead(2000, 3000);*/ //Specific type needs to be specified.
	}
	
	public final String outname;
	private final ConcurrentReadInputStream cris;
	private final ArrayBlockingQueue<ArrayList<Read>> REGEN_PIPE=new ArrayBlockingQueue<ArrayList<Read>>(40);
	public long merged=0;
	public long merged2=0;
	public long removedSingletonDupe=0;
	public long removedLQ=0;
	public long removedShort=0;
	public long processed=0;
	public long basesInitiallyMapped=0;
	public long basesOverlapping=0;
	public long basesMapped=0;
	public long basesRemoved=0;
//	public long numSwapped=0;
	private long readsWritten;
	private long basesWritten;
	
	private boolean asymmetricReads=false;
	
	private final HashMap<String, Block> table=new HashMap<String, Block>(4096);
	
	public final boolean paired;
	public final int blocksize;

	public static boolean USE_CRIS=true; //Similar speed either way.  "true" may be better with many threads.
	public static boolean MOVE_SINGLETONS_TO_END=false;

	public static long READ_LIMIT=-1; //Max number of reads to process
	public static final int WRITE_BUFFER=8000; //Bigger number uses more memory, for less frequent writes.
	public static final int MAX_BLOCKSIZE_TO_SORT=8000000;
	public static boolean overwrite=false;
	public static boolean append=false;
	
	public static final boolean DONT_COMPRESS_TEMP_FILES=false;
	public static boolean MERGE_DUPLICATES=false;
	public static final boolean KILL_BAD_PAIRS=true;
	public static boolean SAME_STRAND_PAIRS=false;
	public static boolean REQUIRE_CORRECT_STRANDS_PAIRS=true;
	public static boolean REMOVE_SINGLETON_DUPLICATES_OF_PAIRS=true;
	public static boolean USE_STRICT_MERGE=false;
	
	public static boolean SWAP_READ1_TO_PLUS=false;
	public static boolean MERGE_OPPOSITE_STRAND_DUPLICATES=false; //Requires SWAP_READ1_TO_PLUS=true
	
	public static final boolean UNLOAD_CHROMS_WHEN_DONE=true;
	
	public static boolean FIX_SHORT_PAIRED_READS=false;
	
	public static boolean TRIM_LOW_QUALITY_TAILS=false;
	public static byte TRIM_QUALITY=7;
	public static byte TRIM_WINDOW=3;
	
	public static boolean REGENERATE_MATCH_STRING=false;
	public static int REGEN_THREADS=Shared.threads();
	
	private final ReadStreamWriter blockwriter1;
	private final ReadStreamWriter blockwriter2;
	
//	private final TranslateColorspaceRead tcr2=new TranslateColorspaceRead(200, 2400);
}

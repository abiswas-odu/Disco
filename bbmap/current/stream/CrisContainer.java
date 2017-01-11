package stream;

import java.util.ArrayList;
import java.util.Comparator;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import structures.ListNum;

public class CrisContainer implements Comparable<CrisContainer> {
	
	public CrisContainer(String fname, Comparator<Read> comparator_){
		comparator=comparator_;
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, false, true);
		cris=ConcurrentReadInputStream.getReadInputStream(-1, true, ff, null, null, null);
		cris.start();
		fetch();
	}
	
	public CrisContainer(ConcurrentReadInputStream cris_, Comparator<Read> comparator_){
		comparator=comparator_;
		cris=cris_;
		fetch();
	}
	
	public ArrayList<Read> fetch(){
		final ArrayList<Read> old=list;
		ListNum<Read> ln=cris.nextList();
		list=(ln==null ? null : ln.list);
		if(list.size()<1){list=null;}
		read=(list==null ? null : list.get(0));
		if(lastNum>=0){cris.returnList(lastNum, list==null);}
		if(ln!=null){lastNum=ln.id;}
//		if(list==null){close();}
		
		assert((read==null)==(list==null || list.size()==0));
		
		return old;
	}
	
	public boolean close(){
		return ReadWrite.closeStream(cris);
	}
	
	public Read peek(){return read;}
	
	@Override
	public int compareTo(CrisContainer other) {
		assert(read!=null);
		assert(other.read!=null);
		return comparator.compare(read, other.read);
	}
	
	public int compareTo(Read other) {
		return comparator.compare(read, other);
	}
	
	public boolean hasMore(){
		return read!=null;
	}
	
	public ConcurrentReadInputStream cris(){return cris;}
	
	final ConcurrentReadInputStream cris;
	private Read read;
	private long lastNum=-1;
	private ArrayList<Read> list;
	private final Comparator<Read> comparator;
	
}

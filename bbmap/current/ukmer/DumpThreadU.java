package ukmer;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import kmer.DumpThread;
import shared.Shared;
import shared.Tools;
import stream.ByteBuilder;

import fileIO.ByteStreamWriter;

/**
 * @author Brian Bushnell
 * @date Nov 16, 2015
 *
 */
public class DumpThreadU extends Thread{
	
	public static boolean dump(final int k, final int mincount, final AbstractKmerTableU[] tables, final ByteStreamWriter bsw){
		final int threads=DumpThread.NUM_THREADS>0 ? DumpThread.NUM_THREADS : Tools.min(tables.length, (Tools.mid(1, Shared.threads()-1, 6)));
		final AtomicInteger lock=new AtomicInteger(0);
		final ArrayList<DumpThreadU> list=new ArrayList<DumpThreadU>(threads);
		for(int i=0; i<threads; i++){
			list.add(new DumpThreadU(k, mincount, lock, tables, bsw));
		}
		for(DumpThreadU t : list){t.start();}
		boolean success=true;
		for(DumpThreadU t : list){
			while(t.getState()!=Thread.State.TERMINATED){
				try {
					t.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			success&=t.success;
		}
		return success;
	}
	
	public DumpThreadU(final int k_, final int mincount_, final AtomicInteger nextTable_, final AbstractKmerTableU[] tables_, final ByteStreamWriter bsw_){
		k=k_;
		mincount=mincount_;
		nextTable=nextTable_;
		tables=tables_;
		bsw=bsw_;
	}
	
	@Override
	public void run(){
		final ByteBuilder bb=new ByteBuilder(16300);
		for(int i=nextTable.getAndIncrement(); i<tables.length; i=nextTable.getAndIncrement()){
			AbstractKmerTableU t=tables[i];
			t.dumpKmersAsBytes_MT(bsw, bb, k, mincount);
		}
		if(bb.length()>0){
			synchronized(bsw){bsw.addJob(bb);}
		}
		success=true;
	}
	
	final int k;
	final int mincount;
	final AtomicInteger nextTable; 
	final AbstractKmerTableU[] tables;
	final ByteStreamWriter bsw;
	boolean success=false;
	
}

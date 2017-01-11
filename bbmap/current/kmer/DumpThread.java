package kmer;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import stream.ByteBuilder;

import fileIO.ByteStreamWriter;
import shared.Shared;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Nov 16, 2015
 *
 */
public class DumpThread extends Thread{
	
	public static boolean dump(final int k, final int mincount, final AbstractKmerTable[] tables, final ByteStreamWriter bsw){
		final int threads=NUM_THREADS>0 ? NUM_THREADS : Tools.min(tables.length, (Tools.mid(1, Shared.threads()-1, 6)));
		final AtomicInteger lock=new AtomicInteger(0);
		final ArrayList<DumpThread> list=new ArrayList<DumpThread>(threads);
		for(int i=0; i<threads; i++){
			list.add(new DumpThread(k, mincount, lock, tables, bsw));
		}
		for(DumpThread t : list){t.start();}
		boolean success=true;
		for(DumpThread t : list){
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
	
	public DumpThread(final int k_, final int mincount_, final AtomicInteger nextTable_, final AbstractKmerTable[] tables_, final ByteStreamWriter bsw_){
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
			AbstractKmerTable t=tables[i];
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
	final AbstractKmerTable[] tables;
	final ByteStreamWriter bsw;
	boolean success=false;
	
	public static int NUM_THREADS=-1;
	
}

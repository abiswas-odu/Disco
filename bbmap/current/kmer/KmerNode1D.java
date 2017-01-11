package kmer;

import stream.ByteBuilder;
import fileIO.ByteStreamWriter;

/**
 * @author Brian Bushnell
 * @date Oct 22, 2013
 *
 */
public class KmerNode1D extends KmerNode {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public KmerNode1D(long pivot_){
		super(pivot_);
	}
	
	public KmerNode1D(long pivot_, int value_){
		super(pivot_);
		value=value_;
	}
	
	public final KmerNode makeNode(long pivot_, int value_){
		return new KmerNode1D(pivot_, value_);
	}
	
	public final KmerNode makeNode(long pivot_, int[] values_){
		throw new RuntimeException("Unimplemented");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final int set(long kmer, int[] vals) {
		throw new RuntimeException("Unimplemented.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	protected int value(){return value;}
	
	protected int[] values(int[] singleton){
		assert(singleton.length==1);
		singleton[0]=value;
		return singleton;
	}
	
	public int set(int value_){return value=value_;}
	
	protected int set(int[] values_){
		throw new RuntimeException("Unimplemented");
	}
	
	int numValues(){return value<1 ? 0 : 1;}
	
	/*--------------------------------------------------------------*/
	/*----------------       Private Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	boolean canResize() {
		return false;
	}
	
	@Override
	public boolean canRebalance() {
		return true;
	}
	
	@Deprecated
	@Override
	public int arrayLength() {
		throw new RuntimeException("Unsupported.");
	}
	
	@Deprecated
	@Override
	void resize() {
		throw new RuntimeException("Unsupported.");
	}
	
	@Deprecated
	@Override
	public void rebalance() {
		throw new RuntimeException("Please call rebalance(ArrayList<KmerNode>) instead, with an empty list.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount){
		if(value<1){return true;}
		if(value>=mincount){bsw.printlnKmer(pivot, value, k);}
		if(left!=null){left.dumpKmersAsBytes(bsw, k, mincount);}
		if(right!=null){right.dumpKmersAsBytes(bsw, k, mincount);}
		return true;
	}
	
	@Override
	public final boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount){
		if(value<1){return true;}
		if(value>=mincount){
			toBytes(pivot, value, k, bb);
			bb.append('\n');
			if(bb.length()>=16000){
				ByteBuilder bb2=new ByteBuilder(bb);
				synchronized(bsw){bsw.addJob(bb2);}
				bb.clear();
			}
		}
		if(left!=null){left.dumpKmersAsBytes_MT(bsw, bb, k, mincount);}
		if(right!=null){right.dumpKmersAsBytes_MT(bsw, bb, k, mincount);}
		return true;
	}
	
	@Override
	protected final StringBuilder dumpKmersAsText(StringBuilder sb, int k, int mincount){
		if(value<1){return sb;}
		if(sb==null){sb=new StringBuilder(32);}
		if(value>=mincount){sb.append(AbstractKmerTable.toText(pivot, value, k)).append('\n');}
		if(left!=null){left.dumpKmersAsText(sb, k, mincount);}
		if(right!=null){right.dumpKmersAsText(sb, k, mincount);}
		return sb;
	}
	
	@Override
	protected final ByteBuilder dumpKmersAsText(ByteBuilder bb, int k, int mincount){
		if(value<1){return bb;}
		if(bb==null){bb=new ByteBuilder(32);}
		if(value>=mincount){bb.append(AbstractKmerTable.toBytes(pivot, value, k)).append('\n');}
		if(left!=null){left.dumpKmersAsText(bb, k, mincount);}
		if(right!=null){right.dumpKmersAsText(bb, k, mincount);}
		return bb;
	}
	
	final boolean TWOD(){return false;}
	
	/*--------------------------------------------------------------*/
	/*----------------       Invalid Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	int value;
	
}

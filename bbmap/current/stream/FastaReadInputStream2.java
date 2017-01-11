package stream;

import java.util.ArrayList;
import java.util.Arrays;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;

public class FastaReadInputStream2 extends ReadInputStream {
	
	public static void main(String[] args){
		
		FastaReadInputStream2 fris=new FastaReadInputStream2(args[0], true);
		
		Read r=fris.next();
		int i=0;
		while(r!=null){
			System.out.println(r.toText(false));
			r=fris.next();
			if(i++>3){break;}
		}
		
	}
	
	public FastaReadInputStream2(String fname, boolean allowSubprocess_){
		
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, allowSubprocess_, false);
		
		if(!ff.fasta() && !ff.stdio()){
			System.err.println("Warning: Did not find expected fasta file extension for filename "+fname);
		}

		tf=ByteFile.makeByteFile(ff, false);
		interleaved=false;

	}

	@Override
	public void start() {
//		if(cris!=null){cris.start();}
	}
	
	
	@Override
	public boolean hasMore() {
		if(buffer==null || next>=buffer.size()){
			if(tf.isOpen()){
				fillBuffer();
			}else{
				assert(generated>0) : "Was the file empty?";
			}
		}
		return (buffer!=null && next<buffer.size());
	}

	@Override
	public Read next() {
		if(!hasMore()){
			if(verbose){System.err.println("hasMore() returned false;  buffer="+(buffer==null ? null : buffer.size())+", next="+next+", consumed="+consumed);}
			return null;
		}
		Read r=buffer.get(next);
		buffer.set(next, null);
		next++;
		consumed++;
		return r;
	}
	
	@Override
	public synchronized ArrayList<Read> nextList() {
		if(next!=0){throw new RuntimeException("'next' should not be used when doing blockwise access.");}
		if(buffer==null || next>=buffer.size()){fillBuffer();}
		ArrayList<Read> r=buffer;
		buffer=null;
		if(r!=null && r.size()==0){r=null;}
		consumed+=(r==null ? 0 : r.size());
//		System.err.println(hashCode()+" produced "+r[0].numericID);
		return r;
	}
	
	private synchronized void fillBuffer(){
		if(verbose){System.err.println("Filling buffer.  buffer="+(buffer==null ? null : buffer.size()));}
		assert(buffer==null || next>=buffer.size());
		
		buffer=null;
		next=0;
		
		buffer=toReadList(tf, BUF_LEN, nextReadID, interleaved, headerA);

		if(verbose){System.err.println("Filled buffer.  buffer="+(buffer==null ? null : buffer.size()));}
		
		nextReadID+=buffer.size();
		if(buffer.size()<BUF_LEN){
			if(verbose){System.err.println("Closing tf");}
			tf.close();
		}
		
		generated+=buffer.size();
		if(verbose){System.err.println("generated="+generated);}
	}
	
	private ArrayList<Read> toReadList(ByteFile tf, int maxReadsToReturn, long numericID, boolean interleaved, String[] headerA){
		if(finished){return null;}
		if(verbose){System.err.println("FastaRIS fetching a list.");}
		
		ArrayList<Read> list=new ArrayList<Read>(Data.min(16384, maxReadsToReturn));
		
		int added=0;
		
		Read prev=null;
		if(verbose){System.err.println("added="+added+", max="+maxReadsToReturn);}
		while(added<maxReadsToReturn){
			Read r=makeNextRead(tf, maxReadsToReturn, numericID, headerA);
			if(verbose){System.err.println("Made "+r);}
			if(r==null){
				if(verbose){System.err.println("makeNextRead returned null.");}
				break;
			}
			if(interleaved){
				if(prev==null){prev=r;}
				else{
					prev.mate=r;
					r.mate=prev;
					list.add(prev);
					added++;
					numericID++;
					prev=null;
				}
			}else{
				list.add(r);
				added++;
				numericID++;
			}
		}
		
		assert(list.size()<=maxReadsToReturn);
		if(verbose){System.err.println("FastaRIS returning a list.  Size="+list.size());}
		return list;
	}
	
	private Read makeNextRead(ByteFile tf, int maxReadsToReturn, long numericID, String[] headerA){
		if(finished){
			if(verbose){System.err.println("Returning null because finished.");}
			return null;
		}
		assert(currentLine==null || currentLine[0]!=carrot);
		
		while(currentLine==null){
			currentLine=tf.nextLine();
			if(currentLine==null){
				if(verbose){System.err.println("Returning null because tf.nextLine()==null: A");}
				return null;
			}
			if(currentLine[0]==carrot){
				headerA[0]=new String(currentLine);
				currentLoc=0;
				currentSection=0;
				currentLine=null;
			}
		}
		if(verbose){System.err.println("currentLine="+new String(currentLine));}
		
		assert(currentLine==null || currentLine[0]!=carrot);
		
		StringBuilder sb=new StringBuilder();
		Read r=null;
		while(r==null){
			if(verbose){System.err.println("r==null, looping; current="+new String(currentLine)+"\nsb="+sb);}
			if(!SPLIT_READS || (currentLoc==0 && (currentLine.length<=(TARGET_READ_LEN-sb.length())))){
//				sb.append(currentLine);
				for(byte b : currentLine){sb.append((char)b);}
				currentLoc=currentLine.length;
			}else{
				while(sb.length()<TARGET_READ_LEN && currentLoc<currentLine.length){
					sb.append((char)currentLine[currentLoc]);
					currentLoc++;
				}
			}
			assert(currentLine==null || currentLine[0]!=carrot);
			assert(sb.length()<=TARGET_READ_LEN) : sb.length()+"\n"+sb;
			
			if(sb.length()==TARGET_READ_LEN){
				assert(currentLine==null || currentLine[0]!=carrot);
				r=makeRead(sb, numericID);
				currentSection++;
				return r;
			}else{
				assert(currentLine==null || currentLine[0]!=carrot);
				assert(currentLoc>=currentLine.length) : currentLoc+", "+currentLine.length+", "+
						TARGET_READ_LEN+", "+sb.length()+"\n"+currentLine+"\n"+sb;
				currentLine=null;
				currentLoc=0;
				while(currentLine==null){
					currentLine=tf.nextLine();
					if(currentLine==null || currentLine[0]==carrot){
						if(sb.length()>=MIN_READ_LEN){
							if(verbose){System.err.println("Made read of length "+sb.length());}
							r=makeRead(sb, numericID);
						}else{
							if(verbose){System.err.println("Read was too short at length "+sb.length()+"\n"+sb);}
							sb.setLength(0);
						}
						if(verbose){System.err.println("headerA was "+headerA[0]);}
						headerA[0]=(currentLine==null ? null : new String(currentLine));
						currentLoc=0;
						currentSection=0;
//						assert(false) : "'"+new String(currentLine)+"', "+headerA[0];
						currentLine=null;
						if(r!=null){
							if(verbose){System.err.println("Returning read "+r);}
							return r;
						}
						if(headerA[0]==null){
							if(verbose){System.err.println("Returning null because tf.nextLine()==null: B");}
							return null;
						}
					}
					assert(currentLine==null || currentLine[0]!=carrot);
				}
				assert(currentLine==null || currentLine[0]!=carrot);
			}
			assert(currentLine==null || currentLine[0]!=carrot);
		}
		assert(currentLine==null || currentLine[0]!=carrot);
		if(verbose){System.err.println("Returning null because loop exited (should be unreachable).");}
		return null;
	}
	
	private Read makeRead(StringBuilder sb, long numericID){
		byte[] quals=null;
		byte[] bases=new byte[sb.length()];
		if(FAKE_QUALITY){
			quals=new byte[sb.length()];
			Arrays.fill(quals, (byte)(30));
		}
		for(int i=0; i<bases.length; i++){
			bases[i]=(byte)Character.toUpperCase(sb.charAt(i));
		}
		assert(bases[0]!=carrot) : new String(bases)+"\n"+numericID+"\n"+headerA[0];
		String hd=(currentSection>0 ? headerA[0].substring(1)+"_"+currentSection : new String(headerA[0].substring(1)));
//		assert(currentSection==0);
		Read r=new Read(bases, (byte)0, (byte)0, 0, 0, hd, quals, numericID);
		return r;
	}
	
	public boolean close(){
		return tf.close();
	}

	@Override
	public synchronized void restart() {
		generated=0;
		consumed=0;
		next=0;
		nextReadID=0;
		buffer=null;
		
		currentLine=null;
		currentLoc=0;
		currentSection=0;
		finished=false;
		
		tf.reset();
	}

	@Override
	public boolean paired() {
		return interleaved;
	}

	private ArrayList<Read> buffer=null;
	private int next=0;
	
	private final ByteFile tf;
	private final boolean interleaved;

	private final int BUF_LEN=Shared.READ_BUFFER_LENGTH;

	public long generated=0;
	public long consumed=0;
	private long nextReadID=0;
	
	private final String[] headerA=new String[1];
	
	public static boolean SPLIT_READS=true;
	public static int TARGET_READ_LEN=500;
	public static int MIN_READ_LEN=40;
	public static int DEFAULT_WRAP=100;
	
	public static boolean verbose=false;
	public static boolean FAKE_QUALITY=false;
	
	private byte[] currentLine=null;
	private int currentLoc=0;
	private int currentSection=0;
	private boolean finished=false;
	private final byte carrot='>';
	
}

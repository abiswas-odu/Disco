package fileIO;
import java.io.File;
import java.io.InputStream;
import java.util.ArrayList;

import dna.Data;
import shared.Shared;


public abstract class ByteFile {
	
//	public static final ByteFile makeByteFile(String fname){
//		return makeByteFile(fname, false, true); 
//	}
	
	public static final ByteFile makeByteFile(String fname, boolean tryAllExtensions, boolean allowSubprocess){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess, false);
		return makeByteFile(ff, tryAllExtensions);
	}
	
	public static final ByteFile makeByteFile(FileFormat ff, boolean tryAllExtensions){
		if(!Shared.LOW_MEMORY && (FORCE_MODE_BF2 || (!FORCE_MODE_BF1 && Shared.threads()>4/* && (ReadWrite.isCompressed(fname) || ReadWrite.isSam(fname))*/))){
//			if(allowSubprocess && ((ReadWrite.USE_UNPIGZ || ReadWrite.USE_GUNZIP) && (fname.endsWith(".gz") || fname.endsWith(".gzip")))){}
			return new ByteFile2(ff, tryAllExtensions);
		}
		return new ByteFile1(ff, tryAllExtensions); 
	}
	
//	protected ByteFile(String fname, boolean tryAllExtensions, boolean allowSubprocess_){
//		allowSubprocess=allowSubprocess_;
//		fname=fname.replace('\\', '/');
//		File f=new File(fname);
//
//		if(tryAllExtensions && !fname.startsWith("jar:") && !f.exists()){
//			name=ReadWrite.findFileExtension(fname);
//			f=new File(name);
//		}else{
//			name=fname;
//		}
//	}
	
	protected ByteFile(FileFormat ff_, boolean tryAllExtensions){
		ff=ff_;
		assert(ff.read()) : ff;
	}
	
	public final ArrayList<byte[]> toByteLines(){
		
		byte[] s=null;
		ArrayList<byte[]> list=new ArrayList<byte[]>(4096);
		
		for(s=nextLine(); s!=null; s=nextLine()){
			list.add(s);
		}
		
		return list;
	}
	
	public final long countLines(){
		byte[] s=null;
		long count=0;
		for(s=nextLine(); s!=null; s=nextLine()){count++;}
		reset();
		
		return count;
	}
	
	public abstract void reset();
	
	public final boolean exists(){
		return name().equals("stdin") || name().startsWith("stdin.") || name().startsWith("jar:") || new File(name()).exists(); //TODO Ugly and unsafe hack for files in jars
	}

	public abstract InputStream is();
	public abstract long lineNum();
	
	/** Returns true if there was an error */
	public abstract boolean close();
	
	public abstract byte[] nextLine();
	
	public final void pushBack(byte[] line){
		assert(pushBack==null);
		pushBack=line;
	}
	
	public abstract boolean isOpen();
	
	public final String name(){return ff.name();}
	public final boolean allowSubprocess(){return ff.allowSubprocess();}
	
	public final FileFormat ff;

	public static boolean FORCE_MODE_BF1=!(Data.GENEPOOL || Data.WINDOWS);
	public static boolean FORCE_MODE_BF2=false;
	
	protected final static byte slashr='\r', slashn='\n', carrot='>', plus='+', at='@';//, tab='\t';
	
	byte[] pushBack=null;
	
}

package jgi;

import java.io.PrintStream;
import java.util.Arrays;

import stream.SamLine;

import dna.Gene;
import fileIO.FileFormat;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Timer;

/**
 * @author Brian Bushnell
 * @date Jul 23, 2013
 *
 */
public class SplitSam4Way {
	
	public static void main(String[] args){
		new SplitSam4Way(args);
	}
	
	private void printOptions(){
		outstream.println("Syntax:\n");
		outstream.println("java -ea -Xmx128m -cp <path> jgi.SplitSam4Way <input> <out plus> <out minus> <out chimeric> <out unmapped>");
		outstream.println("If you do not want one of the output files, use the word 'null'.\n");
	}
	
	public SplitSam4Way(String[] args){
		if(args==null || args.length!=5){
			printOptions();
			System.exit(0);
		}

		for(String s : args){if(s.startsWith("out=standardout") || s.startsWith("out=stdout")){outstream=System.err;}}
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		Timer t=new Timer();
		long reads=0, bases=0;
		long preads=0, mreads=0, creads=0, ureads=0;
		
		String fin=args[0];
		String fplus=args[1];
		String fminus=args[2];
		String fchimeric=args[3];
		String funmapped=args[4];
		
		TextFile tf=new TextFile(fin, true, false);
		TextStreamWriter plus=("null".equalsIgnoreCase(fplus) ? null : new TextStreamWriter(fplus, true, false, true, FileFormat.SAM));
		TextStreamWriter minus=("null".equalsIgnoreCase(fminus) ? null : new TextStreamWriter(fminus, true, false, true, FileFormat.SAM));
		TextStreamWriter chimeric=("null".equalsIgnoreCase(fchimeric) ? null : new TextStreamWriter(fchimeric, true, false, true, FileFormat.SAM));
		TextStreamWriter unmapped=("null".equalsIgnoreCase(funmapped) ? null : new TextStreamWriter(funmapped, true, false, true, FileFormat.SAM));

		if(plus!=null){plus.start();}
		if(minus!=null){minus.start();}
		if(chimeric!=null){chimeric.start();}
		if(unmapped!=null){unmapped.start();}
		
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.charAt(0)=='@'){
				if(plus!=null){plus.println(line);}
				if(minus!=null){minus.println(line);}
				if(chimeric!=null){chimeric.println(line);}
				if(unmapped!=null){unmapped.println(line);}
			}else{
				SamLine sl=new SamLine(line);
				reads++;
//				bases+=sl.seq.length();
				bases+=sl.seq.length;
				
				if(!sl.mapped() || !sl.nextMapped() || !sl.hasMate() || !sl.primary()){
					if(unmapped!=null){unmapped.println(line);}
					ureads++;
//					System.out.println("unmapped: "+sl.mapped()+", "+sl.nextMapped()+", "+sl.hasMate()+", "+!sl.primary());
				}else if(!sl.pairedOnSameChrom() || sl.strand()==sl.nextStrand()){
					if(chimeric!=null){chimeric.println(line);}
					creads++;
//					System.out.println("chimeric: "+sl.pairedOnSameChrom()+", "+(sl.strand()==sl.nextStrand())+", "+sl.strand()+", "+sl.nextStrand()+", "+new String(sl.rname())+", "+new String(sl.rnext()));
				}else if((sl.firstFragment() ? sl.strand() : sl.nextStrand())==Gene.PLUS){
					if(plus!=null){plus.println(line);}
					preads++;
				}else if((sl.firstFragment() ? sl.strand() : sl.nextStrand())==Gene.MINUS){
					if(minus!=null){minus.println(line);}
					mreads++;
				}else{
					throw new RuntimeException("Unhandled case: "+sl.firstFragment()+", "+sl.lastFragment()+", "+sl.strand()+", "+sl.nextStrand()+"\n"+sl+"\n");
				}
			}
		}
		
		if(plus!=null){plus.poisonAndWait();}
		if(minus!=null){minus.poisonAndWait();}
		if(chimeric!=null){chimeric.poisonAndWait();}
		if(unmapped!=null){unmapped.poisonAndWait();}
		t.stop();
		
		
		double rpnano=reads/(double)(t.elapsed);
		double bpnano=bases/(double)(t.elapsed);

		String rpstring=(reads<100000 ? ""+reads : reads<100000000 ? (reads/1000)+"k" : (reads/1000000)+"m");
		String bpstring=(bases<100000 ? ""+bases : bases<100000000 ? (bases/1000)+"k" : (bases/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}

		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		outstream.println("Plus Reads:         "+preads);
		outstream.println("Minus Reads:        "+mreads);
		outstream.println("Chimeric Reads:     "+creads);
		outstream.println("Unmapped Reads:     "+ureads);
		
		
	}
	
	private PrintStream outstream=System.err;
	
}

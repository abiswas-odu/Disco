package jgi;

import java.util.Random;

import dna.AminoAcid;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;

/**
 * @author Brian Bushnell
 * @date Jan 3, 2013
 *
 */
public class RandomGenome {
	
	public static void main(String[] args){
		ReadWrite.ZIPLEVEL=2;
		Random randy=new Random();
		int chroms=Integer.parseInt(args[0]);
		int len=Integer.parseInt(args[1]);
		
		String fname=args[2];
		TextStreamWriter tsw=new TextStreamWriter(fname, false, false, true);
		tsw.start();
		
		for(int chrom=1; chrom<=chroms; chrom++){
			tsw.println(">"+chrom);
			StringBuilder sb=new StringBuilder(101);
			for(int i=0, j=0; i<len; i++, j++){
				char c;
				if((i/10000)%4==3){
					c='N';
				}else{
					c=(char)AminoAcid.numberToBase[randy.nextInt(4)];
				}
				sb.append(c);
				if(j==100){
					sb.append('\n');
					tsw.print(sb);
					sb=new StringBuilder(101);
					j=0;
				}
			}
			if(sb.length()>0){
				sb.append('\n');
				tsw.print(sb);
			}
		}
		tsw.poison();
		tsw.waitForFinish();
		
	}
	
}

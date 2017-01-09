package align2;

import dna.Data;
import dna.Gene;
import fileIO.TextFile;

public class Evaluate {
	
	public static void main(String[] args){
		TextFile tf=new TextFile(args[0], false, false);
		
		String[] lines=tf.toStringLines();
		tf.close();
		
		int trials=lines.length;
		if(args.length>0){
			if(args.length>1){
				trials=Integer.parseInt(args[1]);
			}
		}
		
		int correct=0;
		
		for(String s : lines){
			if(isCorrect(s)){correct++;}
		}
		
		int incorrect=trials-correct;
		int falsePositive=lines.length-correct;
		int falseNegative=trials-lines.length;

		Data.sysout.println("Trials:        \t"+trials);
		Data.sysout.println("Correct:       \t"+correct+"\t"+String.format("%.4f",correct*100f/trials)+"%");
		Data.sysout.println("Incorrect:     \t"+incorrect+"\t"+String.format("%.4f",incorrect*100f/trials)+"%");
		Data.sysout.println("False Positive:\t"+falsePositive+"\t"+String.format("%.4f",falsePositive*100f/trials)+"%");
		Data.sysout.println("False Negative:\t"+falseNegative+"\t"+String.format("%.4f",falseNegative*100f/trials)+"%");
		
	}
	
	public static boolean isCorrect(String s){
		String[] line=s.split("\t");
		
		String[] answer=line[0].split("_");
		int trueChrom=Gene.toChromosome(answer[1]);
		byte trueStrand=Byte.parseByte(answer[2]);
		int trueLoc=Integer.parseInt(answer[3]);

		int calledChrom=Gene.toChromosome(line[2]);
		byte calledStrand=Gene.toStrand(line[1]);
		int calledLoc=Integer.parseInt(line[3]);
		
		return (trueChrom==calledChrom && trueStrand==calledStrand && Math.abs(trueLoc-calledLoc)<1000);
		
	}
	
}

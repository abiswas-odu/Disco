package dna;
import java.util.HashSet;

import var.VarLine;


public class Coverage{

	public Coverage(Gene gg){
		g=gg;
	}

	public final Gene g;
	public HashSet<VarLine> varSet; //TODO:  Could change these to arrays and sort them.
	public int min=Integer.MAX_VALUE;
	public int max=0;
	public int covered=0;
	public int uncovered=0;
	public long sum=0;
	public float avg;
	public float covRatio;

	public int[] missingChromRelative;
	public int[] missingGeneRelative;

}
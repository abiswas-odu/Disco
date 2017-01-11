package var2;

import shared.Tools;

public class VarFilter {
	

	public boolean parse(String a, String b, String arg){
		if(a.equals("minreads")){
			minReads=Integer.parseInt(b);
		}else if(a.equals("minqualitymax") || a.equals("minmaxquality")){
			minMaxQuality=Integer.parseInt(b);
		}else if(a.equals("minedistmax") || a.equals("minmaxedist")){
			minMaxEdist=Integer.parseInt(b);
		}else if(a.equals("minmapqmax") || a.equals("minmaxmapq")){
			minMaxMapq=Integer.parseInt(b);
		}else if(a.equals("minidmax") || a.equals("minmaxid")){
			minMaxIdentity=Float.parseFloat(b);
			if(minMaxIdentity>1){minMaxIdentity/=100;}
		}

		else if(a.equals("minpairingrate") || a.equals("minpairrate")){
			minPairingRate=Float.parseFloat(b);
		}else if(a.equals("minstrandratio")){
			minStrandRatio=Float.parseFloat(b);
		}else if(a.equals("minscore")){
			minScore=Float.parseFloat(b);
		}else if(a.equals("minquality") || a.equals("minavgquality") || a.equals("maq")){
			minAvgQuality=Float.parseFloat(b);
		}else if(a.equals("minedist") || a.equals("minavgedist") || a.equals("mae")){
			minAvgEdist=Float.parseFloat(b);
		}else if(a.equals("minmapq") || a.equals("minavgmapq")){
			minAvgMapq=Float.parseFloat(b);
		}else if(a.equals("minallelefraction") || a.equals("minallelefrequency") || a.equals("maf")){
			minAlleleFraction=Float.parseFloat(b);
		}else if(a.equals("minidentity") || a.equals("mid") || a.equals("minid")){
			minIdentity=Float.parseFloat(b);
			if(minIdentity>1){minIdentity/=100;}
		}else if(a.equals("lowcoveragepenalty") || a.equals("lowcovpenalty") || a.equals("covpenalty")){
			Var.lowCoveragePenalty=Float.parseFloat(b);
			assert(Var.lowCoveragePenalty>=0) : "Low coverage penalty must be at least 0.";
		}
		
		else if(a.equals("rarity")){
			rarity=Float.parseFloat(b);
			assert(rarity>=0 && rarity<=1);
			minAlleleFraction=Tools.min(minAlleleFraction, rarity);
		}
		
		else if(a.equals("clearfilters")){
			if(Tools.parseBoolean(b)){clear();}
		}else{
			return false;
		}
		return true;
	}
	
	public void clear(){
		minReads=0;
		minMaxQuality=0;
		minMaxEdist=0;
		minMaxMapq=0;
		minMaxIdentity=0;

		minPairingRate=0;
		minStrandRatio=0;
		minScore=0;
		minAvgQuality=0;
		minAvgEdist=0;
		minAvgMapq=0;
		minAlleleFraction=0;
		minIdentity=0;
	}
	

	public void setFrom(VarFilter filter) {
		minReads=filter.minReads;
		minMaxQuality=filter.minMaxQuality;
		minMaxEdist=filter.minMaxEdist;
		minMaxMapq=filter.minMaxMapq;
		minMaxIdentity=filter.minMaxIdentity;

		minPairingRate=filter.minPairingRate;
		minStrandRatio=filter.minStrandRatio;
		minScore=filter.minScore;
		minAvgQuality=filter.minAvgQuality;
		minAvgEdist=filter.minAvgEdist;
		minAvgMapq=filter.minAvgMapq;
		minAlleleFraction=filter.minAlleleFraction;
		minIdentity=filter.minIdentity;
	}
	
	public boolean passesFast(Var v){
		if(v.count()<minReads){return false;}
		if(v.baseQMax<minMaxQuality){return false;}
		if(v.endDistMax<minMaxEdist){return false;}
		if(v.mapQMax<minMaxMapq){return false;}
		return true;
	}
	
	public boolean passesFilter(Var v, float pairingRate, float totalQualityAvg, float totalMapqAvg, int ploidy, ScafMap map){
		final int count=v.count();
		if(count<minReads){return false;}
		if(v.baseQMax<minMaxQuality){return false;}
		if(v.endDistMax<minMaxEdist){return false;}
		if(v.mapQMax<minMaxMapq){return false;}
		if(v.idMax*0.001f<minMaxIdentity){return false;}

		//Slower, uses division.
//		if(pairingRate>0 && minPairingRate>0 && v.pairingRate()<minPairingRate){return false;}
//		if(minStrandRatio>0 && v.strandRatio()<minStrandRatio){return false;}
//		if(minAvgQuality>0 && v.baseQAvg()<minAvgQuality){return false;}
//		if(minAvgEdist>0 && v.edistAvg()<minAvgEdist){return false;}
//		if(minAvgMapq>0 && v.mapQAvg()<minAvgMapq){return false;}

		if(pairingRate>0 && minPairingRate>0 && count*minPairingRate>v.properPairCount){return false;}
		if(minAvgQuality>0 && count*minAvgQuality>v.baseQSum){return false;}
		if(minAvgEdist>0 && count*minAvgEdist>v.endDistSum){return false;}
		if(minAvgMapq>0 && count*minAvgMapq>v.mapQSum){return false;}
		if(minIdentity>0 && count*minIdentity*1000>v.idSum){return false;}
		
		if(minStrandRatio>0 && v.strandRatio()<minStrandRatio){return false;}
		
		//assert(v.coverage()>0);
		if(v.coverage()>0 && minAlleleFraction>0 && v.alleleFraction()<minAlleleFraction){return false;}
		
		if(minScore>0 && v.phredScore(pairingRate, totalQualityAvg, totalMapqAvg, rarity, ploidy, map)<minScore){return false;}
		
		return true;
	}
	
	public String toString(float pairingRate, int ploidy){
		StringBuilder sb=new StringBuilder();
		
		sb.append("pairingRate=").append(pairingRate).append("\n");
		sb.append("ploidy=").append(ploidy).append("\n");
		
		sb.append("minReads=").append(minReads).append("\n");
		sb.append("minMaxQuality=").append(minMaxQuality).append("\n");
		sb.append("minMaxEdist=").append(minMaxEdist).append("\n");
		sb.append("minMaxMapq=").append(minMaxMapq).append("\n");
		sb.append("minMaxIdentity=").append(minMaxIdentity).append("\n");

		sb.append("minPairingRate=").append(minPairingRate).append("\n");
		sb.append("minStrandRatio=").append(minStrandRatio).append("\n");
		sb.append("minScore=").append(minScore).append("\n");
		sb.append("minAvgQuality=").append(minAvgQuality).append("\n");
		sb.append("minAvgEdist=").append(minAvgEdist).append("\n");
		sb.append("minAvgMapq=").append(minAvgMapq).append("\n");
		sb.append("minAlleleFraction=").append(minAlleleFraction);
		sb.append("minIdentity=").append(minIdentity);
		
		return sb.toString();
	}
	
	public int minReads=2;
	public int minMaxQuality=15;
	public int minMaxEdist=20;
	public int minMaxMapq=0;
	public float minMaxIdentity=0f;
	
	public float minPairingRate=0.1f;
	public float minStrandRatio=0.1f;
	public float minScore=20;
	public float minAvgQuality=12;
	public float minAvgEdist=10;
	public float minAvgMapq=0;
	public float minAlleleFraction=0.1f;
	public float minIdentity=0f;
	public float rarity=1f;
	
}

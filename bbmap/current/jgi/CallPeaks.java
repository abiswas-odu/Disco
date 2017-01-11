package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import stream.ByteBuilder;
import structures.LongList;
import dna.Parser;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import shared.ReadStats;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Oct 28, 2014
 *
 */
public class CallPeaks {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		Timer t=new Timer();
		CallPeaks cp=new CallPeaks(args);
		cp.process(t);
	}
	
	public CallPeaks(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		if(printClass){outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");}
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens

			if(Parser.isJavaFlag(arg)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Tools.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("in")){
				in=b;
			}else if(a.equals("out")){
				out=b;
			}else if(a.equals("minheight") || a.equals("h")){
				minHeight=Long.parseLong(b);
			}else if(a.equals("minvolume") || a.equals("v")){
				minVolume=Long.parseLong(b);
			}else if(a.equals("minwidth") || a.equals("w")){
				minWidth=Integer.parseInt(b);
			}else if(a.equals("minpeak") || a.equals("minp")){
				minPeak=Integer.parseInt(b);
			}else if(a.equals("maxpeak") || a.equals("maxp")){
				maxPeak=Integer.parseInt(b);
			}else if(a.equals("maxpeakcount") || a.equals("maxpc") || a.equals("maxpeaks")){
				maxPeakCount=Integer.parseInt(b);
			}else if(a.equals("smoothradius")){
				smoothRadius=Integer.parseInt(b);
			}else if(a.equals("smoothprogressive")){
				smoothProgressiveFlag=Tools.parseBoolean(b);
			}else if(a.equals("maxradius")){
				maxRadius=Integer.parseInt(b);
			}else if(a.equals("progressivemult")){
				progressiveMult=Float.parseFloat(b);
			}else if(a.equals("ploidy")){
				ploidyClaimed=Integer.parseInt(b);
			}else if(a.equals("column") || a.equals("col") || a.equals("countcolumn")){
				countColumn=Integer.parseInt(b);
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(in==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				in=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		
		if(out==null){out="stdout.txt";}
		
		ffout=FileFormat.testOutput(out, FileFormat.TEXT, null, true, overwrite, append, false);
		ffin=FileFormat.testInput(in, FileFormat.TEXT, null, true, false);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void process(Timer t){
		LongList hist=loadHistogram(ffin);
		ArrayList<Peak> peaks=callPeaks(hist);
		long sum=Tools.sum(hist.array);
		hist=null;
		printPeaks(peaks, k, sum);
		t.stop();
		System.err.println("\nFound "+peaks.size()+" peaks in "+t);
		
		peaks=null;
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	public static boolean printPeaks(long[] array, String fname, boolean ow, long minHeight, long minVolume, int minWidth, int minPeak, int maxPeak, int maxPeakCount,
			int k, int ploidy, ArrayList<String> list){
		if(list==null){list=new ArrayList<String>();}
		list.add("out="+fname);
		list.add("ow="+ow);
		list.add("minheight="+minHeight);
		list.add("minvolume="+minVolume);
		list.add("minwidth="+minWidth);
		list.add("minpeak="+minPeak);
		list.add("maxpeak="+maxPeak);
		list.add("maxpeaks="+maxPeakCount);
		list.add("k="+(k<1 ? 31 : k));
		if(ploidy>0){list.add("ploidy="+ploidy);}
		CallPeaks cp=new CallPeaks(list.toArray(new String[0]));
		ArrayList<Peak> peaks=cp.callPeaks(array, array.length);
		cp.printPeaks(peaks, k, Tools.sum(array));
		return cp.errorState;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static LongList loadHistogram(FileFormat ff){
		LongList list=new LongList(8000);
		TextFile tf=new TextFile(ff);
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.startsWith("#")){
				//ignore
			}else{
				String[] split=line.split("\\s+");
				if(split.length==1){
					list.add(Long.parseLong(split[0]));
				}else{
					list.set(Integer.parseInt(split[0]), Long.parseLong(split[countColumn]));
				}
			}
		}
		boolean errorState_=tf.close();
		assert(!errorState_) : "Encountered an error when reading "+ff.name()+".\n" +
				"To skip this error message, run with the '-da' flag.";
		
		return list;
	}
	
	private static ArrayList<Peak> condense(ArrayList<Peak> in, int maxCount){
		if(in==null || in.isEmpty()){return in;}
		maxCount=Tools.max(Tools.min(in.size(), maxCount), 1);
		ArrayList<Peak> out=new ArrayList<Peak>(Tools.min(maxCount, in.size()));
		
		final long hlimit, vlimit;
		
		{
			long[] heights=new long[in.size()];
			for(int i=0; i<in.size(); i++){
				Peak p=in.get(i);
				heights[i]=(callByRawCount ? p.centerHeight2() : p.centerHeight);
			}
			Arrays.sort(heights);
			hlimit=heights[heights.length-maxCount];
		}
		
		{
			int mc2=(maxCount+1)/2;
			long[] volumes=new long[in.size()];
			for(int i=0; i<in.size(); i++){
				Peak p=in.get(i);
				volumes[i]=(callByRawCount ? p.volume2 : p.volume);
			}
			Arrays.sort(volumes);
			vlimit=volumes[volumes.length-mc2];
		}
		
		for(Peak p : in){
			final long height=(callByRawCount ? p.centerHeight2() : p.centerHeight);
			final long volume=(callByRawCount ? p.volume2 : p.volume);
			if(volume>=vlimit || height>=hlimit){out.add(p);}
		}
		
		for(Peak p : in){
			final long height=(callByRawCount ? p.centerHeight2() : p.centerHeight);
			final long volume=(callByRawCount ? p.volume2 : p.volume);
			if(volume<vlimit && height<hlimit){
				Peak p2=out.get(0);
				for(Peak temp : out){
					if(Tools.absdif(p.center, temp.center)<Tools.absdif(p.center, p2.center)){
						p2=temp;
					}
				}
				if(p2.compatibleWith(p)){
					p2.absorb(p);
				}//else discard
			}
		}
		
		return out;
	}
	
	private void capWidth(ArrayList<Peak> peaks, float maxWidthMult, long[] counts){
		float mult=1/maxWidthMult;
		for(Peak p : peaks){
			p.start=(int)Math.round(Tools.max(p.start, p.center*mult));
			p.stop=(int)Math.round(Tools.min(p.stop, p.center*maxWidthMult));
			p.recalculate(counts);
		}
		
//		for(int i=0; i<peaks.)
	}
	
	private void printPeaks(ArrayList<Peak> peaks, int k, long uniqueKmers){
		if(ffout==null){return;}
		ByteStreamWriter bsw=new ByteStreamWriter(ffout);
		bsw.start();
		
		if(peaks.size()>0){
			try {
				final Peak p0=peaks.get(0);
				final int center0=p0.center;
				final int ploidyEstimate=calcPloidy(peaks);
				final int ploidy=ploidyClaimed>0 ? ploidyClaimed : ploidyEstimate;
				final long genomeSize=genomeSize(peaks);
				final long repeatSize=repeatSize(peaks, ploidy);
				final long haploidSize=genomeSize/ploidy;
				final long hetLocs=calcHetLocations(peaks, ploidy, k);
				final double hetRate=hetLocs/(double)haploidSize;
				final double repeatRate=repeatSize*1.0/genomeSize;
				
				Peak ploidyPeak=p0, mainPeak=p0;
				int target=center0*ploidy, haploidCov;
				for(Peak p : peaks){
					if(p.volume>mainPeak.volume){
						mainPeak=p;
					}
					if(Tools.absdif(p.center, target)<Tools.absdif(ploidyPeak.center, target)){
						ploidyPeak=p;
					}
				}
				if(Tools.max(target,ploidyPeak.center)/(float)Tools.min(target,ploidyPeak.center)<1.3f){
					haploidCov=ploidyPeak.center;
				}else{
					haploidCov=target;
				}
				
//				System.err.println("ploidyEstimate="+ploidyEstimate);
//				System.err.println("genomeSize="+genomeSize);
//				System.err.println("repeatSize="+repeatSize);
//				System.err.println("haploidSize="+haploidSize);
//				System.err.println("hetLocs="+hetLocs);
//				System.err.println("biggestPeak="+biggestPeak(peaks));
//				System.err.println("homozygousPeak("+ploidy+")="+homozygousPeak(peaks, ploidy));
//				System.err.println("homozygousPeak("+ploidyEstimate+")="+homozygousPeak(peaks, ploidyEstimate));
				
				if(ploidy!=ploidyEstimate){
					System.err.println("Warning - ploidy detected at "+ploidyEstimate+" differs from stated ploidy of "+ploidyClaimed);
				}

				if(k>0){bsw.println("#k\t"+k);}
				bsw.println("#unique_kmers\t"+uniqueKmers);
				bsw.println("#main_peak\t"+mainPeak.center);
				bsw.println("#genome_size\t"+genomeSize);
				if(ploidy>1 || true){bsw.println("#haploid_genome_size\t"+(genomeSize/ploidy));}
				bsw.println("#fold_coverage\t"+center0);
				if(ploidy>1 || true){bsw.println("#haploid_fold_coverage\t"+haploidCov);}
				bsw.println("#ploidy\t"+ploidy);
				if(ploidy!=ploidyEstimate){bsw.println("#ploidy_detected\t"+ploidyEstimate);}
				if(ploidy>1){bsw.println("#het_rate\t"+String.format("%.5f", hetRate));}
				bsw.println("#percent_repeat\t"+String.format("%.3f", (100*repeatRate)));
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		bsw.println("#start\tcenter\tstop\tmax\tvolume");
		ByteBuilder bb=new ByteBuilder(200);
		for(Peak p : peaks){
			p.toBytes(bb);
			bsw.println(bb);
			bb.setLength(0);
		}
		errorState|=bsw.poisonAndWait();
	}
	
	private long genomeSize(ArrayList<Peak> peaks){
		if(peaks.size()<1){return 0;}
		
		long sizeSum=0;
		final Peak p0=peaks.get(0);
		final int center0=p0.center;
		final double mult=1.0/(Tools.max(1, center0));
		for(Peak p : peaks){
			long size=p.volume*(Math.round(p.center*mult));
			sizeSum+=size;
		}
		return sizeSum;
	}
	
	private long repeatSize(ArrayList<Peak> peaks, int ploidy){
		if(peaks.size()<2){return 0;}
		assert(ploidy>0) : ploidy;
		final int homozygousLoc=homozygousPeak(peaks, ploidy);
		final Peak p0=peaks.get(0);
		final int center0=p0.center;
		final double mult=1.0/(Tools.max(1, center0));
		
		long sizeSum=0;
		for(int i=homozygousLoc+1; i<peaks.size(); i++){
			Peak p=peaks.get(i);
			long size=p.volume*(Math.round(p.center*mult));
			sizeSum+=size;
		}
		return sizeSum;
	}
	
	private int biggestPeak(ArrayList<Peak> peaks){
		if(peaks.size()<2){return peaks.size()-1;}

		final Peak p0=peaks.get(0);
		Peak biggest=p0;
		int loc=0;
		for(int i=1; i<peaks.size(); i++){
			Peak p=peaks.get(i);
			if(p.volume>biggest.volume){
				loc=i;
				biggest=p;
			}
		}
		return loc;
	}
	
	private int secondBiggestPeak(ArrayList<Peak> peaks){
		if(peaks.size()<2){return peaks.size()-1;}
		
		Peak biggest=peaks.get(0);
		Peak second=peaks.get(1);
		int bloc=0;
		int sloc=1;
		if(second.volume>biggest.volume){
			Peak temp=second;
			second=biggest;
			biggest=temp;
			bloc=1;
			sloc=0;
		}
		for(int i=2; i<peaks.size(); i++){
			Peak p=peaks.get(i);
			if(p.volume>second.volume){
				sloc=i;
				second=p;
				if(second.volume>biggest.volume){
					Peak temp=second;
					second=biggest;
					biggest=temp;
					sloc=bloc;
					bloc=i;
				}
			}
		}
		return sloc;
	}
	
	private int homozygousPeak(ArrayList<Peak> peaks, final int ploidy){
		if(peaks.size()<2){return peaks.size()-1;}
		assert(ploidy>0) : ploidy;

		final Peak p0=peaks.get(0);
		final int target=p0.center*ploidy;
//		System.err.println("target="+target);
		int bestDif=Integer.MAX_VALUE;
		int loc=0;
		for(int i=0; i<peaks.size(); i++){
			Peak p=peaks.get(i);
			int dif=Tools.absdif(target, p.center);
//			System.err.println("dif="+dif+" for peak "+i+", center "+p.center);
			if(dif<bestDif){
				bestDif=dif;
				loc=i;
//				System.err.println("New best at loc "+i);
			}
		}
		return loc;
	}
	
	private int calcPloidy(ArrayList<Peak> peaks){
		if(peaks.size()<2){return 1;}

		final Peak p0=peaks.get(0);
		final Peak biggest=peaks.get(biggestPeak(peaks));
		
		if(biggest!=p0){
			int ratio=(int)(biggest.center/(float)p0.center);
			return ratio;
		}else{//p0 is biggest.
			final Peak second=peaks.get(secondBiggestPeak(peaks));
			if(second.volume*4<biggest.volume){return 1;} //Probably second is a repeat peak.
			int ratio=(int)(second.center/(float)biggest.center);
			return ratio;
		}
	}
	
	private long calcHetLocations(ArrayList<Peak> peaks, final int ploidy, final int k){
		if(peaks.size()<2){return 0;}
		assert(ploidy>0) : ploidy;
		final int homozygousLoc=homozygousPeak(peaks, ploidy);
		final Peak homoPeak=peaks.get(homozygousLoc);
		long sum=0;
		final int lim=ploidy/2;
		for(int i=0; i<homozygousLoc; i++){
			final Peak p=peaks.get(i);
			final int copyCount=Math.round((p.center*ploidy)/(float)homoPeak.center);
//			System.err.println("lim="+lim+". For peak "+i+", copyCount="+copyCount+", volume="+p.volume);
			if(copyCount>lim){break;}
//			double peakLocs=(p.volume/(double)k);
			sum=sum+p.volume;
		}
		return sum/k;
	}
	
	public ArrayList<Peak> callPeaks(LongList list){
		return callPeaks(list.array, list.size);
	}
	
	public ArrayList<Peak> callPeaks(long[] original, int length){
		
		final long[] array;
		if(smoothRadius>0){
			if(smoothProgressiveFlag){
				array=smoothProgressive(original, smoothRadius);
			}else{
				array=smooth(original, smoothRadius);
			}
		}else{
			array=original;
		}
		
		ArrayList<Peak> peaks=new ArrayList<Peak>();
		
		int dip0=-1;
		for(int i=1; i<length; i++){
			if(array[i-1]<array[i]){
				dip0=i-1;
				break;
			}
		}
		if(dip0<0){return peaks;}
//		assert(false) : dip0+", "+array[dip0);
		
		final int UP=0, DOWN=1;
		int mode=UP;
		int start=dip0, center=-1;
		long prev=array[dip0];
		long sum=prev;
		long sum2=prev*dip0;
		for(int i=dip0+1; i<length; i++){
			final long x=array[i];
			
//			if(i<16){System.err.println("i="+i+", x="+x+", mode="+mode+", center="+center+", start="+start+", dip0="+dip0);}
			
			if(mode==UP){
				if(x<prev){
					mode=DOWN;
					center=i-1;
				}
			}else{
				if(x>prev){
					mode=UP;
					int stop=i-1;
					long max=array[center];
					if(center>=minPeak && center<=maxPeak && max>=minHeight && (stop-start)>=minWidth && sum>=minVolume){
						for(int j=center-1; j>=0; j--){//find middle of mesas
							if(array[j]!=max){
								center=(center+j+2)/2;
								break;
							}
						}
						{
							long valley=array[stop];
							for(int j=stop; j>=0; j--){//find middle of valleys
								if(array[j]!=valley){
									if(valley==0){stop=j+1;}
									else{stop=(stop+j+2)/2;}
									break;
								}
							}
						}
						
						Peak p=new Peak(center, start, stop, max, array[start], array[stop], sum, sum2);
						peaks.add(p);
					}else{
//						Peak p=new Peak(center, start, stop, max, sum);
//						System.err.println("*"+p);
					}
					start=stop;
					stop=-1;
					sum=sum2=0;
					center=-1;
					if(i>maxPeak){break;}
					while(i<array.length && array[i]==0){i++;}//Skip zero regions
				}
			}
			
			sum+=x;
			sum2+=(x*i);
			prev=x;
		}
		
		if(mode==DOWN){
			int stop=length;
			long max=array[center];
			for(int j=center-1; j>=0; j--){//find middle of mesas
				if(array[j]!=max){
					center=(center+j+2)/2;
					break;
				}
			}
			{
				long valley=array[stop-1];
				for(int j=stop-1; j>=0; j--){//find middle of valleys
					if(array[j]!=valley){
						if(valley==0){stop=j+1;}
						else{stop=(stop+j+2)/2;}
						break;
					}
				}
			}
			if(center>=minPeak && center<=maxPeak && max>=minHeight && (stop-start)>=minWidth && sum>=minVolume){
				Peak p=new Peak(center, start, stop, max, array[start], array[Tools.min(stop, length-1)], sum, sum2);
				peaks.add(p);
			}else{
//				Peak p=new Peak(center, start, stop, max, sum);
//				System.err.println("*"+p);
			}
		}
		
		capWidth(peaks, maxWidthMult, array);
		
		if(maxPeakCount<peaks.size()){
			peaks=condense(peaks, maxPeakCount);
		}
		
		capWidth(peaks, maxWidthMult, array);
		
		if(peaks.size()>1){
			Peak biggest=peaks.get(biggestPeak(peaks));
			while(peaks.size()>1 && peaks.get(0).volume<0.0001*biggest.volume){
				peaks.remove(0);
			}
		}
		
		if(array!=original){
			recalculate(peaks, original);
		}
		
		return peaks;
	}
	
	private static void recalculate(ArrayList<Peak> peaks, long[] array){
		for(Peak p : peaks){
			p.recalculate(array);
		}
	}
	
	/**
	 * Display usage information.
	 */
	private static void printOptions(){
		System.err.println("Please consult the shellscript for usage information.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Smoothing          ----------------*/
	/*--------------------------------------------------------------*/
	
	public static long[] smoothProgressive(final long[] data, int radius0){
		int radius=radius0;
		long div=radius*radius;
		double mult=1.0/div;
		long[] smoothed=new long[data.length];
		for(int i=0, next=5; i<data.length; i++){
			long sum=sumPoint(data, i, radius);
			double product=sum*mult;
//			if(data[i]>=product){smoothed[i]=(long)Math.ceil(product);}
//			else{smoothed[i]=(long)product;}
			smoothed[i]=Math.round(product);
			if(i>next){
				next=(int)Math.ceil(1+next*progressiveMult);
				radius+=1;
				div=radius*radius;
				mult=1.0/div;
				if(radius>maxRadius){next=Integer.MAX_VALUE;}
//				System.err.println(radius+", "+div);
			}
		}
		return smoothed;
	}
	
	public static long[] smooth(final long[] data, int radius){
		final long div=radius*radius;
		final double mult=1.0/div;
		long[] smoothed=new long[data.length];
		for(int i=0; i<data.length; i++){
			long sum=sumPoint(data, i, radius);
			double product=sum*mult;
//			if(data[i]>=product){smoothed[i]=(long)Math.ceil(product);}
//			else{smoothed[i]=(long)product;}
			smoothed[i]=Math.round(product);
		}
		return smoothed;
	}
	
	private static long sumPoint(long[] data, int loc, int radius){
		long sum=0;
		int start=loc-radius+1;
		int stop=loc+radius-1;
		for(int i=start, x=1; i<loc; i++, x++){
			int i2=Tools.max(i, 0);
			sum+=data[i2]*x;
		}
		for(int i=loc, x=radius, max=data.length-1; i<=stop; i++, x--){
			int i2=Tools.min(i, max);
			sum+=data[i2]*x;
		}
		return sum;
	}
	
	public static int[] smoothProgressive(final int[] data, int radius0){
		int radius=radius0;
		long div=radius*radius;
		double mult=1.0/div;
		int[] smoothed=new int[data.length];
		for(int i=0, next=5; i<data.length; i++){
			long sum=sumPoint(data, i, radius);
			double product=sum*mult;
//			if(data[i]>=product){smoothed[i]=(long)Math.ceil(product);}
//			else{smoothed[i]=(long)product;}
			smoothed[i]=(int)Math.round(product);
			if(i>next){
				next=(int)Math.ceil(next*2);
				radius+=1;
				div=radius*radius;
				mult=1.0/div;
				if(radius>10){next=Integer.MAX_VALUE;}
//				System.err.println(radius+", "+div);
			}
		}
		return smoothed;
	}
	
	public static int[] smooth(final int[] data, int radius){
		final long div=radius*radius;
		final double mult=1.0/div;
		int[] smoothed=new int[data.length];
		for(int i=0; i<data.length; i++){
			long sum=sumPoint(data, i, radius);
			double product=sum*mult;
//			if(data[i]>=product){smoothed[i]=(int)Math.ceil(product);}
//			else{smoothed[i]=(int)product;}
			smoothed[i]=(int)Math.round(product);
		}
		return smoothed;
	}
	
	private static long sumPoint(int[] data, int loc, int radius){
		long sum=0;
		int start=loc-radius+1;
		int stop=loc+radius-1;
		for(int i=start, x=1; i<loc; i++, x++){
			int i2=Tools.max(i, 0);
			sum+=data[i2]*x;
		}
		for(int i=loc, x=radius, max=data.length-1; i<=stop; i++, x--){
			int i2=Tools.min(i, max);
			sum+=data[i2]*x;
		}
		return sum;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class Peak{
		
		Peak(int center_, int start_, int stop_, long centerHeight_, long startHeight_, long stopHeight_, long volume_, long volume2_){
			
			center=center_;
			start=start_;
			stop=stop_;
			
			centerHeight=centerHeight_;
			startHeight=startHeight_;
			stopHeight=stopHeight_;
			volume=volume_;
			volume2=volume2_;
			
			assert(center>=0) : this;
			assert(start<center) : this;
			assert(stop>center) : this;
		}
		
		public boolean compatibleWith(Peak p) {
			int min=Tools.min(center, p.stop);
			int max=Tools.max(stop, p.center);
//			assert(min*maxWidthMult>=max) : this+", "+p+", "+(min*maxWidthMult)+", "+max;
			return min*maxWidthMult>=max;
		}

		/**
		 * @param array
		 */
		public void recalculate(long[] array) {
			centerHeight=array[center];
			startHeight=array[start];
			stopHeight=array[stop];
			volume=0;
			volume2=0;
			for(int i=start; i<stop; i++){
				long x=array[i];
				volume+=x;
				volume2+=(x*i);
			}
		}

		/**
		 * @param p
		 */
		public void absorb(Peak p) {
			
			if(center>p.center){
				assert(p.stop<stop) : "\n"+this+"\n"+p+"\n";
				if(start>p.start){
					start=p.start;
					startHeight=p.startHeight;
				}
			}else{
				assert(p.start>start) : "\n"+this+"\n"+p+"\n";
				if(stop<p.stop){
					stop=p.stop;
					stopHeight=p.stopHeight;
				}
			}
			
			long c1=callByRawCount ? centerHeight2() : centerHeight;
			long c2=callByRawCount ? p.centerHeight2() : p.centerHeight;
			
			if(c1<c2){
				center=p.center;
				centerHeight=p.centerHeight;
			}
			
			volume+=p.volume;
			volume2+=p.volume2;
			
		}

		int width(){return stop-start;}
		
		@Override
		public String toString(){
			return start+"\t"+center+"\t"+stop+"\t"+centerHeight+"\t"+volume;
		}
		
		public ByteBuilder toBytes(ByteBuilder bb){
			if(bb==null){bb=new ByteBuilder();}
			bb.append(start);
			bb.append('\t');
			bb.append(center);
			bb.append('\t');
			bb.append(stop);
			bb.append('\t');
			bb.append(centerHeight);
			bb.append('\t');
			bb.append(volume);
			bb.append('\t');
			return bb;
		}

		/** Inclusive */
		public int start;
		public int center;
		/** Exclusive */
		public int stop;

		//Unique counts
		public long startHeight;
		public long centerHeight;
		public long stopHeight;
		public long volume;

		public long volume2;

		//Raw counts
		public long startHeight2(){return startHeight*start;}
		public long centerHeight2(){return centerHeight*center;}
		public long stopHeight2(){return stopHeight*stop;}
		
		
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private long minHeight=2;
	private long minVolume=2;
	private int minWidth=3;
	private int minPeak=2;
	private int maxPeak=Integer.MAX_VALUE;
	private int maxPeakCount=10;
	private float maxWidthMult=2.5f;
	private int smoothRadius=0;
	private boolean smoothProgressiveFlag=false;
	private int k=31;
	
	private int ploidyClaimed=-1;
	
	private String in;
	private String out;
	
	private final FileFormat ffin;
	private final FileFormat ffout;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static int maxRadius=10;
	public static float progressiveMult=2;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private static int countColumn=1;
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	public static boolean printClass=true;
	
	public static boolean callByRawCount=true;
	
}

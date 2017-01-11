#!/bin/bash
#callvariants in=<file.sam> out=<file.vcf>

usage(){
echo "
Written by Brian Bushnell
Last modified January 4, 2017

Description:  Calls variants from sam or bam input.
In default mode, all input files are combined and treated as a single sample.
In multisample mode, each file is treated as an individual sample,
and gets its own column in the VCF file.  Unless overridden, input file
names are used as sample names.

Usage:  callvariants.sh in=<file,file,...> vcf=<file> ref=<file>

Input may be sorted or unsorted.
The reference should be fasta.


I/O parameters:
in=<file>       Input; may be one file or multiple comma-delimited files.
list=<file>     Optional text file containing one input file per line.
                Use list or in, but not both.
out=<file>      Output variant list in native format.  If the name ends
                with .vcf then it will be vcf format.
vcf=<file>      Output variant list in vcf format.
ref=<file>      Reference fasta.  Required to display ref alleles.
                Variant calling wil be more accurate with the reference.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
extended=f      Print additional variant statistics columns.
sample=         Optional comma-delimited list of sample names.
multisample=f   (multi) Set to true if there are multiple sam/bam files,
                and each should be tracked as an individual sample.
vcf0=           Optional comma-delimited list of per-sample outputs.
                Only used in multisample mode.

Processing Parameters:
prefilter=f     Use a Bloom filter to exclude variants seen fewer than
                minreads times.  Doubles the runtime but greatly reduces
                memory usage.  The results are identical.
samstreamer=t   (ss) Load reads multithreaded to increase speed.
coverage=t      (cc) Calculate coverage, to better call variants.
ploidy=1        Set the organism's ploidy.
rarity=1.0      Penalize the quality of variants with allele fraction 
                lower than this.  For example, if you are interested in
                4% frequency variants, you could set both rarity and
                minallelefraction to 0.04.  This is affected by ploidy - 
                a variant with frequency indicating at least one copy
                is never penalized.
covpenalty=0.8  (lowcoveragepenalty) A lower penalty will increase the 
                scores of low-coverage variants, and is useful for 
                low-coverage datasets.
border=5        Ignore this many bases on either end of reads.
useidentity=t   Include average read identity in score calculation.
usepairing=t    Include pairing rate in score calculation.
usebias=t       Include strand bias in score calculation.
useedist=t      Include read-end distance in score calculation.
homopolymer=t   Penalize scores of substitutions matching adjacent bases.
nscan=t         Consider the distance of a variant from contig ends when 
                calculating strand bias.
32bit=f         Set to true to allow coverage tracking over depth 65545.
strandedcov=t   Tracks per-strand ref coverage to print the DP4 field.
                Requires more memory when enabled.

Realignment parameters:
realign=f       Realign all reads with more than a couple mismatches.
                Decreases speed.  Recommended for aligners other than BBMap.
unclip=f        Convert clip symbols from exceeding the ends of the
                realignment zone into matches and substitutitions.
repadding=70    Pad alignment by this much on each end.  Typically,
                longer is more accurate for long indels, but greatly
                reduces speed.
rerows=602      Use this many rows maximum for realignment.  Reads longer
                than this cannot be realigned.
recols=2000     Reads may not be aligned to reference seqments longer 
                than this.  Needs to be at least read length plus
                max deletion length plus twice padding.
msa=            Select the aligner.  Options:
                   MultiStateAligner11ts:     Default.
                   MultiStateAligner9PacBio:  Use for PacBio reads, or for
                   Illumina reads mapped to PacBio/Nanopore reads.

Variant-Calling Cutoffs:
minreads=2              Ignore variants seen in fewer reads.
minqualitymax=15        Ignore variants with lower max base quality.
minedistmax=20          Ignore variants with lower max distance from read ends.
minmapqmax=0            Ignore variants with lower max mapq.
minidmax=0              Ignore variants with lower max read identity.
minpairingrate=0.1      Ignore variants with lower pairing rate.
minstrandratio=0.1      Ignore variants with lower plus/minus strand ratio.
minquality=12.0         Ignore variants with lower average base quality.
minedist=10.0           Ignore variants with lower average distance from ends.
minmapq=0.0             Ignore variants with lower average mapq.
minallelefraction=0.1   Ignore variants with lower allele fraction.  This
                        should be adjusted for high ploidies.
minid=0                 Ignore variants with lower average read identity.
minscore=20.0           Ignore variants with lower Phred-scaled score.
clearfilters            Reset all filters to zero.


Java Parameters:
-Xmx            This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx4g"
z2="-Xms4g"
EA="-ea"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	parseXmx "$@"
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 4000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

callvariants() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
		module load samtools/1.3.1
	fi
	local CMD="java $EA $z $z2 -cp $CP var2.CallVariants $@"
	echo $CMD >&2
	eval $CMD
}

callvariants "$@"

#!/bin/bash
#merge in=<infile> out=<outfile>

function usage(){
echo "
BBMergeGapped v2.0
Written by Brian Bushnell
Last modified February 17, 2015

***** DEPRECATED *****
Should still work, but is no longer maintained.
Please use bbmerge.sh with the 'extend' flag instead, for nonoverlapping reads.

Description:  Merges paired reads into single reads by overlap detection.
With sufficient coverage, can also merge nonoverlapping reads using gapped kmers.

Usage:	bbmerge.sh in=<input> out=<merged reads> outbad=<unmerged reads>

Input may be fasta or fastq, compressed or uncompressed.



Input parameters:
in2=null            Second input file for paired reads.
extra=null	      Additional files to use for input (generating hash table) but not for output.
interleaved=auto    May be set to true or false to force the input read file to override autodetection of the input file as paired interleaved.
reads=-1            Only process this number of reads, then quit (-1 means all).

Output parameters:
out=<file>          File for merged reads.
outbad=<file>       File for unmerged reads.
outinsert=<file>    File list of read names and their insert sizes.
hist=null           Insert length histogram output file.
ziplevel=2          Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster.

Trimming parameters:
qtrim=f             Trim read ends to remove bases with quality below minq.  Performed BEFORE merging.
                    Values: t (trim both ends), f (neither end), r (right end only), l (left end only).
trimq=10            Trim quality threshold.
minlength=20        (ml) Reads shorter than this after trimming (before merging) will be discarded.  Pairs will be discarded only if both are shorter.
trimonfailure=t     (tof) If detecting insert size by overlap fails, the reads will be trimmed and this will be re-attempted.

Other parameters:
join=t              Create merged reads.  If set to false, you can still generate an insert histogram.
useoverlap=t        Attempt merge based on paired ead overlap.
minoverlapbases=12  Minimum number of overlapping bases to merge reads.
minoverlapinsert=16 Do not look for insert sizes shorter than this.
mininsert=0	      Reads with insert sizes less than this (after merging) will be discarded.
gap=null            Sets gap size for merging via gapped kmers.
                    'gap=50;50,100' would run one pass with a gap size of 50 and another with both 50 and 100.
                    This script sets memory appropriately for ungapped merging only, though.

Java Parameters:
-Xmx                This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
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

z="-Xmx1g"
z2="-Xms1g"
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
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

function merge() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
		module load samtools
	fi
	local CMD="java $EA $z $z2 -cp $CP jgi.MateReadsMT $@"
	echo $CMD >&2
	eval $CMD
}

merge "$@"

#!/bin/bash
#kmercountexact in=<infile>

usage(){
echo "
Written by Brian Bushnell
Last modified November 2, 2016

Description:  Counts the number of unique kmers in a file.
Generates a kmer frequency histogram and genome size estimate (in peaks output),
and prints a file containing all kmers and their counts.
Supports K=1 to infinity, though not all values are allowed.
SEE ALSO: bbnorm.sh/khist.sh, which have similar functionality.

Usage:   kmercountexact.sh in=<file> khist=<file> peaks=<file> out=<file>

Input may be fasta or fastq, compressed or uncompressed.
Output may be stdout or a file.  out, khist, and peaks are optional.


Input parameters:
in=<file>           Primary input file.
in2=<file>          Second input file for paired reads.

Output parameters:
out=<file>          Print kmers and their counts.
fastadump=t         Print kmers and counts as fasta versus 2-column tsv.
mincount=1          Only print kmers with at least this depth.
reads=-1            Only process this number of reads, then quit (-1 means all).
dumpthreads=-1      Use this number of threads for dumping kmers (-1 means auto).

Hashing parameters:
k=31                Kmer length (1-31 is fastest).
prealloc=t          Pre-allocate memory rather than dynamically growing; faster and more memory-efficient.  A float fraction (0-1) may be specified, default 1.
prefilter=0         If set to a positive integer, use a countmin sketch to ignore kmers with depth of that value or lower.
prehashes=2         Number of hashes for prefilter.
prefiltersize=0.2   Fraction of memory to use for prefilter.
minq=6              Ignore kmers containing bases with quality below this. (TODO)
minprob=0.0         Ignore kmers with overall probability of correctness below this.
threads=X           Spawn X hashing threads (default is number of logical processors).
onepass=f           If true, prefilter will be generated in same pass as kmer counts.  Much faster but counts will be lower, by up to prefilter's depth limit.
rcomp=t             Store and count each kmer together and its reverse-complement.

Histogram parameters:
khist=<file>        Print kmer frequency histogram.
histcolumns=2       2 columns: (depth, count).  3 columns: (depth, rawCount, count).
histmax=100000      Maximum depth to print in histogram output.
histheader=t        Set true to print a header line.
nzo=t               (nonzeroonly) Only print lines for depths with a nonzero kmer count.
smooth=f            Smooth the kmer histogram.
smoothradius=1      Initial radius of progressive smoothing function.
maxradius=10        Maximum radius of progressive smoothing function.
progressivemult=2   Increment radius each time depth increases by this factor.
gchist=f            Add an extra histogram column with the average GC%.

Peak calling parameters:
peaks=<file>        Write the peaks to this file.  Default is stdout. 
                    Also contains the genome size estimate in bp.
minHeight=2         (h) Ignore peaks shorter than this.
minVolume=2         (v) Ignore peaks with less area than this.
minWidth=2          (w) Ignore peaks narrower than this.
minPeak=2           (minp) Ignore peaks with an X-value below this.
maxPeak=BIG         (maxp) Ignore peaks with an X-value above this.
maxPeakCount=12     (maxpc) Print up to this many peaks (prioritizing height).
ploidy=-1           Specify ploidy; otherwise it will be autodetected.

Sketch parameters (for making a MinHashSketch):
sketch=<file>       Write a minhash sketch to this file.
sketchlen=10000     Output the top 10000 kmers.  Only kmers with at least mincount are included.
sketchname=         Name of output sketch.
sketchid=           taxID of output sketch.

Quality parameters:
qtrim=f             Trim read ends to remove bases with quality below minq.
                    Values: t (trim both ends), f (neither end), r (right end only), l (left end only).
trimq=4             Trim quality threshold.
minavgquality=0     (maq) Reads with average quality (before trimming) below this will be discarded.

Overlap parameters (for overlapping paired-end reads only):
merge=f             Attempt to merge reads before counting kmers.
ecco=f              Error correct via overlap, but do not merge reads.   

Java Parameters:
-Xmx                This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.
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

kmercountexact() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z $z2 -cp $CP jgi.KmerCountExact $@"
	echo $CMD >&2
	eval $CMD
}

kmercountexact "$@"

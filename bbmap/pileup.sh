#!/bin/bash
#pileup in=<infile> out=<outfile>

usage(){
echo "
Written by Brian Bushnell
Last modified December 5, 2016

Description:  Calculates per-scaffold coverage information from an unsorted sam or bam file.

Usage:        pileup.sh in=<input> out=<output>


Input Parameters:
in=<file>           The input sam file; this is the only required parameter.
ref=<file>          Scans a reference fasta for per-scaffold GC counts, not otherwise needed.
fastaorf=<file>     An optional fasta file with ORF header information in PRODIGAL's output format.  Must also specify 'outorf'.
unpigz=t            Decompress with pigz for faster decompression.

Output Parameters:
out=<file>          (covstats) Per-scaffold coverage info.
rpkm=<file>         Per-scaffold RPKM/FPKM counts.
twocolumn=f         Change to true to print only ID and Avg_fold instead of all 6 columns.
countgc=t           Enable/disable counting of read GC content.
outorf=<file>       Per-orf coverage info to this file (only if 'fastaorf' is specified).
outsam=<file>       Print the input sam stream to this file (or stdout).  Useful for piping data.
hist=<file>         Histogram of # occurrences of each depth level.
basecov=<file>      Coverage per base location.
bincov=<file>       Binned coverage per location (one line per X bases).
binsize=1000        Binsize for binned coverage output.
keepshortbins=t     (ksb) Keep residual bins shorter than binsize.
normcov=<file>      Normalized coverage by normalized location (X lines per scaffold).
normcovo=<file>     Overall normalized coverage by normalized location (X lines for the entire assembly).
normb=-1            If positive, use a fixed number of bins per scaffold; affects 'normcov' and 'normcovo'.
normc=f             Normalize coverage to fraction of max per scaffold; affects 'normcov' and 'normcovo'.
delta=f             Only print base coverage lines when the coverage differs from the previous base.
nzo=f               Only print scaffolds with nonzero coverage.
concise=f           Write 'basecov' in a more concise format.
header=t            (hdr) Include headers in output files.
headerpound=t       (#) Prepend header lines with '#' symbol.
stdev=t             Calculate coverage standard deviation.
covminscaf=0        (minscaf) Don't print coverage for scaffolds shorter than this.
covwindow=0         Calculate how many bases are in windows of this size with
                    low average coverage.  Produces an extra stats column.
covwindowavg=5      Average coverage below this will be classified as low.

Processing Parameters:
strandedcov=f       Track coverage for plus and minus strand independently.
startcov=f          Only track start positions of reads.
secondary=t         Use secondary alignments, if present.
softclip=f          Include soft-clipped bases in coverage.
minmapq=0           (minq) Ignore alignments with mapq below this.
physical=f          (physcov) Calculate physical coverage for paired reads.  This includes the unsequenced bases.
tlen=t              Track physical coverage from the tlen field rather than recalculating it.
arrays=auto         Set to t/f to manually force the use of coverage arrays.  Arrays and bitsets are mutually exclusive.
bitsets=auto        Set to t/f to manually force the use of coverage bitsets.
32bit=f             Set to true if you need per-base coverage over 64k; does not affect per-scaffold coverage precision.
                    This option will double RAM usage (when calculating per-base coverage).
delcoverage=t       (delcov) Count bases covered by deletions as covered.
                    True is faster than false.
samstreamer=t       (ss) Load reads multithreaded to increase speed.

Java Parameters:
-Xmx                This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.

Output format (tab-delimited):
ID, Avg_fold, Length, Ref_GC, Covered_percent, Covered_bases, Plus_reads, Minus_reads, Read_GC, Median_fold, Std_Dev

ID:                Scaffold ID
Length:            Scaffold length
Ref_GC:            GC ratio of reference
Avg_fold:          Average fold coverage of this scaffold
Covered_percent:   Percent of scaffold with any coverage (only if arrays or bitsets are used)
Covered_bases:     Number of bases with any coverage (only if arrays or bitsets are used)
Plus_reads:        Number of reads mapped to plus strand
Minus_reads:       Number of reads mapped to minus strand
Read_GC:           Average GC ratio of reads mapped to this scaffold
Median_fold:       Median fold coverage of this scaffold (only if arrays are used)
Std_Dev:           Standard deviation of coverage (only if arrays are used)

Notes:

Only supports SAM format for reads and FASTA for reference (though either may be gzipped).
Sorting is not needed, so output may be streamed directly from a mapping program.
Requires approximately 1 bit per reference base plus 100 bytes per scaffold (even if no reference is specified).
This script will attempt to autodetect and correctly specify the -Xmx parameter to use all memory on the target node.
If this fails with a message including 'Error: Could not create the Java Virtual Machine.', then...
Please decrease the -Xmx parameter.  It should be set to around 85% of the available memory.
For example, -Xmx20g needs around 23 GB of virtual (and physical) memory when qsubbed.
If the program fails with a message including 'java.lang.OutOfMemoryError:', then...
-Xmx needs to be increased, which probably also means it needs to be qsubbed with a higher memory allocation.

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

pileup() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
		module load samtools
	fi
	local CMD="java $EA $z -cp $CP jgi.CoveragePileup $@"
	echo $CMD >&2
	eval $CMD
}

pileup "$@"

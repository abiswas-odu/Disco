#!/bin/bash
#bbcountunique in=<infile> out=<outfile>

usage(){
echo "
Written by Brian Bushnell
Last modified August 1, 2016

Description:  Generates a kmer uniqueness histogram, binned by file position.
There are 3 columns for single reads, 6 columns for paired:
count        number of reads or pairs processed
r1_first     percent unique 1st kmer of read 1
r1_rand      percent unique random kmer of read 1
r2_first     percent unique 1st kmer of read 2
r2_rand      percent unique random kmer of read 2
pair         percent unique concatenated kmer from read 1 and 2

Usage:	bbcountunique.sh in=<input> out=<output>


Optional parameters (and their defaults)

Input parameters:
in2=null            Second input file for paired reads
interleaved=auto    Set true/false to override autodetection of the input file as paired interleaved.
samplerate=1        Set to below 1 to sample a fraction of input reads.
reads=-1            Only process this number of reads, then quit (-1 means all)

Output parameters:
out=<file>          File for output stats

Processing parameters:
k=25                Kmer length (range 1-31).
interval=25000      Print one line to the histogram per this many reads.
cumulative=f        Show cumulative numbers rather than per-interval numbers.
percent=t           Show percentages of unique reads.
count=f             Show raw counts of unique reads.
printlastbin=f      (plb) Print a line for the final undersized bin.
minprob=0           Ignore kmers with a probability of correctness below this (based on q-scores).


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

bbcountunique() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z $z2 -cp $CP jgi.CalcUniqueness $@"
	echo $CMD >&2
	eval $CMD
}

bbcountunique "$@"

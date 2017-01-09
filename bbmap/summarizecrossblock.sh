#!/bin/bash
#summarizecrossblock in=<infile> out=<outfile>

usage(){
echo "
Written by Brian Bushnell
Last modified June 10, 2016

Description:  Does nothing.  Should be fast.

Usage:  summarizecrossblock.sh in=<input file> out=<output file>

Input may be fasta or fastq, compressed or uncompressed.


Standard parameters:
in=<file>       A text file of files, or a comma-delimited list of files.
                Each is a path to results.txt output from Crossblock.
out=<file>      Output file for the summary.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Processing parameters:
None yet!

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

z="-Xmx200m"
EA="-ea"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	parseXmx "$@"
}
calcXmx "$@"

summarizecrossblock() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z -cp $CP driver.SummarizeCrossblock $@"
	echo $CMD >&2
	eval $CMD
}

summarizecrossblock "$@"

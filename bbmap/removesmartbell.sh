#!/bin/bash
#removesmartbell in=<infile> out=<outfile>

usage(){
echo "
Written by Brian Bushnell
Last modified February 17, 2015

Description:  Remove Smart Bell adapters from PacBio reads.

Usage:        removesmartbell in=<input> out=<output> split=t

Input may be fasta or fastq, compressed or uncompressed (not H5 files).


Parameters:
in=file         Specify the input file, or stdin.
out=file        Specify the output file, or stdout.
split=f         'split=t' splits reads at adapters; split=f masks adapters.

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

z="-Xmx400m"
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

removesmartbell() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
	fi
	local CMD="java $EA $z -cp $CP pacbio.RemoveAdapters2 $@"
	echo $CMD >&2
	eval $CMD
}

removesmartbell "$@"

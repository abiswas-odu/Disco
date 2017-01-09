#!/bin/bash
#commonkmers in=<infile> out=<outfile>

function usage(){
echo "
Written by Brian Bushnell
Last modified February 20, 2015

Description:  Prints the most common kmers in each sequence.
This is intended for short kmers only!

Usage:  commonkmers.sh in=<file> out=<file>


Parameters:
k=2                     Kmer length, 0-12.
display=3               Print this many kmers per sequence.
count=f                 Print the kmer counts as well.

ow=f                    (overwrite) Overwrites files that already exist.
app=f                   (append) Append to files that already exist.
zl=4                    (ziplevel) Set compression level, 1 (low) to 9 (max).
qin=auto                ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.

Java Parameters:
-Xmx                    This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
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

z="-Xmx800m"
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

function commonkmers() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z -cp $CP jgi.SmallKmerFrequency $@"
	echo $CMD >&2
	eval $CMD
}

commonkmers "$@"

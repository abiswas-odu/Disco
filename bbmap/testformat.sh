#!/bin/bash
#testformat in=<infile>

usage(){
echo "
Written by Brian Bushnell
Last modified August 4, 2016

Description:  Tests file extensions and contents to determine format, quality, compression, interleaving, and read length.

Usage:        testformat.sh <file>

More than one file may be specified.
Note that ASCII-33 (sanger) and ASCII-64 (illumina) cannot always be differentiated.

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

z="-Xmx120m"
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

testformat() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
	fi
	local CMD="java $EA $z -cp $CP fileIO.FileFormat $@"
#	echo $CMD >&2
	eval $CMD
}

testformat "$@"

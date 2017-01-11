#!/bin/bash

function usage(){
echo "
Written by Brian Bushnell
Last modified February 17, 2015

Description:  Selects reads with designated numeric IDs.

Usage:  getreads.sh in=<file> id=<number,number,number...> out=<file>

The first read (or pair) has ID 0, the second read (or pair) has ID 1, etc.

Parameters:
in=<file>       Specify the input file, or stdin.
out=<file>      Specify the output file, or stdout.
id=             Comma delimited list of numbers or ranges, in any order.
                For example:  id=5,93,17-31,8,0,12-13

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

function tf() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
	fi
	local CMD="java $EA $z -cp $CP jgi.GetReads $@"
	echo $CMD >&2
	eval $CMD
}

tf "$@"
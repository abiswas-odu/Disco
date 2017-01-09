#!/bin/bash
#dedupe in=<infile> out=<outfile>

usage(){
echo "
Written by Brian Bushnell and Jonathan Rood
Last modified September 15, 2015

Dedupe2 is identical to Dedupe except it supports hashing unlimited kmer
prefixes and suffixes per sequence.  Dedupe supports at most 2 of each,
but uses slightly more memory.  You can manually set the number of kmers to
hash per read with the numaffixmaps (nam) flag.  Dedupe will automatically
call Dedupe2 if necessary (if nam=3 or higher).

For documentation, please consult dedupe.sh; syntax is identical.
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
NATIVELIBDIR="$DIR""jni/"

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

dedupe() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java -Djava.library.path=$NATIVELIBDIR $EA $z $z2 -cp $CP jgi.Dedupe2 $@"
	echo $CMD >&2
	eval $CMD
}

dedupe "$@"

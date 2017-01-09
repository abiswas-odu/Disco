#!/bin/bash
#ecc.sh in=<infile> out=<outfile>

usage(){
echo "
Description:  Corrects substitution errors in reads using kmer depth information.
Can also normalize and/or bin reads by kmer depth.

Usage:	ecc.sh in=<input> out=<reads to keep> outt=<reads to toss> hist=<histogram output>

Please see bbnorm.sh for more information.  
All the flags are the same, only the parameters (near the bottom of this file) differ.
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

z="-Xmx31g"
z2="-Xms31g"
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
	freeRam 31000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

correct() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z $z2 -cp $CP jgi.KmerNormalize bits=16 ecc=t passes=1 keepall dr=f prefilter $@"
	echo $CMD >&2
	eval $CMD
}

correct "$@"

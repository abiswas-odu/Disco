#!/bin/bash
#merge in=<infile> out=<outfile>

function usage(){
echo "
bbmerge-auto.sh is a wrapper for BBMerge that attempts to use all available
memory, instead of a fixed amount.  This is for use with the Tadpole options
of error-correction (ecct) and extension, which require more memory.
For merging by overlap only, please use bbmerge.sh.  If you set memory
manually with the -Xmx flag, bbmerge.sh and bbmerge-auto.sh are equivalent.

For information about usage and parameters, please run bbmerge.sh.
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

z="-Xmx14g"
z2="-Xms14g"
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
	freeRam 15000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

function merge() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java -Djava.library.path=$NATIVELIBDIR $EA $z $z2 -cp $CP jgi.BBMerge $@"
	echo $CMD >&2
	eval $CMD
}

merge "$@"

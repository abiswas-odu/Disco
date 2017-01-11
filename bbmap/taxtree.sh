#!/bin/bash

usage(){
echo "
Written by Brian Bushnell.
Last modified December 15, 2015

Description:  Creates tree.taxtree from names.dmp and nodes.dmp.
These are in taxdmp.zip available at ftp://ftp.ncbi.nih.gov/pub/taxonomy/
The taxtree file is needed for programs that can deal with taxonomy, 
like Seal and SortByTaxa.

Usage:  taxtree.sh names.dmp nodes.dmp tree.taxtree.gz

Java Parameters:
-Xmx    This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
        -Xmx20g will specify 20 gigs of RAM.  The max is typically 85% of physical memory.

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
NATIVELIBDIR="$DIR""jni/"

z="-Xmx2g"
z2="-Xms2g"
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
	freeRam 2000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"


taxtree() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z $z2 -cp $CP tax.TaxTree $@"
	echo $CMD >&2
	eval $CMD
}

taxtree "$@"

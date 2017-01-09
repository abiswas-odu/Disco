#!/bin/bash
#tadwrapper in=<infile>

usage(){
echo "
Written by Brian Bushnell
Last modified August 18, 2016

Description:  Generates multiple assemblies with Tadpole
to estimate the optimal kmer length.

Usage:
tadwrapper.sh in=reads.fq out=contigs%.fa k=31,62,93

Parameters:
out=<file>          Output file name.  Must contain a % symbol.
outfinal=<file>     Optional.  If set, the best assembly file 
                    will be renamed to this.
k=31                Comma-delimited list of kmer lengths.
delete=f            Delete assemblies before terminating.
quitearly=f         Quit once metrics stop improving with longer kmers.
bisect=f            Recursively assemble with the kmer midway between
                    the two best kmers until improvement halts.
expand=f            Recursively assemble with kmers shorter or longer
                    than the current best until improvement halts.

All other parameters are passed to Tadpole.
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

tadwrapper() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z $z2 -cp $CP assemble.TadpoleWrapper $@"
	echo $CMD >&2
	eval $CMD
}

tadwrapper "$@"

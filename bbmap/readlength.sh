#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified December 19, 2018
Description:  Generates a length histogram of input reads.

Usage:	readlength.sh in=<input file>

in=<file>    	The 'in=' flag is needed only if the input file is not the first parameter.  'in=stdin.fq' will pipe from standard in.
in2=<file>   	Use this if 2nd read of pairs are in a different file.
out=<file>   	Write the histogram to this file.  Default is stdout.
bin=10       	Set the histogram bin size.
max=80000    	Set the max read length to track.
round=f      	Places reads in the closest bin, rather than the highest bin of at least readlength.
nzo=t        	(nonzeroonly) Do not print empty bins.
reads=-1     	If nonnegative, stop after this many reads.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

#This block allows symlinked shellscripts to correctly set classpath.
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
z2="-Xms400m"
EA="-ea"
EOOM=""
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

stats() {
	if [[ $SHIFTER_RUNTIME == 1 ]]; then
		#Ignore NERSC_HOST
		shifter=1
	elif [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_144_64bit
		module load pigz
	elif [[ $NERSC_HOST == denovo ]]; then
		module unload java
		module load java/1.8.0_144
		module load pigz
	elif [[ $NERSC_HOST == cori ]]; then
		module use /global/common/software/m342/nersc-builds/denovo/Modules/jgi
		module use /global/common/software/m342/nersc-builds/denovo/Modules/usg
		module unload java
		module load java/1.8.0_144
		module load pigz
	fi
	local CMD="java $EA $EOOM $z -cp $CP jgi.MakeLengthHistogram $@"
#	echo $CMD >&2
	eval $CMD
}

stats "$@"

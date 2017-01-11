#!/bin/bash
#cutprimers in=<file> out=<file>

usage(){
echo "
Written by Brian Bushnell
Last modified February 17, 2015

Description:  Cuts out sequences corresponding to primers identified in sam files.

Usage:	cutprimers.sh in=<file> out=<file> sam1=<file> sam2=<file>

Parameters:

in=<file>           File containing reads. in=stdin.fa will pipe from stdin.
out=<file>          Output sequences. out=stdout will pipe to stdout.
sam1=<file>         Sam file containing mapped locations of primer sequence 1.
sam2=<file>         Sam file containing mapped locations of primer sequence 2.
fake=t              Output 1bp 'N' reads in cases where there is no primer.
include=f           Include the flanking primer sequences in output.

Java Parameters:

-Xmx                This will be passed to Java to set memory usage, overriding 
                    the program's automatic memory detection. -Xmx20g will specify 
                    20 gigs of RAM, and -Xmx200m will specify 200 megs.  
                    The max is typically 85% of physical memory.

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
	freeRam 2000m 42
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

cutprimers() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z -cp $CP jgi.CutPrimers $@"
	echo $CMD >&2
	eval $CMD
}

cutprimers "$@"

#!/bin/bash
#gi2taxid in=<infile> out=<outfile>

usage(){
echo "
Written by Brian Bushnell.
Last modified May 11, 2016

Description:  Renames fasta sequences with gi numbers to NCBI taxa IDs.

Usage:  gi2taxid.sh in=<file> out=<file> gi=<file>

Parameters:
in=<file>           Input sequences; required parameter.  Must be fasta.
out=<file>          Destination for renamed sequences.
invalid=<file>      Destination for headers with no taxid.
gi=<file>           2-column tsv with gi and taxid numbers.
prefix=f            Append the taxid as a prefix to the old header.
ziplevel=2          (zl) Compression level for gzip output.
pigz=f              Spawn a pigz (parallel gzip) process for faster 
                    compression than Java.  Requires pigz to be installed.

Java Parameters:
-Xmx                This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx800m will specify 800 megs.  The max is typically 85% of physical memory.

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

z="-Xmx7g"
z2="-Xms7g"
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
	freeRam 7000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"


gi2taxid() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z $z2 -cp $CP tax.RenameGiToNcbi $@"
	echo $CMD >&2
	eval $CMD
}

gi2taxid "$@"

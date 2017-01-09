#!/bin/bash
#gitable gi_taxid_nucl.dmp.gz gitable.int1d.gz

usage(){
echo "
Written by Brian Bushnell.
Last modified May 11, 2016

Description:  Creates gitable.int1d from gi_taxid_nucl.dmp and/or gi_taxid_prot.dmp.
gitable.int1d is a much more efficient representation,
allowing easy translation of gi numbers to ncbi taxids.
Dump files are at ftp://ftp.ncbi.nih.gov/pub/taxonomy/

Usage:  gitable.sh gi_taxid_nucl.dmp.gz,gi_taxid_prot.dmp.gz gitable.int1d.gz


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

z="-Xmx12g"
z2="-Xms12g"
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
	freeRam 12000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"


gitable() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z $z2 -cp $CP tax.GiToNcbi $@"
	echo $CMD >&2
	eval $CMD
}

gitable "$@"

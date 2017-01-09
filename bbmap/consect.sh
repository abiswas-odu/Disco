#!/bin/bash
#consect in=<infile> out=<outfile>

usage(){
echo "
Written by Brian Bushnell
Last modified October 25, 2016

Description:  Generates the conservative consensus of multiple
error-correction tools.  Corrections will be accepted only
if all tools agree.  This tool is designed for substitutions only,
not indel corrections.

Usage:  consect.sh in=<file,file,file,...> out=<file>

Input may be fasta or fastq, compressed or uncompressed.


Standard parameters:
in=             A comma-delimited list of files; minimum of 3.
                All files must have reads in the same order.
                The first file must contain the uncorrected reads.
                All additional files must contain corrected reads.
out=<file>      Output of consensus reads.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Processing parameters:
None yet!

Java Parameters:
-Xmx            This will be passed to Java to set memory usage, overriding 
                the program's automatic memory detection.  -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.

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

z="-Xmx2g"
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

consect() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z -cp $CP jgi.Consect $@"
	echo $CMD >&2
	eval $CMD
}

consect "$@"

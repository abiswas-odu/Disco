#!/bin/bash
#mutate in=<infile> out=<outfile> id=<identity>

usage(){
echo "
Written by Brian Bushnell
Last modified June 8, 2016

Description:  Creates a mutant version of a genome.

Usage:  mutate.sh in=<input file> out=<output file> id=<identity>

Input may be fasta or fastq, compressed or uncompressed.


Standard parameters:
in=<file>       Input genome.
out=<file>      Output genome.
subrate=0       Substitution rate, 0 to 1.     
indelrate=0     Indel rate, 0 to 1.
id=1            Target identity, 0 to 1; 1 means 100%.
                If this is used it will override subrate and indelrate;
                99% of the mutations will be substitutions, and 1% indels.
prefix=         Set this flag to rename the new contigs with this prefix
                and a number.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
showspeed=t     (ss) Set to 'f' to suppress display of processing speed.
ziplevel=2      (zl) Set to 1 (lowest) through 9 (max) to change compression
                level; lower compression is faster.

Processing parameters:
None yet!

Java Parameters:
-Xmx            This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.

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

z="-Xmx4g"
z2="-Xms4g"
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
	freeRam 4000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

mutate() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z $z2 -cp $CP jgi.MutateGenome $@"
	echo $CMD >&2
	eval $CMD
}

mutate "$@"

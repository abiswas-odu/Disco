#!/bin/bash
#sortbyname in=<infile> out=<outfile> groups=<number>

usage(){
echo "
Written by Brian Bushnell
Last modified November 9, 2016

Description:  Sorts reads by name.
Writes temp files if memory is exceeded.

Usage:   sortbyname.sh in=<file> out=<file>

Input may be fasta, fastq, or sam, compressed or uncompressed.

Optional parameters (and their defaults)

in=<file>           Input file.
out=<file>          Output file.
delete=t            Delete temporary files.
name=t              Sort reads by name.
length=f            Sort reads by length.
quality=f           Sort reads by quality.
ascending=t         Sort ascending.

Note: name, length, and quality are mutually exclusive.
Sorting by quality actually sorts by average expected error rate,
so ascending will place the highest-quality reads first.

Java Parameters:
-Xmx                This will be passed to Java to set memory usage, overriding
                    the program's automatic memory detection.  -Xmx20g will 
                    specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
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

sortbyname() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z $z2 -cp $CP sort.SortByName $@"
	echo $CMD >&2
	eval $CMD
}

sortbyname "$@"

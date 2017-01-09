#!/bin/bash
#rename in=<infile> out=<outfile>

function usage(){
echo "
Written by Brian Bushnell
Last modified July 31, 2015

Description:  Reduces Silva entries down to one entry per taxa.

Usage:  reducesilva.sh in=<file> out=<file> column=<1>


Parameters:
column              The taxonomic level.  0=species, 1=genus, etc.
ow=f                (overwrite) Overwrites files that already exist.
zl=4                (ziplevel) Set compression level, 1 (low) to 9 (max).
fastawrap=70        Length of lines in fasta output.


Sampling parameters:
reads=-1            Set to a positive number to only process this many INPUT sequences, then quit.

Java Parameters:
-Xmx                This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.

To read from stdin, set 'in=stdin'.  The format should be specified with an extension, like 'in=stdin.fq.gz'
To write to stdout, set 'out=stdout'.  The format should be specified with an extension, like 'out=stdout.fasta'

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

function reducesilva() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z -cp $CP driver.ReduceSilva $@"
	echo $CMD >&2
	eval $CMD
}

reducesilva "$@"

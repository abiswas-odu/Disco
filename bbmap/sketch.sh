#!/bin/bash
#sketch in=<infile> out=<outfile>

usage(){
echo "
Written by Brian Bushnell
Last modified December 22, 2016

Description:  Creates one or more MinHashSketches from a fasta file,
optionally annotated with taxonomic information.


Usage:  sketch.sh in=<fasta file> out=<sketch file>


Standard parameters:
in=<file>           A fasta file containing one or more sequences.
out=<file>          Output filename.  If multiple files are desired it must
                    contain the # symbol.
files=1             Number of parallel files, for faster loading.  This is
                    independent of the number of sketches produced; sketches
                    will be randomly distributed between file.
size=10000          Desired size of sketches.
minsize=100         Do not generate sketches for genomes smaller than this.
k=31                Kmer length.
rcomp=t             Look at reverse-complement kmers also.
mode=single         Possible modes:
                       single: Write one sketch.
                       sequence: Write one sketch per sequence.
                       taxa: Write one sketch per taxonomic unit.
                          Requires more memory, and taxonomic annotation.
                       img: Write one sketch per IMG id.
tree=               Specify a taxtree file.  On Genepool, use 'auto'.
gi=                 Specify a gitable file.  On Genepool, use 'auto'.
imgdump=            Specify an IMG dump file containing NCBI taxIDs,
                    for IMG mode.
taxid=              Set the sketch taxid, in single-sketch mode.
name=               Set the sketch name, in single-sketch mode.
delta=t             Delta-compress sketches.
a48=t               Encode sketches as ASCII-48 rather than hex.
prefilter=f         For huge datasets full of junk like nt, this flag
                    will save memory by ignoring taxa smaller than minsize.
                    Requires taxonomic information (tree and gi).

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

sketch() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
	fi
	local CMD="java $EA $z -cp $CP sketch.SketchMaker $@"
	echo $CMD >&2
	eval $CMD
}

sketch "$@"

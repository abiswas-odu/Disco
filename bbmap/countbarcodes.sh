#!/bin/bash
#filterbarcodes in=<infile> out=<outfile>

usage(){
echo "
Written by Brian Bushnell
Last modified October 16, 2015

Description: Counts the number of reads with each barcode.

Usage:   countbarcodes.sh in=<file> counts=<file>

Input may be stdin or a fasta or fastq file, raw or gzipped.
If you pipe via stdin/stdout, please include the file type; e.g. for gzipped fasta input, set in=stdin.fa.gz


Optional parameters (and their defaults)

Input parameters:
in=<file>           Input reads, whose names end in a colon then barcode.
counts=<file>       Output of counts.
interleaved=auto    (int) If true, forces fastq input to be paired and interleaved.
qin=auto            ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.
unpigz=t            Use pigz to decompress.
expected=           Comma-delimited list of expected bar codes.
valid=              Comma-delimited list of valid bar codes.
countundefined=t    Count barcodes that contain non-ACGT symbols.
printheader=t       Print a header.
maxrows=-1          Optionally limit the number of rows printed.

Output parameters:
out=<file>          Write bar codes and counts here.  'out=stdout' will pipe to standard out.

Java Parameters:
-Xmx                This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
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

z="-Xmx200m"
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

countbarcodes() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z -cp $CP jgi.CountBarcodes $@"
	echo $CMD >&2
	eval $CMD
}

countbarcodes "$@"

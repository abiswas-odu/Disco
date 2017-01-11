#!/bin/bash
#translate6frames in=<infile> out=<outfile>

usage(){
echo "
Written by Brian Bushnell
Last modified February 17, 2015

Description:  Translates nucleotide sequences to all 6 amino acid frames,
or amino acids to a canonical nucleotide representation.

Usage:        translate6frames.sh in=<input file> out=<output file>

Input may be fasta or fastq, compressed or uncompressed.


Optional parameters (and their defaults)

Input parameters:
in=<file>           Main input. in=stdin.fa will pipe from stdin.
in2=<file>          Input for 2nd read of pairs in a different file.
interleaved=auto    (int) t/f overrides interleaved autodetection.
qin=auto            Input quality offset: 33 (Sanger), 64, or auto.
aain=f              False if input is nucleotides, true for amino acids.
reads=-1            If positive, quit after processing X reads or pairs.

Output parameters:
out=<file>          Write output here.  'out=stdout.fa' goes to standard out.
out2=<file>         Use this to write 2nd read of pairs to a different file.
overwrite=t         (ow) Grant permission to overwrite files.
append=f            Append to existing files.
ziplevel=2          (zl) Compression level; 1 (min) through 9 (max).
fastawrap=80        Length of lines in fasta output.
qout=auto           Output quality offset: 33 (Sanger), 64, or auto.
aaout=t             False to output nucleotides, true for amino acids.
tag=t               Tag read id with the frame, adding e.g. ' fr1'
frames=6            Only print this many frames.  
                    If you already know the sense, set 'frames=3'.

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
	freeRam 2000m 42
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

translate6frames() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z -cp $CP jgi.TranslateSixFrames $@"
	echo $CMD >&2
	eval $CMD
}

translate6frames "$@"

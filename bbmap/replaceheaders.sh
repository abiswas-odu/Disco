#!/bin/bash
#reheader in=<infile> out=<outfile>

function usage(){
echo "
Written by Brian Bushnell
Last modified May 23, 2016

Description:  Replaces read names with names from another file.
The other file can either be sequences or simply names, with
one name per line (and no > or @ symbols).  If you use one name
per line, please give the file a .header extension.

Usage:  replaceheaders.sh in=<file> hin=<headers file> out=<out file>


Parameters:
in=                 Input sequences.  Use in2 for a second paired file.
in=                 Header input sequences.  Use hin2 for a second paired file.
out=                Output sequences.  Use out2 for a second paired file.
qfin=<.qual file>   Optional input quality file.
qfout=<.qual file>  Optional output quality file.
ow=f                (overwrite) Overwrites files that already exist.
zl=4                (ziplevel) Set compression level, 1 (low) to 9 (max).
int=f               (interleaved) Determines whether INPUT file is considered interleaved.
fastawrap=70        Length of lines in fasta output.
qin=auto            ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.
qout=auto           ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input).

Renaming modes (if not default):
addprefix=f         Rename the read by prepending the new name to the existing name.

Sampling parameters:
reads=-1            Set to a positive number to only process this many INPUT reads (or pairs), then quit.

Java Parameters:
-Xmx                This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.

Supported input formats are fastq, fasta, fast+qual, scarf, and bread (BBMap's native format)
Supported output formats are fastq, fasta, fast+qual, bread, sam, and bam (bam only if samtools is installed)
Supported compression formats are gz, zip, and bz2
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

function reheader() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z -cp $CP jgi.ReplaceHeaders $@"
	echo $CMD >&2
	eval $CMD
}

reheader "$@"

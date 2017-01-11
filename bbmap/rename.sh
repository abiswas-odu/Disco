#!/bin/bash
#rename in=<infile> out=<outfile>

function usage(){
echo "
Written by Brian Bushnell
Last modified June 24, 2015

Description:  Renames reads to <prefix>_<number> where you specify the prefix and the numbers are ordered.

Usage:  rename.sh in=<file> in2=<file2> out=<outfile> out2=<outfile2> prefix=<>

in2 and out2 are for paired reads and are optional.
If input is paired and there is only one output file, it will be written interleaved.


Parameters:
ow=f                (overwrite) Overwrites files that already exist.
zl=4                (ziplevel) Set compression level, 1 (low) to 9 (max).
int=f               (interleaved) Determines whether INPUT file is considered interleaved.
fastawrap=70        Length of lines in fasta output.
minscaf=1           Ignore fasta reads shorter than this.
qin=auto            ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.
qout=auto           ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input).
qfin=<.qual file>   Read qualities from this qual file, for the reads coming from 'in=<fasta file>'
qfin2=<.qual file>  Read qualities from this qual file, for the reads coming from 'in2=<fasta file>'
qfout=<.qual file>  Write qualities from this qual file, for the reads going to 'out=<fasta file>'
qfout2=<.qual file> Write qualities from this qual file, for the reads coming from 'out2=<fasta file>'
ignorebadquality=f  (ibq) Fix out-of-range quality values instead of crashing with a warning.

Renaming modes (if not default):
renamebyinsert=f    Rename the read to indicate its correct insert size.
renamebymapping=f   Rename the read to indicate its correct mapping coordinates.
renamebytrim=f      Rename the read to indicate its correct post-trimming length.
addprefix=f         Rename the read by prepending the prefix to the existing name.
prefixonly=f        Only use the prefix; don't add _<number>

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

function rename() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z -cp $CP jgi.RenameReads $@"
	echo $CMD >&2
	eval $CMD
}

rename "$@"

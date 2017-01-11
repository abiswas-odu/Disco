#!/bin/bash
#demuxbyname in=<infile> out=<outfile>

function usage(){
echo "
Written by Brian Bushnell
Last modified December 12, 2016

Description:  Demultiplexes reads into multiple files based on their names.
Opposite of muxbyname.

Usage:
demuxbyname.sh in=<file> in2=<file2> out=<outfile> out2=<outfile2> names=<string,string,string...>

Alternately:
demuxbyname.sh in=<file> out=<outfile> delimiter=whitespace prefixmode=f
This will demultiplex by the substring after the last whitespace.

demuxbyname.sh in=<file> out=<outfile> length=8 prefixmode=t
This will demultiplex by the first 8 characters of read names.

in2 and out2 are for paired reads and are optional.
If input is paired and there is only one output file, it will be written interleaved.

Parameters and their defaults:

in=<file>           Input file.
out=<file>          Output files for reads with matched headers (must contain % symbol).
outu=<file>         Output file for reads with unmatched headers.
prefixmode=t        (pm) Match prefix of read header.  If false, match suffix of read header.
substringmode=f     (substring) Names can be substrings of read headers.
names=              List of strings (or files containing strings) to parse from read names.
length=0            If positive, use a suffix or prefix of this length from read name instead of or in addition to the list of names.
                    For example, you could create files based on the first 8 characters of read names.
delimiter=          For prefix or suffix mode, specifying a delimiter will allow exact matches.
                    This allows demultiplexing based on names that are found without specifying a list of names.
                    Special delimiters are space, tab, and whitespace.  Otherwise, the delimiter will be used
                    as a literal string (for example, ':' or 'HISEQ').  This is similar to the length flag but allows variable length.

Common parameters:

ow=t                (overwrite) Overwrites files that already exist.
app=f               (append) Append to files that already exist.
zl=4                (ziplevel) Set compression level, 1 (low) to 9 (max).
int=auto            (interleaved) Determines whether INPUT file is considered interleaved.
qin=auto            ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.
qout=auto           ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input).
                    

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

z="-Xmx1200m"
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

function demuxbyname() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
		module load samtools
	fi
	local CMD="java $EA $z -cp $CP jgi.DemuxByName $@"
	echo $CMD >&2
	eval $CMD
}

demuxbyname "$@"

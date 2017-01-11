#!/bin/bash
#stats in=<infile>

usage(){
echo "
Written by Brian Bushnell
Last modified September 6, 2016

Description:  Generates basic assembly statistics such as scaffold count, 
              N50, L50, GC content, gap percent, etc.  For multiple files,
              please use statswrapper.sh.  Works with fasta and fastq only
              (gzipped is fine).

Usage:        stats.sh in=<file>


Parameters:
in=file         Specify the input fasta file, or stdin.
gc=file         Writes ACGTN content per scaffold to a file.
gchist=file     Filename to output scaffold gc content histogram.
shist=file      Filename to output cumulative scaffold length histogram.
gcbins=200      Number of bins for gc histogram.
n=10            Number of contiguous Ns to signify a break between contigs.
k=13            Estimate memory usage of BBMap with this kmer length.
minscaf=0       Ignore scaffolds shorter than this.
phs=f           (printheaderstats) Set to true to print total size of headers.
n90=t           (printn90) Print the N/L90 metrics.
extended=f      Print additional metrics such as N/L90 and log sum.
logoffset=1000  Minimum length for calculating log sum.
logbase=2       Log base for calculating log sum.
pdl=f           (printduplicatelines) Set to true to print lines in the 
                scaffold size table where the counts did not change.
n_=t            This flag will prefix the terms 'contigs' and 'scaffolds'
                with 'n_' in formats 3-6.
addname=f       Adds a column for input file name, for formats 3-6.

format=<0-7>    Format of the stats information; default 1.
	format=0 prints no assembly stats.
	format=1 uses variable units like MB and KB, and is designed for compatibility with existing tools.
	format=2 uses only whole numbers of bases, with no commas in numbers, and is designed for machine parsing.
	format=3 outputs stats in 2 rows of tab-delimited columns: a header row and a data row.
	format=4 is like 3 but with scaffold data only.
	format=5 is like 3 but with contig data only.
	format=6 is like 3 but the header starts with a #.
	format=7 is like 1 but only prints contig info.

gcformat=<0-4>  Select GC output format; default 1.
	gcformat=0:	(no base content info printed)
	gcformat=1:	name	length	A	C	G	T	N	GC
	gcformat=2:	name	GC
	gcformat=4:	name	length	GC
	Note that in gcformat 1, A+C+G+T=1 even when N is nonzero.

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

z="-Xmx120m"
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

stats() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z -cp $CP jgi.AssemblyStats2 $@"
#	echo $CMD >&2
	eval $CMD
}

stats "$@"

#!/bin/bash
#filterbyname in=<infile> out=<outfile>

function usage(){
echo "
Written by Brian Bushnell
Last modified December 18, 2015

Description:  Filters sequences by exact sequence matches.

Usage:  filterbysequence.sh in=<file> out=<file> ref=<file> include=<t/f>


I/O Parameters:
in=                 Primary input. 'in2' will specify a second file.
out=                Primary out. 'out2' will specify a second file.
ref=                A reference file or comma-delimited list of files.
literal=            A literal sequence or comma-delimited list of sequences.
ow=t                (overwrite) Overwrites files that already exist.
zl=2                (ziplevel) Set compression level, 1 (low) to 9 (max).

Processing Parameters:
include=f           Set to 'true' to include the filtered sequences rather
                    than excluding them.
rcomp=t             Match reverse complements as well.
casesensitive=f     (case) Require matching case.
storebases=t        (sb) Store ref bases.  Requires more memory.  If false,
                    case-sensitive matching cannot be done, and the matching
                    will be probabilistic based 128-bit hashcodes.
threads=auto        (t) Specify the number of worker threads.


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

z="-Xmx800m"
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
	freeRam 800m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

function filterbysequence() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
		module load samtools
	fi
	local CMD="java $EA $z -cp $CP jgi.FilterBySequence $@"
	echo $CMD >&2
	eval $CMD
}

filterbysequence "$@"

#!/bin/bash
#muxbyname in=<file,file,file...> out=<outfile>

usage(){
echo "
Written by Brian Bushnell
Last modified June 22, 2016

Description:  Multiplexes reads from multiple files after renaming them based on their initial file.
Opposite of demuxbyname.

Usage:  muxbyname.sh in=<file,file,file...> out=<output file>
Input files may also be given without an in= prefix, so that you can use wildcards:

muxbyname.sh *.fastq out=muxed.fastq


Standard parameters:
in=<file,file>  A list of input files.
in2=<file,file> Read 2 input if reads are in paired files.
out=<file>      Primary output, or read 1 output.
out2=<file>     Read 2 output if reads are in paired files.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
showspeed=t     (ss) Set to 'f' to suppress display of processing speed.
ziplevel=2      (zl) Set to 1 (lowest) through 9 (max) to change compression
                level; lower compression is faster.

Processing parameters:
None yet!

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

z="-Xmx400m"
z2="-Xms400m"
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

muxbyname() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z -cp $CP driver.RenameAndMux $@"
	echo $CMD >&2
	eval $CMD
}

muxbyname "$@"

#!/bin/bash
#taxonomy in=<infile> out=<outfile>

usage(){
echo "
Written by Brian Bushnell
Last modified September 6, 2016

Description:   Prints the full taxonomy of a string.
String may be a gi number, NCBI taxID, or Latin name.
An NCBI identifier should just be a number or ncbi|number.
A gi number should be gi|number.

Usage:  taxonomy.sh tree=<tree file> <identifier>
Alternate usage: taxonomy.sh tree=<tree file> in=<file>

Usage examples:
taxonomy.sh tree=tree.taxtree.gz homo_sapiens canis_lupus 9606
taxonomy.sh tree=tree.taxtree.gz gi=gitable.int1.d.gz in=refseq.fasta

Processing parameters:
in=<file>       A file containing named sequences, or just the names.
out=<file>      Output file.  If blank, use stdout.
tree=<file>     A taxonomic tree made by TaxTree, such as tree.taxtree.gz.
                On Genepool, you can set table=auto to use the default.
table=<file>    A table translating gi numbers to NCBI taxIDs.
                Only needed if gi numbers will be used.
                On Genepool, you can set table=auto to use the default.
level=null      Set to a taxonomic level like phylum to just print that level.
minlevel=-1     For multi-level printing, do not print levels below this.
maxlevel=life   For multi-level printing, do not print levels above this.
column=-1       If set to a non-negative integer, parse the taxonomy
                information from this column in a tab-delimited file.
                Example if column=1:
                read1 TAB gi|944259871|gb|KQL24128.1| TAB score:42
                becomes
                read1 TAB  k__Viridiplantae;p__Streptophyta;... TAB score:42

Parameters without an '=' symbol will be considered organism identifiers.

* Note *
Tree and table files are in /global/projectb/sandbox/gaag/bbtools/tax
For non-Genepool users, or to make new ones, use taxtree.sh and gitable.sh

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

z="-Xmx8g"
z2="-Xms8g"
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
	freeRam 2000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

taxonomy() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z -cp $CP tax.PrintTaxonomy $@"
	echo $CMD >&2
	eval $CMD
}

taxonomy "$@"

#!/bin/bash
#clumpify in=<infile> out=<outfile> groups=<number>

usage(){
echo "
Written by Brian Bushnell
Last modified January 5, 2017

Description:  Sorts sequences to put similar reads near each other.
Can be used for increased compression or error correction.
Note - Paired reads should be interleaved for both input and output,
not split into two files.

Usage:   clumpify.sh in=<file> out=<file> groups=<number>

Input may be fasta or fastq, compressed or uncompressed.

Optional parameters (and their defaults)

in=<file>           Input file.
out=<file>          Output file.  May not be standard out.
groups=auto         Use this many intermediate files (to save memory).
                    1 group is fastest.  Auto will estimate the number
                    of groups needed based on the file size.
rcomp=f             Give read clumps the same orientation to increase 
                    compression.  Should be disabled for paired reads.
overwrite=f         (ow) Set to false to force the program to abort rather 
                    than overwrite an existing file.
qin=auto            Auto-detect input quality encoding.  May be set to:
                       33:  ASCII-33 (Sanger) encoding.
                       64:  ASCII-64 (old Illumina) encoding.
                    All modern sequence is encoded as ASCII-33.
qout=auto           Use input quality encoding as output quality encoding.
changequality=f     (cq) If true, fix broken quality scores such as Ns with
                    Q>0.  Default is false to ensure lossless compression.

Compression parameters:
ziplevel=6          (zl) Gzip compression level (1-11).  Higher is slower.
                    Level 11 is only available if pigz is installed and is
                    extremely slow to compress, but faster to decompress.
                    Naming the output file to *.bz2 will use bzip2 instead of
                    gzip for ~95 additional compression, which requires
                    bzip2 or pbzip2 in the path.
blocksize=128       Size of blocks for pigz, in kb.  Higher gives slightly
                    better compression.
shortname=f         Make the names as short as possible.  'shortname=shrink'
                    will shorten the names where possible, but retain the 
                    flowcell and barcode information.
reorder=f           Reorder clumps for additional compression.  Only valid
                    when groups=1, passes=1, and ecc=f.  Possible modes:
                       f:  Do not reorder clumps.
                       c:  Reorder using consensus reads.  Uses additional
                           time and memory.
                       p:  Reorder using pair information.  Requires paired
                           reads.  Yields the highest compression.
                       a:  Automatically choose between 'c' and 'p'.  The
                           flag reorder with no argument will set 'reorder=a'.
quantize=f          Bin the quality scores, like NextSeq.  This greatly
                    increases compression, but information is lost.

Temp file parameters:
compresstemp=auto   (ct) Gzip temporary files.  By default temp files will be
                    compressed if the output file is compressed.
delete=t            Delete temporary files.
usetmpdir=f         Use tmpdir for temp files.
tmpdir=             By default, this is the environment variable TMPDIR.

Hashing parameters:
k=31                Use kmers of this length (1-31).  Shorter kmers may
                    increase compression, but 31 is recommended for error
                    correction.
mincount=0          Don't use pivot kmers with count less than this.
                    Setting mincount=2 can increase compression.
                    Increases time and memory usage.
seed=1              Random number generator seed for hashing.  
                    Set to a negative number to use a random seed.
hashes=4            Use this many masks when hashing.  0 uses raw kmers.
                    Often hashes=0 increases compression, but it should
                    not be used with error-correction.
border=1            Do not use kmers within this many bases of read ends.

Deduplication parameters:
dedupe=f            Remove duplicate reads.  For pairs, both must match.
markduplicates=f    Don't remove; just append ' duplicate' to the name.
allduplicates=f     Mark or remove all copies of duplicates, instead of
                    keeping the highest-quality copy.
dupesubs=2          (subs) Maximum substitutions allowed between duplicates.
optical=f           If true, mark or remove optical duplicates only.
                    This means they are Illumina reads within a certain
                    distance on the flowcell.  Normal Illumina names needed.
dupedist=40         (dist) Max distance to consider for optical duplicates.
                    Higher removes more duplicates but is more likely to
                    remove PCR rather than optical duplicates.
                    This is platform-specific; recommendations:
                       NextSeq      28
                       HiSeq 1T     40
                       HiSeq 2500   40
spantiles=t         Allow reads to be considered optical duplicates if they
                    are on different tiles, but are within dupedist in
                    either the x or y axis.
scanlimit=5         (scan) Continue for this many reads after encountering a
                    nonduplicate.  Improves detection of inexact duplicates.

Pairing/ordering parameters (for use with error-correction):
unpair=f            For paired reads, clump all of them rather than just
                    read 1.  Destroys pairing.  Without this flag, for paired
                    reads, only read 1 will be error-corrected.
repair=f            After clumping and error-correction, restore pairing.
                    If groups>1 this will sort by name which will destroy
                    clump ordering; with a single group, clumping will
                    be retained.

Error-correction parameters:
ecc=f               Error-correct reads.  Requires multiple passes for
                    complete correction.
ecco=f              Error-correct paired reads via overlap before clumping.
passes=1            Use this many error-correction passes.  6 passes are 
                    suggested.
consensus=f         Output consensus sequence instead of clumps.

Advanced error-correction parameters:
mincc=4             (mincountcorrect) Do not correct to alleles occuring less
                    often than this.
minss=4             (minsizesplit) Do not split into new clumps smaller than 
                    this.
minsfs=0.17         (minsizefractionsplit) Do not split on pivot alleles in
                    areas with local depth less than this fraction of clump size.
minsfc=0.20         (minsizefractioncorrect) Do not correct in areas with local
                    depth less than this.
minr=30.0           (minratio) Correct to the consensus if the ratio of the
                    consensus allele to second-most-common allele is >=minr,
                    for high depth.  Actual ratio used is:
                    min(minr, minro+minorCount*minrm+quality*minrqm).
minro=1.9           (minratiooffset) Base ratio.
minrm=1.8           (minratiomult) Ratio multiplier for secondary allele count.
minrqm=0.08         (minratioqmult) Ratio multiplier for base quality.
minqr=2.8           (minqratio) Do not correct bases when cq*minqr>rqsum.
minaqr=0.70         (minaqratio) Do not correct bases when cq*minaqr>5+rqavg.
minid=0.97          (minidentity) Do not correct reads with identity to 
                    consensus less than this.
maxqadjust=0        Adjust quality scores by at most maxqadjust per pass.
maxqi=-1            (maxqualityincorrect) Do not correct bases with quality 
                    above this (if positive).
maxci=-1            (maxcountincorrect) Do not correct alleles with count 
                    above this (if positive).
findcorrelations=t  Look for correlated SNPs in clumps to split into alleles.
maxcorrelations=12  Maximum number of eligible SNPs per clump to consider for
                    correlations.  Increasing this number can reduce false-
                    positive corrections at the possible expense of speed.

Java Parameters:
-Xmx                This will be passed to Java to set memory usage, overriding
                    the program's automatic memory detection.  -Xmx20g will 
                    specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.

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
	#Note that this uses slighly less (82%) of memory than normal to account for multiple pigz instances.
	freeRam 2000m 82
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

clumpify() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		#module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	java -version
	local CMD="java $EA $z $z2 -cp $CP clump.Clumpify $@"
	echo $CMD >&2
	eval $CMD
}

clumpify "$@"

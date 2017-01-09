#!/bin/bash
#bbduk2 in=<file> out=<file> fref=<file>

usage(){
echo "
Written by Brian Bushnell
Last modified June 27, 2016

BBDuk2 is like BBDuk but can kfilter, kmask, and ktrim in a single pass.
It does not replace BBDuk, and is only provided to allow maximally efficient
pipeline integration when multiple steps will be performed.  The syntax is 
slightly different.

Description:  Compares reads to the kmers in a reference dataset, optionally 
allowing an edit distance. Splits the reads into two outputs - those that 
match the reference, and those that don't. Can also trim (remove) the matching 
parts of the reads rather than binning the reads.

Usage:	bbduk2.sh in=<input file> out=<output file> fref=<contaminant files>

Input may be stdin or a fasta or fastq file, compressed or uncompressed.
If you pipe via stdin/stdout, please include the file type; e.g. for gzipped 
fasta input, set in=stdin.fa.gz


Input parameters:
in=<file>           Main input. in=stdin.fq will pipe from stdin.
in2=<file>          Input for 2nd read of pairs in a different file.
fref=<file,file>    Comma-delimited list of fasta reference files for filtering.
rref=<file,file>    Comma-delimited list of fasta reference files for right-trimming.
lref=<file,file>    Comma-delimited list of fasta reference files for left-trimming.
mref=<file,file>    Comma-delimited list of fasta reference files for masking.
fliteral=<seq,seq>  Comma-delimited list of literal sequences for filtering.
rliteral=<seq,seq>  Comma-delimited list of literal sequences for right-trimming.
lliteral=<seq,seq>  Comma-delimited list of literal sequences for left-trimming.
mliteral=<seq,seq>  Comma-delimited list of literal sequences for masking.
touppercase=f       (tuc) Change all bases upper-case.
interleaved=auto    (int) t/f overrides interleaved autodetection.
qin=auto            Input quality offset: 33 (Sanger), 64, or auto.
reads=-1            If positive, quit after processing X reads or pairs.
copyundefined=f     (cu) Process non-AGCT IUPAC reference bases by making all
                    possible unambiguous copies.  Intended for short motifs
                    or adapter barcodes, as time/memory use is exponential.

Output parameters:
out=<file>          (outnonmatch) Write reads here that do not contain 
                    kmers matching the database.  'out=stdout.fq' will pipe 
                    to standard out.
out2=<file>         (outnonmatch2) Use this to write 2nd read of pairs to a 
                    different file.
outm=<file>         (outmatch) Write reads here that contain kmers matching
                    the database.
outm2=<file>        (outmatch2) Use this to write 2nd read of pairs to a 
                    different file.
outs=<file>         (outsingle) Use this to write singleton reads whose mate 
                    was trimmed shorter than minlen.
stats=<file>        Write statistics about which contamininants were detected.
refstats=<file>     Write statistics on a per-reference-file basis.
rpkm=<file>         Write RPKM for each reference sequence (for RNA-seq).
dump=<file>         Dump kmer tables to a file, in fasta format.
nzo=t               Only write statistics about ref sequences with nonzero hits.
overwrite=t         (ow) Grant permission to overwrite files.
showspeed=t         (ss) 'f' suppresses display of processing speed.
ziplevel=2          (zl) Compression level; 1 (min) through 9 (max).
fastawrap=80        Length of lines in fasta output.
qout=auto           Output quality offset: 33 (Sanger), 64, or auto.
statscolumns=3      (cols) Number of columns for stats output, 3 or 5.
                    5 includes base counts.
rename=f            Rename reads to indicate which sequences they matched.
refnames=f          Use names of reference files rather than scaffold IDs.
trd=f               Truncate read and ref names at the first whitespace.
ordered=f           Set to true to output reads in same order as input.

Histogram output parameters:
bhist=<file>        Base composition histogram by position.
qhist=<file>        Quality histogram by position.
qchist=<file>       Count of bases with each quality value.
aqhist=<file>       Histogram of average read quality.
bqhist=<file>       Quality histogram designed for box plots.
lhist=<file>        Read length histogram.
gchist=<file>       Read GC content histogram.
gcbins=100          Number gchist bins.  Set to 'auto' to use read length.

Histograms for sam files only (requires sam format 1.4 or higher):

ehist=<file>        Errors-per-read histogram.
qahist=<file>       Quality accuracy histogram of error rates versus quality 
                    score.
indelhist=<file>    Indel length histogram.
mhist=<file>        Histogram of match, sub, del, and ins rates by read location.
idhist=<file>       Histogram of read count versus percent identity.
idbins=100          Number idhist bins.  Set to 'auto' to use read length.

Processing parameters:
k=27                Kmer length used for finding contaminants.  Contaminants 
                    shorter than k will not be found.  k must be at least 1.
rcomp=t             Look for reverse-complements of kmers in addition to 
                    forward kmers.
maskmiddle=t        (mm) Treat the middle base of a kmer as a wildcard, to 
                    increase sensitivity in the presence of errors.
minkmerhits=1       (mkh) Reads need at least this many matching kmers 
                    to be considered as matching the reference.
hammingdistance=0   (hdist) Maximum Hamming distance for ref kmers (subs only).
                    Memory use is proportional to (3*K)^hdist.
qhdist=0            Hamming distance for query kmers; impacts speed, not memory.
editdistance=0      (edist) Maximum edit distance from ref kmers (subs 
                    and indels).  Memory use is proportional to (8*K)^edist.
hammingdistance2=0  (hdist2) Sets hdist for short kmers, when using mink.
qhdist2=0           Sets qhdist for short kmers, when using mink.
editdistance2=0     (edist2) Sets edist for short kmers, when using mink.
forbidn=f           (fn) Forbids matching of read kmers containing N.
                    By default, these will match a reference 'A' if 
                    hdist>0 or edist>0, to increase sensitivity.
removeifeitherbad=t (rieb) Paired reads get sent to 'outmatch' if either is 
                    match (or either is trimmed shorter than minlen).  
                    Set to false to require both.
findbestmatch=f     (fbm) If multiple matches, associate read with sequence 
                    sharing most kmers.  Reduces speed.
skipr1=f            Don't do kmer-based operations on read 1.
skipr2=f            Don't do kmer-based operations on read 2.
ecco=f              For overlapping paired reads only.  Performs error-
                    correction with BBMerge prior to kmer operations.
recalibrate=f       (recal) Recalibrate quality scores.  Requires calibration
                    matrices generated by CalcTrueQuality.
sam=<file,file>     If recalibration is desired, and matrices have not already
                    been generated, BBDuk will create them from the sam file.

Speed and Memory parameters:
threads=auto        (t) Set number of threads to use; default is number of 
                    logical processors.
prealloc=f          Preallocate memory in table.  Allows faster table loading 
                    and more efficient memory usage, for a large reference.
monitor=f           Kill this process if it crashes.  monitor=600,0.01 would 
                    kill after 600 seconds under 1% usage.
minrskip=1          (mns) Force minimal skip interval when indexing reference 
                    kmers.  1 means use all, 2 means use every other kmer, etc.
maxrskip=1          (mxs) Restrict maximal skip interval when indexing 
                    reference kmers. Normally all are used for scaffolds<100kb, 
                    but with longer scaffolds, up to maxrskip-1 are skipped.
rskip=              Set both minrskip and maxrskip to the same value.
                    If not set, rskip will vary based on sequence length.
qskip=1             Skip query kmers to increase speed.  1 means use all.
speed=0             Ignore this fraction of kmer space (0-15 out of 16) in both
                    reads and reference.  Increases speed and reduces memory.
Note: Do not use more than one of 'speed', 'qskip', and 'rskip'.

Trimming/Filtering/Masking parameters:
Note - for BBDuk2, kmer filtering, trimming, and masking are independent,
and all can be performed at the same time.

ktrim=f             Trim reads to remove bases matching reference kmers.
                    Values: 
                            f (don't trim), 
                            r (trim to the right), 
                            l (trim to the left)
kmask=f             Replace bases matching ref kmers with another symbol.
                    Allows any non-whitespace character other than t or f,
                    and processes short kmers on both ends.  'kmask=lc' will
                    convert masked bases to lowercase.
mink=0              Look for shorter kmers at read tips down to this length, 
                    when k-trimming or masking.  0 means disabled.  Enabling
                    this will disable maskmiddle.
qtrim=f             Trim read ends to remove bases with quality below trimq.
                    Performed AFTER looking for kmers.
                    Values: 
                            rl (trim both ends), 
                            f (neither end), 
                            r (right end only), 
                            l (left end only),
                            w (sliding window)
trimq=6             Regions with average quality BELOW this will be trimmed.
minlength=10        (ml) Reads shorter than this after trimming will be 
                    discarded.  Pairs will be discarded if both are shorter.
mlf=0               (minlengthfraction) Reads shorter than this fraction of 
                    original length after trimming will be discarded.
maxlength=          Reads longer than this after trimming will be discarded.
                    Pairs will be discarded only if both are longer.
minavgquality=0     (maq) Reads with average quality (after trimming) below 
                    this will be discarded.
maqb=0              If positive, calculate maq from this many initial bases.
chastityfilter=f    (cf) Discard reads with id containing ' 1:Y:' or ' 2:Y:'.
barcodefilter=f     Remove reads with unexpected barcodes if barcodes is set,
                    or barcodes containing 'N' otherwise.  A barcode must be
                    the last part of the read header.
barcodes=           Comma-delimited list of barcodes or files of barcodes.
maxns=-1            If non-negative, reads with more Ns than this 
                    (after trimming) will be discarded.
mcb=0               (minconsecutivebases) Discard reads without at least 
                    this many consecutive called bases.
ottm=f              (outputtrimmedtomatch) Output reads trimmed to shorter 
                    than minlength to outm rather than discarding.
tp=0                (trimpad) Trim this much extra around matching kmers.
tbo=f               (trimbyoverlap) Trim adapters based on where paired 
                    reads overlap.
strictoverlap=t     Adjust sensitivity for trimbyoverlap mode.
minoverlap=14       Require this many bases of overlap for detection.
mininsert=50        Require insert size of at least this for overlap. 
                    Should be reduced to 16 for small RNA sequencing.
tpe=f               (trimpairsevenly) When kmer right-trimming, trim both 
                    reads to the minimum length of either.
forcetrimleft=0     (ftl) If positive, trim bases to the left of this position
                    (exclusive, 0-based).
forcetrimright=0    (ftr) If positive, trim bases to the right of this position
                    (exclusive, 0-based).
forcetrimright2=0   (ftr2) If positive, trim this many bases on the right end.
forcetrimmod=0      (ftm) If positive, right-trim length to be equal to zero,
                    modulo this number.
restrictleft=0      If positive, only look for kmer matches in the 
                    leftmost X bases.
restrictright=0     If positive, only look for kmer matches in the 
                    rightmost X bases.
mingc=0             Discard reads with GC content below this.
maxgc=1             Discard reads with GC content above this.
gcpairs=t           Use average GC of paired reads.
                    Also affects gchist.

Entropy/Complexity parameters:
entropy=-1          Set between 0 and 1 to filter reads with entropy below
                    that value.  Higher is more stringent.
entropywindow=50    Calculate entropy using a sliding window of this length.
entropyk=5          Calculate entropy using kmers of this length.
minbasefrequency=0  Discard reads with a minimum base frequency below this.

Cardinality estimation:
cardinality=f           (loglog) Count unique kmers using the LogLog algorithm.
loglogk=31              Use this kmer length for counting.
loglogbuckets=1999      Use this many buckets for counting.

Java Parameters:

-Xmx                This will be passed to Java to set memory usage, overriding 
                    the program's automatic memory detection. -Xmx20g will 
                    specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  
                    The max is typically 85% of physical memory.

There is a changelog at /bbmap/docs/changelog_bbduk.txt
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
NATIVELIBDIR="$DIR""jni/"

z="-Xmx1g"
z2="-Xms1g"
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
	freeRam 1400m 42
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

bbduk2() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
		module load samtools
	fi
	local CMD="java -Djava.library.path=$NATIVELIBDIR $EA $z $z2 -cp $CP jgi.BBDuk2 $@"
	echo $CMD >&2
	eval $CMD
}

bbduk2 "$@"

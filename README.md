# omega3

Omeag3 is a multi threaded memory efficient overlap-layout-consensus (OLC) metagenome assembler - **omega3**. 

### Setup and Installation
1. Download the tarball with compiled executables for Linux or the source code at: [https://github.com/abiswas-odu/Omega3](https://github.com/abiswas-odu/Omega3). The code has been tested on both Linux and MacOS systems, but not under Windows.
2. If you decide to sownload the executable, type `make` to build.
3. If compiled successfully, the executabled required will be generated. The assembler is invoked through a run script `./runOmega3.sh`. Use `./runOmega3.sh -h` for help information.

### Features

* Quick summary
**omega3** is a massively improved multi-threaded memory efficient version of [omega](http://bioinformatics.oxfordjournals.org/content/early/2014/07/06/bioinformatics.btu395.short), its unique capabilities include:

    1. Modularization of contained and duplicate reads removal, and initial graph construction after transitive edge reduction step: Big data set can be processed in chunks so that memory limitation problem is solved.

    2. More effective and accurate graph reduction steps result in shorter running time and much higher assembly quality.

    3. The assembler includes a paired end read based scaffolding step. 
    
 ### Current Version
* v3.0


### Assembly of Metagenomic sequencing reads from raw Illumina data

#### Preprocessing of the Illumina data

Since omega3 works best with reads without errors, preprocessing plays an important role in deciding the quality of the assembly results. We have 3 basc pre processing steps. Trimming, filtering, and eror correction.

##### Trimming, filtering, (merging), and eror correction

We have tested Brian Bushnell's suite of tools [BBTools](http://sourceforge.net/projects/bbmap/files/) extensively on Illumina data and have obtained good results. Suppose the Illumina reads data set is called `$reads`, the steps we recommend are following:

```
#!sh

# Use bbduk.sh to quality and length trim the Illumina reads and remove adapter sequences
# 1. ftm = 5, right trim read length to a multiple of 5
# 2. k = 23, Kmer length used for finding contaminants
# 3. ktrim=r, Trim reads to remove bases matching reference kmers to the right
# 4. mink=11, look for shorter kmers at read tips down to 11 bps
# 5. qhdist=1, hamming distance for query kmers
# 6. tbo, trim adapters based on where paired reads overlap
# 7. tpe, when kmer right-trimming, trim both reads to the minimum length of either
# 8. qtrim=r, trim read right ends to remove bases with low quality
# 9. trimq=10, regions with average quality below 10 will be trimmed.
# 10. minlength=70, reads shorter than 70bps after trimming will be discarded.
# 11. ref=$adapters, adapters shipped with bbnorm tools
# 12. –Xmx8g, use 8G memory
# 13. 1>trim.o 2>&1, redirect stderr to stdout, and save both to file *trim.o*
adapters= bbmap_dir/resources/adapters.fa
phiX_adapters= bbmap_dir/resources/phix174_ill.ref.fa.gz
bbduk.sh in=$reads out=trim.fq.gz ftm=5 k=23 ktrim=r mink=11 qhdist=1 tbo tpe qtrim=r trimq=10 minlength=70 ref=$adapters -Xmx8g 1>trim.o 2>&1
bbduk.sh in=trim.fq.gz out=filter.fq.gz ref=$phiX_adapters hdist=1 k=31 threads=48
```

##### Error correction with Tadpole

Tarpole is a memory efficient error correction tool that runs within reasonable time. We suggest using this for error correction. 

```
#!bash
# 1. mode=correct, use tadpole for correction
# 2. ecc=t, error correct via kmer counts
# 3. shave, remove dead ends in the kmer graph
# 4. rinse, remove bubbles
# 5. -Xmx120g, use 120G memory
# 6. markbadbases=2, change bases fully covered by less than 2 kmers to N
tadpole.sh in=merged.fq.gz out=merged_EC.fq.gz mode=correct markbadbases=2 ecc=t shave rinse -Xmx120g &> ecc.out
# 1. maxns=0, reads more than 0 Ns after trimming wil be discarded.
reformat.sh in=merged_EC.fq.gz out=merged_EC_reformat.fq.gz maxns=0 -Xmx120g 1>> ecc.out 2>&1
```

### Assembly of Error Corrected Data

The omega3 assembler is invoked through the run script `./runOmega3.sh`. The basic quick start commands with default parameters are as follows. The default parameters are based on empherical tests on real metagenomic datasets.     

```
#!/bin/bash

# Seperated paired end reads
runOmega3.sh -d ${output_directory} -in1 {read_1.fastq}  -in2 ${read2_2.fastq} -n ${num_threads} -m {max_mem_usage} -o ${64gen} 

# Interleaved paired end reads
runOmega3.sh -d ${output_directory} -inP {read_P.fastq} -n ${num_threads} -m {max_mem_usage} -o ${64gen}

# Single end reads
runOmega3.sh -d ${output_directory} -inS {read.fastq} -n ${num_threads} -m {max_mem_usage} -o ${64gen} 
```
For all the options of omega3, use `./runOmega3.sh -h`

#### Setting assembly parameters

The assembly parameters can be modified to attempt better assembly. This can be done through a parameter file passed using the `-p` parameter to the run script. 

The configurable parameters include the following:

```
##########################################
### Assembly configuration file for Omega3
##########################################

# Minimum overlap length to create an edge
minOvl=40

# minimum reads count in an edge to be not dead end edge (default: 20)
minReadsCountInEdgeToBeNotDeadEnd = 20
# minimum edge length to be not dead end edge (default: 1000)
minEdgeLengthToBeNotDeadEnd = 1000

# minimum reads count for an edge to be kept even if it has 0 flow (default: 15)
minReadsCountToHave0Flow = 15
# minimum edge length or an edge to be kept even if it has 0 flow (default: 1500)
minEdgeLengthToHave0Flow = 1500

# minimum reads count in an edge to be 1 minimum flow (default: 10)
minReadsCountInEdgeToBe1MinFlow = 10
# minimum edge length to be 1 minimum flow (default: 1000)
minEdgeLengthToBe1MinFlow = 1000

# minimum overlap length difference to clip branches (default: 15)
minOvlDiffToClip = 15

# minimum fold difference to consider branches to be short (default: 5)
minFoldToBeShortBranch = 5

# minumum unique mate pair support to join edge (default: 3)
minUinqSupport=3
# minumum non-unique mate pair support to join edge (default: 0)
minNonUniqSupport=0

# minimum contig length to be reported (default: 1000)
minContigLengthTobeReported = 1000
```
#### Omega3 Assembler Output

Please see the OUTPUT.md file for description of the output files.  

### Dependencies

* C++11, gcc4.9+

### Questions?

* [Abhishek Biswas](mailto:ab.prof@gmail.com)
* [Chongle Pan](mailto:chongle.pan@gmail.com)

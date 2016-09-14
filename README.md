# Omega3

Omega3 is a multi threaded and multiprocess distributed memory overlap-layout-consensus (OLC) metagenome assembler - **Omega3**. 

### Setup and Installation

### Dependencies

1. GNU GCC with C++11 support i.e. gcc4.9+ or above
2. MPI Library with MPI-3 support i.e. OpenMPI 1.8 and above or cray-mpich/7.4.0 and above. By default the mpic++ wrapper is needed. If you are on a Cray cluster and the wrapper is "CC". You will need to edit the compiler.mk file. Uncomment the line "CC := CC" and comment out "CC := mpic++".   
 
### Installation Steps
1. Download the tarball with compiled executables for Linux or the source code at: [https://github.com/abiswas-odu/Omega3](https://github.com/abiswas-odu/Omega3). The code has been tested on both Linux and MacOS systems, but not under Windows.
2. If you decide to download the executable, type `make` to build.
3. If compiled successfully, the required executables will be built. 

### Quickly Running The Assembler

There are two basic versions of the assembler one for running on a single machine and another for running with MPI on a cluster.  

* __Single Machine Version:__ This version of the assembler should be used if you are going to run the assembler on a single machine with one or more cores. The assembler is invoked through a run script `./runOmega3.sh`. Make sure the RAM on the machine is more than the disk space size of the reads. The quick start command as shown below will be used in a batch job submission script or directly typed on the commandline terminal.   

```
#!/bin/bash

# Seperated paired end reads
runOmega3.sh -d ${output_directory} -in1 {read_1.fastq}  -in2 ${read2_2.fastq} -n ${num_threads} -o ${OUTPUT_DIR} 

```
Use `./runOmega3.sh -h` for help information.

* __MPI Version:__ This version of the assembler should be used if you are going to run the assembler with MPI support on a cluster. The run script to invoke the assembler depends on the cluster management and job scheduling system.

	1. If you have ORTE i.e. __mpirun__ is avilable, invoke the assembler using the run script `runOmega3-MPI.sh`. 
	2. If you have SLRUM i.e. __srun__ is available, invoke the assembler using the run script `runOmega3-MPI-SLRUM.sh`.
	3. If you have ALPS i.e. __aprun__ is available, invoke the assembler using the run script `runOmega3-MPI-ALPS.sh`.
 
For the basic MPI version make sure the RAM on the nodes is more than the disk space size of the reads. If you have a large dataset, then use the Remote Memory Access (RMA) version. The RMA version of the assembler will equally distribute about 70% of the memory usage across all the MPI nodes. The quick start commands are:
```
#!/bin/bash

### MPI Verion 
### Seperated paired end reads
runOmega3-MPI.sh -d ${output_directory} -in1 {read_1.fastq}  -in2 ${read2_2.fastq} -n ${num_threads} -o ${OUTPUT_DIR} 

### MPI Remote Memory Access(RMA) Verion 
### Seperated paired end reads
runOmega3-MPI.sh -d ${output_directory} -in1 {read_1.fastq}  -in2 ${read2_2.fastq} -n ${num_threads} -o ${OUTPUT_DIR} -rma 

```
Use `runOmega3-MPI.sh -h` for help information.

### Features

* Quick summary
**omega3** is a massively improved multi-threaded, multi-process distributed memory version of [omega](http://bioinformatics.oxfordjournals.org/content/early/2014/07/06/bioinformatics.btu395.short), its unique capabilities include:

    1. Modularization of contained and duplicate reads removal, and initial graph construction after transitive edge reduction step: Big data set can be processed in chunks so that memory limitation problem is solved.

    2. More effective and accurate graph reduction steps result in shorter running time and much higher assembly quality.

    3. The assembler includes a paired end read based scaffolding step. 
    
 ### Current Version
* v3.0.1

### Assembly of Metagenomic sequencing reads from raw Illumina data

#### Preprocessing of the Illumina data

Since Omega3 works best with reads without errors, preprocessing plays an important role in deciding the quality of the assembly results. The 3 basc pre-processing steps are trimming, filtering and eror correction.

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

#### Assembly on a Single Node

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

In case the program crashes due to exceeding wall clock time, the assembler can be restarted with the same command. 

#### Controlling memory usage

The memory usage of Omega can be controlled using the `-m` option to the run script as shown above. The default memory usage is to take all the system resources. In case that has to be avoided or the program crashes ot is too slow due to memory page swapping, the user can set a ubber bound on the memory. The minumum memory to assemble a dataset is:

```
Min Required Memory (GB) = (Disk Space of Reads) + (1GB * num_threads)
``` 
The program will run faster if more memory is made available.

#### Restarting Omega3 for repeat assembly and handling assembly crashes

Omega3 assembler can be restarted with changed assembly and scaffolding parameters using the `-osg` option. Setting this option while invoking `runOmega3.sh` will reuse the overlap graph constructed earlier and only perform the graph simplification step. This will significantly reduce executime time of assemblies on the same dataset with different parameters.    

Omega3 assembler can also be restarted after a crash caused due to exceeding wall clock time or out of memory errors. The job must be restarted with the same command as before and Omega3 will attempt to continue the assembly. Do not set the `-osg` option in this case.   

#### Setting assembly parameters

The assembly parameters can be modified to attempt better assembly. This can be done through a parameter file passed using the `-p` parameter to the run script. 

The configurable parameters include the following:

```
##############################################################
###   Assembly and scaffolding configurations for Omega3   ###
##############################################################

#### Parameters for building an overlap graph ####

# Minimum overlap length (bp) required to insert an edge between two reads during graph construction 
# Increase this to reduce N50 and mis-assemblies
MinOverlap4BuildGraph = 40  


###################################################


#### Parameters for simplifying an overlap graph ####
# You can run graph simplification using different settings below on the same assembly graph without re-doing graph construction.

# Parameters for Omega output

# Print contigs or not. Options are "false" (default) or "true". Printing takes a non-trivial amount of wall-clock time.
PrintContigs = false
# Print scaffolds or not. Options are "true" (default) or "false"
PrintScaffolds = true
# Minimum length of contigs or scaffolds to be printed (default: 1000 bp)
minSequenceLengthTobePrinted = 1000


# Minimum overlap length (bp) required to keep an edge between two reads during graph simplification
# This minimum overlap length must be equal to (Default) or larger than the MinOverlap4BuildGraph above
# This allows you to try different minimum overlap lengths for assembly without re-doing assembly graph construction
# Edges with shorter overlap length than this parameter will be ignored during graph simplification
# Increase this to reduce N50 and mis-assemblies
MinOverlap4SimplifyGraph = 40

# Minimum overlap length difference (bp) to clip branches (default: 25 bp)
# If a read has multiple edges, Omega clips the branches with overlap lengths less than the largest overlap of this read by this difference or more.
# Increase this to reduce N50 and mis-assemblies
minOverlapDifference4ClipBranches = 25 


# Parameters for joins edges or scaffolding edges using paired-end information 
# Minumum number of paired-end reads that provide unique support to merge two edges (default: 3)
# Increase this to reduce N50 and mis-assemblies
minUniquePEsupport=3  
# Minumum number of paired-end reads that provide non-unique support to merge two edges (default: 0)
minNonUniquePEsupport=0


# Parameters for dead-end edge removal

# Minimum number of reads in an edge to be not dead-end edge (default: 10)
minReadsCountInEdgeToBeNotDeadEnd = 10
# Minimum edge length (bp) to be not dead-end edge (default: 1000)
minEdgeLengthToBeNotDeadEnd = 1000
# Minimum fold difference between two branches' lengths to consider a branch to be short (default: 5)
minFoldToBeShortBranch = 5

# Parameters for flow analysis

# Minimum number of reads for an edge to be kept even if it has 0 flow (default: 15)
minReadsCountToHave0Flow = 15
# Minimum edge length for an edge to be kept even if it has 0 flow (default: 1500)
minEdgeLengthToHave0Flow = 1500
# Minimum number of reads in an edge to be assigned with 1 minimum flow (default: 20)
minReadsCountInEdgeToBe1MinFlow = 20
# Minimum edge length to be assigned with 1 minimum flow (default: 2000 bp)
minEdgeLengthToBe1MinFlow = 2000

```
#### Omega3 Assembler Output

Please see the OUTPUT.md file for description of the output files.  

### Questions?

* [Abhishek Biswas](mailto:ab.prof@gmail.com)
* [Chongle Pan](mailto:chongle.pan@gmail.com)

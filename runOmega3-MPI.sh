#!/bin/bash 

#Get path of run script
if [ -L $0 ] ; then
    exePath=$(dirname $(readlink -f $0)) ;
else
    exePath=$(dirname $0) ;
fi;

#Basic parameters defaults
numThreads=`nproc`
numProcs=3
readFile1=""
readFile2=""
readFileS=""
readFileP=""
asmParaFileP="${exePath}/omega3.cfg"
asmParaFileP2="${exePath}/omega3_2.cfg"
asmParaFileP3="${exePath}/omega3_3.cfg"
constructGraph="Y"
simplifyGraph="Y"
phymem=`grep MemTotal /proc/meminfo | awk '{print $2}'`
maxMem=$((phymem/1048576))
while [[ $# > 0 ]]
do
key="$1"
case $key in
    -h|--help)              # help output
    echo -e "Usage:\n"
    echo -e "   runOmega3.sh [OPTION]...<PARAM>...\n\n"
    echo -e "<PARAMS>\n"
    echo -e "   -inS\t single read filename.\n"
    echo -e "   -in1\t forward paired read filename.\n"
    echo -e "   -in2\t reverse paired read filename.\n"
    echo -e "   -inP\t interleaved paired read filename.\n"
    echo -e "   -d\t output directory path.\n"
    echo -e "   -o\t output filename prefix.\n"
    echo -e "   -p\t assembly parameter file.\n"
    echo -e "<OPTIONS>\n"
    echo -e "   -h\t help.\n"
    echo -e "   -m\t maximum memory to be used (DEFAULT: $maxMem GB).\n"
    echo -e "   -np\t number of MPI processes (DEFAULT: $numProcs).\n"
    echo -e "   -n\t number of threads (DEFAULT: $numThreads).\n"
    echo -e "   -obg\t only build overlap graph (DEFAULT: False).\n"
    echo -e "   -osg\t only simplify existing overlap graph (DEFAULT: False).\n"
    exit 1
    ;;
    -p|--parameterFile)     # overlap length to consider
    asmParaFileP="$2"
    shift # past argument
    ;;
    -p2|--parameterFile2)     # overlap length to consider
    asmParaFileP2="$2"
    shift # past argument
    ;;
    -p3|--parameterFile3)     # overlap length to consider
    asmParaFileP3="$2"
    shift # past argument
    ;;
    -d|--outdir)			# output directory
    dataOutPath="$2"
    shift # past argument
    ;;
    -in1)					# Forward paired end read file
    readFile1="$2"
    shift # past argument
    ;;
    -in2)					# Reverse paired end read file
    readFile2="$2"
    shift # past argument
    ;;
    -inS)					# Single read file
    readFileS="$2"
    shift # past argument
    ;;
    -inP)					# Interleaved paired end read file
    readFileP="$2"
    shift # past argument
    ;;
    -o|--outprefix)			# Output file prefix
    outPrefix="$2"
    shift # past argument
    ;;
    -n|--numthreads)		# Threads to use
    numThreads="$2"
    shift # past argument
    ;;
    -np|--numprocs)		# Threads to use
    numProcs="$2"
    shift # past argument
    ;;
    -m|--maxmem)			# Maximum memory limit
    maxMem="$2"
    shift # past argument
    ;;
    -obg|--OnlyBuildGraph)                        # Maximum memory limit
    simplifyGraph=""
    ;;
    -osg|--OnlySimplifyGraph)                        # Maximum memory limit
    constructGraph=""
    ;;
    *)
    echo "ERROR: Unidentified user variable $key"
    exit 1        				# unknown option
    ;;
esac
shift # past argument or value
done

if [ -z "$dataOutPath" ] && [ -z "$outPrefix" ] && [ -z "$readFile1" ] && [ -z "$readFile2" ] && [ -z "$readFileS" ] && [ -z "$readFileP" ] ; then
   echo "Input incomplete. Not all required parameters specified. Run with -h for help. Exiting..."
   exit 1
fi

#Check required exe files are in the same directory as the script
BUILDGEXE="${exePath}/buildG"
if [ -f $BUILDGEXE ] ; then
   echo "Graph construction module $BUILDGEXE exists."
else
   echo "Graph construction module $BUILDGEXE does not exist in script directory."
fi

PARSIMPLIFYEXE="${exePath}/parsimplify"
if [ -f $PARSIMPLIFYEXE ] ; then
   echo "Partial graph simplification module $PARSIMPLIFYEXE exists."
else
   echo "Partial graph simplification module $PARSIMPLIFYEXE does not exist in script directory."
fi

SIMPLIFYEXE="${exePath}/fullsimplify"
if [ -f $SIMPLIFYEXE ] ; then
   echo "Graph simplification module $PARSIMPLIFYEXE exists."
else
   echo "Graph simplification module $PARSIMPLIFYEXE does not exist in script directory."
fi
#All exe files in place


#check if output directory exists, if not create it
if [ -d $dataOutPath ] ; then
   echo "Output directory exists."
else
   echo "Cresting output directory: $dataOutPath"
   `mkdir $dataOutPath`
fi
 
outGraphPrefix="${dataOutPath}/graph/${outPrefix}"
outSimplifyPrefix="${dataOutPath}/assembly/${outPrefix}"

logFile="${dataOutPath}/${outPrefix}.log"
#Create edge file list
i=0
while [ $i -lt $numThreads ]
do
   edgeFiles="${edgeFiles}${dataOutPath}/graph/${outPrefix}_${i}_parGraph.txt,"
   i=$(( $i+1 ))
done
edgeFiles=${edgeFiles%?}

#Create contained read file list
i=0
while [ $i -lt $numThreads ]
do
   containedReads="${containedReads}${dataOutPath}/graph/${outPrefix}_${i}_containedReads.txt,"
   i=$(( $i+1 ))
done
containedReads=${containedReads%?}


echo Starting Time is $(date)

#create output directories
if [ -d "${dataOutPath}/graph" ] ; then
	echo "Previous result directory \"graph\" exists. Will keep contained read data and partial graphs from previous run."
else
	`mkdir ${dataOutPath}/graph`
fi

if [ -d "${dataOutPath}/assembly" ] ; then
   echo "Previous result directory \"assembly\" exists. Will nuke previous run..."
else
   `mkdir ${dataOutPath}/assembly`
fi

#Build graph
if [ -z "$readFile1" ] && [ -z "$readFile2" ] && [ -z "$readFileS" ] && [ -z "$readFileP" ] ; then
   echo "No input files specified. Exiting..."
   exit 1
elif [ -z "$readFile1" ] && [ -z "$readFile2" ] && [ -z "$readFileS" ] ; then 
   if [ "$constructGraph" = "Y" ] ; then
      mpirun -np $numProcs ${exePath}/buildG -pe 1 ${readFileP} -f $outGraphPrefix -p ${asmParaFileP} -t ${numThreads} -m ${maxMem}  > ${logFile}
   fi
   if [ "$simplifyGraph" = "Y" ] ; then
      ${exePath}/fullsimplify -fpi ${readFileP} -e ${edgeFiles} -crd ${containedReads} -simPth ${exePath} -p ${asmParaFileP} -p2 ${asmParaFileP2} -p3 ${asmParaFileP3} -o $outSimplifyPrefix -t ${numThreads} -log DEBUG4 >> ${logFile} 2>&1
      `cat ${outSimplifyPrefix}_contigsFinal_1.fasta ${outSimplifyPrefix}_contigsFinal_2.fasta ${outSimplifyPrefix}_contigsFinal_3.fasta > ${outSimplifyPrefix}_contigsFinalCombined.fasta`
      `cat ${outSimplifyPrefix}_scaffoldsFinal_1.fasta ${outSimplifyPrefix}_scaffoldsFinal_2.fasta ${outSimplifyPrefix}_scaffoldsFinal_3.fasta > ${outSimplifyPrefix}_scaffoldsFinalCombined.fasta`
      `cp ${outSimplifyPrefix}_contigsFinalCombined.fasta ${dataOutPath}`
      `cp ${outSimplifyPrefix}_scaffoldsFinalCombined.fasta ${dataOutPath}`
   fi
elif [ -z "$readFile1" ] && [ -z "$readFile2" ] && [ -z "$readFileP" ] ; then
   if [ "$constructGraph" = "Y" ] ; then
      mpirun -np $numProcs ${exePath}/buildG -se 1 ${readFileS} -f $outGraphPrefix -p ${asmParaFileP} -t ${numThreads} -m ${maxMem}  &> ${logFile}
   fi
   if [ "$simplifyGraph" = "Y" ] ; then
      ${exePath}/fullsimplify -fs ${readFileS} -e ${edgeFiles} -crd ${containedReads} -simPth ${exePath} -p ${asmParaFileP} -p2 ${asmParaFileP2} -p3 ${asmParaFileP3} -o $outSimplifyPrefix -t ${numThreads} -log DEBUG4 >> ${logFile} 2>&1
      `cat ${outSimplifyPrefix}_contigsFinal_1.fasta ${outSimplifyPrefix}_contigsFinal_2.fasta ${outSimplifyPrefix}_contigsFinal_3.fasta > ${outSimplifyPrefix}_contigsFinalCombined.fasta`
      `cat ${outSimplifyPrefix}_scaffoldsFinal_1.fasta ${outSimplifyPrefix}_scaffoldsFinal_2.fasta ${outSimplifyPrefix}_scaffoldsFinal_3.fasta > ${outSimplifyPrefix}_scaffoldsFinalCombined.fasta`
      `cp ${outSimplifyPrefix}_contigsFinalCombined.fasta ${dataOutPath}`
      `cp ${outSimplifyPrefix}_scaffoldsFinalCombined.fasta ${dataOutPath}`
   fi
elif [ -z "$readFileS" ] && [ -z "$readFileP" ] ; then
   if [ "$constructGraph" = "Y" ] ; then
      mpirun -np $numProcs ${exePath}/buildG -pe 2 ${readFile1} ${readFile2} -f $outGraphPrefix -p ${asmParaFileP} -t ${numThreads} -m ${maxMem}  &> ${logFile}
   fi
   if [ "$simplifyGraph" = "Y" ] ; then
      ${exePath}/fullsimplify -fp ${readFile1},${readFile2} -e ${edgeFiles} -crd ${containedReads} -simPth ${exePath} -p ${asmParaFileP} -p2 ${asmParaFileP2} -p3 ${asmParaFileP3} -o $outSimplifyPrefix -t ${numThreads} -log DEBUG4 >> ${logFile} 2>&1
      `cat ${outSimplifyPrefix}_contigsFinal_1.fasta ${outSimplifyPrefix}_contigsFinal_2.fasta ${outSimplifyPrefix}_contigsFinal_3.fasta > ${outSimplifyPrefix}_contigsFinalCombined.fasta`
      `cat ${outSimplifyPrefix}_scaffoldsFinal_1.fasta ${outSimplifyPrefix}_scaffoldsFinal_2.fasta ${outSimplifyPrefix}_scaffoldsFinal_3.fasta > ${outSimplifyPrefix}_scaffoldsFinalCombined.fasta`
      `cp ${outSimplifyPrefix}_contigsFinalCombined.fasta ${dataOutPath}`
      `cp ${outSimplifyPrefix}_scaffoldsFinalCombined.fasta ${dataOutPath}`
   fi
elif [ -z "$readFile1" ] && [ -z "$readFile2" ] ; then
   if [ "$constructGraph" = "Y" ] ; then
      mpirun -np $numProcs ${exePath}/buildG -pe 1 ${readFileP} -se 1 ${readFileS} -f $outGraphPrefix -p ${asmParaFileP} -t ${numThreads} -m ${maxMem}  &> ${logFile}
   fi
   if [ "$simplifyGraph" = "Y" ] ; then
      ${exePath}/fullsimplify -fpi ${readFileP} -fs ${readFileS} -e ${edgeFiles} -crd ${containedReads} -simPth ${exePath} -p ${asmParaFileP} -p2 ${asmParaFileP2} -p3 ${asmParaFileP3} -o $outSimplifyPrefix -t ${numThreads} -log DEBUG4 >> ${logFile} 2>&1
      `cat ${outSimplifyPrefix}_contigsFinal_1.fasta ${outSimplifyPrefix}_contigsFinal_2.fasta ${outSimplifyPrefix}_contigsFinal_3.fasta > ${outSimplifyPrefix}_contigsFinalCombined.fasta`
      `cat ${outSimplifyPrefix}_scaffoldsFinal_1.fasta ${outSimplifyPrefix}_scaffoldsFinal_2.fasta ${outSimplifyPrefix}_scaffoldsFinal_3.fasta > ${outSimplifyPrefix}_scaffoldsFinalCombined.fasta`
      `cp ${outSimplifyPrefix}_contigsFinalCombined.fasta ${dataOutPath}`
      `cp ${outSimplifyPrefix}_scaffoldsFinalCombined.fasta ${dataOutPath}`
   fi
elif [ -z "$readFileP" ] ; then
   if [ "$constructGraph" = "Y" ] ; then
      mpirun -np $numProcs ${exePath}/buildG -pe 2 ${readFile1} ${readFile2} -se 1 ${readFileS} -f $outGraphPrefix -p ${asmParaFileP} -t ${numThreads} -m ${maxMem}  &> ${logFile}
   fi
   if [ "$simplifyGraph" = "Y" ] ; then
      ${exePath}/fullsimplify -fp ${readFile1},${readFile2} -fs ${readFileS} -e ${edgeFiles} -crd ${containedReads} -simPth ${exePath} -p ${asmParaFileP} -p2 ${asmParaFileP2} -p3 ${asmParaFileP3} -o $outSimplifyPrefix -t ${numThreads} -log DEBUG4 >> ${logFile} 2>&1
      `cat ${outSimplifyPrefix}_contigsFinal_1.fasta ${outSimplifyPrefix}_contigsFinal_2.fasta ${outSimplifyPrefix}_contigsFinal_3.fasta > ${outSimplifyPrefix}_contigsFinalCombined.fasta`
      `cat ${outSimplifyPrefix}_scaffoldsFinal_1.fasta ${outSimplifyPrefix}_scaffoldsFinal_2.fasta ${outSimplifyPrefix}_scaffoldsFinal_3.fasta > ${outSimplifyPrefix}_scaffoldsFinalCombined.fasta`
      `cp ${outSimplifyPrefix}_contigsFinalCombined.fasta ${dataOutPath}`
      `cp ${outSimplifyPrefix}_scaffoldsFinalCombined.fasta ${dataOutPath}`
   fi
else
   echo "Invalid combination of input files. You can specify oe single end file with one interleaved paired file or speprate paired file. Exiting..."
   exit 1
fi

echo Ending Time is $(date)

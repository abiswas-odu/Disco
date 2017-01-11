#!/bin/bash 

#Get path of run script
if [ -L $0 ] ; then
    exePath=$(dirname $(readlink -f $0)) ;
else
    exePath=$(dirname $0) ;
fi;

#Basic parameters defaults
numThreads=`nproc`
numProcs=$PBS_NUM_NODES
readFile1=""
readFile2=""
readFileS=""
readFileP=""
outPrefix="disco"
asmParaFileP="${exePath}/disco.cfg"
asmParaFileP2="${exePath}/disco_2.cfg"
asmParaFileP3="${exePath}/disco_3.cfg"
constructGraph="Y"
simplifyGraph="Y"
useMPIRMA=""
dataOutPath="."
phymem=`grep MemTotal /proc/meminfo | awk '{print $2}'`
maxMem=$((phymem/1048576))
while [[ $# > 0 ]]
do
key="$1"
case $key in
    -h|--help)              # help output
    echo -e "Usage: For clusters using PBS with OpenMPI's mpirun\n"
    echo -e "   runDisco.sh [OPTION]...<PARAM>...\n\n"
    echo -e "<PARAMS>\n"
    echo -e "   -inS\t single read filenames (comma seperated fasta/fastq/fastq.gz files).\n"
    echo -e "   -in1\t forward paired read filename (single fasta/fastq/fastq.gz file).\n"
    echo -e "   -in2\t reverse paired read filename (single fasta/fastq/fastq.gz file).\n"
    echo -e "   -inP\t interleaved paired read filenames (comma seperated fasta/fastq/fastq.gz files).\n"
    echo -e "   -d\t output directory path.\n"
    echo -e "   -o\t output filename prefix.\n"
    echo -e "<OPTIONS>\n"
    echo -e "   -h\t help.\n"
    echo -e "   -m\t maximum memory to be used (DEFAULT: $maxMem GB).\n"
    echo -e "   -np\t number of MPI processes (DEFAULT: $numProcs).\n"
    echo -e "   -n\t number of threads (DEFAULT: $numThreads).\n"
    echo -e "   -obg\t only build overlap graph (DEFAULT: False).\n"
    echo -e "   -osg\t only simplify existing overlap graph (DEFAULT: False).\n"
    echo -e "   -rma\t enable distributed memory with remote memory acess (DEFAULT: False).\n"
    echo -e "   -p\t assembly parameter file for 1st assembly iteration.\n"
    echo -e "   -p2\t assembly parameter file for 2nd assembly iteration.\n"
    echo -e "   -p3\t assembly parameter file for 3rd assembly iteration.\n"
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
    -np|--numprocs)		# MPI processes to use
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
    -rma|--MPIRMA)                        # Maximum memory limit
    useMPIRMA="Y"
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
if [ "$useMPIRMA" = "Y" ] ; then
   BUILDGEXE="${exePath}/buildG-MPIRMA"
else
   BUILDGEXE="${exePath}/buildG-MPI"
fi

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
   echo "Graph simplification module $SIMPLIFYEXE exists."
else
   echo "Graph simplification module $SIMPLIFYEXE does not exist in script directory."
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
j=0
while [ $j -lt $numProcs ]
do
   i=1
   while [ $i -lt $numThreads ]
   do
      edgeFiles="${edgeFiles}${dataOutPath}/graph/${outPrefix}_${j}_${i}_parGraph.txt,"
      i=$(( $i+1 ))
   done
   j=$(( $j+1 ))
done
edgeFiles=${edgeFiles%?}

#Create contained read file list
j=0
while [ $j -lt $numProcs ]
do
   i=1
   while [ $i -lt $numThreads ]
   do
      containedReads="${containedReads}${dataOutPath}/graph/${outPrefix}_${j}_${i}_containedReads.txt,"
      i=$(( $i+1 ))
   done
   j=$(( $j+1 ))
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
elif [ -z "$readFile1" ] && [ -z "$readFile2" ] && [ -z "$readFileS" ] ; then     #Multiple interleaved PE file as input
   #BBTools Preprocessing
   OLDIFS=$IFS
   IFS=',' 
   array=($readFileP)
   IFS=$OLDIFS
   trimFtlEccOutput="" #File list for error corrected output
   for element in "${array[@]}"
   do
      fullName="$(echo -e "${element}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"	#Remove any leading or trailing spaces
      fName=`basename $fullName`
      trimFtlEccOutput="${trimFtlEccOutput},tecc.ftl.trm.${fName}" 
   done
   trimFtlEccOutput=${trimFtlEccOutput#?}
   if [ "$constructGraph" = "Y" ] ; then
<<<<<<< HEAD
      mpirun -np $numProcs --map-by ppr:1:node:pe=${numThreads} ${BUILDGEXE} -pe ${trimFtlEccOutput} -f $outGraphPrefix -s ${outSimplifyPrefix} -simPth ${exePath} -p ${asmParaFileP} -t ${numThreads} -m ${maxMem}  > ${logFile}
   fi
   if [ "$simplifyGraph" = "Y" ] ; then
      ${exePath}/fullsimplify -fpi ${trimFtlEccOutput} -e ${edgeFiles} -crd ${containedReads} -simPth ${exePath} -p ${asmParaFileP} -p2 ${asmParaFileP2} -p3 ${asmParaFileP3} -o $outSimplifyPrefix -t ${numThreads} -log DEBUG4 >> ${logFile} 2>&1
   fi
elif [ -z "$readFile1" ] && [ -z "$readFile2" ] && [ -z "$readFileP" ] ; then              #Multiple SE file as input
   #BBTools Preprocessing
   OLDIFS=$IFS
   IFS=',' 
   array=($readFileS)
   IFS=$OLDIFS
   trimFtlEccOutput="" #File list for error corrected output
   for element in "${array[@]}"
   do
      fullName="$(echo -e "${element}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"	#Remove any leading or trailing spaces
      fName=`basename $fullName`
      trimFtlEccOutput="${trimFtlEccOutput},tecc.ftl.trm.${fName}" 
   done
   trimFtlEccOutput=${trimFtlEccOutput#?}
   if [ "$constructGraph" = "Y" ] ; then
       mpirun -np $numProcs --map-by ppr:1:node:pe=${numThreads} ${BUILDGEXE} -se ${trimFtlEccOutput} -f $outGraphPrefix -s ${outSimplifyPrefix} -simPth ${exePath} -p ${asmParaFileP} -t ${numThreads} -m ${maxMem}  &> ${logFile}
   fi
   if [ "$simplifyGraph" = "Y" ] ; then
      ${exePath}/fullsimplify -fs ${trimFtlEccOutput} -e ${edgeFiles} -crd ${containedReads} -simPth ${exePath} -p ${asmParaFileP} -p2 ${asmParaFileP2} -p3 ${asmParaFileP3} -o $outSimplifyPrefix -t ${numThreads} -log DEBUG4 >> ${logFile} 2>&1
   fi
elif [ -z "$readFileS" ] && [ -z "$readFileP" ] ; then		#Multiple P1/P2 files as input
    #BBTools Preprocessing
   # Process pair R1
   OLDIFS=$IFS
   IFS=',' 
   array1=($readFile1)
   array2=($readFile2)
   IFS=$OLDIFS
   trimFtlEccOutput="" #File list for error corrected output
   i=0
   for element in "${array1[@]}"
   do
      trimFtlEccOutput="${trimFtlEccOutput},int.tecc.ftl.trm.${i}.${extension1}"      #File list for error corrected output
      i=$(( $i+1 ))
   done
   trimFtlEccOutput=${trimFtlEccOutput#?}
   if [ "$constructGraph" = "Y" ] ; then
       mpirun -np $numProcs --map-by ppr:1:node:pe=${numThreads} ${BUILDGEXE} -pe ${trimFtlEccOutput} -f $outGraphPrefix -s ${outSimplifyPrefix} -simPth ${exePath} -p ${asmParaFileP} -t ${numThreads} -m ${maxMem}  &> ${logFile}
   fi
   if [ "$simplifyGraph" = "Y" ] ; then
      ${exePath}/fullsimplify -fp ${trimFtlEccOutput} -e ${edgeFiles} -crd ${containedReads} -simPth ${exePath} -p ${asmParaFileP} -p2 ${asmParaFileP2} -p3 ${asmParaFileP3} -o $outSimplifyPrefix -t ${numThreads} -log DEBUG4 >> ${logFile} 2>&1
   fi
elif [ -z "$readFile1" ] && [ -z "$readFile2" ] ; then		#Multiple Interleaved PE file and a Multiple SE file as input
   #BBTools Preprocessing
   #Process SE files
   OLDIFS=$IFS
   IFS=',' 
   arrayS=($readFileS)
   arrayP=($readFileP)
   IFS=$OLDIFS
   trimFtlEccOutputS="" #File list for error corrected output
   for element in "${arrayS[@]}"
   do
      fullName="$(echo -e "${element}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"	#Remove any leading or trailing spaces
      fName=`basename $fullName`
      trimFtlEccOutputS="${trimFtlEccOutputS},tecc.ftl.trm.${fName}" 
   done
   trimFtlEccOutputS=${trimFtlEccOutputS#?}
   #Process PE files
   trimFtlEccOutputP="" #File list for error corrected output
   for element in "${arrayP[@]}"
   do
      fullName="$(echo -e "${element}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"	#Remove any leading or trailing spaces
      fName=`basename $fullName`
      trimFtlEccOutputP="${trimFtlEccOutputP},tecc.ftl.trm.${fName}" 
   done
   trimFtlEccOutputP=${trimFtlEccOutputP#?}
   if [ "$constructGraph" = "Y" ] ; then
       mpirun -np $numProcs --map-by ppr:1:node:pe=${numThreads} ${BUILDGEXE} -pe ${trimFtlEccOutputP} -se ${trimFtlEccOutputS} -f $outGraphPrefix -s ${outSimplifyPrefix} -simPth ${exePath} -p ${asmParaFileP} -t ${numThreads} -m ${maxMem}  &> ${logFile}
=======
      mpirun -np $numProcs --map-by ppr:1:node:pe=${numThreads} ${BUILDGEXE} -pe ${trimFtlEccOutput} -f $outGraphPrefix -p ${asmParaFileP} -t ${numThreads} -m ${maxMem}  > ${logFile}
   fi
   if [ "$simplifyGraph" = "Y" ] ; then
      ${exePath}/fullsimplify -fpi ${trimFtlEccOutput} -e ${edgeFiles} -crd ${containedReads} -simPth ${exePath} -p ${asmParaFileP} -p2 ${asmParaFileP2} -p3 ${asmParaFileP3} -o $outSimplifyPrefix -t ${numThreads} -log DEBUG4 >> ${logFile} 2>&1
   fi
elif [ -z "$readFile1" ] && [ -z "$readFile2" ] && [ -z "$readFileP" ] ; then              #Multiple SE file as input
   #BBTools Preprocessing
   OLDIFS=$IFS
   IFS=',' 
   array=($readFileS)
   IFS=$OLDIFS
   trimFtlEccOutput="" #File list for error corrected output
   for element in "${array[@]}"
   do
      fullName="$(echo -e "${element}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"	#Remove any leading or trailing spaces
      fName=`basename $fullName`
      trimFtlEccOutput="${trimFtlEccOutput},tecc.ftl.trm.${fName}" 
   done
   trimFtlEccOutput=${trimFtlEccOutput#?}
   if [ "$constructGraph" = "Y" ] ; then
       mpirun -np $numProcs --map-by ppr:1:node:pe=${numThreads} ${BUILDGEXE} -se ${trimFtlEccOutput} -f $outGraphPrefix -s ${outSimplifyPrefix} -simPth ${exePath} -p ${asmParaFileP} -t ${numThreads} -m ${maxMem}  &> ${logFile}
   fi
   if [ "$simplifyGraph" = "Y" ] ; then
      ${exePath}/fullsimplify -fs ${trimFtlEccOutput} -e ${edgeFiles} -crd ${containedReads} -simPth ${exePath} -p ${asmParaFileP} -p2 ${asmParaFileP2} -p3 ${asmParaFileP3} -o $outSimplifyPrefix -t ${numThreads} -log DEBUG4 >> ${logFile} 2>&1
   fi
elif [ -z "$readFileS" ] && [ -z "$readFileP" ] ; then		#Multiple P1/P2 files as input
    #BBTools Preprocessing
   # Process pair R1
   OLDIFS=$IFS
   IFS=',' 
   array1=($readFile1)
   array2=($readFile2)
   IFS=$OLDIFS
   trimFtlEccOutput="" #File list for error corrected output
   i=0
   for element in "${array1[@]}"
   do
      trimFtlEccOutput="${trimFtlEccOutput},int.tecc.ftl.trm.${i}.${extension1}"      #File list for error corrected output
      i=$(( $i+1 ))
   done
   trimFtlEccOutput=${trimFtlEccOutput#?}
   if [ "$constructGraph" = "Y" ] ; then
       mpirun -np $numProcs --map-by ppr:1:node:pe=${numThreads} ${BUILDGEXE} -pe ${trimFtlEccOutput} -f $outGraphPrefix -s ${outSimplifyPrefix} -simPth ${exePath} -p ${asmParaFileP} -t ${numThreads} -m ${maxMem}  &> ${logFile}
   fi
   if [ "$simplifyGraph" = "Y" ] ; then
      ${exePath}/fullsimplify -fp ${trimFtlEccOutput} -e ${edgeFiles} -crd ${containedReads} -simPth ${exePath} -p ${asmParaFileP} -p2 ${asmParaFileP2} -p3 ${asmParaFileP3} -o $outSimplifyPrefix -t ${numThreads} -log DEBUG4 >> ${logFile} 2>&1
   fi
elif [ -z "$readFile1" ] && [ -z "$readFile2" ] ; then		#Multiple Interleaved PE file and a Multiple SE file as input
   #BBTools Preprocessing
   #Process SE files
   OLDIFS=$IFS
   IFS=',' 
   arrayS=($readFileS)
   arrayP=($readFileP)
   IFS=$OLDIFS
   trimFtlEccOutputS="" #File list for error corrected output
   for element in "${arrayS[@]}"
   do
      fullName="$(echo -e "${element}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"	#Remove any leading or trailing spaces
      fName=`basename $fullName`
      trimFtlEccOutputS="${trimFtlEccOutputS},tecc.ftl.trm.${fName}" 
   done
   trimFtlEccOutputS=${trimFtlEccOutputS#?}
   #Process PE files
   trimFtlEccOutputP="" #File list for error corrected output
   for element in "${arrayP[@]}"
   do
      fullName="$(echo -e "${element}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"	#Remove any leading or trailing spaces
      fName=`basename $fullName`
      trimFtlEccOutputP="${trimFtlEccOutputP},tecc.ftl.trm.${fName}" 
   done
   trimFtlEccOutputP=${trimFtlEccOutputP#?}
   if [ "$constructGraph" = "Y" ] ; then
       mpirun -np $numProcs --map-by ppr:1:node:pe=${numThreads} ${BUILDGEXE} -pe ${trimFtlEccOutputP} -se ${trimFtlEccOutputS} -f $outGraphPrefix -s ${outSimplifyPrefix} -simPth ${exePath} -p ${asmParaFileP} -t ${numThreads} -m ${maxMem}  &> ${logFile}
   fi
   if [ "$simplifyGraph" = "Y" ] ; then
      ${exePath}/fullsimplify -fpi ${trimFtlEccOutputP} -fs ${trimFtlEccOutputS} -e ${edgeFiles} -crd ${containedReads} -simPth ${exePath} -p ${asmParaFileP} -p2 ${asmParaFileP2} -p3 ${asmParaFileP3} -o $outSimplifyPrefix -t ${numThreads} -log DEBUG4 >> ${logFile} 2>&1
   fi
else
   echo "Invalid combination of input files. You can specify either a set of comma seperated interleaved paired file or two separate paired files not both. Any number of comma seperated single end files can be provided. Exiting..."
   exit 1
fi
`cat ${outSimplifyPrefix}_contigsFinal_*.fasta > ${outSimplifyPrefix}_contigsFinalCombined.fasta`
`cat ${outSimplifyPrefix}_scaffoldsFinal_*.fasta > ${outSimplifyPrefix}_scaffoldsFinalCombined.fasta`
`cp ${outSimplifyPrefix}_contigsFinalCombined.fasta ${dataOutPath}`
`cp ${outSimplifyPrefix}_scaffoldsFinalCombined.fasta ${dataOutPath}`
`python ${exePath}/assemblyStats.py denovo -i ${dataOutPath}/${outPrefix}_contigsFinalCombined.fasta`
`python ${exePath}/assemblyStats.py denovo -i ${dataOutPath}/${outPrefix}_scaffoldsFinalCombined.fasta`
echo Ending Time is $(date)

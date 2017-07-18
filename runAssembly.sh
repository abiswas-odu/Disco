#!/bin/bash 

#Get path of run script
if [ -L $0 ] ; then
    exePath=$(dirname $(readlink -f $0)) ;
else
    exePath=$(dirname $0) ;
fi;

#Basic parameters defaults
numThreads=`nproc`
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
dataOutPath="."
phymem=`grep MemTotal /proc/meminfo | awk '{print $2}'`
maxMem=$((phymem/1048576))
while [[ $# > 0 ]]
do
key="$1"
case $key in
    -h|--help)              # help output
    echo -e "Usage:\n"
    echo -e "   runDisco.sh [OPTION]...<PARAM>...\n\n"
    echo -e "<PARAMS>\n"
    echo -e "   -inS\t single read filenames (comma separated fasta/fastq/fastq.gz files).\n"
    echo -e "   -in1\t forward paired read filename (single fasta/fastq/fastq.gz file).\n"
    echo -e "   -in2\t reverse paired read filename (single fasta/fastq/fastq.gz file).\n"
    echo -e "   -inP\t interleaved paired read filenames (comma separated fasta/fastq/fastq.gz files).\n"
    echo -e "   -d\t output directory path.(DEFAULT: current directory)\n"
    echo -e "   -o\t output filename prefix.(DEFAULT: disco)\n"
    echo -e "<OPTIONS>\n"
    echo -e "   -h\t help.\n"
    echo -e "   -m\t maximum memory to be used (DEFAULT: $maxMem GB).\n"
    echo -e "   -n\t number of threads (DEFAULT: $numThreads).\n"
    echo -e "   -obg\t only build overlap graph (DEFAULT: False).\n"
    echo -e "   -osg\t only simplify existing overlap graph (DEFAULT: False).\n"
    echo -e "   -p\t assembly parameter file for 1st assembly iteration.\n"
    echo -e "   -p2\t assembly parameter file for 2nd assembly iteration.\n"
    echo -e "   -p3\t assembly parameter file for 3rd assembly iteration.\n"
    exit 1
    ;;
    -p|--parameterFile)     # parameter file for first interation of assembly
    asmParaFileP="$2"
    shift # past argument
    ;;
    -p2|--parameterFile2)     # parameter file for second interation of assembly
    asmParaFileP2="$2"
    shift # past argument
    ;;
    -p3|--parameterFile3)     # parameter file for subsequebt interation of assembly
    asmParaFileP3="$2"
    shift # past argument
    ;;
    -d|--outdir)			# output directory
    dataOutPath="$2"
    shift # past argument
    ;;
    -in1)					# Forward paired end read file -- single file
    readFile1="$2"
    shift # past argument
    ;;
    -in2)					# Reverse paired end read file -- single file
    readFile2="$2"
    shift # past argument
    ;;
    -inS)					# Single read file  -- multiple files
    readFileS="$2"
    shift # past argument
    ;;
    -inP)					# Interleaved paired end read file -- multiple files
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
   echo "Graph simplification module $SIMPLIFYEXE exists."
else
   echo "Graph simplification module $SIMPLIFYEXE does not exist in script directory."
fi
#All exe files in place


#check if output directory exists, if not create it
if [ -d $dataOutPath ] ; then
   echo "Output directory exists."
else
   echo "Creating output directory: $dataOutPath"
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

#Start Assembly
if [ -z "$readFile1" ] && [ -z "$readFile2" ] && [ -z "$readFileS" ] && [ -z "$readFileP" ] ; then
   echo "No input files specified. Exiting..."
   exit 1
elif [ -z "$readFile1" ] && [ -z "$readFile2" ] && [ -z "$readFileS" ] ; then		#Multiple interleaved PE file as input
   #BBTools Preprocessing
   OLDIFS=$IFS
   IFS=',' 
   array=($readFileP)
   IFS=$OLDIFS
   trimFtlOutput=""       #File list for filtered output
   rmTrimOutput=""		#File list to delete trimming output
   rmTrimFtlOutput=""     #File list for deleting filtered output
   trimFtlEccOutput="" #File list for error corrected output
   for element in "${array[@]}"
   do
      fullName="$(echo -e "${element}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"	#Remove any leading or trailing spaces
      fName=`basename $fullName`
      ${exePath}/bbmap/bbduk.sh in=${fullName} out=trm.${fName} ktrim=r k=11 mink=7 hdist=1 tpe tbo ref=${exePath}/bbmap/resources/adapters.fa ftm=5 qtrim=r trimq=15
      ${exePath}/bbmap/bbduk.sh in=trm.${fName} out=ftl.trm.${fName} k=23 hdist=1 ref=${exePath}/bbmap/resources/sequencing_artifacts.fa.gz,${exePath}/bbmap/resources/phix174_ill.ref.fa.gz
      trimFtlOutput="${trimFtlOutput},ftl.trm.${fName}"       
      rmTrimOutput="${rmTrimOutput} trm.${fName}"		
      rmTrimFtlOutput="${rmTrimFtlOutput} ftl.trm.${fName}"     
      trimFtlEccOutput="${trimFtlEccOutput},tecc.ftl.trm.${fName}" 
   done
   trimFtlOutput=${trimFtlOutput#?}
   trimFtlEccOutput=${trimFtlEccOutput#?}
   rmTrimOutput=${rmTrimOutput#?}
   rmTrimFtlOutput=${rmTrimFtlOutput#?}
   ${exePath}/bbmap/tadpole.sh in=${trimFtlOutput} out=${trimFtlEccOutput} ecc k=23 prealloc prefilter=2 tossjunk
   #Delete intermediate files after pre-processing
   rm ${rmTrimOutput} ${rmTrimFtlOutput}
   if [ "$constructGraph" = "Y" ] ; then
      ${exePath}/buildG -pe ${trimFtlEccOutput} -f $outGraphPrefix -p ${asmParaFileP} -t ${numThreads} -m ${maxMem}  > ${logFile}
   fi
   if [ "$simplifyGraph" = "Y" ] ; then
      ${exePath}/fullsimplify -fpi ${trimFtlEccOutput} -e ${edgeFiles} -crd ${containedReads} -simPth ${exePath} -p ${asmParaFileP} -p2 ${asmParaFileP2} -p3 ${asmParaFileP3} -o $outSimplifyPrefix -t ${numThreads} -log DEBUG4 >> ${logFile} 2>&1
   fi
elif [ -z "$readFile1" ] && [ -z "$readFile2" ] && [ -z "$readFileP" ] ; then		#Multiple SE file as input
   #BBTools Preprocessing
   OLDIFS=$IFS
   IFS=',' 
   array=($readFileS)
   IFS=$OLDIFS
   trimFtlOutput=""       #File list for filtered output
   rmTrimOutput=""		#File list to delete trimming output
   rmTrimFtlOutput=""     #File list for deleting filtered output
   trimFtlEccOutput="" #File list for error corrected output
   for element in "${array[@]}"
   do
      fullName="$(echo -e "${element}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"	#Remove any leading or trailing spaces
      fName=`basename $fullName`
      ${exePath}/bbmap/bbduk.sh in=${fullName} out=trm.${fName} ktrim=r k=11 mink=07 hdist=1 tpe tbo ref=${exePath}/bbmap/resources/adapters.fa ftm=5 qtrim=r trimq=15
      ${exePath}/bbmap/bbduk.sh in=trm.${fName} out=ftl.trm.${fName} k=23 hdist=1 ref=${exePath}/bbmap/resources/sequencing_artifacts.fa.gz,${exePath}/bbmap/resources/phix174_ill.ref.fa.gz
      trimFtlOutput="${trimFtlOutput},ftl.trm.${fName}"       
      rmTrimOutput="${rmTrimOutput} trm.${fName}"		
      rmTrimFtlOutput="${rmTrimFtlOutput} ftl.trm.${fName}"     
      trimFtlEccOutput="${trimFtlEccOutput},tecc.ftl.trm.${fName}" 
   done
   trimFtlOutput=${trimFtlOutput#?}
   trimFtlEccOutput=${trimFtlEccOutput#?}
   rmTrimOutput=${rmTrimOutput#?}
   rmTrimFtlOutput=${rmTrimFtlOutput#?}
   ${exePath}/bbmap/tadpole.sh in=${trimFtlOutput} out=${trimFtlEccOutput} ecc k=23 prealloc prefilter=2 tossjunk
   #Delete intermediate files after pre-processing
   rm ${rmTrimOutput} ${rmTrimFtlOutput}
   if [ "$constructGraph" = "Y" ] ; then
      ${exePath}/buildG -se ${trimFtlEccOutput} -f $outGraphPrefix -p ${asmParaFileP} -t ${numThreads} -m ${maxMem}  &> ${logFile}
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
   trimFtlOutput1=""       #File list for filtered output
   rmTrimOutput1=""		#File list to delete trimming output
   rmTrimFtlOutput1=""     #File list for deleting filtered output
   trimFtlOutput2=""       #File list for filtered output
   rmTrimOutput2=""		#File list to delete trimming output
   rmTrimFtlOutput2=""     #File list for deleting filtered output
   trimFtlEccOutput="" #File list for error corrected output
   i=0
   for element in "${array1[@]}"
   do
      fullName1="$(echo -e "${element}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"	#Remove any leading or trailing spaces
      fullName2="$(echo -e "${array2[$i]}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"	#Remove any leading or trailing spaces
      fName1=`basename $fullName1`
      fName2=`basename $fullName2`
       ${exePath}/bbmap/bbduk.sh in=${fullName1} in2=${fullName2} out=trm.${fName1} out2=trm.${fName2} ktrim=r k=23 mink=11 hdist=1 tpe tbo ref=${exePath}/bbmap/resources/adapters.fa ftm=5 qtrim=r trimq=10
   ${exePath}/bbmap/bbduk.sh in=trm.${fName1} in2=trm.${fName2} out=ftl.trm.${fName1} out2=ftl.trm.${fName2} k=31 hdist=1 ref=${exePath}/bbmap/resources/sequencing_artifacts.fa.gz,${exePath}/bbmap/resources/phix174_ill.ref.fa.gz
      trimFtlOutput1="${trimFtlOutput1},ftl.trm.${fName1}"
      trimFtlOutput2="${trimFtlOutput2},ftl.trm.${fName2}"       
      rmTrimOutput1="${rmTrimOutput1} trm.${fName1}"		
      rmTrimOutput2="${rmTrimOutput2} trm.${fName2}"		
      rmTrimFtlOutput1="${rmTrimFtlOutput1} ftl.trm.${fName1}"     
      rmTrimFtlOutput2="${rmTrimFtlOutput2} ftl.trm.${fName2}"     
      extension1="${fName1#*.}"
      trimFtlEccOutput="${trimFtlEccOutput},int.tecc.ftl.trm.${i}.${extension1}"      #File list for error corrected output
      i=$(( $i+1 ))
   done
   trimFtlOutput1=${trimFtlOutput1#?}
   trimFtlOutput2=${trimFtlOutput2#?}
   rmTrimOutput1=${rmTrimOutput1#?}
   rmTrimOutput2=${rmTrimOutput2#?}
   rmTrimFtlOutput1=${rmTrimFtlOutput1#?}
   rmTrimFtlOutput2=${rmTrimFtlOutput2#?}
   trimFtlEccOutput=${trimFtlEccOutput#?}
   #Error Correction
   ${exePath}/bbmap/tadpole.sh in=${trimFtlOutput1} in2=${trimFtlOutput2} out=${trimFtlEccOutput} ecc k=31 prealloc prefilter=2 tossjunk
   #Delete intermediate files after pre-processing
   rm ${rmTrimOutput1} ${rmTrimFtlOutput1} ${rmTrimOutput2} ${rmTrimFtlOutput2}
   if [ "$constructGraph" = "Y" ] ; then
      ${exePath}/buildG -pe ${trimFtlEccOutput} -f $outGraphPrefix -p ${asmParaFileP} -t ${numThreads} -m ${maxMem}  &> ${logFile}
   fi
   if [ "$simplifyGraph" = "Y" ] ; then
      ${exePath}/fullsimplify -fpi ${trimFtlEccOutput} -e ${edgeFiles} -crd ${containedReads} -simPth ${exePath} -p ${asmParaFileP} -p2 ${asmParaFileP2} -p3 ${asmParaFileP3} -o $outSimplifyPrefix -t ${numThreads} -log DEBUG4 >> ${logFile} 2>&1
   fi
elif [ -z "$readFile1" ] && [ -z "$readFile2" ] ; then		#Multiple Interleaved PE file and a Multiple SE file as input
   #BBTools Preprocessing
   #Process SE files
   OLDIFS=$IFS
   IFS=',' 
   arrayS=($readFileS)
   arrayP=($readFileP)
   IFS=$OLDIFS
   trimFtlOutputS=""       #File list for filtered output
   rmTrimOutputS=""		#File list to delete trimming output
   rmTrimFtlOutputS=""     #File list for deleting filtered output
   trimFtlEccOutputS="" #File list for error corrected output
   for element in "${arrayS[@]}"
   do
      fullName="$(echo -e "${element}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"	#Remove any leading or trailing spaces
      fName=`basename $fullName`
      ${exePath}/bbmap/bbduk.sh in=${fullName} out=trm.${fName} ktrim=r k=23 mink=11 hdist=1 tpe tbo ref=${exePath}/bbmap/resources/adapters.fa ftm=5 qtrim=r trimq=10
      ${exePath}/bbmap/bbduk.sh in=trm.${fName} out=ftl.trm.${fName} k=31 hdist=1 ref=${exePath}/bbmap/resources/sequencing_artifacts.fa.gz,${exePath}/bbmap/resources/phix174_ill.ref.fa.gz
      trimFtlOutputS="${trimFtlOutputS},ftl.trm.${fName}"       
      rmTrimOutputS="${rmTrimOutputS} trm.${fName}"		
      rmTrimFtlOutputS="${rmTrimFtlOutputS} ftl.trm.${fName}"     
      trimFtlEccOutputS="${trimFtlEccOutputS},tecc.ftl.trm.${fName}" 
   done
   trimFtlOutputS=${trimFtlOutputS#?}
   trimFtlEccOutputS=${trimFtlEccOutputS#?}
   rmTrimOutputS=${rmTrimOutputS#?}
   rmTrimFtlOutputS=${rmTrimFtlOutputS#?}
  
   #Process PE files
   trimFtlOutputP=""       #File list for filtered output
   rmTrimOutputP=""		#File list to delete trimming output
   rmTrimFtlOutputP=""     #File list for deleting filtered output
   trimFtlEccOutputP="" #File list for error corrected output
   for element in "${arrayP[@]}"
   do
      fullName="$(echo -e "${element}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"	#Remove any leading or trailing spaces
      fName=`basename $fullName`
      ${exePath}/bbmap/bbduk.sh in=${fullName} out=trm.${fName} ktrim=r k=23 mink=11 hdist=1 tpe tbo ref=${exePath}/bbmap/resources/adapters.fa ftm=5 qtrim=r trimq=10
      ${exePath}/bbmap/bbduk.sh in=trm.${fName} out=ftl.trm.${fName} k=31 hdist=1 ref=${exePath}/bbmap/resources/sequencing_artifacts.fa.gz,${exePath}/bbmap/resources/phix174_ill.ref.fa.gz
      trimFtlOutputP="${trimFtlOutputP},ftl.trm.${fName}"       
      rmTrimOutputP="${rmTrimOutputP} trm.${fName}"		
      rmTrimFtlOutputP="${rmTrimFtlOutputP} ftl.trm.${fName}"     
      trimFtlEccOutputP="${trimFtlEccOutputP},tecc.ftl.trm.${fName}" 
   done
   trimFtlOutputP=${trimFtlOutputP#?}
   trimFtlEccOutputP=${trimFtlEccOutputP#?}
   rmTrimOutputP=${rmTrimOutputP#?}
   rmTrimFtlOutputP=${rmTrimFtlOutputP#?}
   #Error Correction
   ${exePath}/bbmap/tadpole.sh in=${trimFtlOutputP},${trimFtlOutputS} out=${trimFtlEccOutputP},${trimFtlEccOutputS} ecc k=31 prealloc prefilter=2 tossjunk
   #Delete intermediate files after pre-processing
   rm ${rmTrimOutputS} ${rmTrimFtlOutputS} ${rmTrimOutputP} ${rmTrimFtlOutputP}
   if [ "$constructGraph" = "Y" ] ; then
      ${exePath}/buildG -pe ${trimFtlEccOutputP} -se ${trimFtlEccOutputS} -f $outGraphPrefix -p ${asmParaFileP} -t ${numThreads} -m ${maxMem}  &> ${logFile}
   fi
   if [ "$simplifyGraph" = "Y" ] ; then
      ${exePath}/fullsimplify -fpi ${trimFtlEccOutputP} -fs ${trimFtlEccOutputS} -e ${edgeFiles} -crd ${containedReads} -simPth ${exePath} -p ${asmParaFileP} -p2 ${asmParaFileP2} -p3 ${asmParaFileP3} -o $outSimplifyPrefix -t ${numThreads} -log DEBUG4 >> ${logFile} 2>&1
   fi
elif [ -z "$readFileP" ] ; then			#Multiple sepatate P1/P2 files and a Multiple SE file as input
   #BBTools Preprocessing
   # Process pair reads
   OLDIFS=$IFS
   IFS=',' 
   array1=($readFile1)
   array2=($readFile2)
   arrayS=($readFileS)
   IFS=$OLDIFS
   trimFtlOutput1=""       #File list for filtered output
   rmTrimOutput1=""		#File list to delete trimming output
   rmTrimFtlOutput1=""     #File list for deleting filtered output
   trimFtlOutput2=""       #File list for filtered output
   rmTrimOutput2=""		#File list to delete trimming output
   rmTrimFtlOutput2=""     #File list for deleting filtered output
   trimFtlEccOutput="" #File list for error corrected output
   i=0
   for element in "${array1[@]}"
   do
      fullName1="$(echo -e "${element}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"	#Remove any leading or trailing spaces
      fullName2="$(echo -e "${array2[$i]}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"	#Remove any leading or trailing spaces
      fName1=`basename $fullName1`
      fName2=`basename $fullName2`
       ${exePath}/bbmap/bbduk.sh in=${fullName1} in2=${fullName2} out=trm.${fName1} out2=trm.${fName2} ktrim=r k=23 mink=11 hdist=1 tpe tbo ref=${exePath}/bbmap/resources/adapters.fa ftm=5 qtrim=r trimq=10
   ${exePath}/bbmap/bbduk.sh in=trm.${fName1} in2=trm.${fName2} out=int.ftl.trm.${fName1} k=31 hdist=1 ref=${exePath}/bbmap/resources/sequencing_artifacts.fa.gz,${exePath}/bbmap/resources/phix174_ill.ref.fa.gz
      trimFtlOutput="${trimFtlOutput},int.ftl.trm.${fName1}"
      rmTrimOutput1="${rmTrimOutput1} trm.${fName1}"		
      rmTrimOutput2="${rmTrimOutput2} trm.${fName2}"		
      rmTrimFtlOutput="${rmTrimFtlOutput} int.ftl.trm.${fName1}"     
      extension1="${fName1#*.}"
      trimFtlEccOutput="${trimFtlEccOutput},tecc.int.ftl.trm.${fName1}"      #File list for error corrected output
      i=$(( $i+1 ))
   done
   trimFtlOutput=${trimFtlOutput#?}
   rmTrimOutput1=${rmTrimOutput1#?}
   rmTrimOutput2=${rmTrimOutput2#?}
   rmTrimFtlOutput=${rmTrimFtlOutput#?}
   trimFtlEccOutput=${trimFtlEccOutput#?}

   trimFtlOutputS=""       #File list for filtered output
   rmTrimOutputS=""		#File list to delete trimming output
   rmTrimFtlOutputS=""     #File list for deleting filtered output
   trimFtlEccOutputS="" #File list for error corrected output
   for element in "${arrayS[@]}"
   do
      fullName="$(echo -e "${element}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"	#Remove any leading or trailing spaces
      fName=`basename $fullName`
      ${exePath}/bbmap/bbduk.sh in=${fullName} out=trm.${fName} ktrim=r k=23 mink=11 hdist=1 tpe tbo ref=${exePath}/bbmap/resources/adapters.fa ftm=5 qtrim=r trimq=10
      ${exePath}/bbmap/bbduk.sh in=trm.${fName} out=ftl.trm.${fName} k=31 hdist=1 ref=${exePath}/bbmap/resources/sequencing_artifacts.fa.gz,${exePath}/bbmap/resources/phix174_ill.ref.fa.gz
      trimFtlOutputS="${trimFtlOutputS},ftl.trm.${fName}"       
      rmTrimOutputS="${rmTrimOutputS} trm.${fName}"		
      rmTrimFtlOutputS="${rmTrimFtlOutputS} ftl.trm.${fName}"     
      trimFtlEccOutputS="${trimFtlEccOutputS},tecc.ftl.trm.${fName}" 
   done
   trimFtlOutputS=${trimFtlOutputS#?}
   trimFtlEccOutputS=${trimFtlEccOutputS#?}
   rmTrimOutputS=${rmTrimOutputS#?}
   rmTrimFtlOutputS=${rmTrimFtlOutputS#?}
   #Error Correction
   ${exePath}/bbmap/tadpole.sh in=${trimFtlOutput},${trimFtlOutputS} out=${trimFtlEccOutput},${trimFtlEccOutputS} ecc k=31 prealloc prefilter=2 tossjunk
   #Delete intermediate files after pre-processing
   rm ${rmTrimOutput1} ${rmTrimFtlOutput} ${rmTrimOutput2} ${rmTrimFtlOutputS}
   if [ "$constructGraph" = "Y" ] ; then
      
      ${exePath}/buildG -pe ${trimFtlEccOutput} -se ${trimFtlEccOutputS} -f $outGraphPrefix -p ${asmParaFileP} -t ${numThreads} -m ${maxMem}  &> ${logFile}
   fi
   if [ "$simplifyGraph" = "Y" ] ; then
      ${exePath}/fullsimplify -fp ${trimFtlEccOutput} -fs ${trimFtlEccOutputS} -e ${edgeFiles} -crd ${containedReads} -simPth ${exePath} -p ${asmParaFileP} -p2 ${asmParaFileP2} -p3 ${asmParaFileP3} -o $outSimplifyPrefix -t ${numThreads} -log DEBUG4 >> ${logFile} 2>&1
      `cat ${outSimplifyPrefix}_contigsFinal_*.fasta > ${outSimplifyPrefix}_contigsFinalCombined.fasta`
      `cat ${outSimplifyPrefix}_scaffoldsFinal_*.fasta > ${outSimplifyPrefix}_scaffoldsFinalCombined.fasta`
      `cp ${outSimplifyPrefix}_contigsFinalCombined.fasta ${dataOutPath}`
      `cp ${outSimplifyPrefix}_scaffoldsFinalCombined.fasta ${dataOutPath}`
   fi
else
   echo "Invalid combination of input files. You can specify either a set of comma separated interleaved paired file or two separate paired files not both. Any number of comma separated single end files can be provided. Exiting..."
   exit 1
fi
`cat ${outSimplifyPrefix}_contigsFinal_*.fasta > ${outSimplifyPrefix}_contigsFinalCombined.fasta`
`cat ${outSimplifyPrefix}_scaffoldsFinal_*.fasta > ${outSimplifyPrefix}_scaffoldsFinalCombined.fasta`
`cp ${outSimplifyPrefix}_contigsFinalCombined.fasta ${dataOutPath}`
`cp ${outSimplifyPrefix}_scaffoldsFinalCombined.fasta ${dataOutPath}`
`python ${exePath}/assemblyStats.py denovo -i ${dataOutPath}/${outPrefix}_contigsFinalCombined.fasta`
`python ${exePath}/assemblyStats.py denovo -i ${dataOutPath}/${outPrefix}_scaffoldsFinalCombined.fasta`
echo Ending Time is $(date)

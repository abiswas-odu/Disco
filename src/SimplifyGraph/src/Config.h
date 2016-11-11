#ifndef CONFIG_H
#define CONFIG_H

//============================================================================
// Name        : Config.h
// Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey, Abhishek Biswas
// Version     : v3.0
// Copyright   : 2017 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : Config header file
//============================================================================

//multi-thread library OPENMP
#include <omp.h>

// C headers:
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <assert.h> // for assert()
#include <spawn.h>
#include "Utils.h"

// C++ headers:
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <string>
#include <vector>
#include <utility>
#include <limits>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <streambuf>
#include <map>
#include <typeinfo>
#include <cstring>
#include <set>
#include <sys/wait.h>
#include <zlib.h>

#define SSTR( x ) dynamic_cast< std::ostringstream & >(( std::ostringstream() << std::dec << x ) ).str()

// Local headers
#include "logcpp/log.h"
using namespace std;

#define DISK_GRAPH_UPDATE 18000.0

#ifndef BUFFER_SIZE
#define BUFFER_SIZE 100
#endif

#define DEBUG


//============================================================================
// global variables with default value
//============================================================================
// Minimum overlap length to create an edge
extern unsigned int minOvl;

// minimum reads count in an edge to be not dead end edge (default: 10)
extern unsigned int minReadsCountInEdgeToBeNotDeadEnd;	
// minimum edge length to be not dead end edge (default: 1000)
extern unsigned int minEdgeLengthToBeNotDeadEnd;

// minimum reads count for an edge to be kept even if it has 0 flow (default: 15)
extern unsigned int minReadsCountToHave0Flow;
// minimum edge length or an edge to be kept even if it has 0 flow (default: 1500)
extern unsigned int minEdgeLengthToHave0Flow;

// minimum reads count in an edge to be 1 minimum flow (default: 10)
extern unsigned int minReadsCountInEdgeToBe1MinFlow;	
// minimum edge length to be 1 minimum flow (default: 1000)
extern unsigned int minEdgeLengthToBe1MinFlow;			

// Minimum overlap length difference to clip branches (default: 15)
extern unsigned int minOvlDiffToClip;	
// Minimum fold difference to consider branches to be short (default: 5)
extern unsigned int minFoldToBeShortBranch;	
// Minimum size to consider branches to be short (default: 5)
extern unsigned int minSizeToBeShortBranch;

//Minumum unique mate pair support to join edge (default: 3)
extern unsigned int minUinqSupport;
//Minumum non-unique mate pair support to join edge (default: 0)
extern unsigned int minNonUniqSupport;

// minimum contig length to be reported (default: 1000)
extern unsigned int minContigLengthTobeReported;

// minimum reads in a contig to be reported (default: 5)
extern unsigned int minNumberofReadsTobePrinted;

// Fraction of reads used before iteration
extern double maxReadsUsed;

// print contigs or not
extern bool printContigs;

// print scaffolds or not
extern bool printScaffolds;

 // print unused reads or not
extern bool printUnused;

// print GFA/GFA2 graph
extern bool printGFA;
extern bool printGFA2;

//============================================================================


//============================================================================
// variables typedef
//============================================================================
typedef unsigned char UINT8;
typedef signed char INT8;
typedef unsigned short UINT16;
typedef short INT16;
typedef unsigned long UINT32;
typedef long INT32;
typedef unsigned long long UINT64;
typedef long long INT64;


//============================================================================
//	Exit code that displays the place of exit and message.
//============================================================================
#define MYEXIT(a) { cerr << endl << "Exit from File: " << __FILE__ << " Line: " << __LINE__ << " Function: " << __FUNCTION__ << "()" << endl << "Message: " << a << endl; exit(0);}


//============================================================================
// Clock logging for debugging
//============================================================================
#define CLOCKSTART double begin = omp_get_wtime(); INT64 mem_start = checkMemoryUsage(); \
		FILE_LOG(logDEBUG) << endl << ">>> Function start: "<< __FUNCTION__ << "()" << endl;
#define CLOCKSTOP double end = omp_get_wtime();  INT64 mem_end = checkMemoryUsage(); \
		FILE_LOG(logINFO) << "<<< Function stop: " << __FUNCTION__ << "(), Elapsed time: " << double(end - begin)<< " seconds, Memory usage: " \
		<< mem_end << " - " <<  mem_start << " = "<< mem_end - mem_start << " MB."<< endl << "----" << endl;
#define PRINTMEM FILE_LOG(logDEBUG) << checkMemoryUsage() << " MB used" << endl;

//#define CLOCKSTART double begin = omp_get_wtime(); cout<<"Currently in file: " << __FILE__ << " Function: "<< __FUNCTION__ << "()" << endl;
//#define CLOCKSTOP double end = omp_get_wtime(); cout << "Function " << __FUNCTION__ << "() finished in " << double(end - begin) << " Seconds." << endl << endl;

// Get the memory usage with a Linux kernel.
inline unsigned int checkMemoryUsage()
{
    // get KB memory into count
    unsigned int count=0;

    #if defined(__linux__)
    ifstream f("/proc/self/status"); // read the linux file
    while(!f.eof()){
        string key;
        f>>key;
        if(key=="VmData:"){     // size of data
            f>>count;
        break;
        }
    }
    f.close();
    #endif

    // return MBs memory (size of data)
    return (count/1024);
}


//============================================================================
// Class of Config
//=============================================================================
class Config
{

public:
    // Default constructor
    Config();

    // Default destructor
    ~Config();

    // Variables
    static vector<string> readSingleFilenamesList;
    static vector<string> readPairedFilenamesList;
    static vector<string> readInterPairedFilenamesList;
    static vector<string> edgeFilenamesList;
    static string outFilenamePrefix;
    static vector<string> containedReadsFile;
    static UINT64 parallelThreadPoolSize;
    static string simplifyGraphPath;
    static string paramFileName;
    static string paramFileName2;
    static string paramFileName3;

    // Get options
    static bool setConfig(int argc, char **argv);

    //Read and set parameter file
    static void setParameters(string pFile);

    // Print help
    static void printHelp();

    // Get read filenames
    static vector<string> getSingleReadFilenames() {return readSingleFilenamesList;}

    // Get paired-end read filenames
    static vector<string> getPairedReadFilenames() {return readPairedFilenamesList;}

    // Get paired interleaved read filenames
    static vector<string> getInterPairedReadFilenames() {return readInterPairedFilenamesList;}

    // Get edge filenames
    static vector<string> getEdgeFilenames() {return edgeFilenamesList;}

    // Get contained read filenames
    static vector<string> getContainedReadsFile() {return containedReadsFile;}

    // Get output prefix name
    static string getOutputFilenamePrefix() {return outFilenamePrefix;}

    // Get verbosity level of log messages
    static std::string getLogLevel(){return FILELog::ToString(FILELog::ReportingLevel());}

    // Get maximum allowed threads
    static UINT64 getThreadPoolSize() {return parallelThreadPoolSize;}

    // Get output prefix name
    static string getSimplifyGraphPath() {return simplifyGraphPath;}

    //Get parameter file name
    static string getParamFileName() {return paramFileName;}

    //Get parameter file name
    static string getParamFileName2() {return paramFileName2;}

    //Get parameter file name
    static string getParamFileName3() {return paramFileName3;}

    //
};


#endif

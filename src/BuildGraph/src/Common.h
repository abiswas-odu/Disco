/*
 * Common.h
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider, Abhishek Biswas
 */


#ifndef COMMON_H_
#define COMMON_H_

#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <time.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <streambuf>
#include <sys/stat.h>
#include <map>
#include <string>
#include <cstring>
#include <queue>
#include <omp.h>
#include <cmath>
#include <ctime>
#include <functional>
#include <memory>
#include <zlib.h>

using namespace std;

typedef unsigned char UINT8;
typedef unsigned short UINT16;
typedef short INT16;
typedef unsigned long UINT32;
typedef long INT32;
typedef unsigned long long UINT64;
typedef long long INT64;

//Multi-thread parallel options
#define MEGA_PAR_GRAPH_SIZE 80000
#define MAX_PAR_GRAPH_SIZE 40000
#define MID_PAR_GRAPH_SIZE 20000
#define MIN_PAR_GRAPH_SIZE 1000

#define MEGA_THD_MEMORY_SIZE 20000
#define MAX_THD_MEMORY_SIZE 10000
#define MID_THD_MEMORY_SIZE 5000
#define MIN_THD_MEMORY_SIZE 1000


//	Exit code that displays the place of exit and message.
#define MYEXIT(a) { cout << endl << "Exit from File: " << __FILE__ << " Line: " << __LINE__ << " Function: " << __FUNCTION__ << "()" << endl << "Message: " << a << endl; exit(0);}
// Print which function is currently executing. Only for functions that take long time


#define SSTR( x ) dynamic_cast< std::ostringstream & >(( std::ostringstream() << std::dec << x ) ).str()

// To keep time information of functions.
#define CLOCKSTART INT64 mem_start = checkMemoryUsage(); double begin = omp_get_wtime(); cout<<"Currently in file: " << __FILE__ << " Function: "<< __FUNCTION__ << "()" << endl;
#define CLOCKSTOP INT64 mem_end = checkMemoryUsage(); double end = omp_get_wtime(); cout << "Function " << __FUNCTION__ << "() finished in " <<(end - begin)<< " Seconds." << endl << "Memory used: " << mem_end << " - " <<  mem_start << " = "<< mem_end - mem_start << " MB."<< endl << endl;

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
};

// Get the maximum memory of the machine.
inline unsigned int getMaxMemory()
{
    // get KB memory into count
    unsigned int count=0;

    #if defined(__linux__)
    ifstream f("/proc/meminfo"); // read the linux file
    while(!f.eof()){
        string key;
        f>>key;
        if(key=="MemTotal:"){     // size of data
            f>>count;
        break;
        }
    }
    f.close();
    #endif

    // return GBs memory (size of data)
    return (count/(1024*1024));
};


inline void str_reverse( char *str ) {
    char *str_end = strchr( str, 0 );
    std::reverse( str, str_end );
}

inline void split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}
inline std::vector<std::string> splitTok(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
    return elems;
}

// trim from start (in place)
inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
}

// trim from end (in place)
inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

// trim from both ends (in place)
inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

// trim from start (copying)
inline std::string ltrimmed(std::string s) {
    ltrim(s);
    return s;
}

// trim from end (copying)
inline std::string rtrimmed(std::string s) {
    rtrim(s);
    return s;
}

// trim from both ends (copying)
inline std::string trimmed(std::string s) {
    trim(s);
    return s;
}

#endif /* COMMON_H_ */

#ifndef UTILS_HPP
#define UTILS_HPP

/********************************************************************
 * utils.hpp
 *
 * utils.cpp includes useful common uility functions.
 * 
 * Created by JJ Chai  on 10/17/2013 based on Ted's utils.hpp.
 * Copyright (c) 2013 JJ Chai (ORNL). Allrights reserved.
********************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <sys/stat.h>
#include <stdio.h>
#include <string>

#include "Config.h"
// for getting current working directory
#ifdef WINDOWS
	#include <direct.h>
	#define GetCurrentDir _getcwd
#else
	#include <unistd.h>
	#define GetCurrentDir getcwd
#endif

#define	BUFFER_SIZE 100			/* buffer size for reading/writing a stream */

/**
 * utils: utility functions for file (directory) path, filename,
 *        base filename without extension, and format changing.
 */
namespace Utils
{
    /** get path separator **/
    std::string getPathSeparator();

    /** check file exist **/
    bool isFileExist(const std::string & filename);

    /** if file exists, then remove the file **/
    void ifFileExistRemove(const std::string & filepath);
    
    /** create directory if it did not exist **/
    void mkdirIfNonExist(const std::string & dirname);
    
    /** check if file is empty **/
    bool isFileEmpty(const std::string & filename);
    
    /** if file is empty, remove the file **/
    void ifFileEmptyRemove(const std::string & filename);

    /** get base name from program path, string after the last path separator **/
    std::string getBaseName(const std::string & path_str);

    /** get directory from path, string until the last path separator **/
    std::string getDir(const std::string & path_str);

    /** if error, then print error and exit program **/
    void exitWithError(const std::string & error);

    /** method to convert string to int **/
    int stringToInt(const std::string & input_string);

    /** method to convert string to unsigned int **/
    unsigned int stringToUnsignedInt(const std::string & input_string);

     /** method to convert string to float **/
    float stringToFloat( const std::string& input_string );

    /** method to convert string to double **/
    double stringToDouble(const std::string & input_string);

    /** method to convert int to string **/
    std::string intToString(const int & input_value);

    /** method to convert unsigned int to string **/
    std::string unsignedIntToString(const unsigned int & input_value);

    /** method to convert double to string **/
    std::string doubleToString(const double & input_value);

    /** check directory **/
    bool isDirectory(const std::string & directory);

    /** check fasta format **/
    bool isFastaFormat(const std::string & input_string);

    /** check fastq format **/
    bool isFastqFormat(const std::string & input_string);
    
    /** get filebase from filename without extension and path **/
    std::string getFilename(const std::string & filename);

    /** get filebase from filename without extension **/
    std::string getFilebase(const std::string & filename);

    /** get filebase from filename without extension **/
    std::string getFilebase2(const std::string & filename);

    /** Get current date/time, format is YYYY-MM-DD.HH:mm:ss **/
    std::string currentDateTime();

    /** string vector to delemiter separated string **/
    std::string vectorToString(const std::vector<std::string> & input_vector, const std::string delemit);
    
    /** delimited string to string vector **/
    std::vector<std::string> StringToVector(const std::string & str, const char & delimit);
    
    /** determine if this filename matches this pattern **/
    bool matchPattern(  const std::string & pattern, const std::string & filename );
    
    /** find all files with specified extensions in a directory **/
    //void getFiles(const std::string & directory, const std::string & extension, std::vector<std::string> & filenames);
    
    /** find the last non-empty line of a file **/
    std::string last_line(const std::string & filename);

    /* Get ref sequence names and save them in a vector of strings, from stream */
    bool getRefNames(FILE * stream, std::vector<std::string> & refNames, std::vector<unsigned short> & refLens);

    /* Get current working directory */
    std::string get_cwd();

    /* Read each line in a file and store in a vector of strings */
    bool saveLinesToVec(const std::string & fileName, std::vector<std::string> & stringVector, std::vector<unsigned short> & lengthVector);

    /* Get length of pacbio read from its name */
    unsigned short getPacBioLength(const std::string & readName);

    /* Generate reverse complement of a DNA sequence*/
    std::string reverseComplement(const std::string & seq);

    /* Split a string by character delimiter */
    std::vector<std::string> split(const std::string &s, char delim);

    // trim from start (in place)
    inline void ltrim(std::string &s);
    // trim from end (in place)
    inline void rtrim(std::string &s);
    // trim from both ends (in place)
    inline void trim(std::string &s);
    // trim from start (copying)
    inline std::string ltrimmed(std::string s);
    // trim from end (copying)
    inline std::string rtrimmed(std::string s);
    // trim from both ends (copying)
    std::string trimmed(std::string s);

    void writeCheckPointFile(std::string allFileNamePrefix, std::string message);

    struct compare {
        bool operator()(const std::string& first, const std::string& second) {
            return first.size() < second.size();
        }
    };
}


#endif //UTILS.H

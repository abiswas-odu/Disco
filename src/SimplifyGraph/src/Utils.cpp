/********************************************************************
 * utils.cpp
 *
 * utils.cpp includes useful common uility functions.
 * 
 * Created by JJ Chai on 10/17/2013 based on Ted's script.
 * Copyright (c) 2013 JJ Chai (ORNL). Allrights reserved.
********************************************************************/

#include "Utils.h"
#include <unistd.h>

namespace Utils
{

    /** get path separator **/
    std::string getPathSeparator()
    {
    #if _WIN32
        return "\\" ;
    #else
        return "/" ;
    #endif
    }

    /** check file exist **/
    bool isFileExist(const std::string & filename)
    {
        // infile
        std::ifstream infile(filename.c_str());

        // return true if none of the stream's error state flags (eofbit, failbit and badbit) is set
        return infile.good();
    }
   
    /** if file exists, then remove the file **/
    void ifFileExistRemove(const std::string & filename)
    {
        // if file exists,
        if (isFileExist(filename))
        {
            // then remove if there is no error
            if (remove(filename.c_str()) != 0)
                exitWithError("*** Error: Failed to remove: " + filename);
        }
    }

    /** create directory if it did not exist **/
    void mkdirIfNonExist(const std::string & dirname)
    {
        if (!Utils::isDirectory(dirname)) {
		std::string cmd = "mkdir -p "+ dirname;
            const int res = system(cmd.c_str());
            if (res != 0 )
                exitWithError("*** Error: Failed command: " + cmd);
            
        }
    }
    
    /** check if file is empty **/
    bool isFileEmpty(const std::string & filename)
    {

        
        //check if file exists first
        if (!isFileExist(filename)) {
            //std::cout << "*** file " << filename << "does not exist .... " << std::endl ;
            return true;
        }
        else
        {
            //infile
            std::ifstream infile(filename.c_str());
            
            //check file size
            infile.seekg(0,std::ios::end);
            size_t size = infile.tellg();
            infile.close();
            
            if( size == 0 )
                return true;
            else
                return false;
        }
        
    }
    
    /** if file is empty, remove the file **/
    void ifFileEmptyRemove(const std::string & filename)
    {
        // if file is empty
        if (isFileEmpty(filename))
        {
            // then remove if there is no error
            if (remove(filename.c_str()) != 0)
                exitWithError("*** Error: Failed to remove: " + filename);
        }
    }
    
    /** get base name from path, string after the last path separator **/
    std::string getBaseName(const std::string & path_str)
    {
        std::string baseName;

        // position
        size_t last_separator_pos = path_str.find_last_of(getPathSeparator());

        // find path separator
        if (last_separator_pos != std::string::npos)
            baseName.assign(path_str.begin()+last_separator_pos+1,path_str.end());
        
        // couldn't find path separator
        else
            baseName = path_str;

        // return 
        return baseName;
    }

    /** get directory from path, string until the last path separator **/
    std::string getDir(const std::string & path_str)
    {
        std::string dirName;

        // position
        size_t last_separator_pos = path_str.find_last_of(getPathSeparator());

        // find path separator
        if (last_separator_pos != std::string::npos)
            dirName.assign(path_str.begin(),path_str.begin()+last_separator_pos+1);
        // couldn't find path separator
        else
            dirName = "." + getPathSeparator();

        // return 
        return dirName;
    }

    /** exit program when error happens **/
    void exitWithError(const std::string &error) 
    {
        std::cerr << error;
        exit(EXIT_FAILURE);
    }

    /** method to convert string to int **/
    int stringToInt( const std::string &input_string ) 
    {
        std::istringstream input_string_stream(input_string);
        int x;
        if (!(input_string_stream >> x)) {
            std::cerr << input_string;
            exitWithError("*** Error: string was not converted to int correctly!\n");
        }   

        // return 
        return x;
    } 

    /** method to convert string to unsigned int **/
    unsigned int stringToUnsignedInt( const std::string &input_string ) 
    {
        std::istringstream input_string_stream(input_string);
        unsigned int x;
        if (!(input_string_stream >> x)) {
            exitWithError("*** Error: string was not converted to unsigned int correctly!\n");
        }   

        // return 
        return x;
    } 

    /** method to convert string to float **/
    float stringToFloat( const std::string& input_string ) 
    {
        std::istringstream input_string_stream(input_string);
        float x;
        if (!(input_string_stream >> x)) {
            std::cerr << "input string is " << input_string << "\n";
            exitWithError("*** Error: string was not converted to float correctly!");
        }   

        // return 
        return x;
    } 
    /** method to convert string to double **/
    double stringToDouble( const std::string& input_string ) 
    {
        std::istringstream input_string_stream(input_string);
        double x;
        if (!(input_string_stream >> x)) {
            std::cerr << "input string is " << input_string << "\n";
            exitWithError("*** Error: string was not converted to double correctly!");
        }   

        // return 
        return x;
    } 

    /** method to convert int to string **/
    std::string intToString( const int & input_value )
    {
        std::string input_value_str;
        std::ostringstream convert;
        if (!(convert << input_value)) {
            exitWithError("*** Error: unsigned int was not converted to string correctly!");
        }   
        input_value_str = convert.str();

        // return 
        return input_value_str;
    }

    /** method to convert unsigned int to string **/
    std::string unsignedIntToString( const unsigned int & input_value )
    {
        std::string input_value_str;
        std::ostringstream convert;
        if (!(convert << input_value)) {
            exitWithError("*** Error: unsigned int was not converted to string correctly!");
        }   
        input_value_str = convert.str();

        // return 
        return input_value_str;
    }

    /** method to convert double to string **/
    std::string doubleToString( const double & input_value )
    {
        std::string input_value_str;
        std::ostringstream convert;
        if (!(convert << input_value)) {
            exitWithError("*** Error: unsigned int was not converted to string correctly!");
        }   
        input_value_str = convert.str();

        // return 
        return input_value_str;
    }

    /** check directory **/
    bool isDirectory( const std::string & directory )
    {
        struct stat sb;
        if (stat(directory.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode) )
            return true;
        else
            return false;
    }

    /** check fasta format **/
    bool isFastaFormat(const std::string & input_string)
    {

        // list of fasta format
        std::vector<std::string> file_format;
        file_format.push_back("fasta");
        file_format.push_back("fa");
        file_format.push_back("fas");
        file_format.push_back("fna");
        file_format.push_back("ffn");

        // get file extension
        std::string file_ext = input_string.substr(input_string.find_last_of(".") + 1);

        // transform to lower characters
        std::transform(file_ext.begin(), file_ext.end(), file_ext.begin(), ::tolower);

        // loop list
        for(unsigned int i=0; i < file_format.size(); i++) {
            if (file_ext == file_format[i]) 
                return true;
        }
        // return 
        return false;
    }

    /** check fastq format **/
    bool isFastqFormat(const std::string & input_string)
    {

        // list of fasta format
        std::vector<std::string> file_format;
        file_format.push_back("fastq");
        file_format.push_back("fq");
        file_format.push_back("faq");

        // get file extension
        std::string file_ext = input_string.substr(input_string.find_last_of(".") + 1);
        std::transform(file_ext.begin(), file_ext.end(), file_ext.begin(), ::tolower);

        // loop list
        for(unsigned int i=0; i < file_format.size(); i++) {
            if (file_ext == file_format[i]) 
                // return 
                return true;
        }
        // return 
        return false;
    }

    /** get filebase from filename without extension and path **/
    std::string getFilename(const std::string & filename)
    {
        std::string filebase;
        filebase = getBaseName(filename);
        filebase = getFilebase(filebase);
        return filebase;
    }

    /** get filebase from filename without extension **/
    std::string getFilebase(const std::string & filename)
    {
        std::string filebase;

        size_t pos = filename.find_last_of(".");
        if (pos != std::string::npos)
            filebase.assign(filename.begin(),filename.begin()+pos);
        else
            filebase = filename;

        // return 
        return filebase;
    }

    /** get filebase from filename i.e., AAA.CC.X -> AAA **/
    std::string getFilebase2(const std::string & filename)
    {
        std::string filebase;

        // parse base filename
        std::istringstream filename_stream(filename);
        std::vector<std::string> fields_vector;
        for (std::string field; getline(filename_stream, field, '.'); ) { 
           fields_vector.push_back(field);
        }   
        if (fields_vector.size() > 2) {
            filebase = ""; 
            for(unsigned int i = 0; i < fields_vector.size() - 2; i++) {
                filebase += fields_vector[i] + '.';
            }   
            filebase = filebase.substr(0, filebase.size()-1);
        }   
        else {
            filebase = filename;
        }

        // return 
        return filebase;
    }

    /** Get current date/time, format is YYYY-MM-DD.HH:mm:ss **/
    std::string currentDateTime() 
    {
        time_t     now = time(0);
        struct tm  tstruct;
        char       buf[80];
        tstruct = *localtime(&now);
        strftime(buf, sizeof(buf), "[%Y-%m-%d.%X]", &tstruct);

        // return 
        return buf;
    }

    /** string vector to delemiter separated string **/
    std::string vectorToString(const std::vector<std::string> & input_vector, const std::string delemit)
    {
        // variables
        std::stringstream ss;
        std::string comma_string;

        // loop list
        for (unsigned int i=0; i<input_vector.size(); i++ ) {
            if (i != 0)
                ss << delemit;
            ss << input_vector[i];
        }
  
        // comma_string
        comma_string = ss.str();

        // return 
        return comma_string;
    }
    
    /** delimiter separated string to string vector **/
    std::vector<std::string> StringToVector(const std::string & str, const char & delimit)
    {
        std::vector<std::string> strVec;
        std::size_t pos = 0;
        std::size_t delimiter_pos = str.find(delimit,pos);	// position of the delimiter
	std::string field;	// extract field
        while (delimiter_pos != std::string::npos) {
		field = str.substr(pos, delimiter_pos-pos);
            strVec.push_back(field);
            pos = delimiter_pos + 1;
            delimiter_pos = str.find(delimit,pos);
        }
        field.assign(str.begin()+pos,str.end());
        strVec.push_back(field);
        return strVec;
    }
    
    /** determine if this filename matches this pattern **/
    bool matchPattern(  const std::string & pattern, const std::string & filename )
    {
        // if filename is ".", "..", or "lost+found", return false
        if ( (filename.compare(".") == 0) ||
             (filename.compare("..") == 0) ||
             (filename.compare("lost+found") == 0 ) ||
             (filename.length() <= pattern.size())
            ){
            return false;
        }
        
        // if pattern is empty or "*", match all files
        if ((pattern.compare("") == 0) || (pattern.compare("*") == 0)) {
            return true;
        }
        // get file extension
        std::string file_ext = filename.substr(filename.find_last_of(".") + 1);
        
        // compare extension with the pattern
        if (file_ext.compare(pattern) == 0)
            return true;
        else
            return false;
    }
    
    /** find all files with specified extensions in a directory **/
    /*
    void getFiles(const std::string & directory, const std::string & extension,
                  std::vector<std::string> & filenames)
    {
        DirectoryStructure searchDir(directory);
        
        searchDir.setPattern(extension);
        searchDir.getFiles(filenames);
        
        int i, fileNum;
        fileNum = (int) searchDir.getFileCount();
        vector<string>::reverse_iterator rit = filenames.rbegin();
        // full path of the input files and the output files
        for (i = 0; i < fileNum; rit++,i++) {
            //cout << *rit << " -> ";
            *rit = directory + getPathSeparator() + *rit;
            //cout << *rit << endl;
        }

    }
    */
    
    /** find the last non-empty line of a file **/
    std::string last_line(const std::string & filename)
    {
        std::ifstream fin(filename.c_str());
        std::string lastline;
        
        if(fin)
        {
            char ch;
            
            // go to one spot before the EOF
            fin.seekg(-1,std::ios_base::end);
            
            fin.get(ch);
            
            // skip the last EOF
            if (ch =='\n')
            {
                //std::cout << "file ended with EOF \n";
                fin.seekg(-2,std::ios_base::cur);
            }
            else
            {
                fin.unget();
                fin.seekg(-2,std::ios_base::cur);
            }
            
            // indicator of end of looping
            bool keepLooping = true;
            
            while (keepLooping) {
                
                // get current byte's data
                fin.get(ch);

                // if the data was at or before the beginning
                // first line is the last line
                // stop here
                if ((int)fin.tellg() <= 1) {
                    //std::cout << "fin.tellg = " << (int)fin.tellg() << std::endl;
                    fin.seekg(0);
                    keepLooping = false;
                    //std::cout << "now at beginning \n";
                }
                // if the data is a new line, stop
                else if (ch == '\n')
                {
                    keepLooping = false;
                    //std::cout << "now at a new line \n";
                }
                
                // keep moving until the first 2 cases
                else
                {
                    fin.seekg(-2,std::ios_base::cur);
                }
            }
            
            //std::cout << "looping ended" << std::endl;
            
            getline(fin,lastline);
            
            fin.close();
        }
        
        
        return lastline;
    }


    /* Get ref sequence names and length, save them in a vector of strings and vector of unsigned integers, from stream */
    bool getRefNames(FILE * stream, std::vector<std::string> & refNames, std::vector<unsigned short> & refLens)
    {
	    std::string line = "";
	    char buffer[BUFFER_SIZE];
	    while(1)
	    {
		    if(fgets(buffer, BUFFER_SIZE, stream)==NULL) /* eof is reached before any characters could be read, or a read error occurs */
		    {
			    //perror("In function GetRefNames, Encountered in function fgets");
			    break;
		    }
		    
		    line += buffer;
		   
		    if (line.at(line.size() - 1) == '\n')        /* A whole line has been read, processing it */
		    {
			    line.erase(line.size()-1, 1);        /* Delete last character from line, which is EOL */
			    if ( line.substr(0,3).compare("@SQ") ==0) /* SQ line */
			    {
				    size_t SNPos = line.find("SN:");
				    size_t nextTabPos = line.find("\t", SNPos+3);
				    refNames.push_back(line.substr(SNPos+3,nextTabPos-SNPos-3)); /* length of the PacBio read */
				    size_t LNPos = line.find("LN:");
				    nextTabPos = line.find("\t", LNPos+3);
				    refLens.push_back(std::stoi(line.substr(LNPos+3,nextTabPos-LNPos-3))); /* length of the PacBio read */
			    }
			    line = "";              /* Set line to empty again */
		    }
	    }
	    return true;
    }

    /* Get current working directory */
    std::string get_cwd()
    {
	char path[FILENAME_MAX];
	return ( GetCurrentDir(path, sizeof(path)) ? std::string(path) : std::string("") );
    }

    /* Read each line in a file and store in a vector of strings */
    bool saveLinesToVec(const std::string & fileName, std::vector<std::string> & stringVector, std::vector<unsigned short> & lengthVector)
    {
	    stringVector.clear();
	    lengthVector.clear();
	    std::ifstream fin(fileName.c_str());
	    size_t delimit_pos;
	    if(fin.is_open())
	    {
		    std::string line;
		    while(std::getline(fin, line))
		    {
			    if (line.length() > 0)
			    {
				    delimit_pos = line.find("\t");
				    stringVector.push_back(line.substr(0, delimit_pos));
				    //lengthVector.push_back(std::stoi(line.substr(delimit_pos+1, std::string::npos)));
			    }
		    }
		    fin.close();
	    }

	    return true;
    }

    /* Get length of pacbio read from its name */
    unsigned short getPacBioLength(const std::string & readName)
    {
	    size_t lastSlash = readName.find_last_of("/");
	    size_t lastUnderscore = readName.find_last_of("_");
	    unsigned short startPos = std::stoi(readName.substr(lastSlash+1, lastUnderscore-lastSlash));
	    unsigned short endPos = std::stoi(readName.substr(lastUnderscore+1, std::string::npos));
	    return (endPos - startPos);
    }
    /*
     * ===  FUNCTION  ======================================================================
     *         Name:  reverseComplement
     *  Description:  Returns the reverse complement of a read.
     * =====================================================================================
     */
    std::string reverseComplement(const std::string & seq)
    {
    	uint64_t sLength = seq.length();
    	std::string reverse(sLength,'0');
    	for(uint64_t i = 0;i < sLength; i++)
    	{
    		if( seq[i] & 0X02 ) // C or G
    			reverse.at(sLength -  i - 1 ) = seq[i] ^ 0X04;
    		else // A <==> T
    			reverse.at(sLength -  i - 1 ) = seq[i] ^ 0X15;
    	}
    	return reverse; // return the reverse complement as a string
    }

    void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    	std::stringstream ss(s);
    	std::string item;
        while(std::getline(ss, item, delim)) {
            elems.push_back(item);
        }
    }
    std::vector<std::string> split(const std::string &s, char delim) {
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
    std::string trimmed(std::string s) {
        trim(s);
        return s;
    }

    void writeCheckPointFile(std::string allFileNamePrefix, std::string message)
    {
    	//Write checkpoint file and set graph construction complete
		std::ofstream filePointer;
		std::string fileName = allFileNamePrefix+"_SimplificationCheckpointInfo.txt";
		filePointer.open(fileName.c_str(), std::ios_base::app);
		if(!filePointer)
			MYEXIT("Unable to open file: "+fileName);
		filePointer<<message<<std::endl;
		filePointer.close();
    }
    void populateThresh(std::map<UINT64,UINT64> *ref)
    {
    	ref->insert(std::pair<UINT64,UINT64>(22286068,60000));
    	ref->insert(std::pair<UINT64,UINT64>(107718722,62300));
    	ref->insert(std::pair<UINT64,UINT64>(106998276,62300));
    	ref->insert(std::pair<UINT64,UINT64>(770370712,102100));
    	ref->insert(std::pair<UINT64,UINT64>(146001066,2800));
    }
}

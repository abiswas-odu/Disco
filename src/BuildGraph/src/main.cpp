/*
 * main.cpp
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider, Abhishek Biswas
 * Version: 3.o
 */


#include "../src/Common.h"
#include "../src/Dataset.h"
#include "../src/Edge.h"
#include "../src/HashTable.h"
#include "../src/OverlapGraph.h"
#include "../src/Read.h"

void parseArguments(int argc, char **argv, vector<string> & pairedEndFileNames, vector<string> & singleEndFileNames,
		string & allFileName,UINT64 & maxthreads, UINT64 & writeGraphSize, UINT64 &maxMemSizeGB,string &parameterFile);

UINT64 readOverlapParameter(string parameterFile);

int main(int argc, char **argv)
{
	CLOCKSTART;

	vector<string> pairedEndFileNames, singleEndFileNames;
	string allFileName;
	string parameterFile;
	UINT64 maxThreads = omp_get_max_threads();
	UINT64 writeGraphSize = MID_PAR_GRAPH_SIZE;
	UINT64 maxMemSizeGB = getMaxMemory();
	cout<<"Max available memory: "<<maxMemSizeGB<< " GB"<<endl;
	parseArguments(argc, argv, pairedEndFileNames, singleEndFileNames, allFileName,
			maxThreads, writeGraphSize,maxMemSizeGB,parameterFile);
	UINT64 minimumOverlapLength = readOverlapParameter(parameterFile);
	cout<<"Max usable memory: "<<maxMemSizeGB<< " GB"<<endl;

	Dataset *dataSet = new Dataset(pairedEndFileNames, singleEndFileNames, allFileName, minimumOverlapLength);
	OverlapGraph *overlapGraph;
	HashTable *hashTable=new HashTable();
	hashTable->insertDataset(dataSet, minimumOverlapLength,maxThreads);
	overlapGraph=new OverlapGraph(hashTable,maxThreads,writeGraphSize,maxMemSizeGB,allFileName); //hashTable deleted by this function after building the graph also writes graph
	delete dataSet;
	delete overlapGraph;
	CLOCKSTOP;
}

/**********************************************************************************************************************
	Parse the input arguments
**********************************************************************************************************************/

void parseArguments(int argc, char **argv, vector<string> & pairedEndFileNames, vector<string> & singleEndFileNames,
		string & allFileName,UINT64 & maxthreads, UINT64 & writeGraphSize, UINT64 &maxMemSizeGB,string &parameterFile)
{
	allFileName = "";
	vector<string> argumentsList;
	cout << "PRINTING ARGUMENTS" << endl;
	for(int i = 0; i < argc; i++)
	{
		cout << argv[i] << ' ';
	}
	cout << endl;
	while(argc--)
			argumentsList.push_back(*argv++);

	if(argumentsList.size() == 1)
	{
		cerr << endl << "Usage: buildG [OPTION]...[PRARAM]..." << endl;
		cerr << "  -pe\tnumber of files and paired-end file names" <<endl; 			// Paired-end file name in fasta/fastq format. mate pairs should be one after another in the file.
		cerr << "  -se\tnumber of files and single-end file names" <<endl; 			// Single-end file name in fasta/fastq format.
		cerr << "  -f\tAll file name prefix" <<endl; 			// all output file with have this name with different extensions.
		cerr << "  -t\tmaximum threads used" << endl; 	// Maximum OMP threads used
		cerr << "  -m\tmaximum memory usage allowed (default: max available; info: use atleast disk space size of reads + 1GB per thread specified by -t)" << endl; 	// Maximum memory to be used before written to disk
		exit(0);
	}

	for(UINT64 i = 1; i <= argumentsList.size()-1; i++)
	{
		if(argumentsList[i] == "-pe")
		{
			UINT64 numberOfPairedEndDatasets = atoi(argumentsList[++i].c_str());
			for(UINT64 j = 0; j < numberOfPairedEndDatasets; j++)
			{
				pairedEndFileNames.push_back(argumentsList[++i]);
			}
		}
		else if(argumentsList[i] == "-se")
		{
			UINT64 numberOfSingleEndDatasets = atoi(argumentsList[++i].c_str());
			for(UINT64 j = 0; j < numberOfSingleEndDatasets; j++)
			{
				singleEndFileNames.push_back(argumentsList[++i]);
			}
		}
		else if (argumentsList[i] == "-f")
			allFileName = argumentsList[++i];
		else if (argumentsList[i] == "-t")
			maxthreads = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-w")
			writeGraphSize = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-m")
			maxMemSizeGB = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-p")
			parameterFile = argumentsList[++i];
		else
		{
			cerr << endl << "Usage: buildG [OPTION]...[PRARAM]..." << endl;
			cerr << "  -pe\tnumber of files and paired-end file names" <<endl; 			// Paired-end file name in fasta/fastq format. mate pairs should be one after another in the file.
			cerr << "  -se\tnumber of files and single-end file names" <<endl; 			// Single-end file name in fasta/fastq format.
			cerr << "  -f\tAll file name prefix" <<endl; 			// all output file with have this name with different extensions.
			cerr << "  -t\tmaximum threads used" << endl; 	// Maximum OMP threads used
			cerr << "  -m\tmaximum memory usage allowed (default: max available; info: use at least disk space size of reads + 1GB per thread specified by -t)" << endl; 	// Maximum memory to be used before written to disk
			if (argumentsList[i] == "-h" || argumentsList[i] == "--help")
				exit(0);
			else
			{
				cerr << "Unknown option: " << argumentsList[i] << endl << endl;
				exit(1);
			}
		}
	}
}

UINT64 readOverlapParameter(string parameterFile)
{
	ifstream filePointer;
	filePointer.open(parameterFile.c_str());
	UINT64 minimumOverlapLength=30;
	if(!filePointer.is_open()) {
		cerr<<"Unable to open parameter file: "<<parameterFile<<endl;
		exit(1);
	}
	else {
		string par_text="";
		while(getline(filePointer,par_text)) {

			if(par_text.find("=") != std::string::npos)
			{
				vector<string> tok = splitTok(par_text,'=');
				string parName = trimmed(tok[0]);
				string parVal = trimmed(tok[1]);
				if(parName=="MinOverlap4BuildGraph")
					minimumOverlapLength=stoi(parVal);
			}
		}
	}
	return minimumOverlapLength;
}

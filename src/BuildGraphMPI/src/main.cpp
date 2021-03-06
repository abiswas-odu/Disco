/*
 * main.cpp
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider, Abhishek Biswas
 * Version: 3.0
 */


#include "Common.h"
#include "Dataset.h"
#include "Edge.h"
#include "HashTable.h"
#include "OverlapGraph.h"
#include "Read.h"

void parseArguments(int argc, char **argv, vector<string> & pairedEndFileNames, vector<string> & singleEndFileNames,
		string & graphFileName,string &fnamePrefixSimplify, string &simplifyPartialPath,
		UINT64 & maxthreads, UINT64 & writeGraphSize, UINT64 &maxMemSizeGB,string &parameterFile);

UINT64 readOverlapParameter(string parameterFile);

int main(int argc, char **argv)
{
	int numprocs=1, myid=0, len, provided=0;
	double start, end;
	char name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
	if (provided < MPI_THREAD_FUNNELED)
	{
	   printf("Error: the MPI library doesn't provide the required thread level\n");
	   MPI_Abort(MPI_COMM_WORLD, 0);
	}
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
	MPI_Get_processor_name(name, &len);
	start = MPI_Wtime();

	cout<<"Software: Disco Assembler (Distributed Computation) [November 2016]"<<endl;
	cout<<"Version : 1.0"<<endl;
	cout<<"Developed by: Biswas, Abhishek; Pan, Chongle et.al."<<endl;
	cout<<"Affiliation: Oak Ridge National Lab / University of Tennessee"<<endl;


	printf("Rank %d running on %s with %d threads.\n", myid, name, omp_get_max_threads());
	vector<string> pairedEndFileNames, singleEndFileNames;
	string graphFileNames;
	string parameterFile;
	UINT64 maxThreads = DEF_THREAD_COUNT;
	UINT64 writeGraphSize = MID_PAR_GRAPH_SIZE;
	UINT64 maxMemSizeGB = getMaxMemory();
	string fnamePrefixSimplify;
	string simplifyPartialPath;
	cout<<"Max available memory: "<<maxMemSizeGB<< " GB"<<endl;
	parseArguments(argc, argv, pairedEndFileNames, singleEndFileNames, graphFileNames,fnamePrefixSimplify,simplifyPartialPath,
				maxThreads, writeGraphSize,maxMemSizeGB,parameterFile);
	UINT64 minimumOverlapLength = readOverlapParameter(parameterFile);
	cout<<"Max usable memory: "<<maxMemSizeGB<< " GB"<<endl;
	Dataset *dataSet = new Dataset(pairedEndFileNames, singleEndFileNames, graphFileNames, minimumOverlapLength,maxThreads);
	HashTable *hashTable=new HashTable();
	hashTable->insertDataset(dataSet, minimumOverlapLength,maxThreads);
	OverlapGraph *overlapGraph=nullptr;
	overlapGraph=new OverlapGraph(hashTable,maxThreads,writeGraphSize,maxMemSizeGB,graphFileNames,fnamePrefixSimplify,
			simplifyPartialPath, myid,numprocs); //hashTable deleted by this function after building the graph also writes graph
	MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
	end = MPI_Wtime();
	delete hashTable;	//  Do not need the hash table any more.
	delete dataSet;
	delete overlapGraph;
	if (myid == 0) { /* use time on master node */
	  printf("Runtime for %d processes = %f\n", numprocs, end-start);
	}
	MPI_Finalize();
}

/**********************************************************************************************************************
	Parse the input arguments
**********************************************************************************************************************/

void parseArguments(int argc, char **argv, vector<string> & pairedEndFileNames, vector<string> & singleEndFileNames,
		string & graphFileName,string &fnamePrefixSimplify, string &simplifyPartialPath,
		UINT64 & maxthreads, UINT64 & writeGraphSize, UINT64 &maxMemSizeGB,string &parameterFile)
{
	graphFileName = "";
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
		cerr << "  -f\tGraph file name prefix" <<endl; 			// graph output file with have this name with different extensions.
		cerr << "  -s\tSimplified edges file name prefix" <<endl; 			// edges output file with have this name with different extensions.
		cerr << "  -simPth\tPartial simplification executable path" <<endl; 			// path to partial simplification exe
		cerr << "  -t\tmaximum threads used" << endl; 	// Maximum OMP threads used
		cerr << "  -m\tmaximum memory usage allowed (default: max available; info: use atleast disk space size of reads + 1GB per thread specified by -t)" << endl; 	// Maximum memory to be used before written to disk
		exit(0);
	}

	for(UINT64 i = 1; i <= argumentsList.size()-1; i++)
	{
		if(argumentsList[i] == "-pe")
		{
			string fileNames = argumentsList[++i];
			vector<string> files = splitTok(fileNames,',');
			for(UINT64 j = 0; j < files.size(); j++)
			{
				pairedEndFileNames.push_back(files[j]);
			}
		}
		else if(argumentsList[i] == "-se")
		{
			vector<string> files = splitTok(argumentsList[++i].c_str(),',');
			for(UINT64 j = 0; j < files.size(); j++)
			{
				singleEndFileNames.push_back(files[j]);
			}
		}
		else if (argumentsList[i] == "-f")
			graphFileName = argumentsList[++i];
		else if (argumentsList[i] == "-t")
			maxthreads = std::stoull(argumentsList[++i].c_str(),nullptr,0);
		else if (argumentsList[i] == "-w")
			writeGraphSize = std::stoull(argumentsList[++i].c_str(),nullptr,0);
		else if (argumentsList[i] == "-m")
			maxMemSizeGB = std::stoull(argumentsList[++i].c_str(),nullptr,0);
		else if (argumentsList[i] == "-p")
			parameterFile = argumentsList[++i];
		else if (argumentsList[i] == "-s")
			fnamePrefixSimplify = argumentsList[++i];
		else if (argumentsList[i] == "-simPth")
			simplifyPartialPath = argumentsList[++i];
		else
		{
			cerr << endl << "Usage: buildG [OPTION]...[PRARAM]..." << endl;
			cerr << "  -pe\tnumber of files and paired-end file names" <<endl; 			// Paired-end file name in fasta/fastq format. mate pairs should be one after another in the file.
			cerr << "  -se\tnumber of files and single-end file names" <<endl; 			// Single-end file name in fasta/fastq format.
			cerr << "  -f\tGraph file name prefix" <<endl; 			// graph output file with have this name with different extensions.
			cerr << "  -s\tSimplified edges file name prefix" <<endl; 			// edges output file with have this name with different extensions.
			cerr << "  -simPth\tPartial simplification executable path" <<endl; 			// path to partial simplification exe
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

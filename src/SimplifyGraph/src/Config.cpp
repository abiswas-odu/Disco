//============================================================================
// Name        : Config.cpp
// Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey, Abhishek Biswas
// Version     : v3
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : Config cpp file
//============================================================================

#include "Config.h"


//=============================================================================
// Default constructor
//=============================================================================
Config::Config() {
}


//=============================================================================
// Default destructor
//=============================================================================
Config::~Config() {
}

//=============================================================================
//Here is to create and initialize the static members
//=============================================================================
vector<string> Config::readSingleFilenamesList;
vector<string> Config::readPairedFilenamesList;
vector<string> Config::readInterPairedFilenamesList;
vector<string> Config::edgeFilenamesList;
string Config::outFilenamePrefix = "disco";
vector<string> Config::containedReadsFile;
string Config::simplifyGraphPath ="";
UINT64 Config::parallelThreadPoolSize(4);
string Config::paramFileName="disco.cfg";
string Config::paramFileName2="disco_2.cfg";
string Config::paramFileName3="disco_3.cfg";

// global variables with default value
unsigned int minOvl=20;

unsigned int minReadsCountInEdgeToBeNotDeadEnd = 2;
unsigned int minEdgeLengthToBeNotDeadEnd = 500;

unsigned int minReadsCountToHave0Flow = 2;
unsigned int minEdgeLengthToHave0Flow = 200;

unsigned int minReadsCountInEdgeToBe1MinFlow = 5;
unsigned int minEdgeLengthToBe1MinFlow = 500;

unsigned int minOvlToClip = 30;

unsigned int minOvlDiffToClip = 10;
unsigned int minFoldToBeShortBranch = 5;
unsigned int minSizeToBeShortBranch = 200;

unsigned int minUinqSupport=20;
unsigned int minNonUniqSupport=0;

unsigned int minContigLengthTobeReported = 300;

unsigned int minNumberofReadsTobePrinted = 2;

double maxReadsUsed = 0.7;

bool printContigs=false;
bool printScaffolds=true;
bool printUnused=false;
bool printGFA = false;
bool printGFA2 =false;
//=============================================================================
// print help usage
//=============================================================================
void Config::printHelp()
{
	cout << endl
		<< "  Usage:" << endl
		<< "    fullsimplify [OPTION]...<PARAM>..." << endl
		<< endl
		<< "  <PARAM>" << std::endl
		<< "    -fs\t contained single read reduction read filename(s) (comma separated fasta/fastq)" << endl
		<< "    -fp\t contained paired-end read reduction read filename(s) in pairs of 2 (comma separated fasta/fastq)" << endl
		<< "    -fpi\t contained interleaved paired-end read reduction read filename(s) (comma separated fasta/fastq)" << endl
		<< "    -crd\t Contained read file (default: none)" << endl
		<< "    -e\t overlapped edge property graph filename(s) (comma separated edge list)" << endl
		<< "    -p\t\t assembly parameter file (default: parameter.cfg)" << endl
		<< "    -o\t\t all output filename prefix" << endl
		<< "    -simPth\t\t path to partial simplification executable" << endl
		<< endl
		<< "  [OPTION]" << std::endl
		<< "    -h/--help\t only print out the help contents" << endl
		<< "    -ovl\t minimum overlap length (default: 0, use all overlap found in edge property graph files)" << endl
		<< "    -log\t verbosity level of log messages: ERROR, WARNING, INFO (default: INFO)" << endl
		<< endl;
}

void Config::setParameters(string pFile)
{
	ifstream filePointer;
	filePointer.open(pFile.c_str());
	if(!filePointer.is_open()) {
		MYEXIT("Unable to open parameter file: "+pFile);
	}
	else {
		string par_text="";
		while(getline(filePointer,par_text)) {
			string par_text_trm = Utils::trimmed(par_text);
			if(par_text_trm.find("=") != std::string::npos && par_text_trm.at(0)!='#')
			{
				vector<string> tok = Utils::split(par_text_trm,'=');
				string parName = Utils::trimmed(tok[0]);
				string parVal = Utils::trimmed(tok[1]);

				if(parName=="minReadsCountInEdgeToBeNotDeadEnd")
					minReadsCountInEdgeToBeNotDeadEnd=stoi(parVal);
				else if(parName=="minEdgeLengthToBeNotDeadEnd")
					minEdgeLengthToBeNotDeadEnd=stoi(parVal);
				else if(parName=="minReadsCountInEdgeToBe1MinFlow")
					minReadsCountInEdgeToBe1MinFlow=stoi(parVal);
				else if(parName=="minEdgeLengthToBe1MinFlow")
						minEdgeLengthToBe1MinFlow=stoi(parVal);
				else if(parName=="minReadsCountToHave0Flow")
					minReadsCountToHave0Flow=stoi(parVal);
				else if(parName=="minEdgeLengthToHave0Flow")
					minEdgeLengthToHave0Flow=stoi(parVal);
				else if(parName== "minSequenceLengthTobePrinted")
					minContigLengthTobeReported=stoi(parVal);
				else if(parName== "minNumberofReadsTobePrinted")
					minNumberofReadsTobePrinted=stoi(parVal);
				else if(parName== "minOverlapDifference4ClipBranches")
					minOvlDiffToClip=stoi(parVal);
				else if(parName== "minFoldToBeShortBranch")
					minFoldToBeShortBranch=stoi(parVal);
				else if(parName== "MinOverlap4Clip")
					minOvlToClip=stoi(parVal);
				else if(parName== "minUniquePEsupport")
					minUinqSupport=stoi(parVal);
				else if(parName== "minNonUniquePEsupport")
					minNonUniqSupport=stoi(parVal);
				else if(parName== "MinOverlap4SimplifyGraph")
					minOvl=stoi(parVal);
				else if(parName== "minSizeToBeShortBranch")
					minSizeToBeShortBranch=stoi(parVal);
				else if(parName== "MinOverlap4BuildGraph")
					continue;
				else if(parName == "PrintContigs") {
					if(parVal=="true") printContigs = true;
				}
				else if(parName == "PrintUnused") {
					if(parVal=="true") printUnused = true;
				}
				else if(parName == "PrintGFA") {
					if(parVal=="true") printGFA = true;
				}
				else if(parName == "PrintGFA2") {
					if(parVal=="true") printGFA2 = true;
				}
				else if(parName== "PrintScaffolds") {
					if(parVal=="false") printScaffolds = false;
				}
				else if(parName== "maxReadsUsed") {
					maxReadsUsed=stod(parVal);
				}
				else
					MYEXIT("Unknown parameter in parameter file: "+parName);
			}
		}
		filePointer.close();
	}

}

bool Config::setConfig(int argc, char **argv)
{
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
		Config::printHelp();
		return false;
	}

	for(UINT64 i = 1; i <= argumentsList.size()-1; i++)
	{
		// -h/--help: only print out the help contents
		if (argumentsList[i] == "-h" || argumentsList[i] == "--help")
		{
			Config::printHelp();
			exit(0);
		}
		// -fs: contained single read reduction read filename(s) (comma separated fasta/fastq)
		else if (argumentsList[i] == "-fs") {
			string inputFilenames=argumentsList[++i];
			stringstream ss(inputFilenames);
			string item;

			while (getline(ss, item, ','))
			{
				Config::readSingleFilenamesList.push_back(item);
			}
		}
		// -fpi: contained interleaved paired-end read reduction read filename(s) (comma separated fasta/fastq)
		else if (argumentsList[i] == "-fpi") {
			string inputFilenames=argumentsList[++i];
			stringstream ss(inputFilenames);
			string item;
			while (getline(ss, item, ','))
			{
				Config::readInterPairedFilenamesList.push_back(item);
			}
		}
		// -fp: contained paired-end read reduction read filename(s) in pairs of two (comma separated fasta/fastq)
		else if (argumentsList[i] == "-fp") {
			string inputFilenames=argumentsList[++i];
			stringstream ss(inputFilenames);
			string item;

			while (getline(ss, item, ','))
			{
				Config::readPairedFilenamesList.push_back(item);
			}
			if(readPairedFilenamesList.size()%2)
			{
				cout<<"Paired-end files should be in pairs of two (e.g. A1.fq,A2.fq,B1.fq,B2.fq)."<<endl;
				exit(0);
			}
		}
		// -e: overlapped edge property graph filename(s) (comma separated file list)
		else if (argumentsList[i] == "-e") {
			string inputFilenames=argumentsList[++i];
			stringstream ss(inputFilenames);
			string item;

			while (getline(ss, item, ','))
			{
				Config::edgeFilenamesList.push_back(item);
			}
		}
		// -log: verbosity level of log messages
		else if (argumentsList[i] == "-log") {
			FILELog::ReportingLevel() = FILELog::FromString(argumentsList[++i]);
		}
		// -crd: contained read filename(s) (comma separated file list)
		else if (argumentsList[i] == "-crd")
		{
			string inputFilenames=argumentsList[++i];
			stringstream ss(inputFilenames);
			string item;

			while (getline(ss, item, ','))
			{
				Config::containedReadsFile.push_back(item);
			}
		}
		else if (argumentsList[i] == "-p")
		{
			Config::paramFileName = argumentsList[++i];
		}
		else if (argumentsList[i] == "-p2")
		{
			Config::paramFileName2 = argumentsList[++i];
		}
		else if (argumentsList[i] == "-p3")
		{
			Config::paramFileName3 = argumentsList[++i];
		}
		else if (argumentsList[i] == "-simPth")
		{
			Config::simplifyGraphPath = argumentsList[++i];
		}
		// -o: all output filename prefix (default: disco)
		else if (argumentsList[i] == "-o") {
			Config::outFilenamePrefix = argumentsList[++i];
		}
		// -t: max omp thread allowed
		else if (argumentsList[i] == "-t") {
			Config::parallelThreadPoolSize = stoi(argumentsList[++i]);
		}
		else
		{
			Config::printHelp();
			return false;
		}
	}
	return true;
}

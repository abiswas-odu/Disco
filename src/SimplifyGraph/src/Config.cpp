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
string Config::outFilenamePrefix = "new_omega_out";
string Config::containedReadsFile = "";
string Config::simplifyGraphPath ="";
UINT64 Config::parallelThreadPoolSize(4);
string Config::paramFileName="parameter.cfg";

// global variables with default value
unsigned int minOvl=40;

unsigned int minReadsCountInEdgeToBeNotDeadEnd = 5;
unsigned int minEdgeLengthToBeNotDeadEnd = 500;

unsigned int minReadsCountToHave0Flow = 10;
unsigned int minEdgeLengthToHave0Flow = 1000;

unsigned int minReadsCountInEdgeToBe1MinFlow = 10;
unsigned int minEdgeLengthToBe1MinFlow = 1000;

unsigned int minOvlDiffToClip = 25;
unsigned int minFoldToBeShortBranch = 5;

unsigned int minUinqSupport=3;
unsigned int minNonUniqSupport=0;

unsigned int minContigLengthTobeReported = 500;

bool printContigs=false;

bool printScaffolds=true;
//=============================================================================
// print help usage
//=============================================================================
void Config::printHelp()
{
	cout << endl
		<< "  Usage:" << endl
		<< "    omega2 [OPTION]...<PARAM>..." << endl
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

void Config::setParameters()
{
	ifstream filePointer;
	filePointer.open(paramFileName.c_str());
	if(!filePointer.is_open()) {
		MYEXIT("Unable to open parameter file: "+paramFileName);
	}
	else {
		string par_text="";
		while(getline(filePointer,par_text)) {

			if(par_text.find("=") != std::string::npos)
			{
				vector<string> tok = Utils::split(par_text,'=');
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
				else if(parName== "minOverlapDifference4ClipBranches")
					minOvlDiffToClip=stoi(parVal);
				else if(parName== "minFoldToBeShortBranch")
					minFoldToBeShortBranch=stoi(parVal);
				else if(parName== "minUniquePEsupport")
					minUinqSupport=stoi(parVal);
				else if(parName== "minNonUniquePEsupport")
					minNonUniqSupport=stoi(parVal);
				else if(parName== "MinOverlap4SimplifyGraph")
					minOvl=stoi(parVal);
				else if(parName== "MinOverlap4BuildGraph")
					continue;
				else if(parName == "PrintContigs") {
					if(parVal=="true") printContigs = true;
				}
				else if(parName== "PrintScaffolds") {
					if(parVal=="false") printScaffolds = false;
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
		// -e: ovelapped edge property graph filename(s) (comma separated edge list)
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
		else if (argumentsList[i] == "-crd")
		{
			Config::containedReadsFile = argumentsList[++i];
		}
		else if (argumentsList[i] == "-p")
		{
			Config::paramFileName = argumentsList[++i];
		}
		else if (argumentsList[i] == "-simPth")
		{
			Config::simplifyGraphPath = argumentsList[++i];
		}
		// -o: all output filename prefix (default: new_omega_out)
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

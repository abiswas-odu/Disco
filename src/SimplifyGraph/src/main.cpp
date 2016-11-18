//============================================================================
// Name        : main.cpp
// Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey, Abhishek Biswas
// Version     : v3.0
// Copyright   : 2017 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : Main code
//============================================================================

#include "DataSet.h"
#include "OverlapGraph.h"
#include "logcpp/log.h"
#include "Utils.h"

#define FINAL_ITER 4

int OverlapGraph::s_nReads_in_goodEdges = 0;
int OverlapGraph::s_nGoodEdges = 0;
TLogLevel loglevel = logINFO;                   /* verbosity level of logging */
string outputFilenamePrefix = "omega3";

bool SimplifyGraph(const vector<std::string> &read_SingleFiles,const vector<std::string> &read_PairFiles,
		vector<std::string> &read_PairInterFiles, vector<vector<int>> &checkPointParams, const vector<std::string> &edgeFilenameList,
		string simplifyPartialPath, DataSet *dataSet,
		UINT64 minOvl, UINT64 &ctgCount, UINT64 &scfCount, UINT64 parallelThreadPoolSize,UINT64 containedCtr, int interationCount);

void SetParameters(int interationCount);

UINT64 readCheckpointInfo(vector< vector<int> > &checkPointParams, UINT64 &ctgCount,UINT64 &scfCount);

int main(int argc, char **argv) {

	// Parse and print command line options
	if(!Config::setConfig(argc, argv)){
		cerr << "Error: wrong configurations" << endl;
		return false;
	}
	vector<string> readSingleFilenameList 	= Config::getSingleReadFilenames();
	vector<string> readPairedFilenameList 	= Config::getPairedReadFilenames();
	vector<string> readInterPairedFilenameList 	= Config::getInterPairedReadFilenames();
	vector<string> edgeFilenameList 	= Config::getEdgeFilenames();
	outputFilenamePrefix 			= Config::getOutputFilenamePrefix();
	vector<string> containedReadsFileName = Config::getContainedReadsFile();
	string simplifyPartialPath = Config::getSimplifyGraphPath();
	loglevel 				= FILELog::ReportingLevel();
	UINT64 threadPoolSize = Config::getThreadPoolSize();
	UINT64 ctgCount=0, scfCount=0;
	FILE_LOG(logINFO) << "File(s) including reads: ";
	if(loglevel > 1){
		for(vector<std::string>::iterator it = readSingleFilenameList.begin(); it!=readSingleFilenameList.end(); ++it)
			FILE_LOG(logINFO) << *it << "\t";
		for(vector<std::string>::iterator it = readPairedFilenameList.begin(); it!=readPairedFilenameList.end(); ++it)
			FILE_LOG(logINFO) << *it << "\t";
		for(vector<std::string>::iterator it = readInterPairedFilenameList.begin(); it!=readInterPairedFilenameList.end(); ++it)
			FILE_LOG(logINFO) << *it << "\t";
	}
	FILE_LOG(logINFO) << endl;
	FILE_LOG(logINFO) << "File(s) including edges: ";
	if(loglevel > 1){
		for(vector<std::string>::iterator it = edgeFilenameList.begin(); it!=edgeFilenameList.end(); ++it)
			FILE_LOG(logINFO) << *it << "\t";
	}
	FILE_LOG(logINFO) << endl;
	FILE_LOG(logINFO) << "Output file names' prefix is: " << outputFilenamePrefix << endl;

	CLOCKSTART;
	//Load reads before simplification...
	DataSet *dataSet = new DataSet(readSingleFilenameList,readPairedFilenameList, readInterPairedFilenameList); // construct dataset from reads file(s)
	//Load and get count of contained reads
	UINT64 containedCtr = dataSet->storeContainedReadInformation(containedReadsFileName);

	FILE_LOG(logINFO) << "Total number of unique reads loaded from read file(s): "
		<< dataSet->size() << "\n";
	//Read checkpoint file and set checkpoint parameters
	vector< vector<int> > checkPointParams;
	int startItr=readCheckpointInfo(checkPointParams,ctgCount, scfCount);

	//Do simplification iterations
	bool continueSimplification= true;
	for(int i=startItr;i < FINAL_ITER && continueSimplification; i++)
	{
		//Read parameter file and set assembly parameters
		SetParameters(i);
		continueSimplification=SimplifyGraph(readSingleFilenameList,readPairedFilenameList, readInterPairedFilenameList,
				checkPointParams,edgeFilenameList, simplifyPartialPath, dataSet,
					minOvl, ctgCount, scfCount, threadPoolSize,containedCtr, i);

		//Clear edge information stored in the reads before the second iteration
		#pragma omp parallel for schedule(guided) num_threads(threadPoolSize)
		for(UINT64 i = 1; i <= dataSet->size() ; i++) // For each read.
		{
			dataSet->at(i)->ClearEdgeInfo();
			dataSet->at(i)->setUsedRead(false);
		}
	}
	//Print unused reads
	if(printUnused)
	{
		dataSet->writeUnUsedReads(outputFilenamePrefix);
	}
	delete dataSet;
	CLOCKSTOP;
	return 0;
}

bool SimplifyGraph(const vector<std::string> &read_SingleFiles,const vector<std::string> &read_PairFiles,
		vector<std::string> &read_PairInterFiles, vector<vector<int>> &checkPointParams,const vector<std::string> &edgeFilenameList,
		string simplifyPartialPath, DataSet *dataSet, UINT64 minOvl, UINT64 &ctgCount,UINT64 &scfCount,
		UINT64 threadPoolSize, UINT64 containedCtr, int interationCount)
{
	CLOCKSTART;
	FILE_LOG(logINFO) <<"Graph Simplification Iteration: "<<interationCount<<endl;
	UINT64 usedReads = 0;
	Utils::writeCheckPointFile(outputFilenamePrefix,"Iteration="+SSTR(interationCount));
	UINT64 nonContainedReads = dataSet->size()-containedCtr;
	for(int i=1;i < interationCount; i++)	//Load used reads
	{
		string usedReadFileName = outputFilenamePrefix+"_UsedReads_"+SSTR(i)+".txt";
		usedReads += dataSet->LoadUsedReads(usedReadFileName);
	}
	if(usedReads>(maxReadsUsed*nonContainedReads))
	{
		FILE_LOG(logINFO) <<"Graph simplification iteration terminated. Most reads used already. Assembly simplification complete."<<endl;
		return false;
	}
	OverlapGraph *overlapGraph;
	//Check if the simple edges have been processed
	if(checkPointParams[interationCount-1][0]==0)	//Simple edges need to be processed into composite graph
	{
		overlapGraph = new OverlapGraph(edgeFilenameList, simplifyPartialPath, dataSet,
				minOvl, threadPoolSize);
		//Write checkpoint graph
		overlapGraph->printAllEdges(outputFilenamePrefix+"_CurrGraph_.txt");
		Utils::writeCheckPointFile(outputFilenamePrefix,"ParSimplify=1");
	}
	else											//Simple edges already processed into composite graph; load global graph
	{
		string parGlobalGraph = outputFilenamePrefix+"_CurrGraph_.txt";
		overlapGraph = new OverlapGraph(parGlobalGraph, simplifyPartialPath, dataSet,
						minOvl, threadPoolSize);
	}
	//Initial Simplification
	if(checkPointParams[interationCount-1][1]==0)		//Check if initial simplification complete
	{
		overlapGraph->graphPathFindInitial();
		//Write checkpoint graph
		overlapGraph->printAllEdges(outputFilenamePrefix+"_CurrGraph_.txt");
		Utils::writeCheckPointFile(outputFilenamePrefix,"InitialSimplify=1");
	}
	else
		FILE_LOG(logINFO) <<"Skipping initial simplification..."<<endl;
	//ClipBranches and remove similar edges
	if(checkPointParams[interationCount-1][2]==0)		//Check if aggressive simplification complete
	{
		overlapGraph->simplifyGraph();
		//Write checkpoint graph
		overlapGraph->printAllEdges(outputFilenamePrefix+"_CurrGraph_.txt");
		Utils::writeCheckPointFile(outputFilenamePrefix,"AggressiveSimplify=1");
	}
	else
		FILE_LOG(logINFO) <<"Skipping aggressive simplification..."<<endl;
	// Flow analysis
	if(checkPointParams[interationCount-1][3]==0)		//Check if flow analysis complete
	{
		overlapGraph->calculateFlowStream();
		overlapGraph->removeAllEdgesWithoutFlow();
		//Write checkpoint graph
		overlapGraph->printAllEdges(outputFilenamePrefix+"_CurrGraph_.txt");
		Utils::writeCheckPointFile(outputFilenamePrefix,"FlowAnalysis=1");
	}
	else
		FILE_LOG(logINFO) <<"Skipping flow analysis simplification..."<<endl;
	if(checkPointParams[interationCount-1][4]==0)		//Check if post-flow analysis simplification complete
	{
		overlapGraph->simplifyGraph();
		//Write checkpoint graph
		overlapGraph->printAllEdges(outputFilenamePrefix+"_CurrGraph_.txt");
		Utils::writeCheckPointFile(outputFilenamePrefix,"PostFlowAnalysis=1");
	}
	else
		FILE_LOG(logINFO) <<"Skipping post flow analysis simplification..."<<endl;

	//Print contig files before scaffolding
	if(printContigs && checkPointParams[interationCount-1][5]==0)
	{
		string edge_file = outputFilenamePrefix+"_contigEdgesFinal_"+SSTR(interationCount)+".txt";
		string edge_cov_file = outputFilenamePrefix+"_contigEdgeCoverageFinal_"+SSTR(interationCount)+".txt";
		string contig_file = outputFilenamePrefix+"_contigsFinal_"+SSTR(interationCount)+".fasta";
		string usedReadFileName = outputFilenamePrefix+"_UsedReads_"+SSTR(interationCount)+".txt";
		//overlapGraph->printContigs(contig_file, edge_file, edge_cov_file,usedReadFileName,"contig",ctgCount);
		overlapGraph->streamContigs(read_SingleFiles,read_PairFiles, read_PairInterFiles,
				contig_file, edge_file, edge_cov_file,usedReadFileName,"contig",ctgCount);
		//Write checkpoint graph
		overlapGraph->printAllEdges(outputFilenamePrefix+"_CurrGraph_.txt");
		Utils::writeCheckPointFile(outputFilenamePrefix,"PrintCtg="+SSTR(ctgCount));
	}
	else
		FILE_LOG(logINFO) <<"Skipping contig printing..."<<endl;
	//Print GFA file
	if(printGFA)
	{
		string gfa_file = outputFilenamePrefix+"_Graph_"+SSTR(interationCount)+".gfa";
		ofstream gfaFilePointer;
		gfaFilePointer.open(gfa_file.c_str());
		if(!gfaFilePointer)
			MYEXIT("Unable to open file: "+gfa_file);
		overlapGraph->generateGFAOutput(gfaFilePointer);
		gfaFilePointer.close();
	}
	//Print GFA2 file
	if(printGFA2)
	{
		string gfa_file = outputFilenamePrefix+"_Graph_"+SSTR(interationCount)+".gfa2";
		ofstream gfaFilePointer;
		gfaFilePointer.open(gfa_file.c_str());
		if(!gfaFilePointer)
			MYEXIT("Unable to open file: "+gfa_file);
		overlapGraph->generateGFA2Output(gfaFilePointer);
		gfaFilePointer.close();
	}
	//Start scaffolding
	if(checkPointParams[interationCount-1][6]==0)		//Check if scaffolding is complete
	{
		overlapGraph->calculateMeanAndSdOfInnerDistance();
		UINT64 iteration=0, counter = 0;
		do
		{
			// Mate pair paths are used to simplify the graph in this step
			FILE_LOG(logINFO) << endl;
			FILE_LOG(logINFO) << "===============================================================================================================================================" <<endl;
			FILE_LOG(logINFO) << "FIRST LOOP ITERATION " << ++iteration << endl;
			FILE_LOG(logINFO) << "===============================================================================================================================================" <<endl;
			counter = overlapGraph->findSupportByMatepairsAndMerge();
			overlapGraph->simplifyScaffoldGraph();
		} while (counter > 0 && iteration < loopLimit); // To avoid infinite loops

		iteration = 0;
		do
		{
			// Scaffolder
			FILE_LOG(logINFO) << endl;
			FILE_LOG(logINFO) << "===============================================================================================================================================" <<endl;
			FILE_LOG(logINFO) << "SECOND LOOP ITERATION " << ++iteration << endl;
			FILE_LOG(logINFO) << "===============================================================================================================================================" <<endl;
			counter = overlapGraph->scaffolder();
			overlapGraph->simplifyScaffoldGraph();

		} while (counter > 0 && iteration < loopLimit);// To avoid infinite loops

		if(printScaffolds)
		{
			string edge_file = outputFilenamePrefix+"_scaffoldEdgesFinal_"+SSTR(interationCount)+".txt";
			string contig_file = outputFilenamePrefix+"_scaffoldsFinal_"+SSTR(interationCount)+".fasta";
			string edge_cov_file = outputFilenamePrefix+"_scaffoldEdgeCoverageFinal_"+SSTR(interationCount)+".txt";
			string usedReadFileName = outputFilenamePrefix+"_UsedReads_"+SSTR(interationCount)+".txt";
			overlapGraph->streamContigs(read_SingleFiles,read_PairFiles, read_PairInterFiles,
					contig_file, edge_file, edge_cov_file,usedReadFileName,"scaff",scfCount);
			//overlapGraph->printContigs(contig_file, edge_file, edge_cov_file,usedReadFileName,"scaff",scfCount);
		}
		//Write checkpoint graph
		overlapGraph->printAllEdges(outputFilenamePrefix+"_CurrGraph_.txt");
		Utils::writeCheckPointFile(outputFilenamePrefix,"Scaffold="+SSTR(scfCount));
	}
	else
		FILE_LOG(logINFO) <<"Skipping scaffolding printing..."<<endl;
	//Print the total used read count
	usedReads = 0;
	#pragma omp parallel for schedule(guided) reduction(+:usedReads) num_threads(threadPoolSize)
	for(UINT64 i = 1; i <= dataSet->size() ; i++) // For each read.
	{
		if(dataSet->at(i)->isUsedRead())
			usedReads++;
	}
	FILE_LOG(logINFO) <<"Iteration:"<<interationCount<<" Graph simplification has used a total of "<<usedReads<<" reads."<<endl;
	delete overlapGraph;
	if(usedReads>(maxReadsUsed*nonContainedReads))
	{
		FILE_LOG(logINFO) <<"Graph simplification iteration terminated. Most reads used already. Assembly simplification complete."<<endl;
		return false;
	}
	return true;
	CLOCKSTOP;
}

void SetParameters(int interationCount)
{
	//Load parameters based on iteration count...
	if(interationCount==1)
		Config::setParameters(Config::getParamFileName());
	else if(interationCount==2)
		Config::setParameters(Config::getParamFileName2());
	else
		Config::setParameters(Config::getParamFileName3());
	FILE_LOG(logINFO) << "Log level is " << FILELog::ReportingLevel() << ":\t" << Config::getLogLevel() << endl;
	FILE_LOG(logINFO) << "Minimum overlap length is: " << minOvl << endl;
	FILE_LOG(logINFO) << "Maximum read count in dead-end edge is: " << minReadsCountInEdgeToBeNotDeadEnd << endl;
	FILE_LOG(logINFO) << "Maximum edge length in dead-end edge is: " << minEdgeLengthToBeNotDeadEnd << endl;
	FILE_LOG(logINFO) << "Minimum read count in edges with flow is: " << minReadsCountInEdgeToBe1MinFlow << endl;
	FILE_LOG(logINFO) << "Minimum edge length of edges with flow is: " << minEdgeLengthToBe1MinFlow << endl;
	FILE_LOG(logINFO) << "Minimum edge length for edges to be reported is: " << minContigLengthTobeReported << endl;
	FILE_LOG(logINFO) << "Minimum overlap length difference for branches to clip: " << minOvlDiffToClip << endl;
	FILE_LOG(logINFO) << "Minimum fold difference to consider branches to be short: " << minFoldToBeShortBranch << endl;
}

UINT64 readCheckpointInfo(vector< vector<int> > &checkPointParams, UINT64 &ctgCount,UINT64 &scfCount)
{
	//Initialize checkpoint vector
	for(int i=1;i < FINAL_ITER; i++)
	{
		vector<int> paramValues = {0,0,0,0,0,0,0};
		checkPointParams.push_back(paramValues);
	}
	ifstream filePointer;
	string fileName = outputFilenamePrefix+"_SimplificationCheckpointInfo.txt";
	filePointer.open(fileName.c_str());
	int iterCtr=1;
	if(!filePointer.is_open()) {
		return 1;
	}
	else {
		string par_text="";
		while(getline(filePointer,par_text)) {
			if(par_text.find("=") != std::string::npos)
			{
				vector<string> tok = Utils::split(par_text,'=');
				string parName = Utils::trimmed(tok[0]);
				string parVal = Utils::trimmed(tok[1]);
				if(parName=="Iteration")
				{
					iterCtr=atoi(parVal.c_str());
				}
				if(parName=="ParSimplify" && parVal=="1")
					checkPointParams[iterCtr-1][0]=1;
				if(parName=="InitialSimplify" && parVal=="1")
					checkPointParams[iterCtr-1][1]=1;
				if(parName=="AggressiveSimplify" && parVal=="1")
					checkPointParams[iterCtr-1][2]=1;
				if(parName=="FlowAnalysis" && parVal=="1")
					checkPointParams[iterCtr-1][3]=1;
				if(parName=="PostFlowAnalysis" && parVal=="1")
					checkPointParams[iterCtr-1][4]=1;
				if(parName=="PrintCtg")
				{
					checkPointParams[iterCtr-1][5]=1;
					ctgCount=atoi(parVal.c_str());
				}
				if(parName=="Scaffold")
				{
					checkPointParams[iterCtr-1][6]=1;
					scfCount=atoi(parVal.c_str());
				}
			}
		}
	}
	UINT64 startIter=1;
	for(;startIter < FINAL_ITER; startIter++)
	{
		for(int j=0;j < 7; j++)
			if(checkPointParams[startIter-1][j]==0)
				return startIter;
	}
	return 1;
}



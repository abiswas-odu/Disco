//============================================================================
// Name        : main.cpp
// Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey, Abhishek Biswas
// Version     : v3.0
// Copyright   : 2017 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : Main code
//============================================================================

#include "Config.h"
#include "DataSet.h"
#include "OverlapGraph.h"
#include "logcpp/log.h"

#define FINAL_ITER 4
#define MAX_USED 0.6

int OverlapGraph::s_nReads_in_goodEdges = 0;
int OverlapGraph::s_nGoodEdges = 0;
TLogLevel loglevel = logINFO;                   /* verbosity level of logging */
string outputFilenamePrefix = "omega3";

void SimplifyGraph(const vector<std::string> &edgeFilenameList,
		string simplifyPartialPath, DataSet *dataSet,
		UINT64 minOvl,UINT64 &ctgCount, UINT64 parallelThreadPoolSize,UINT64 containedCtr, int interationCount);

void SetParameters(int interationCount);

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
	UINT64 ctgCount=0;
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
	//Read parameter file and set assembly parameters
	SetParameters(1);
	SimplifyGraph(edgeFilenameList, simplifyPartialPath, dataSet,
			minOvl, ctgCount, threadPoolSize,containedCtr, 1);

	//Do more iterations
	for(int i=2;i < FINAL_ITER; i++)
	{
		//Clear edge information stored in the reads before the second iteration
		#pragma omp parallel for schedule(guided) num_threads(threadPoolSize)
		for(UINT64 i = 1; i <= dataSet->size() ; i++) // For each read.
		{
			dataSet->at(i)->ClearEdgeInfo();
			dataSet->at(i)->setUsedRead(false);
		}

		//Read parameter file and set assembly parameters
		SetParameters(i);
		SimplifyGraph(edgeFilenameList, simplifyPartialPath, dataSet,
					minOvl, ctgCount, threadPoolSize,containedCtr, i);
	}
	//Print unused reads
	if(printUnused)
	{
		string unusedReads = outputFilenamePrefix+"_UnusedReads.fasta";
		ofstream unUsedReadsFilePointer;
		unUsedReadsFilePointer.open(unusedReads.c_str());
		if(!unUsedReadsFilePointer)
			MYEXIT("Unable to open file: "+unusedReads);
		for(UINT64 i = 1; i <= dataSet->size() ; i++) // For each read.
		{
			UINT64 mateID = dataSet->getMatePair(i);
			if(mateID!=0)
			{
				if((!dataSet->at(i)->isUsedRead() || !dataSet->at(mateID)->isUsedRead()) && i<mateID)
				{
					unUsedReadsFilePointer<<">"<<i<<".1"<<endl<<dataSet->at(i)->getStringForward()<<endl;
					unUsedReadsFilePointer<<">"<<mateID<<".2"<<endl<<dataSet->at(mateID)->getStringForward()<<endl;
				}
			}
			else
			{
				if(!dataSet->at(i)->isUsedRead())
					unUsedReadsFilePointer<<">"<<i<<".1"<<endl<<dataSet->at(i)->getStringForward()<<endl;
			}
		}
		unUsedReadsFilePointer.close();
	}
	delete dataSet;
	CLOCKSTOP;
	return 0;
}

void SimplifyGraph(const vector<std::string> &edgeFilenameList,
		string simplifyPartialPath, DataSet *dataSet, UINT64 minOvl, UINT64 &ctgCount,
		UINT64 threadPoolSize, UINT64 containedCtr, int interationCount)
{
	CLOCKSTART;
	FILE_LOG(logINFO) <<"Graph Simplification Iteration: "<<interationCount<<endl;

	//Load already used reads
	for(int i=1;i < interationCount; i++)
	{
		string usedReadFileName = outputFilenamePrefix+"_UsedReads_"+SSTR(i)+".txt";
		UINT64 usedReads = dataSet->LoadUsedReads(usedReadFileName);
		UINT64 nonContainedReads = dataSet->size()-containedCtr;
		if(usedReads>(MAX_USED*nonContainedReads))
		{
			FILE_LOG(logINFO) <<"Graph simplification iteration terminated. Most reads used already. Assembly simplification complete."<<endl;
			return;
		}
	}

	OverlapGraph *overlapGraph = new OverlapGraph(edgeFilenameList, simplifyPartialPath, dataSet,
				minOvl, threadPoolSize);

	//Initial Simplification
	overlapGraph->graphPathFindInitial();
	/*std::string graph_file = outputFilenamePrefix+"_graph0.cytoscape";
	ofstream g_out(graph_file.c_str());
	g_out << *overlapGraph;
	g_out.close();*/
	//ClipBranches and remove similar edges
	overlapGraph->simplifyGraph();
	// Flow analysis
	overlapGraph->calculateFlowStream();
	overlapGraph->removeAllEdgesWithoutFlow();
	overlapGraph->simplifyGraph();

	//Print contig files before scaffolding
	if(printContigs)
	{
		string edge_file = outputFilenamePrefix+"_contigEdgesFinal_"+SSTR(interationCount)+".txt";
		string edge_cov_file = outputFilenamePrefix+"_contigEdgeCoverageFinal_"+SSTR(interationCount)+".txt";
		string contig_file = outputFilenamePrefix+"_contigsFinal_"+SSTR(interationCount)+".fasta";
		string usedReadFileName = outputFilenamePrefix+"_UsedReads_"+SSTR(interationCount)+".txt";
		overlapGraph->printContigs(contig_file, edge_file, edge_cov_file,usedReadFileName,"contig",ctgCount);
	}

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
		overlapGraph->printContigs(contig_file, edge_file, edge_cov_file,usedReadFileName,"scaff",ctgCount);
	}

	//Print the total used read count
	UINT64 usedReads = 0;
	#pragma omp parallel for schedule(guided) reduction(+:usedReads) num_threads(threadPoolSize)
	for(UINT64 i = 1; i <= dataSet->size() ; i++) // For each read.
	{
		if(dataSet->at(i)->isUsedRead())
			usedReads++;
	}
	FILE_LOG(logINFO) <<"Iteration:"<<interationCount<<" Graph simplification has used a total of "<<usedReads<<" reads."<<endl;

	delete overlapGraph;
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

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

int OverlapGraph::s_nReads_in_goodEdges = 0;
int OverlapGraph::s_nGoodEdges = 0;
TLogLevel loglevel = logINFO;                   /* verbosity level of logging */
string outputFilenamePrefix = "omega2";

int main(int argc, char **argv) {

	// Parse command line options:0
	if(!Config::setConfig(argc, argv)){
		cerr << "Error: wrong configurations" << endl;
		return false;
	}
	CLOCKSTART;
	//Read parameter file and set assembly parameters
	Config::setParameters();

	vector<string> readSingleFilenameList 	= Config::getSingleReadFilenames();
	vector<string> readPairedFilenameList 	= Config::getPairedReadFilenames();
	vector<string> readInterPairedFilenameList 	= Config::getInterPairedReadFilenames();
	vector<string> edgeFilenameList 	= Config::getEdgeFilenames();
	outputFilenamePrefix 			= Config::getOutputFilenamePrefix();
	string containedReadsFileName = Config::containedReadsFile;
	string simplifyPartialPath = Config::getSimplifyGraphPath();

	loglevel 				= FILELog::ReportingLevel();
	UINT64 threadPoolSize = Config::getThreadPoolSize();

	FILE_LOG(logINFO) << "Log level is " << FILELog::ReportingLevel() << ":\t" << Config::getLogLevel() << endl;
	FILE_LOG(logINFO) << "Minimum overlap length is: " << minOvl << endl;
	FILE_LOG(logINFO) << "File(s) including reads: ";
	if(loglevel > 1){
		for(vector<std::string>::iterator it = readSingleFilenameList.begin(); it!=readSingleFilenameList.end(); ++it)
			FILE_LOG(logINFO) << *it << "\t";
		for(vector<std::string>::iterator it = readPairedFilenameList.begin(); it!=readPairedFilenameList.end(); ++it)
					FILE_LOG(logINFO) << *it << "\t";
	}
	FILE_LOG(logINFO) << endl;
	FILE_LOG(logINFO) << "File(s) including edges: ";
	if(loglevel > 1){
		for(vector<std::string>::iterator it = edgeFilenameList.begin(); it!=edgeFilenameList.end(); ++it)
			FILE_LOG(logINFO) << *it;
	}
	FILE_LOG(logINFO) << endl;
	FILE_LOG(logINFO) << "Output file names' prefix is: " << outputFilenamePrefix << endl;
	FILE_LOG(logINFO) << "Maximum read count in dead-end edge is: " << minReadsCountInEdgeToBeNotDeadEnd << endl;
	FILE_LOG(logINFO) << "Maximum edge length in dead-end edge is: " << minEdgeLengthToBeNotDeadEnd << endl;
	FILE_LOG(logINFO) << "Minimum read count in edges with flow is: " << minReadsCountInEdgeToBe1MinFlow << endl;
	FILE_LOG(logINFO) << "Minimum edge length of edges with flow is: " << minEdgeLengthToBe1MinFlow << endl;
	FILE_LOG(logINFO) << "Minimum edge length for edges to be reported is: " << minContigLengthTobeReported << endl;
	FILE_LOG(logINFO) << "Minimum overlap length difference for branches to clip: " << minOvlDiffToClip << endl;
	FILE_LOG(logINFO) << "Minimum fold difference to consider branches to be short: " << minFoldToBeShortBranch << endl;

	OverlapGraph *overlapGraph = new OverlapGraph(edgeFilenameList, readSingleFilenameList,
			readPairedFilenameList, readInterPairedFilenameList, simplifyPartialPath,
			minOvl, threadPoolSize);

	overlapGraph->graphPathFindInitial(containedReadsFileName);

	string edge_file = outputFilenamePrefix+"_contigEdges_2.txt";
	string edge_cov_file = outputFilenamePrefix+"_contigEdgeCoverage_2.txt";
	string contig_file = outputFilenamePrefix+"_contigsFinal_2.fasta";
	ofstream f_out;
	f_out.open(contig_file.c_str());
	overlapGraph->printContigs(f_out, edge_file, edge_cov_file,"contig", readSingleFilenameList, readPairedFilenameList);
	f_out.close();

	/*std::string graph_file = outputFilenamePrefix+"_graph0.cytoscape";
	ofstream g_out(graph_file.c_str());
	g_out << *overlapGraph;
	g_out.close();*/

	overlapGraph->simplifyGraph();

	edge_file = outputFilenamePrefix+"_contigEdges_3.txt";
	edge_cov_file = outputFilenamePrefix+"_contigEdgeCoverage_3.txt";
	contig_file = outputFilenamePrefix+"_contigsFinal_3.fasta";
	f_out.open(contig_file.c_str());
	overlapGraph->printContigs(f_out, edge_file, edge_cov_file,"contig", readSingleFilenameList, readPairedFilenameList);
	f_out.close();

	// Flow analysis
	overlapGraph->calculateFlowStream();
	overlapGraph->removeAllEdgesWithoutFlow();

	edge_file = outputFilenamePrefix+"_contigEdges_4.txt";
	edge_cov_file = outputFilenamePrefix+"_contigEdgeCoverage_4.txt";
	contig_file = outputFilenamePrefix+"_contigsFinal_4.fasta";
	f_out.open(contig_file.c_str());
	overlapGraph->printContigs(f_out, edge_file, edge_cov_file,"contig", readSingleFilenameList, readPairedFilenameList);
	f_out.close();


	overlapGraph->simplifyGraph();

	//Print contig files before scaffolding
	if(printContigs)
	{
		edge_file = outputFilenamePrefix+"_contigEdgesFinal.txt";
		edge_cov_file = outputFilenamePrefix+"_contigEdgeCoverageFinal.txt";
		contig_file = outputFilenamePrefix+"_contigsFinal.fasta";
		f_out.open(contig_file.c_str());
		overlapGraph->printContigs(f_out, edge_file, edge_cov_file,"contig", readSingleFilenameList, readPairedFilenameList);
		f_out.close();
	}

	overlapGraph->calculateMeanAndSdOfInnerDistance();

	UINT64 iteration=0, counter = 0;
	do
	{
		// Mate pair paths are used to simplify the graph in this step
		cout << endl;
		cout << "===============================================================================================================================================" <<endl;
		cout << "FIRST LOOP ITERATION " << ++iteration << endl;
		cout << "===============================================================================================================================================" <<endl;
		counter = overlapGraph->findSupportByMatepairsAndMerge();
		overlapGraph->simplifyScaffoldGraph();
	} while (counter > 0 && iteration < loopLimit); // To avoid infinite loops

	iteration = 0;
	do
	{
		// Scaffolder
		cout << endl;
		cout << "===============================================================================================================================================" <<endl;
		cout << "SECOND LOOP ITERATION " << ++iteration << endl;
		cout << "===============================================================================================================================================" <<endl;
		overlapGraph->simplifyScaffoldGraph();
		counter = overlapGraph->scaffolder();

	} while (counter > 0 && iteration < loopLimit);// To avoid infinite loops

	if(printScaffolds)
	{
		edge_file = outputFilenamePrefix+"_scaffoldEdgesFinal.txt";
		contig_file = outputFilenamePrefix+"_scaffoldsFinal.fasta";
		edge_cov_file = outputFilenamePrefix+"_scaffoldEdgeCoverageFinal.txt";
		f_out.open(contig_file.c_str());
		overlapGraph->printContigs(f_out, edge_file, edge_cov_file,"scaff", readSingleFilenameList, readPairedFilenameList);
		f_out.close();
	}
	delete overlapGraph;
	CLOCKSTOP;
	return 0;
}

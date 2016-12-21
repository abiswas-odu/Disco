//============================================================================
// Name        : main.cpp
// Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey, Abhishek Biswas
// Version     : v1.2
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : Main code
//============================================================================

#include "OverlapGraphSimple.h"

TLogLevel loglevel = logDEBUG4;                   /* verbosity level of logging */
string outputFilenamePrefix = "disco";

int main(int argc, char **argv) {

	CLOCKSTART;
	string edgeFilename = argv[1];
	string outputFilename = argv[2];
	UINT64 minOvl = stoi(argv[3]);
	UINT64 threadPoolSize = stoi(argv[4]);
	OverlapGraphSimple *overlapGraph = new OverlapGraphSimple(edgeFilename, outputFilename, minOvl, threadPoolSize);
	delete overlapGraph;
	CLOCKSTOP;
	return 0;
}

#ifndef DATASET_H
#define DATASET_H

//============================================================================
// Name        : DataSet.h
// Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey
// Version     : v1.2
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : DataSet header file
//============================================================================

#include "Config.h"
#include "Read.h"
#include "Utils.h"

struct DataSetInfo
{
	UINT64 datasetNumber;			// the ID of the dataset

	bool isPaired;					//0: Not paired-end
									//1: Single

	bool isReadInterleaved;				//0: Not interleaved paired-end reads
										//1: Interleaved paired-end reads

	UINT8 matePairOrientation; 		// Mate pair orientation.
									// 0 = means the forward of this read and the reverse of the second read are paired-end.
										// paired reads are in forward - reverse orientation on different DNA strands, example:
											// TAT------------
											// ------------CAG
									// 1 = means the forward of this read and the forward of the second read are matepairs.
	double avgInnerDistance;			// Inner distance between the reads
	double avgInnerDistanceSD;			// Standard deviation of inner distance between the reads
	UINT64 r1Start;					//Start ReadID of this dataset
	UINT64 r1End;					//Ending ReadID of this dataset
	UINT64 r2Start;					//Start ReadID of this dataset
	UINT64 r2End;					//Ending ReadID of this dataset

	string r1FileName;
	string r2FileName;

};
class DataSet
{
	private:
		/* ====================  DATA MEMBERS  ======================================= */
		vector<Read*> *m_vec_reads;    /* vector pointers to Reads */

		vector<DataSetInfo> *dataSetInfo;

		bool isSequenceLoaded;

		/* ====================  METHODS      ======================================= */
		// Load reads from a read file
		void loadReadsFromReadFile(const std::string &read_file);
		//
		// Load reads from an edge file, in this case, only the length is available
		void loadReadsFromEdgeFile(const std::string &edge_file);

		//Check if read is good or bad
		bool testRead(const string & read);
	public:
		/* ====================  LIFECYCLE     ======================================= */
		DataSet();

		// Load a list of files (either reads or edges)
		DataSet(const vector<std::string> &read_SingleFiles,const vector<std::string> &read_PairFiles,
				vector<std::string> &read_PairInterFiles,string usedReadFileName);

		// Copy constructor
		DataSet(const DataSet &s_dataset);

		~DataSet();

		/* ====================  OPERATORS     ======================================= */
		friend std::ostream& operator<< (std::ostream &out, DataSet & a_data_set);

		DataSet& operator= (const DataSet &s_dataset);

		/* ====================  MUTATORS      ======================================= */
		void addRead(Read *r){m_vec_reads->push_back(r);}

		void rmRead(Read *r);

		void setSequenceLoaded(bool isLoaded) {isSequenceLoaded=isLoaded;}

		// Load contained read information from file
		void storeContainedReadInformation(vector<string> containedReadFile);

		/* ====================  ACCESSORS     ======================================= */ 
		UINT64 size() const{ return m_vec_reads->size();}

		vector<DataSetInfo> * getDataSetInfo() {return dataSetInfo;}

		Read* at(UINT64 ID) const;

		UINT64 getDataSetNumber(UINT64 readID);

		bool getSequenceLoaded() { return isSequenceLoaded;}

		UINT32 getReadCoverage(UINT64 readID, UINT64 indx) const;

		vector<UINT64> getMatePairList(Read *read);

		UINT64 getMatePair(UINT64 r1ID);

		void writeUsedReads(string usedReadfileName);
};

bool compareEdgesByReads (const Edge *edge1, const Edge* edge2);

void printUnorderedMap(const unordered_map<UINT64, UINT64> & readIDMap, ostream & out);
#endif /* DATASET_H */

#ifndef OVERLAPGRAPH_H
#define OVERLAPGRAPH_H

/*
 * ===== CLASS HEADER ========================================================
 * Name        : OverlapGraph.cpp
 * Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey,  Abhishek Biswas
 * Version     : v1.2
 * Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
 * Description : OverlapGraph header file
 *============================================================================
 */

#include "Config.h"
#include "DataSet.h"
#include "Edge.h"

#define loopLimit 15			// Number of time to loop in the main function.
#define allowedUsedReads 3
#define EXPLORE_DEPTH 100
#define insertSizeRangeSD 3 	// 3 means mean +/- 3 SD
#define MAX_INNER_DIST_TRESH 100


extern char **environ;

// This structure is used to store list of pair of edges and their support. Used in two function: 1. when we find path by mate-pairs 2. scaffolder.
// When it's used in the findSupportByMatePairs function, these two edges are adjacent edges supported by all paths between several matepairs
// When it's used in the scaffolder function, these two edges are linked by multiple matepairs and they are not necesarrily adjacent edges, but they can be
struct pairedEdges
{
		/*
		 *  edge 1 is in front of edge 2 in the sequence
		 *  -----edge1----		-----edge2---
		 *  --------------		-------------
		 */
		Edge * edge1;
		Edge * edge2;
		UINT64 uniqSupport;			// number of matepairs uniquely supporting these two edges
		UINT64 nonUniqsupport;			// number of matepairs non-uniquely supporting these two edges
		INT64 distance;		// sum of the end of read1 to the end of edge 1 and the beginning of read2 to the beginning of edge2
		bool isFreed;
		bool operator < (const pairedEdges& rhs) const
		{
		       return uniqSupport > rhs.uniqSupport;
		}

};

struct unitigExt
{
	UINT64 seedSource;		//Source readID of the seed edge
	UINT64 seedDest;        //Destination readID of the seed edge
	UINT64 extDest;         //Destination readID of the extension sequence
	string extSeq;
};

typedef vector<Edge *> t_edge_vec;	// vector of pointers to Edge
class OverlapGraph
{
	private:
		/* ====================  DATA MEMBERS  ======================================= */
		DataSet 		*m_dataset;
		map<UINT64, t_edge_vec*> 	*m_graph;
		UINT64 			m_numberOfNodes;
		UINT64 			m_numberOfEdges;
		UINT64 			m_minOvl;
		bool 			m_flowComputed;
		UINT64 longestMeanOfInsertDistance;
		UINT64          p_ThreadPoolSize;

		/* ====================  METHODS       ======================================= */

		// Get mismatch information (type t_vpair) from mismatch string
		t_vpair* getMismatchFromString(const std::string &mismath_string);

		// Insert an edge in the overlap graph
		// Does not insert its twin edge
		void insertEdge( Edge * edge);   
		void insertFwdEdge( Edge * edge);   

		// Remove an edge from the edge list of the edge's source read, 
		// but do not delete this edge from memory yet.
		void removeEdgeFromSourceRead(Edge *edge);

		// Remove an edge from the overlap graph.
		void removeEdge(Edge *edge);

		void removeFwdEdge(Edge *edge);
		// Contract composite paths in the overlap graph.
		UINT64 contractCompositeEdgesPar(void);

		// Remove dead-ends from the overlap graph.
		UINT64 removeDeadEndNodes(void);

		// Remove dead-end like short branches, if they are much shorter 
		// than other longer edges even though they are longer than
		// the dead-end length threshold
		UINT64 removeShortBranches(void);

		// remove multi-edges with similar strings
		UINT64 removeSimilarEdges(void);

		// Clip branches with much shorter overlap length
		// comparing to other branches.
		UINT64 clipBranches(void);

		// Loops that can be traversed only one way
		UINT64 reduceLoops(void);

		// Calculate bounds and costs of flow for minimum cost flow in the overlap graph.
		void calculateBoundAndCost(const Edge *edge, INT64* FLOWLB, INT64* FLOWUB, INT64* COST);

		// Find an edge from source to destination in the overlap graph.
		Edge *findEdge(const UINT64 & source, const UINT64 & destination);

		// Find edges between two reads
		vector<Edge *> findEdges(const UINT64 & source, const UINT64 & destination);

		// Given lists of edge IDs (positive and negative), merge the corresponding
		// edges in each list, and insert them in the graph.
		// Also delete these edges at the end
		UINT64 merge_edges_with_IDlists(const vector<vector<INT64> > *lists_edgeIDs_to_merge,
				const map<INT64, Edge*> & id_to_edge);

		// locate a read given the p_read (composite numebr - last bit for orientation)
		Edge * locateRead(UINT64 p_read) const;

		// Merge a list of edges (instead of only two edges)
		UINT64 mergeListOfEdges(const vector< t_edge_vec > & list_edges_to_merge);

		//Merge 2 edges with flow, Remove the old edges that doesn't have any flow left, but keep the old edges that still have flow left.
		void merge2Edges(Edge *edge1, Edge *edge2);

		// Sort edges of each read based on ID of the destination read. 
		// This is only for ordering edges for convenience in the output file
		void sortEdgesByLength();

		void sortEdgesByDestID();

		/*Test if read has only 4 bases*/
		bool testRead(const string & read);

	public:
		/* ====================  DATA MEMBERS  ======================================= */
		// Keep track of number of reads contained in the edges with 1 unit of flow assigned 
		static int s_nReads_in_goodEdges; 

		// keep track of number of edges with 1 unit of flow assigned 
		static int s_nGoodEdges;        

		/* ====================  LIFECYCLE     ======================================= */
		OverlapGraph(void);

		OverlapGraph(const vector<std::string> &edge_files, string simplifyPartialPath, DataSet *dataSet,
				const UINT64 minOvl, const UINT64 parallelThreadPoolSize);

		OverlapGraph(const string parGlobalGraph, string simplifyPartialPath, DataSet *dataSet,
				UINT64 minOvl, UINT64 parallelThreadPoolSize);

		~OverlapGraph();

		/* ====================  ACCESSORS     ======================================= */
		UINT64 getNumberOfEdges(void) const {return m_numberOfEdges;}

		UINT64 getNumberOfNodes(void) const {return m_numberOfNodes;}

		bool isUsedEdge(UINT64 lFSize, UINT64 usedReadCtr,UINT64 unUsedMate, Read *source, Read *destination);

		/* ====================  MUTATORS      ======================================= */
		void setMinOvl(const UINT64 & minOvl = 50){m_minOvl = minOvl;}

		// Remove all simple edges without flow
		UINT64 removeAllEdgesWithoutFlow();

		void populate_edge(Edge *edge);

		/* ====================  OPERATORS     ======================================= */
		friend ostream& operator<< (ostream &out, const OverlapGraph & graph);


//		// Load read file again to print contigs
//		UINT64 streamReadFileSequences(const std::string &readFilename, UINT64 &readID, 
//				vector<vector<char> > & contigStrings, 
//				vector<unordered_map<UINT32, vector<char> > > & mismatchMap);


		// Some simple simplification.
		void simplifyGraph(void);

		void graphPathFindInitial();

		// Calculate the minimum cost flow of the overlap graph using file
		void calculateFlowStream(void);

		// Find all the edges in the graph, to be used in print graph and contigs
		void getEdges(t_edge_vec & contigEdges) const;
		
		// Print contigs to file, only the ones longer than the specified printing threshold
		void printContigs(string contig_file, string edge_file,string edge_cov_file,string usedReadFileName, string namePrefix, UINT64 &printed_contigs);

		// Print edges to a file, with an edgeName
		void printEdge(Edge *contigEdges, ostream & filePointer,ostream & fileUsedReadPointer , UINT64 edgeNameID) const;

		// Print all edges in the graph to a file
		void printAllEdges(string edge_file) const;

		// Print an edge
		void printEdge(Edge *contigEdge, ostream & filePointer) const;


		//Checks if an edge already exists
		bool existsEdge(Edge *checkEdge);

		//Functions for parallel graph simplification

		UINT64 contractCompositeEdges(void);
		void readParEdges(string edge_file);

		//Scaffolding functions...

		bool findPathBetweenMatepairs(const Read * read1, const Read * read2,
				UINT8 orient, UINT8 datasetNumber, vector <Edge *> &copyOfPath, vector <UINT64> &copyOfFlags);

		void exploreGraph(Edge* firstEdge, Edge * lastEdge, int distanceOnFirstEdge,
				int distanceOnLastEdge, UINT64 datasetNumber, UINT64 level, vector <Edge *> &firstPath, vector <UINT64> &flags,
				UINT64 &pathFound, vector <Edge *> &listOfEdges, vector <UINT64> &pathLengths);

		UINT64 findSupportByMatepairsAndMerge(void);
		UINT64 scaffolder(void);
		vector<Edge *> * getListOfFeasibleEdges(const Edge *edge);
		bool calculateMeanAndSdOfInnerDistance(void);
		INT64 checkForScaffold(const Edge *edge1, const Edge *edge2, INT64 &averageGapDistance);
		UINT64 findOverlap(const string & string1, const string & string2);
		UINT8 mergedEdgeOrientationDisconnected(const Edge *edge1, const Edge *edge2);
		UINT8 twinEdgeOrientation(UINT8 orientation);
		bool mergeListDisconnected(Edge *edge1, Edge *edge2, UINT64 overlapOffset,
				UINT64 **listReads, UINT32 &lSize);
		bool mergeEdgesDisconnected(Edge *edge1, Edge *edge2, INT64 gapLength);
		void simplifyScaffoldGraph(void);

		void loadSequence(const vector<std::string> &readSingleFilenameList,
				const vector<std::string> &readPairedFilenameList);

		void updateReadsLocations(Edge *edge, EdgeOP operation, DataSet *d);
		void updateEdgeInfo(Read * updateRead, Edge *edge, UINT32 read_index, EdgeOP operation);

		void printEdgeCoverage(Edge *contigEdge, ostream & filePointer, UINT64 edgeNameID) const;

		void generateGFAOutput(ostream & edgeFilePointer);

		void generateGFA2Output(ostream & gfaFilePointer);

		void generateGFA2Edge(ostream & gfaFilePointer, UINT64 edge_id, UINT64 source, string sOri,
				UINT64 destination,string dOri, UINT64 offset);

		void streamContigs(const vector<std::string> &read_SingleFiles,const vector<std::string> &read_PairFiles,
				vector<std::string> &read_PairInterFiles, string contig_file, string edge_file,string edge_cov_file,
				string usedReadFileName, string namePrefix, UINT64 &printed_contigs);

		void populate_read(const UINT64 &readID, const std::string & read_str);
		void loadStringFromReadsFile(const std::string &read_file, UINT64 & readID);

		UINT64 removeLowOvlEdges(void);

		void streamUnitigs(const vector<std::string> &read_SingleFiles,const vector<std::string> &read_PairFiles,
				const vector<std::string> &read_PairInterFiles, string unitig_file, string namePrefix, UINT64 &printed_contigs);
		void getUnitigExtensions(Edge *seedEdge, vector<unitigExt> &unitigEntensions);
		void exploreUnitigExtensions(UINT64 firstSrc, UINT64 nextSrc, vector<unitigExt> &unitigEntensions,string seq);
};

void createRevList(Edge *fwdEdge, UINT64 **returnListReads, UINT64 &lSize);

bool compareEdgesByDestID (const Edge *edge1, const Edge* edge2);
bool compareEdgesByLength (const Edge *edge1, const Edge* edge2);
#endif /* -----  end of class OverlapGraph  ----- */

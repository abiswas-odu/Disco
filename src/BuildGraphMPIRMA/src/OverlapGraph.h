/*
 * OverlapGraph.h
 *
 * Created on: April 22, 2013
 * Author: Abhishek Biswas
 */


#ifndef OVERLAPGRAPH_H_
#define OVERLAPGRAPH_H_
#include "../../BuildGraphMPIRMA/src/Common.h"
#include "../../BuildGraphMPIRMA/src/Dataset.h"
#include "../../BuildGraphMPIRMA/src/Edge.h"
#include "../../BuildGraphMPIRMA/src/HashTable.h"

/**********************************************************************************************************************
	Class to store the overlap graph.
**********************************************************************************************************************/

#define MIN_MARKED 128				//No. of reads to be marked before a send communication is initiated...

enum nodeType {
	UNEXPLORED = 0, // Current node u is not explored yet. Meaning that there is no edge (u,v) in the graph.
	EXPLORED = 1, //  Current node u is explored. Meaning that all edges (u,v) are inserted in the dataset.
	EXPLORED_AND_TRANSITIVE_EDGES_MARKED = 2, // Meaning that all transitive edges (u,v) of current node u is marked and its neighbors transitive edges are also marked. Now it is safe to remove any transitive edge (u,v) from node u.
	EXPLORED_AND_TRANSITIVE_EDGES_REMOVED = 3, // Meaning that all transitive edges (u,v) of current node u has been removed and its neighbors transitive edges are also marked. Now it is safe to remove any transitive edge (u,v) from node u.
	EXPLORED_AND_TRANSITIVE_EDGES_WRITTEN = 4
};

enum markType{
	VACANT = 0,
	INPLAY = 1,
	ELIMINATED = 2
};

class OverlapGraph
{
	private:
		Dataset * dataSet; 											// Pointer to the dataset containing all the reads.
		HashTable * hashTable;										// Pointer to the hash table.
		int *myMarked;												//List of reads already marked locally by this process or lazy globally
		int myProcID;													// Id of the MPI process
		UINT64 numberOfNodes;										// Number of nodes in the overlap graph.
		UINT64 numberOfEdges;										// Number of edges in the overlap graph.
		UINT64 parallelThreadPoolSize;								//No. of OMP threads to spawn
		UINT64 writeParGraphSize;									//No. of vertices to mark before writing graph to memory
		UINT8 twinEdgeOrientation(UINT8 orientation);				// Orientation of the reverse edge.
	public:
		OverlapGraph(HashTable *ht, UINT64 maxThreads,UINT64 maxParGraph,
				UINT64 maxMemSizeGB, string fnamePrefixGraph,string fnamePrefixSimplify,
				string simplifyPartialPath, int myid, int numprocs);								// Another constructor.
		~OverlapGraph();											// Destructor.
		bool markTransitiveEdges(UINT64 readNumber, map<UINT64, vector<Edge*> * > *parGraph); // Mark transitive edges of a read.
		bool buildOverlapGraphFromHashTable(HashTable *ht, string fnamePrefixGraph,string fnamePrefixSimplify,
				string simplifyPartialPath, int numprocs);			// Build the overlap graph using hashtable.
		bool insertEdge(Edge * edge, map<UINT64, vector<Edge*> * > *parGraph); 								// Insert an edge in the partial overlap graph.
		bool insertEdge(Read *read1, Read *read2, UINT64 r1Len, UINT64 r2Len,  UINT8 orient, UINT16 overlapOffset, map<UINT64, vector<Edge*> * > *parGraph); // Insert an edge in the overlap graph.

		bool checkOverlapForContainedRead(string read1, string read2, UINT64 orient, UINT64 start);
		bool checkOverlap(string read1, string read2, UINT64 orient, UINT64 start);

		bool insertAllEdgesOfRead(UINT64 readNumber, map<UINT64,nodeType> * exploredReads, map<UINT64, vector<Edge*> * > *parGraph);	// Insert into the overlap graph all edges of a read.
		bool removeTransitiveEdges(UINT64 readNumber, map<UINT64, vector<Edge*> * > *parGraph);				// Remove all transitive edges from the overlap graph incident to a given read.
		bool saveParGraphToFile(string fileName, map<UINT64,nodeType> * exploredReads, map<UINT64, vector<Edge*> * > *parGraph);   //Save partial graph to file and reduce memory footprint
		void markContainedReads(string fnamePrefix, map<UINT64, UINT64> *fIndxReadIDMap, int numprocs);									// Find superReads for each read and mark them as contained read.
};
#endif /* OVERLAPGRAPH_H_ */

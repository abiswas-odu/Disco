#ifndef OVERLAPGRAPH_H
#define OVERLAPGRAPH_H

/*
 * ===== CLASS HEADER ========================================================
 * Name        : OverlapGraphSimple.cpp
 * Author      : Abhishek Biswas
 * Version     : v3.0
 * Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
 * Description : OverlapGraph header file
 *============================================================================
 */

#include "Config.h"
#include "EdgeSimple.h"

#define SSTR( x ) dynamic_cast< std::ostringstream & >(( std::ostringstream() << std::dec << x ) ).str()


typedef vector<EdgeSimple *> t_edge_vec;	// vector of pointers to EdgeSimple
class OverlapGraphSimple
{
	private:
		/* ====================  DATA MEMBERS  ======================================= */
		UINT64 			m_numberOfNodes;
		UINT64 			m_numberOfEdges;
		UINT64 			m_minOvl;
		UINT64          p_ThreadPoolSize;

		/* ====================  METHODS       ======================================= */

		// Sort edges of each read based on ID of the destination read. 
		// This is only for ordering edges for convenience in the output file
		void sortEdgesByLength();

		void sortEdgesByDestID();

		/*Test if read has only 4 bases*/
		bool testRead(const string & read);

	public:

		/* ====================  LIFECYCLE     ======================================= */

		OverlapGraphSimple(string edge_file, string composite_out_edge_file, UINT64 minOvl, UINT64 parallelThreadPoolSize);

		/* ====================  ACCESSORS     ======================================= */
		UINT64 getNumberOfEdges(void) const {return m_numberOfEdges;}

		UINT64 getNumberOfNodes(void) const {return m_numberOfNodes;}

		/* ====================  OPERATORS     ======================================= */

		//Functions for parallel graph simplification

		void printEdge(EdgeSimple *contigEdge, ostream & filePointer) const;
		void loadParEdgesFromEdgeFile(const std::string &readFilename, map<UINT64, t_edge_vec* > *parGraph, set<UINT64> &markedNodes);
		void insertParFwdEdge( EdgeSimple *edge, map<UINT64, t_edge_vec* > *parGraph);
		void insertParEdge( EdgeSimple *edge, map<UINT64, t_edge_vec* > *parGraph);
		UINT64 contractParCompositeEdges(map<UINT64, t_edge_vec* > *parGraph, set<UINT64> &markedNodes);
		UINT64 removeParDeadEndNodes(map<UINT64, t_edge_vec* > *parGraph, set<UINT64> &markedNodes, vector<UINT64> &nodeList);
		void removeParEdge(EdgeSimple *edge, map<UINT64, t_edge_vec* > *parGraph);
		void removeParFwdEdge(EdgeSimple *edge, map<UINT64, t_edge_vec* > *parGraph);
		void sortParEdgesByDestID(map<UINT64, t_edge_vec* > *parGraph);
		bool existsParEdge(EdgeSimple *checkEdge, map<UINT64, t_edge_vec* > *parGraph);
		void printParEdges(string edge_file, map<UINT64, t_edge_vec* > *parGraph) const;
		void removeParEdgeFromSourceRead(EdgeSimple *edge, map<UINT64, t_edge_vec* > *parGraph);
		UINT8 twinEdgeOrientation(UINT8 orientation);
		UINT64 contractParCompositeEdges_Serial(map<UINT64, t_edge_vec* > *parGraph, set<UINT64> &markedNodes);

};



bool compareParEdgesByDestID (const EdgeSimple *edge1, const EdgeSimple* edge2);
bool compareParEdgesByLength (const EdgeSimple *edge1, const EdgeSimple* edge2);
#endif /* -----  end of class OverlapGraph  ----- */

#ifndef EDGE_H
#define EDGE_H

//============================================================================
// Name        : Edge.h
// Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey, Abhishek Biswas
// Version     : v1.2
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : Edge header file
//============================================================================


#include "Config.h"
#include "DataSet.h"

class Read;

extern TLogLevel loglevel;                      /* verbosity level of logging */
typedef vector<pair<UINT32, UINT32> > t_vpair;
class Edge{
	private:
		/* ====================  DATA MEMBERS  ======================================= */
		Read *m_source;
		Read *m_destination;
		Edge *m_reverseEdge;

		// Packed list of ordered reads in the current edge. NOT including u and v.
		// First 32LSB ReadID
		// First MSB orientation
		/*  Orientations of the ordered reads in the current edge.
		 *  orientation is 0 or 1. 
		 *  0 means read's reverse complement read.getStringReverse(). 
		 *  1 means read's forward string.
		 */
		/*  Middle 31MSB overlap offsets of the ordered reads in the current edge.
		*  The offset of the destination read is calculated as
		*  m_overlapOffset - sum(getOverlapOffset(i)) or it can retrieved from the reverse edge
		*/
		UINT64 *m_listOfReads;

		float m_coverageDepth;           // Average depth of coverage.

		float m_SD;                      // Standard deviation of the base-by-base coverage depth

		// String generated by the edge
		std::string m_string;

		UINT32 m_listSize;				// Size of m_listOfReads


		/* overlap offset in the following example is 6
		*   012345678901234567890123456789
		* u ACTTACGGGATTATACCATCGAGA
		* v       GGGATTATACCATCGAGATTCAAT
		*/
		UINT32 m_overlapOffset;

		/* 0 = u<-----------<v		reverse of u to reverse of v
		 * 1 = u<----------->v		reverse of u to forward of v
		 * 2 = u>-----------<v		forward of u to reverse of v
		 * 3 = u>----------->v		forward of u to forward of v
		 */
		UINT8 m_orient;

		/*  flag of edge that has these bits:
		 *  0x1: not dead end
		 *  0x2: loop
		 */
		INT8 m_flag;

		/*  flag of edge that has these bits:
			 *  0x0: valid edge
			 *  0x1: invalidated edge; to be deleted
		 */
		bool m_invalid;

		/* ====================  METHODS      ======================================= */
		friend Edge*  merge_forward_edges(const Edge & edge1, const Edge & edge2);

		void copyEdge(const Edge &edge);

		void clearEdge();

		UINT64* getReverseListOfReads();

	public:
		/* ====================  DATA MEMBERS  ======================================= */
		UINT32 m_flow;                    // Store the flow in the current edge.
		
		/* ====================  LIFECYCLE     ======================================= */
		// Default constructor.
		Edge(void);								

		// constructor for a non-empty edge
		Edge(Read *source, Read *destination, UINT8 orient, UINT32 overlapOffset);

		Edge(Read *source, Read *destination, UINT8 orient, UINT32 overlapOffset,
				UINT64 *listOfReads, UINT32 listSize);

		Edge(const Edge &edge);
		// Destructor.
		~Edge();								

		/* ====================  OPERATORS     ======================================= */
		Edge& operator=(const Edge &edge);

		friend bool operator== (const Edge &subject, const Edge &query);

		friend Edge* Add( const Edge *edge1, const Edge *edge2);

		/* ====================  MUTATORS      ======================================= */
		// Set the pointer to the reverse edge, only if reverse edge is nullptr

		void setReverseEdge(Edge * edge);

		void markNotDeadEnd() {m_flag |= 1; m_reverseEdge->m_flag |= 1;}

		void markLoop() {m_flag |= (1 << 1);}

		void loadReadString(const std::string read_str, int index);

		void make_reverseEdge();

		void setInvalid() {m_invalid=true;}

		/* ====================  ACCESSORS     ======================================= */

		//Get coverage values of each base pair of the edge
		vector<UINT64> getBaseByBaseCoverageValues(DataSet *d);

		// Get the read object of the source node.
		Read * getSourceRead() const {return m_source;}	

		// Check if an edge is marked as not dead end
		bool isNotDeadEnd() const { return (m_flag & 1);}

		// Check if an edge is a loop
		bool isLoop() const {return (m_flag >> 1) & 1;}

		// Get the read object of the destination node.
		Read * getDestinationRead() const {return m_destination; }	

		// Return the orientation of the edge.
		UINT8 getOrientation() const {return m_orient;}	

		// Return the overlap offset.
		UINT32 getOverlapOffset() const {return m_overlapOffset;}	

		UINT64 getInnerReadInfo(int indx) const
		{
			return m_listOfReads[indx];
		}

		UINT64 getInnerReadID(int indx) const
		{
			return (m_listOfReads[indx] & 0X00000000FFFFFFFF);
		}

		UINT32 getInnerOverlapOffset(int indx) const
		{
			return ((m_listOfReads[indx] >> 32) & 0X000000007FFFFFFF);
		}

		UINT8 getInnerOrientation(int indx) const
		{
			return (m_listOfReads[indx] >> 63);
		}

		bool isListofReads() const { return !(m_listOfReads==nullptr);}

		size_t getListofReadsSize() const { return m_listSize; }

		UINT32 getInnerOverlapSum(size_t start, size_t end) const;

		// Return the pointer to the reverse edge.
		Edge * getReverseEdge() const {return m_reverseEdge;}	

		UINT32 getOverlapLen() const;

		UINT32 getFirstOverlapOffset() const 
		{
			if (!m_listOfReads || m_listSize==0)
				return m_overlapOffset;
			else
				return (getInnerOverlapOffset(0));
		}

		UINT64 getEdgeLength() const {return m_overlapOffset + m_destination->getReadLength();}

		UINT32 getLastOverlapOffset() const;

		std::string getEdgeString() const {return m_string;}

		bool isSmallerEdge() const 
		{
			if(m_source->getReadID() < m_destination->getReadID())
				return true;
			else if (m_source->getReadID() > m_destination->getReadID())
				return false;
			else if(this < m_reverseEdge)
				return true;
			else
				return false;
		}

		// Get base by base coverage
		void updateBaseByBaseCoverageStat(DataSet *d);

		float getCovDepth() const {return m_coverageDepth;}

		float getCovSD() const {return m_SD;}

		INT8 getFlag() const {return m_flag;}

		bool isInvalid() {return m_invalid;}
		/* ====================  OPERATIONS     ======================================= */
		// Break edge at a specific link if the edge is composite, otherwise return empty vector
		// Just do the forward edge
		vector<Edge*> breakForwardEdge(const UINT32 &link, DataSet *d) const;

		// Break edge at specific link, but also do the twin edge accordingly
		vector<Edge*> breakEdge(const UINT32 &link, DataSet *d) const;

};

bool is_mergeable(const Edge *edge1, const Edge *edge2);

UINT8 mergedEdgeOrientation(const UINT8 &orient1, const UINT8 &orient2);

UINT8 get_twin_orient(const UINT8 &orient);

void mergeList(const Edge *edge1, const Edge *edge2, 
		UINT64 **listReads, UINT64 &lSize);

Edge * mergeEdges(const vector<Edge *> & list_edges);

void createFwdList(string readList,UINT64 **returnListReads, UINT64 &lSize);

void createRevList(Edge *fwdEdge, UINT64 **returnListReads, UINT64 &lSize, DataSet *d);

template <typename T>
float get_mean(const vector<T> &numbers);

template <typename T>
float get_sd(const vector<T> &numbers);
#endif /* EDGE_H */

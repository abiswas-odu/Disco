#ifndef EDGE_H
#define EDGE_H

//============================================================================
// Name        : EdgeSimple.h
// Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey, Abhishek Biswas
// Version     : v1.2
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : EdgeSimple header file
//============================================================================


#include "Config.h"
#include "Utils.h"

extern TLogLevel loglevel;                      /* verbosity level of logging */
typedef vector<pair<UINT32, UINT32> > t_vpair;
class EdgeSimple{
	private:
		/* ====================  DATA MEMBERS  ======================================= */
		UINT64 m_source;
		UINT64 m_sourceLen;
		UINT64 m_destination;
		UINT64 m_destinationLen;
		EdgeSimple *m_reverseEdge;

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

		friend EdgeSimple*  merge_forward_edges(const EdgeSimple & edge1, const EdgeSimple & edge2);

		void copyEdge(const EdgeSimple &edge);

		void clearEdge();

	public:

		/* ====================  LIFECYCLE     ======================================= */
		// Default constructor.
		EdgeSimple(void);

		// constructor for a non-empty edge
		EdgeSimple(UINT64 source, UINT64 sourceLen, UINT64 destination, UINT64 destLen, UINT8 orient, UINT32 overlapOffset);

		EdgeSimple(UINT64 source, UINT64 sourceLen, UINT64 destination, UINT64 destLen, UINT8 orient, UINT32 overlapOffset,
				UINT64 *listOfReads, UINT32 listSize);

		EdgeSimple(const EdgeSimple &edge);
		// Destructor.
		~EdgeSimple();

		/* ====================  OPERATORS     ======================================= */
		EdgeSimple& operator=(const EdgeSimple &edge);

		friend bool operator== (const EdgeSimple &subject, const EdgeSimple &query);

		friend EdgeSimple* Add( const EdgeSimple *edge1, const EdgeSimple *edge2);

		/* ====================  MUTATORS      ======================================= */
		// Set the pointer to the reverse edge, only if reverse edge is nullptr
		void setReverseEdge(EdgeSimple * edge);

		void markNotDeadEnd() {m_flag |= 1; m_reverseEdge->m_flag |= 1;}

		void markLoop() {m_flag |= (1 << 1);}

		void loadReadString(const std::string read_str, int index);

		void make_nonComposite_reverseEdge();

		void setInvalid() {m_invalid=true;}

		/* ====================  ACCESSORS     ======================================= */
		// Get the read object of the source node.
		UINT64 getSourceRead() const {return m_source;}

		// Check if an edge is marked as not dead end
		bool isNotDeadEnd() const { return (m_flag & 1);}

		// Check if an edge is a loop
		bool isLoop() const {return (m_flag >> 1) & 1;}

		// Get the read object of the destination node.
		UINT64 getDestinationRead() const {return m_destination; }

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
		EdgeSimple * getReverseEdge() const {return m_reverseEdge;}

		UINT32 getOverlapLen() const;

		UINT32 getFirstOverlapOffset() const 
		{
			if (!m_listOfReads || m_listSize==0)
				return m_overlapOffset;
			else
				return (getInnerOverlapOffset(0));
		}

		UINT64 getEdgeLength() const {return m_overlapOffset + m_destinationLen;}

		UINT32 getLastOverlapOffset() const;

		bool isSmallerEdge() const 
		{
			if(m_source < m_destination)
				return true;
			else if (m_source > m_destination)
				return false;
			else if(this < m_reverseEdge)
				return true;
			else
				return false;
		}

		INT8 getFlag() const {return m_flag;}

		bool isInvalid() {return m_invalid;}

};

bool is_mergeable(const EdgeSimple *edge1, const EdgeSimple *edge2);

UINT8 mergedEdgeOrientation(const UINT8 &orient1, const UINT8 &orient2);

UINT8 get_twin_orient(const UINT8 &orient);

void mergeList(const EdgeSimple *edge1, const EdgeSimple *edge2,
		UINT64 **listReads, UINT64 &lSize);

EdgeSimple * mergeEdges(const vector<EdgeSimple *> & list_edges);

template <typename T>
float get_mean(const vector<T> &numbers);

template <typename T>
float get_sd(const vector<T> &numbers);
#endif /* EDGE_H */

//============================================================================
// Name        : EdgeSimple.cpp
// Author      : Abhishek Biswas
// Version     : v1.2
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : EdgeSimple cpp file
//============================================================================

#include "Config.h"
#include "EdgeSimple.h"
#include <numeric>

EdgeSimple::EdgeSimple(void)
{
	m_orient		= 0;
	m_overlapOffset	= 0;
	m_flag			= 0;
	m_invalid       = false;
	m_source        = 0;
	m_sourceLen     = 0;
	m_destination   = 0;
	m_destinationLen = 0;
	m_listOfReads    = nullptr;
	m_listSize       = 0;
	m_reverseEdge    = nullptr;
}

// Constructor for a (composite) edge, given
EdgeSimple::EdgeSimple(UINT64 source, UINT64 sourceLen, UINT64 destination, UINT64 destLen, UINT8 orient, UINT32 overlapOffset,
		UINT64 *listOfReads, UINT32 listSize)
:m_source(source),m_sourceLen(sourceLen), m_destination(destination), m_destinationLen(destLen),m_reverseEdge(nullptr),
 m_listOfReads(listOfReads), m_listSize(listSize),m_overlapOffset(overlapOffset), m_orient(orient),m_flag(0), m_invalid(false)
{
	// If this edge is a simple edge, make its reverse edge
	// For composite edge, leave the reverse edge to the sum of reverse edges
	if(source == destination)
		m_flag |= (1 << 1);
}

EdgeSimple::EdgeSimple(UINT64 source, UINT64 sourceLen, UINT64 destination, UINT64 destLen, UINT8 orient, UINT32 overlapOffset)
:m_source(source),m_sourceLen(sourceLen), m_destination(destination), m_destinationLen(destLen),m_reverseEdge(nullptr),
 m_listOfReads(nullptr), m_listSize(0), m_overlapOffset(overlapOffset), m_orient(orient), m_flag(0), m_invalid(false)
{
	// If this edge is a simple edge, make its reverse edge
	// For composite edge, leave the reverse edge to the sum of reverse edges
	if(source == destination)
		m_flag |= (1 << 1);
}

void EdgeSimple::copyEdge(const EdgeSimple &edge)
{
	//CLOCKSTART;
//	FILE_LOG(logDEBUG1) << edge << "\n";
	//FILE_LOG(logDEBUG1) << "Copy number variables\n";
	m_source = edge.m_source;
	m_destination = edge.m_destination;
	m_orient = edge.m_orient;
	m_overlapOffset = edge.m_overlapOffset;
	m_flag = edge.m_flag;
	m_invalid  = edge.m_invalid;
	m_listSize=edge.m_listSize;
	//FILE_LOG(logDEBUG1) << "Copy list of reads\n";
	if(edge.m_listOfReads){
		m_listOfReads = new UINT64[edge.m_listSize];
		for(size_t i=0; i < m_listSize; i++)
			m_listOfReads[i]=edge.m_listOfReads[i];
	}
	else 
		m_listOfReads = nullptr;

}

EdgeSimple::EdgeSimple(const EdgeSimple &edge)
{
	copyEdge(edge);
	m_reverseEdge = new EdgeSimple;
	m_reverseEdge->copyEdge(*(edge.m_reverseEdge));
	m_reverseEdge->setReverseEdge(this);
}

EdgeSimple& EdgeSimple::operator=(const EdgeSimple &edge)
{
	if(this == &edge)
		return *this;
	clearEdge();
	//CLOCKSTART;
	copyEdge(edge);
	m_reverseEdge = new EdgeSimple;
	m_reverseEdge->copyEdge(*(edge.m_reverseEdge));
	m_reverseEdge->setReverseEdge(this);
	//CLOCKSTOP;
	return *this;
}

EdgeSimple::~EdgeSimple()
{
	delete[] m_listOfReads;
	m_listOfReads = nullptr;
}
void EdgeSimple::clearEdge()
{
	delete[] m_listOfReads;
	m_listOfReads = nullptr;
	// This function does not delete the reverse edge
}
//Used to make a non composite reverse edge only
void EdgeSimple::make_nonComposite_reverseEdge()
{
	// If reverse edge is already made, do not make it again
	if(m_reverseEdge)
		return;
	// Get the lists of info contained in the composite edges.
	// For simple edges, these will simply return nullptr.
	UINT64 *reverseListOfReads = nullptr;
	// Make the new reverse edge
	// And set its reverse edge to the current edge
	m_reverseEdge = new EdgeSimple(m_destination, m_destinationLen, m_source, m_sourceLen, get_twin_orient(m_orient),
			m_destinationLen + m_overlapOffset - m_sourceLen,reverseListOfReads, m_listSize);

//	FILE_LOG(logDEBUG1) << "EdgeSimple: " << *this << "\n";
//	FILE_LOG(logDEBUG1) << "Reverse edge: " << *m_reverseEdge << "\n";
	m_reverseEdge->setReverseEdge(this);
}

/* Get the overlap length for the first link in the edge */
UINT32 EdgeSimple::getOverlapLen() const
{
	if(!m_listOfReads || m_listSize==0)
		return (m_sourceLen - m_overlapOffset);
	else
		return (m_sourceLen - getInnerOverlapOffset(0));
}

UINT32 EdgeSimple::getInnerOverlapSum(size_t start, size_t end) const
{
	UINT32 overlap_sum=0;
	for(size_t i=start;i<m_listSize && i<end;i++)
		overlap_sum += getInnerOverlapOffset(i);
	return overlap_sum;
}

/* Get the overlap length for the last link in the edge */
UINT32 EdgeSimple::getLastOverlapOffset() const
{
	if(!m_listOfReads || m_listSize==0)
		return m_overlapOffset;
	UINT32 overlap_sum = getInnerOverlapSum(0,m_listSize);
	return (m_overlapOffset - overlap_sum);
}




//=============================================================================
// Merge two edges in the overlap graph.
//=============================================================================
EdgeSimple* Add( const EdgeSimple *edge1, const EdgeSimple *edge2)
{
	assert (is_mergeable(edge1, edge2));
	EdgeSimple *merge_forward = merge_forward_edges(*edge1, *edge2);
	EdgeSimple *merge_reverse = merge_forward_edges(*(edge2->m_reverseEdge), *(edge1->m_reverseEdge));
	merge_forward->setReverseEdge(merge_reverse);
	merge_reverse->setReverseEdge(merge_forward);
	return merge_forward;
}

/* EdgeSimple operator+( const EdgeSimple & edge1, const EdgeSimple & edge2)
 * {
 * 	assert (is_mergeable(&edge1, &edge2));
 * 	CLOCKSTART;
 * 	EdgeSimple merge_forward = merge_forward_edges(edge1, edge2);
 * 	EdgeSimple merge_reverse = merge_forward_edges(*(edge2.m_reverseEdge), *(edge1.m_reverseEdge));
 * 	merge_forward.setReverseEdge(&merge_reverse);
 * 	merge_reverse.setReverseEdge(&merge_forward);
 * 	CLOCKSTOP;
 * 	return merge_forward;
 * }
 */

void EdgeSimple::setReverseEdge(EdgeSimple * edge)
{ 
	if(m_reverseEdge == edge)
		return;
	else if(!m_reverseEdge) {
		m_reverseEdge = edge;
	}
	else{
		m_reverseEdge = edge;
	}
}

// Assisting function for just merging forward edges
EdgeSimple*  merge_forward_edges(const EdgeSimple & edge1, const EdgeSimple & edge2)
{

	// Orientation
	UINT8 orientationForward = mergedEdgeOrientation(edge1.m_orient, edge2.m_orient);

	// overlap offset
	UINT32 overlapOffsetForward = edge1.m_overlapOffset + edge2.m_overlapOffset;

	UINT64 *listReadsForward = nullptr;
	UINT64 lSize=0;

	// Merge the lists from the two edges.
	mergeList(&edge1, &edge2, &listReadsForward, lSize);

	// Make the forward edge
	EdgeSimple *edgeForward = new EdgeSimple(edge1.m_source, edge1.m_sourceLen, edge2.m_destination, edge2.m_destinationLen, orientationForward,
			overlapOffsetForward, listReadsForward, lSize);

	return edgeForward;
}

//=============================================================================
// Merge the list of reads, list of overlap offsets and list of orientations of two edges.
//=============================================================================
void mergeList(const EdgeSimple *edge1, const EdgeSimple *edge2,
		UINT64 **returnListReads, UINT64 &lSize)
{
//	CLOCKSTART;
	lSize=edge1->getListofReadsSize()+edge2->getListofReadsSize()+1;
	size_t lCtr=0;
	UINT64 *listReads = new UINT64[lSize];
	// Copy the list from edge1.
	if (edge1->isListofReads() && !(edge1->getListofReadsSize()==0) ){
		for(size_t i = 0; i < edge1->getListofReadsSize(); ++i){
			listReads[lCtr++]=edge1->getInnerReadInfo(i);
		}
	}
	UINT64 overlapOff = (edge1->getLastOverlapOffset() << 32);
	UINT64 orient = (edge1->getOrientation() & 1);
	orient = orient << 63;
	// Insert the common node of the two edges
	UINT64 rID =  edge1->getDestinationRead() | overlapOff | orient;
	listReads[lCtr++] = rID;

	// Concatenate the list from edge2.
	if (edge2->isListofReads() && !(edge2->getListofReadsSize()==0) ){
		for(size_t i = 0; i <edge2->getListofReadsSize(); ++i){
			listReads[lCtr++]=edge2->getInnerReadInfo(i);
		}
	}
	*returnListReads = listReads;
//	CLOCKSTOP;
}

//=============================================================================
// Check if two edges can be merged into one edge
// For two edges e1(u,v) and e2(v,w), at node v, 
// one of the edges should be an incoming edge and the other should be an outgoing
// edge to match.
//=============================================================================
bool is_mergeable(const EdgeSimple *edge1, const EdgeSimple *edge2)
{
	// First, the destination of edge1 has to be the same as 
	// the source read of edge2
	if (edge1->getDestinationRead() != edge2->getSourceRead()){
		return false;
	}

	// *-----> and >------* *1 and 1*
	// *-----< and <------* *0 and 0*
	else if ((edge1->getOrientation() & 1) == ((edge2->getOrientation() >>1) & 1))
		return true;
	else{
		return false;
	}
}

// Orientation of the edge when two edges are merged.
UINT8 mergedEdgeOrientation(const UINT8 &orient1, const UINT8 &orient2)
{
//	assert(orient1 < 4 && orient2 < 4);
	return ((orient1 & 2) | (orient2 & 1));
}

UINT8 get_twin_orient(const UINT8 &orient)
{
	assert(static_cast<int>(orient) < 4); // Quit if orient is not 0, 1, 2, or 3
	// Exchange the last two bits of the orient, and flip them
	UINT8 twin_orient = ((orient >> 1) ^ 1) | (((orient & 1) ^ 1) << 1) ;
	assert(static_cast<int>(twin_orient) < 4);
	return twin_orient;
}

bool isSameCompositePath(const EdgeSimple &subject, const EdgeSimple &query)
{
	if(subject.isListofReads() && query.isListofReads())
	{
		if(subject.getListofReadsSize()!=query.getListofReadsSize())
			return false;
		else
		{
			for(size_t i=0;i<subject.getListofReadsSize();i++)
			{
				UINT64 sRID = subject.getInnerReadID(i);
				UINT64 dRID = query.getInnerReadID(i);
				if(sRID!=dRID)
					return false;
			}
		}
		return true;
	}
	else if(!subject.isListofReads() && !query.isListofReads())
		return true;
	return false;
}

/*
 * Compares source, destination, orientation, overlap and path of the edges
 *
 */
bool operator==(const EdgeSimple &subject, const EdgeSimple &query)
{
	if(subject.getSourceRead()== query.getSourceRead() &&
			subject.getDestinationRead() == query.getDestinationRead() &&
			subject.m_overlapOffset == query.m_overlapOffset &&
			subject.m_orient == query.m_orient &&
			isSameCompositePath(subject, query))
		return true;
	return false;
}


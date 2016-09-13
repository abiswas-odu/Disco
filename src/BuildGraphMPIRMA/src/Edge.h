/*
 * Edge.h
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 */

#ifndef EDGE_H_
#define EDGE_H_
#include "../../BuildGraphMPIRMA/src/Common.h"
#include "../../BuildGraphMPIRMA/src/Read.h"

/**********************************************************************************************************************
	Class to store an edge.
**********************************************************************************************************************/
class Edge{
	private:
		Read *source; 							// Source read u
		UINT16 srcLen;
		Read *destination; 						// Destination read v
		UINT16 destLen;
		UINT8 overlapOrientation;				// Orientation of overlap
												// 0 = u<-----------<v
												// 1 = u<----------->v
												// 2 = u>-----------<v
												// 3 = u>----------->v
		UINT64 overlapOffset;					// Length of the overlap.
												// overlap offset in the following example is 6
												// 	  012345678901234567890123456789
												//	u ACTTACGGGATTATACCATCGAGA
												//	v       GGGATTATACCATCGAGATTCAAT
		Edge *reverseEdge;						// This points to the reverse edge in the overlap graph.
												// Each edge (u,v) is present twice in the graph.
												// Edge (u,v) is present in the list of edges of node u
												// Edge (v,u) is present in the list of edges of node v.
												// Edge (u,v) and (v,u) are called twin edge and are reverse of one another.
	public:
		bool transitiveRemovalFlag;							// Used to mark transitive edges.
		Edge(void);								// Default constructor.
		Edge(Read *from,UINT16 fromLen, Read *to, UINT16 toLen, UINT64 orient, UINT64 length); 	// Another constructor.
		~Edge();								// Destructor.
		bool makeEdge(Read *from,UINT16 fromLen, Read *to,UINT16 toLen, UINT64 orient, UINT64 length);
		string getStringInEdge(void); 			// return the string in the edge
		bool setReverseEdge(Edge * edge);		// Set the pointer to the reverse edge.
		Read * getSourceRead() {return source;}	// Get the read object of the source node.
		Read * getDestinationRead() {return destination; }	// Get the read object of the destination node.
		UINT8 getOrientation() {return overlapOrientation;}	// Return the orientation of the edge.
		UINT64 getOverlapOffset() {return overlapOffset;}	// Return the overlap offset.
		Edge * getReverseEdge() const {return reverseEdge;}	// Return the pointer to the reverse edge.

		void setSrcLen(UINT16 len) {srcLen=len;};
		void setDestLen(UINT16 len) {destLen=len;};

		UINT16 getSrcLen() { return srcLen;};
		UINT16 getDestLen() { return destLen;};
};

#endif /* EDGE_H_ */

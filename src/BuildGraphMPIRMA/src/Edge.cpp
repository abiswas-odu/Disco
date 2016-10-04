/*
 * Edge.cpp
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 */

#include "Edge.h"

#include "Common.h"


/**********************************************************************************************************************
	Default Constructor
**********************************************************************************************************************/
Edge::Edge(void)
{

	// Initialize the variables.
	overlapOffset = 0;
	overlapOrientation = 10;
	transitiveRemovalFlag = false;
}



/**********************************************************************************************************************
	Another Constructor
**********************************************************************************************************************/
Edge::Edge(Read *from,UINT16 fromLen, Read *to, UINT16 toLen, UINT64 orient, UINT64 length)
{
	makeEdge(from,fromLen, to,toLen, orient, length);
}

/**********************************************************************************************************************
	Default destructor
**********************************************************************************************************************/
Edge::~Edge()
{
	// Free the memory used by the current edge.
}



/**********************************************************************************************************************
	Function to insert a simple edge
**********************************************************************************************************************/
bool Edge::makeEdge(Read *from,UINT16 fromLen, Read *to,UINT16 toLen, UINT64 orient, UINT64 length)
{
	source = from;
	destination = to;
	overlapOrientation = orient;
	overlapOffset = length;
	srcLen=fromLen;
	destLen=toLen;
	// Initialize variables.
	transitiveRemovalFlag = false;
	return true;
}

/**********************************************************************************************************************
	Function to set the pointer to the reverse edge;
**********************************************************************************************************************/
bool Edge::setReverseEdge(Edge * edge)
{
	reverseEdge = edge;
	return true;
}



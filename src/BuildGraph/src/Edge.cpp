/*
 * Edge.cpp
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 */

#include "../src/Edge.h"

#include "../src/Common.h"


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
Edge::Edge(Read *from, Read *to, UINT64 orient, UINT64 length)
{
	makeEdge(from, to, orient, length);
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
bool Edge::makeEdge(Read *from, Read *to, UINT64 orient, UINT64 length)
{
	source = from;
	destination = to;
	overlapOrientation = orient;
	overlapOffset = length;
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



//============================================================================
// Name        : Read.cpp
// Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey, Abhishek Biswas
// Version     : v3.0
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : Read cpp file
//============================================================================


#include "Config.h"
#include "Read.h"

//Expose only is read sequence is being stored
/*Read::Read(const std::string &seq)
{
	m_seq = new dna_bitset(seq);
	m_readID	= 0;
	edgeP = nullptr;
	edgeOriIndex = nullptr;
	noOfEdges=0;
	noOfAllocEdgeMemAvail=0;
	containedReads = nullptr;
	noOfConReads=0;
	containedReadFlag=false;
	usedRead=false;
}*/

Read::Read(const UINT64 seqLen)
{
	//m_seq = new dna_bitset(seqLen);
	readLen=seqLen;
	m_readID	= 0;
	edgeP = nullptr;
	edgeOriIndex = nullptr;
	noOfEdges=0;
	noOfAllocEdgeMemAvail=0;
	containedReads = nullptr;
	noOfConReads=0;
	containedReadFlag=false;
	usedRead=false;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Destructor
 *  Description:  Destructor for Read object
 * =====================================================================================
 */
Read::~Read(void)
{

	//delete m_seq;   //expose only id read sequence is being stored
	delete[] edgeP;
	delete[] edgeOriIndex;
	delete[] containedReads;
}
/*
 * Add a contained read to the contained read list...
 */
void Read::setConRead(UINT64 conReadID, UINT32 conOvlStart, UINT64 orient)
{
	noOfConReads++;
	UINT64 *newConReadArray = new UINT64[noOfConReads];
	for(UINT16 i=0;i<(noOfConReads-1);i++)
	{
		newConReadArray[i]=containedReads[i];
	}

	UINT64 overlapStart = (conOvlStart << 32);
	orient = orient << 62;
	UINT64 newConRead =  conReadID | overlapStart | orient;
	newConReadArray[noOfConReads-1] = newConRead;
	delete[] containedReads;
	containedReads = newConReadArray;
}
void Read::delEdge(Edge *edge, UINT64 readIndx, UINT64 orient)
{
	size_t delPos=noOfEdges;
	orient = orient << 31;
	UINT64 delOriIndx =  readIndx | orient;
	for(size_t i=0;i<noOfEdges;i++)
	{
		if(edge==edgeP[i] && delOriIndx==edgeOriIndex[i]){
			delPos=i;
			break;
		}
	}
	for(UINT16 i=delPos;i<(noOfEdges-1);i++)
	{
		edgeP[i] = edgeP[i+1];
		edgeOriIndex[i] = edgeOriIndex[i+1];
	}
	noOfEdges--;
}
/*
 * Clear all edge information from the Read
 */
void Read::ClearEdgeInfo()
{
	noOfEdges=0;
}
/*std::ostream &operator<<(std::ostream & out, const Read & read)
{
	out << "ID: " << setw(10) << setfill(' ') << read.getReadID()
		<< ", String: \n" << read.getStringForward();
	return out;
}*/

vector<t_edge_loc_pair> Read::getFwdEdges() const
{
	vector<t_edge_loc_pair> edgePairs;
	for(size_t i=0;i<noOfEdges;i++)
	{
		Edge *edgePtr = edgeP[i];
		UINT32 readIndex = getEdgeReadIndx(i);
		UINT64 orient = getEdgeOriIndx(i);
		if(orient==0)
			edgePairs.push_back(std::make_pair(edgePtr,readIndex));
	}
	return edgePairs;
}

vector<t_edge_loc_pair> Read::getBwdEdges() const
{
	vector<t_edge_loc_pair> edgePairs;
	for(size_t i=0;i<noOfEdges;i++)
	{
		Edge *edgePtr = edgeP[i];
		UINT32 readIndex = getEdgeReadIndx(i);
		UINT64 orient = getEdgeOriIndx(i);
		if(orient==1)
			edgePairs.push_back(std::make_pair(edgePtr,readIndex));
	}
	return edgePairs;
}

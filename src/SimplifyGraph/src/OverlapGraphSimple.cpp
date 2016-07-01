/*
 * ===== CLASS IMPLEMENTATION ================================================
 * Name        : OverlapGraphSimple.cpp
 * Author      : Abhishek Biswas
 * Version     : v1.2
 * Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
 * Description : OverlapGraph cpp file
 *============================================================================
 */

#include "OverlapGraphSimple.h"
#include "Utils.h"
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  insertEdge
 *  Description:  Insert an edge in the partial overlap graph passed as argument.
 *                Does not automatically insert its twin edge.
 * =====================================================================================
 */
void OverlapGraphSimple::insertParEdge( EdgeSimple *edge, map<UINT64, t_edge_vec* > *parGraph)
{
	insertParFwdEdge(edge, parGraph);
	insertParFwdEdge(edge->getReverseEdge(),parGraph);
}

void OverlapGraphSimple::insertParFwdEdge(EdgeSimple *edge, map<UINT64, t_edge_vec* > *parGraph)
{
	UINT64 ID = edge->getSourceRead();
	if(parGraph->find(ID) == parGraph->end()){ 			// If there is no edge incident to the node
			t_edge_vec *newList = new t_edge_vec;
			parGraph->insert( std::pair<UINT64, vector<EdgeSimple*> * >(ID, newList));
	}
	parGraph->at(ID)->push_back(edge);						// Insert the edge in the list of edges of ID
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  removeParEdge
 *  Description:  Remove an edge from the partial overlap graph passed as argument,
 *            release its memory if we don't need to save this edge for the contigs.
 *  		  This removes and deletes from memory both the edge and its twin edge
 * =====================================================================================
 */
void OverlapGraphSimple::removeParEdge(EdgeSimple *edge, map<UINT64, t_edge_vec* > *parGraph){
	if(edge == nullptr)
		return;
	removeParEdgeFromSourceRead(edge, parGraph);
	removeParEdgeFromSourceRead(edge->getReverseEdge(), parGraph);
	delete edge->getReverseEdge();
	delete edge;

}
void OverlapGraphSimple::removeParEdgeFromSourceRead(EdgeSimple *edge, map<UINT64, t_edge_vec* > *parGraph)
{
	if (!edge)
		return;
	// all edges incident to source read
	t_edge_vec *fwd_edges = parGraph->at(edge->getSourceRead());
	fwd_edges->erase(std::remove(fwd_edges->begin(), fwd_edges->end(), edge), fwd_edges->end());
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  contractCompositeEdges
 *  Description:  Contract composite paths in the partial overlap graph, changing only marked reads.
 *                u*-------->v>---------*w  => u*--------------------*w
 *                u*--------<v<---------*w  => u*--------------------*w
 * =====================================================================================
 */
UINT64 OverlapGraphSimple::contractParCompositeEdges_Serial(map<UINT64, t_edge_vec* > *parGraph, set<UINT64> &markedNodes)
{
	CLOCKSTART;
	UINT64 counter(0);
	set<UINT64>::iterator it;
	for (it=markedNodes.begin(); it!=markedNodes.end(); ++it)
	{
		UINT64 rid=*it;
		map<UINT64, vector<EdgeSimple*> * >::iterator it=parGraph->find(rid);
		if(it->second->size() == 2) // Check if the node has only two edges.
		{
			// First edge, going into Read
			EdgeSimple *edge1 = it->second->at(0)->getReverseEdge();
			UINT64 source1ID = edge1->getSourceRead();
			// Second edge, going out from Read
			EdgeSimple *edge2 = it->second->at(1);
			UINT64 dest2ID = edge2->getDestinationRead();
			//Check edge destinations are part of the marked set.
			if(markedNodes.find(source1ID) != markedNodes.end() && markedNodes.find(dest2ID) != markedNodes.end())
			{
				// One incoming edge and one outgoing edge.
				// And do not merge if either of the edges is a loop
				if( is_mergeable(edge1, edge2) && !(edge1->isLoop()) && !(edge2->isLoop()) )
				{
					EdgeSimple *new_edge = Add(edge1, edge2);
					insertParEdge(new_edge, parGraph);
					if(edge2 != edge1->getReverseEdge()){
						removeParEdge(edge2, parGraph);
					}
					removeParEdge(edge1, parGraph);
					++counter;	// Counter how many edges merged.
				}
			}
		}
	}
	if(counter > 0){
		FILE_LOG(logINFO) << setw(10) << counter << " composite Edges merged." << "\n";
	}
	CLOCKSTOP;
	return counter;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  contractCompositeEdges
 *  Description:  Contract composite paths in the partial overlap graph, changing only marked reads in parallel.
 *                u*-------->v>---------*w  => u*--------------------*w
 *                u*--------<v<---------*w  => u*--------------------*w
 * =====================================================================================
 */
UINT64 OverlapGraphSimple::contractParCompositeEdges(map<UINT64, t_edge_vec* > *parGraph, set<UINT64> &markedNodes)
{
	CLOCKSTART;
	UINT64 counter(0), delEdges(0), nodeCount = parGraph->size();
	map<UINT64, size_t> readIdMap;
	UINT64 indx=0;
	for (map<UINT64, vector<EdgeSimple*> * >::iterator it=parGraph->begin(); it!=parGraph->end();it++)
		readIdMap.insert(make_pair(it->first,indx++));
	vector<EdgeSimple*> addEdgeList;									//Track composite edges added
	bool *allMarked = new bool[nodeCount];			//Track nodes already processed edges added
	std::fill_n(allMarked, nodeCount, false);

	#pragma omp parallel num_threads(p_ThreadPoolSize)
	{
		UINT64 startReadID=0,prevReadID=0;
		UINT64 startIndx=0;
		#pragma omp critical(assignStart)    //Set initial start points...
		{
			for (map<UINT64, vector<EdgeSimple*> * >::iterator it=parGraph->begin(); it!=parGraph->end();it++)
			{
				if(allMarked[startIndx]==false)
				{
					startReadID=prevReadID=it->first;
					allMarked[startIndx]=true;
					break;
				}
				startIndx++;
			}
		}
		while(startReadID!=0) // Loop till all nodes marked
		{
			if(parGraph->at(startReadID)->size() == 2 && markedNodes.find(startReadID) != markedNodes.end()) // Check if the node has only two edges.
			{
				// First edge, going into Read
				EdgeSimple *edge1 = parGraph->at(startReadID)->at(0)->getReverseEdge();
				// Second edge, going out from Read
				EdgeSimple *edge2 = parGraph->at(startReadID)->at(1);
				// One incoming edge and one outgoing edge.
				// And do not merge if either of the edges is a loop
				//All three nodes must be marked...
				if(is_mergeable(edge1, edge2) && !(edge1->isLoop()) && !(edge2->isLoop())
						&& markedNodes.find(edge1->getSourceRead()) != markedNodes.end()
						&& markedNodes.find(edge2->getDestinationRead()) != markedNodes.end())
				{
					//Invalidate edges to be combined
					edge1->setInvalid();
					edge1->getReverseEdge()->setInvalid();
					edge2->setInvalid();
					edge2->getReverseEdge()->setInvalid();

					//Forward search till we can...
					bool isExtendedable=true;
					vector<UINT64> visitedNodes;
					EdgeSimple *currFwdEdge = new EdgeSimple(*edge2);
					visitedNodes.push_back(edge2->getSourceRead());
					while(isExtendedable)
					{
						UINT64 nextFwdRead = currFwdEdge->getDestinationRead();
						if(parGraph->at(nextFwdRead)->size() == 2) // Check if the node has only two edges.
						{
							EdgeSimple *nextFwdEdge = parGraph->at(nextFwdRead)->at(1);
							UINT64 nextFwdEdgeDest = nextFwdEdge->getDestinationRead();
							if(find(visitedNodes.begin(),visitedNodes.end(),nextFwdEdgeDest) != visitedNodes.end())  //check if this node is already visited or not.
								nextFwdEdge = parGraph->at(nextFwdRead)->at(0);
							nextFwdEdgeDest = nextFwdEdge->getDestinationRead();
							if(find(visitedNodes.begin(),visitedNodes.end(),nextFwdEdgeDest) != visitedNodes.end())   //both destinations visited, forms a loop, so break
								break;
							if(is_mergeable(currFwdEdge, nextFwdEdge) && !(nextFwdEdge->isLoop()) && !(currFwdEdge->isLoop())
									&& markedNodes.find(nextFwdEdgeDest)!= markedNodes.end())
							{
								indx = readIdMap.find(nextFwdRead)->second;
								allMarked[indx]=true;
								visitedNodes.push_back(nextFwdRead);
								//Invalidate edge to be combined
								nextFwdEdge->setInvalid();
								nextFwdEdge->getReverseEdge()->setInvalid();
								EdgeSimple *temp_edge = Add(currFwdEdge, nextFwdEdge);
								//Delete previous composite path created
								delete currFwdEdge->getReverseEdge();
								delete currFwdEdge;
								currFwdEdge = temp_edge;
							}
							else
								isExtendedable=false;
						}
						else
							isExtendedable=false;
					}
					//Reverse search till we can
					isExtendedable=true;
					EdgeSimple *currRevEdge = new EdgeSimple(*edge1);
					while(isExtendedable)
					{
						UINT64 nextRevRead = currRevEdge->getSourceRead();
						if(parGraph->at(nextRevRead)->size() == 2) // Check if the node has only two edges.
						{
							EdgeSimple *nextRevEdge = parGraph->at(nextRevRead)->at(0)->getReverseEdge();
							UINT64 nextRevEdgeSrc = nextRevEdge->getSourceRead();
							if(find(visitedNodes.begin(),visitedNodes.end(),nextRevEdgeSrc) != visitedNodes.end())   //check if this node is already visited or not.
								nextRevEdge = parGraph->at(nextRevRead)->at(1)->getReverseEdge();
							nextRevEdgeSrc = nextRevEdge->getSourceRead();
							if(find(visitedNodes.begin(),visitedNodes.end(),nextRevEdgeSrc) != visitedNodes.end())    //both destinations visited, forms a loop, so break
								break;
							if(is_mergeable(nextRevEdge, currRevEdge) && !(nextRevEdge->isLoop()) && !(currRevEdge->isLoop())
									&& markedNodes.find(nextRevEdgeSrc)!= markedNodes.end())
							{
								indx = readIdMap.find(nextRevRead)->second;
								allMarked[indx]=true;
								visitedNodes.push_back(nextRevRead);
								nextRevEdge->setInvalid();
								nextRevEdge->getReverseEdge()->setInvalid();
								EdgeSimple *temp_edge = Add(nextRevEdge, currRevEdge);
								//Delete previous temporary composite path created
								delete currRevEdge->getReverseEdge();
								delete currRevEdge;
								currRevEdge = temp_edge;
							}
							else
								isExtendedable=false;
						}
						else
							isExtendedable=false;
					}
					//Create complete composite path
					EdgeSimple *new_edge = Add(currRevEdge, currFwdEdge);
					//Delete composite sub-paths
					delete currRevEdge->getReverseEdge();
					delete currRevEdge;
					delete currFwdEdge->getReverseEdge();
					delete currFwdEdge;
					#pragma omp critical(addNewEdge)
					{
						addEdgeList.push_back(new_edge);
					}
				}
			}
			//search next node to operate on
			startReadID=0;
			for (map<UINT64, vector<EdgeSimple*> * >::iterator it=parGraph->find(prevReadID); it!=parGraph->end() && startIndx < nodeCount;it++)
			{
				if(allMarked[startIndx]==false)
				{
					startReadID=prevReadID=it->first;
					allMarked[startIndx]=true;
					break;
				}
				startIndx++;
			}
		}//End of while
	}//end of multi-threading
	delete allMarked;
	//Delete all edges marked as invalid
	for (map<UINT64, vector<EdgeSimple*> * >::iterator it=parGraph->begin(); it!=parGraph->end();it++)
	{
		UINT64 j=0;
		while(j < it->second->size())
		{
			EdgeSimple *e = it->second->at(j);
			if(e->isInvalid())
			{
				removeParEdge(e, parGraph);
				delEdges++;
			}
			else
				j++;
		}
	}
	//Add all the combined edges
	for(size_t i=0;i<addEdgeList.size();i++)
	{
		if(!existsParEdge(addEdgeList[i], parGraph))
		{
			insertParEdge(addEdgeList[i], parGraph);
			counter++;
		}
		else
		{
			delete addEdgeList[i]->getReverseEdge();
			delete addEdgeList[i];
		}
	}
	FILE_LOG(logINFO) << setw(10) << delEdges << " edges deleted." << "\n";
	FILE_LOG(logINFO) << setw(10) << addEdgeList.size() << " composite edges created." << "\n";
	FILE_LOG(logINFO) << setw(10) << counter << " composite Edges added." << "\n";
	FILE_LOG(logINFO) << setw(10) << delEdges-counter << " composite Edges merged." << "\n";
	CLOCKSTOP;
	return counter;
}

//Checks if an edge already exists in partial graph
bool OverlapGraphSimple::existsParEdge(EdgeSimple *checkEdge, map<UINT64, t_edge_vec* > *parGraph)
{
	UINT64 sourceID = checkEdge->getSourceRead();
	map<UINT64, vector<EdgeSimple*> * >::iterator it=parGraph->find(sourceID);
	if(it != parGraph->end())
	{
		for(size_t i=0;i<it->second->size();i++)
		{
			EdgeSimple *e = it->second->at(i);
			if(*e==*checkEdge)
				return true;
		}
	}
	return false;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  removeDeadEndNodes
 *  Description:  remove dead end nodes and their incident edges
 * =====================================================================================
 */
UINT64 OverlapGraphSimple::removeParDeadEndNodes(map<UINT64, t_edge_vec* > *parGraph, set<UINT64> &markedNodes, vector<UINT64> &nodeList)
{
	CLOCKSTART;
	vector<UINT64> nodes_to_remove;
	#pragma omp parallel for schedule(dynamic) num_threads(p_ThreadPoolSize)
	for(size_t i=0; i<nodeList.size(); i++)
	{
		map<UINT64, vector<EdgeSimple*> * >::iterator it=parGraph->find(nodeList[i]);
		if(!it->second->empty())	// If the read has some edges.
		{
			bool isDeadEnd = true;	// flag for dead end edge
			UINT64 inEdge = 0; 	// number of incoming edges to this node
			UINT64 outEdge = 0; 	// number of outgoing edges from this node

			// Find number of in- and out- edges
			for(UINT64 j=0; j < it->second->size(); j++)
			{
				EdgeSimple * edge = it->second->at(j);

				//If edge points to a node not marked in this partition donot assume dead end delete it.
				UINT64 destID = edge->getDestinationRead();
				if(markedNodes.find(destID) == markedNodes.end())
				{
					isDeadEnd = false;
					break;
				}
				/* Break case:
				 * 0. edge already marked as not dead end
				 * 1. composite edge with more than minReadsCountInEdgeToBeNotDeadEnd (deafult: 10)
				 * 2. composite edge longer than minEdgeLengthToBeNotDeadEnd (default: 500)
				 * 3. the edge is loop for the current node
				 * Then flag=1 and exit the loop
				 */

				if (edge->isNotDeadEnd()){
					isDeadEnd = false;
					break;
				}
				if(edge->isListofReads() && edge->getListofReadsSize() >= minReadsCountInEdgeToBeNotDeadEnd) {
					edge->markNotDeadEnd();
					isDeadEnd = false;
					break;
				}
				if(edge->getEdgeLength() >= minEdgeLengthToBeNotDeadEnd) {
					edge->markNotDeadEnd();
					isDeadEnd = false;
					break;
				}
				if(edge->isLoop())
				{
					edge->markNotDeadEnd();
					isDeadEnd = false;
					break;
				}

				if((edge->getOrientation() >> 1) & 1)
					++outEdge;
				else
					++inEdge;
			}
			// no good edges incident to the node and only in-edges or out-edges
			if( isDeadEnd && inEdge*outEdge == 0 && inEdge + outEdge > 0){
				#pragma omp critical(updateNodeList)
				{
					nodes_to_remove.push_back(it->first);
				}
			}
		}

	}
	FILE_LOG(logINFO) << "number of dead end nodes found: " << nodes_to_remove.size() << "\n";

	UINT64 deleted_edges(0);
	// Now delete the edges incident to these nodes
	for(auto it = nodes_to_remove.cbegin(); it != nodes_to_remove.cend(); ++it){
		UINT64 nodeID = *it;
		while(!(parGraph->at(nodeID)->empty()))
		{
			removeParEdge(parGraph->at(nodeID)->front(),parGraph);
			++deleted_edges;
		}
	}
	FILE_LOG(logINFO) << "number of edges deleted: " << deleted_edges << "\n";

	CLOCKSTOP;
	return deleted_edges;
}

OverlapGraphSimple::OverlapGraphSimple(string edge_file, string composite_out_edge_file,
		UINT64 minOvl = 0, UINT64 parallelThreadPoolSize = 0)
	: m_minOvl(minOvl), p_ThreadPoolSize(parallelThreadPoolSize)
{
	CLOCKSTART;
	m_numberOfNodes	= 0;
	m_numberOfEdges	= 0;

	// loop edgeFilenameList
	// Composite edge contraction with remove dead end nodes

	map<UINT64, t_edge_vec* > *parGraph = new map<UINT64, t_edge_vec* >;			//Partial graph
	set<UINT64> markedNodes;
	vector<UINT64> nodeList;
	loadParEdgesFromEdgeFile(edge_file, parGraph, markedNodes);
	sortParEdgesByDestID(parGraph);
	nodeList.assign( markedNodes.begin(), markedNodes.end() );
	// Composite edge contraction with remove dead end nodes from partial graph
	//Do only one round in parallel. Others can be done serial for speed...
	contractParCompositeEdges(parGraph, markedNodes);
	removeParDeadEndNodes(parGraph, markedNodes, nodeList);
	UINT64 counter(0);
	do {
		counter = contractParCompositeEdges_Serial(parGraph, markedNodes);
	//	counter += removeParDeadEndNodes(parGraph, markedNodes, nodeList);
	} while (counter > 0);

	//Write partial simplified graph to file
	//Print all edges file
	printParEdges(composite_out_edge_file, parGraph);
	//Delete partial graph structure
	for (map<UINT64, vector<EdgeSimple*> * >::iterator it=parGraph->begin(); it!=parGraph->end();it++)
	{
		UINT64 readID = it->first, j=0;
		UINT64 edgeVecSize = parGraph->at(readID)->size();
		while(j < edgeVecSize)
		{
			EdgeSimple *e = parGraph->at(readID)->at(j);
			delete e;
			j++;
		}
		delete parGraph->at(readID);
	}
	delete parGraph;

	CLOCKSTOP;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  sortEdgesByLength
 *  Description:  For each node, sort its incident edges by their lengths in increasing order
 * =====================================================================================
 */
void OverlapGraphSimple::sortParEdgesByDestID(map<UINT64, t_edge_vec* > *parGraph)
{

	for (map<UINT64, vector<EdgeSimple*> * >::iterator it=parGraph->begin(); it!=parGraph->end();it++)
	{
		sort(it->second->begin(), it->second->end(), compareParEdgesByDestID);
	}
}

//=============================================================================
// Function to compare two edges by the destination read number. 
// The read number was determined by the read std::string lexicographically.
// Used for sorting.
//=============================================================================
bool compareParEdgesByDestID (const EdgeSimple *edge1, const EdgeSimple* edge2)
{
	if (edge1->getDestinationRead() < edge2->getDestinationRead())
	{
		return true;
	}
	else if (edge1->getDestinationRead() == edge2->getDestinationRead()){
		return (compareParEdgesByLength(edge1, edge2));
	}
	else
		return false;
}


//=============================================================================
// Function to compare two edges by overlap offset. 
// Used for transitive edge removal (which is not included any more)
//=============================================================================
bool compareParEdgesByLength (const EdgeSimple *edge1, const EdgeSimple* edge2)
{
	return (edge1->getEdgeLength() < edge2->getEdgeLength());
}



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  compareEdgesByReads
 *  Description:  Compare edges by number of reads contained in the them
 * =====================================================================================
 */
bool compareEdgesByReads (const EdgeSimple *edge1, const EdgeSimple* edge2)
{
	UINT64 num_reads1 = (edge1->isListofReads() ? edge1->getListofReadsSize() : 0);
	UINT64 num_reads2 = (edge2->isListofReads() ? edge2->getListofReadsSize() : 0);
	return (num_reads1 > num_reads2);
}

void OverlapGraphSimple::loadParEdgesFromEdgeFile(const std::string &edgeFilename, map<UINT64, t_edge_vec* > *parGraph,set<UINT64> &markedNodes)
{
	CLOCKSTART;
	FILE_LOG(logINFO) << "Load edge file: " << edgeFilename << "\n";

	// Open file
	ifstream filePointer;
	filePointer.open(edgeFilename.c_str());
	if(!filePointer.is_open() )
		MYEXIT("Unable to open file: "+edgeFilename);

	UINT64 edgeNumber(0);
	// Read file
	std::string line_text;
	while(getline(filePointer,line_text)) {
		// FILE_LOG(logDEBUG1) << "EdgeSimple " << line_text << "\n";
		// 15	3610522	1,169,0,0,370,201,369,414,0,168,NA
		vector<std::string> line_elements;
		std::stringstream text_ss(line_text);
		std::string element;

		// get the three parts separated by tabs in each line
		while (getline(text_ss, element, '\t')) {
			line_elements.push_back(element);
		}

		// Get source readID in the edge file (corresponding to the ID in reads file),
		// then use the readIDMap to find the ID in the graph
		UINT64 source;
		istringstream t1(line_elements.at(0)); 
		t1 >> source;
		// Do the same for the destination ID
		UINT64 destination;	// Tab delimited Column2: destination read ID
		istringstream t2(line_elements.at(1)); t2 >> destination;
		// Find destination readID from readIDMap
		// Properties
		std::string properties = line_elements.at(2);	// Tab delimited Column3: edge attributes

		// For properties list
		vector<std::string> propertiesList;
		std::stringstream properties_ss(properties);
		while (getline(properties_ss, element, ',')) {
			propertiesList.push_back(element);
		}

		// 15	3610522	1,169,0,0,370,201,369,414,0,168,NA
		// properties
		UINT8 orientation;	// Property Col1: orientation
		int dec_orient;
		std::istringstream p1(propertiesList.at(0)); p1 >> dec_orient;
		orientation = static_cast<UINT8>(dec_orient);

		UINT32 overlapLength;	// Property Col2: overlap length
		std::istringstream p2(propertiesList.at(1)); p2 >> overlapLength;

		UINT64 substitutions;	// Property Col3: substitutions
		std::istringstream p3(propertiesList.at(2)); p3 >> substitutions;

		UINT64 edits;	// Property Col4: edit distance
		std::istringstream p4(propertiesList.at(3)); p4 >> edits;

		// If the edge has overlap that satisfies the requirements
		if (overlapLength >= m_minOvl){
			UINT32 length1;	// Property Col5: read1 length
			std::istringstream p5(propertiesList.at(4)); p5 >> length1;

			UINT32 start1;	// Property Col6: read1 overlap start
			std::istringstream p6(propertiesList.at(5)); p6 >> start1;

//			UINT32 stop1;	// Property Col7: read1 overlap end
//			std::istringstream p7(propertiesList.at(6)); p7 >> stop1;

			UINT32 length2;	// Property Col8: read2 length
			std::istringstream p8(propertiesList.at(7)); p8 >> length2;

//			UINT32 start2;	// Property Col9: read2 overlap start
//			std::istringstream p9(propertiesList.at(8)); p9 >> start2;
//
//			UINT32 stop2;	// Property Col10: read2 overlap length
//			std::istringstream p10(propertiesList.at(9)); p10 >> stop2;

			//Denotes if source/destination is marked in this partition of the graph...
			// 0: Only source is marked
			// 1: Only destination is marked
			// 2: Both source and destination are marked
			// Only relevant if you have a graph partitioned
			UINT32 markFlag=2;
			if(propertiesList.size()>11)
			{
				std::istringstream p11(propertiesList.at(11)); p11 >> markFlag;
			}
			/*  get overlap offset
			 *  Example edge list:
			 *  0 = u<-----------<v		reverse of u to reverse of v
			 *  1 = u<----------->v		reverse of u to forward of v
			 *  2 = u>-----------<v		forward of u to reverse of v
			 *  3 = u>----------->v		forward of u to forward of v
			 * 1	2   2,34,0,0,35,1,34,35,0,33,NA
			 * 2	3   0,33,0,0,35,2,34,35,0,32,NA
			 * 3	4   3,34,1,1,35,1,34,35,0,33,1CA
			 */

			// overlap offset
			UINT32 overlapOffset = start1;	// correct, but not ready yet JJ: why not ready yet?
			EdgeSimple *an_edge = new EdgeSimple(source,length1,destination,length2,orientation,overlapOffset);
			an_edge->make_nonComposite_reverseEdge();
			insertParEdge(an_edge, parGraph);	// insert edge

			//Mark node. This node is candidate for local reduction as it was source in the partial graph and was marked in this subgraph
			if(markFlag==0)
				markedNodes.insert(source);
			else if(markFlag==1)
				markedNodes.insert(destination);
			else
			{
				markedNodes.insert(source);
				markedNodes.insert(destination);
			}
			++edgeNumber;
			if(edgeNumber % 1000000 == 0){
				FILE_LOG(logINFO) << setw(10) << (edgeNumber / 1000000) << ",000,000"  
					<< " edges loaded to memory, "
					<< setw(7) << checkMemoryUsage() << " MB" << "\n";
			}
		}
	}
	filePointer.close();
	FILE_LOG(logINFO) << setw(10) << edgeNumber << " edges loaded to memory, "<< "\n";
	CLOCKSTOP;
}

void OverlapGraphSimple::printEdge(EdgeSimple *contigEdge, ostream & filePointer) const
{
	UINT64 source = contigEdge->getSourceRead();
	UINT64 destination = contigEdge->getDestinationRead();
	UINT64 orientation=contigEdge->getOrientation();
	if(source < destination || (source == destination && contigEdge < contigEdge->getReverseEdge()))
	{
		UINT64 offsetSum = 0;
		if(contigEdge->isListofReads())
		{
			for(size_t j=0;j<contigEdge->getListofReadsSize();j++)
				offsetSum+=contigEdge->getInnerOverlapOffset(j);
		}
		filePointer<<source<<"\t";	// store the edge information first
		filePointer<<destination<<"\t";
		filePointer<<orientation<<",";
		filePointer<<contigEdge->getOverlapOffset()<<",";  //total overlap offset
		filePointer<<contigEdge->getEdgeLength()<<",";  //edge length
		filePointer<<0<<",";					//no substitutions
		filePointer<<0<<"\t";					//no edits
		if(contigEdge->isListofReads())
		{
			for(size_t j=0;j<contigEdge->getListofReadsSize();j++)
			{
				orientation=contigEdge->getInnerOrientation(j);
				filePointer<<"("<<contigEdge->getInnerReadID(j)<<","<<orientation<<","
						<<contigEdge->getInnerOverlapOffset(j)<<")";

			}
		}
		filePointer<<endl;
	}
}


/*
 * ===  FUNCTION  ======================================================================
 *         Name:  printContigs
 *  Description:  Print all the edges in the graph
 * =====================================================================================
 */
void OverlapGraphSimple::printParEdges(string edge_file, map<UINT64, t_edge_vec* > *parGraph) const
{
	CLOCKSTART;
	ofstream filePointer;
	filePointer.open(edge_file.c_str());

	for (map<UINT64, vector<EdgeSimple*> * >::iterator it=parGraph->begin(); it!=parGraph->end();it++)
	{
		UINT64 readID = it->first, j=0;
		UINT64 edgeVecSize = parGraph->at(readID)->size();
		while(j < edgeVecSize)
		{
			EdgeSimple *e = parGraph->at(readID)->at(j);
			if(e->isSmallerEdge()){
				printEdge(e,filePointer);
			}
			j++;
		}
	}
	filePointer.close();
	CLOCKSTOP;
}


/**********************************************************************************************************************
	Orientation of a reverse edge;
	Twin edge of Orientation 0 = <-------< is Orientation 3 = >------->
	Twin edge of Orientation 1 = <-------> is Orientation 1 = <------->
	Twin edge of Orientation 2 = >-------< is Orientation 2 = >-------<
	Twin edge of Orientation 3 = >-------> is Orientation 0 = <-------<
**********************************************************************************************************************/
UINT8 OverlapGraphSimple::twinEdgeOrientation(UINT8 orientation)
{
	UINT8 returnValue;
	if(orientation == 0)
		returnValue = 3;
	else if(orientation == 1)
		returnValue = 1;
	else if(orientation == 2)
		returnValue = 2;
	else if(orientation == 3)
		returnValue = 0;
	else
		MYEXIT("Unsupported edge orientation.")
	return returnValue;
}




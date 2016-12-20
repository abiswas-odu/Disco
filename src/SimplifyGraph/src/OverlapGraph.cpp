/*
 * ===== CLASS IMPLEMENTATION ================================================
 * Name        : OverlapGraph.cpp
 * Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey, Abhishek Biswas
 * Version     : v1.2
 * Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
 * Description : OverlapGraph C++ file
 *============================================================================
 */

#include "OverlapGraph.h"
#include "CS2_stream/cs2.h"
#ifdef INCLUDE_READGZ
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)
#endif
extern TLogLevel loglevel;                      /* verbosity level of logging */
extern std::string outputFilenamePrefix;

/**********************************************************************************************************************
	Check if two edges match.
	e1(u,v) and e2(v,w). At node v, one of the edges should be an incoming edge and the other should be an outgoing
	edge to match.
**********************************************************************************************************************/

bool matchEdgeType(const Edge *edge1, const Edge *edge2)
{
	if     ( (edge1->getOrientation() == 1 || edge1->getOrientation() == 3) && (edge2->getOrientation() == 2 || edge2->getOrientation() == 3) ) // *-----> and >------*
		return true;
	else if( (edge1->getOrientation() == 0 || edge1->getOrientation() == 2) && (edge2->getOrientation() == 0 || edge2->getOrientation() == 1) ) // *------< and <------*
		return true;
	return false;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  getMismatchFromString
 *  Description:  From the mismatch string in the edge file, get the mismatch into t_vpair
 * =====================================================================================
 */
t_vpair* OverlapGraph::getMismatchFromString(const std::string &mismatch_string)
{
	if (mismatch_string == "NA")
		return nullptr;
	else{
	// mismatches to list
		// loop mismatches
		t_vpair *mismatches = new t_vpair;
		std::stringstream mismatches_ss(mismatch_string);
		// JJ: The source and destination characters are not used?ff
		// Can I remove them, and only get the position, like below?
		std::string element;
		getline(mismatches_ss, element, ':');
		std::istringstream element_ss(element);
		UINT32 mismatchPosition;
		element_ss >> mismatchPosition;
		mismatches->push_back(make_pair(0, mismatchPosition));
		return mismatches;
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  insertEdge
 *  Description:  Insert an edge in the global overlap graph.
 *                Does not automatically insert its twin edge.
 * =====================================================================================
 */
void OverlapGraph::insertEdge( Edge *edge)
{
	insertFwdEdge(edge);
	insertFwdEdge(edge->getReverseEdge());
}

void OverlapGraph::insertFwdEdge( Edge *edge)
{
	UINT64 ID = edge->getSourceRead()->getReadID();
	if(m_graph->find(ID) == m_graph->end())            /* A new node in the graph */
	{
		t_edge_vec *newList = new t_edge_vec;
		m_graph->insert( std::pair<UINT64, t_edge_vec* >(ID, newList));
	}
	if(m_graph->at(ID)->empty())
		m_numberOfNodes++;
	m_graph->at(ID)->push_back(edge);
	++m_numberOfEdges;
	updateReadsLocations(edge, INSERTION, m_dataset);
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  removeEdgeFromSourceRead
 *  Description:  Remove an edge from the edge list of the edge's source read, 
 *  		  but do not delete this edge from memory yet.
 * =====================================================================================
 */
void OverlapGraph::removeEdgeFromSourceRead(Edge *edge)
{
	if (!edge)
		return;
	// all edges incident to source read
	t_edge_vec *fwd_edges = m_graph->at(edge->getSourceRead()->getReadID());
	fwd_edges->erase(std::remove(fwd_edges->begin(), fwd_edges->end(), edge), fwd_edges->end());
	if (fwd_edges->empty())
	{
		--m_numberOfNodes;
		//Very carefully consider deleting elements from the map. This will mess-up for loop iterators...
		//delete fwd_edges;
		//m_graph->erase(edge->getSourceRead()->getReadID());
	}
	--m_numberOfEdges;
}
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  removeEdge
 *  Description:  Remove an edge from the overlap graph, release its memory if we don't
 *  		  need to save this edge for the contigs.
 *  		  This removes and deletes from memory both the edge and its twin edge
 * =====================================================================================
 */
void OverlapGraph::removeEdge(Edge *edge){
	if(edge == nullptr)
		return;
	removeFwdEdge(edge->getReverseEdge());
	removeFwdEdge(edge);
	delete edge->getReverseEdge();
	delete edge;

}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  removeEdge
 *  Description:  Remove an edge from the overlap graph, release its memory if we don't 
 *  		  need to save this edge for the contigs.
 *  		  This removes and deletes from memory both the edge and its twin edge
 * =====================================================================================
 */
void OverlapGraph::removeFwdEdge(Edge *edge)
{
	// If edge points to NULL, there is nothing to to
	if(edge == nullptr) {
		return;
	}
	// If the current edge contains some reads. We have to update their location formation.
	updateReadsLocations(edge, DELETION, m_dataset);

	removeEdgeFromSourceRead(edge);

}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  contractCompositeEdges
 *  Description:  Contract composite paths in the overlap graph.
 *                u*-------->v>---------*w  => u*--------------------*w
 *                u*--------<v<---------*w  => u*--------------------*w
 * =====================================================================================
 */
UINT64 OverlapGraph::contractCompositeEdgesPar(void)
{
	CLOCKSTART;
	UINT64 counter(0), delEdges(0), nodeCount = m_graph->size();
	map<UINT64, size_t> readIdMap;
	UINT64 indx=0;
	for (map<UINT64, t_edge_vec* >::iterator it=m_graph->begin(); it!=m_graph->end();it++)
		readIdMap.insert(make_pair(it->first,indx++));
	vector<Edge*> addEdgeList;									//Track composite edges added
	bool *allMarked = new bool[nodeCount];			//Track nodes already processed edges added
	std::fill_n(allMarked, nodeCount, false);

	double beginPSec = omp_get_wtime();
	FILE_LOG(logINFO)<<"Composite edge parallel section start: Threads "<<p_ThreadPoolSize<<endl;
	#pragma omp parallel num_threads(4)
	{
		UINT64 startReadID=0,prevReadID=0;
		UINT64 startIndx=0;
		#pragma omp critical(assignStart)    //Set initial start points...
		{
			for (map<UINT64, t_edge_vec* >::iterator it=m_graph->begin(); it!=m_graph->end();it++)
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
			t_edge_vec *startEdgeList = m_graph->at(startReadID);
			if(startEdgeList->size() == 2) // Check if the node has only two edges.
			{
				// First edge, going into Read
				Edge *edge1 = startEdgeList->at(0)->getReverseEdge();
				// Second edge, going out from Read
				Edge *edge2 = startEdgeList->at(1);
				// One incoming edge and one outgoing edge.
				// And do not merge if either of the edges is a loop
				if( is_mergeable(edge1, edge2) && !(edge1->isLoop()) && !(edge2->isLoop()) )
				{
					//Invalidate edges to be combined
					edge1->setInvalid();
					edge1->getReverseEdge()->setInvalid();
					edge2->setInvalid();
					edge2->getReverseEdge()->setInvalid();

					//Forward search till we can...
					bool isExtendedable=true;
					vector<UINT64> visitedNodes;
					visitedNodes.push_back(edge2->getSourceRead()->getReadID());
					Edge *currFwdEdge = new Edge(*edge2);
					while(isExtendedable)
					{
						UINT64 nextFwdRead = currFwdEdge->getDestinationRead()->getReadID();
						t_edge_vec *edgeList = m_graph->at(nextFwdRead);
						if(edgeList->size() == 2) // Check if the node has only two edges.
						{
							Edge *nextFwdEdge = edgeList->at(1);
							UINT64 nextFwdEdgeDest = nextFwdEdge->getDestinationRead()->getReadID();
							if(find(visitedNodes.begin(),visitedNodes.end(),nextFwdEdgeDest) != visitedNodes.end())  //check if this node is already visited or not.
								nextFwdEdge = edgeList->at(0);
							nextFwdEdgeDest = nextFwdEdge->getDestinationRead()->getReadID();
							if(find(visitedNodes.begin(),visitedNodes.end(),nextFwdEdgeDest) != visitedNodes.end())   //both destinations visited, forms a loop, so break
								break;
							if(is_mergeable(currFwdEdge, nextFwdEdge) && !(nextFwdEdge->isLoop()) && !(currFwdEdge->isLoop()) )
							{
								indx = readIdMap.find(nextFwdRead)->second;
								allMarked[indx]=true;
								visitedNodes.push_back(nextFwdRead);
								//Invalidate edge to be combined
								nextFwdEdge->setInvalid();
								nextFwdEdge->getReverseEdge()->setInvalid();
								//currFwdEdge->updateEdge(nextFwdEdge);
								Edge *temp_edge = Add(currFwdEdge, nextFwdEdge);
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
					Edge *currRevEdge = new Edge(*edge1);
					while(isExtendedable)
					{
						UINT64 nextRevRead = currRevEdge->getSourceRead()->getReadID();
						t_edge_vec *edgeList = m_graph->at(nextRevRead);
						if(edgeList->size() == 2) // Check if the node has only two edges.
						{
							Edge *nextRevEdge = edgeList->at(0)->getReverseEdge();
							UINT64 nextRevEdgeSrc = nextRevEdge->getSourceRead()->getReadID();
							if(find(visitedNodes.begin(),visitedNodes.end(),nextRevEdgeSrc) != visitedNodes.end())   //check if this node is already visited or not.
								nextRevEdge = edgeList->at(1)->getReverseEdge();
							nextRevEdgeSrc = nextRevEdge->getSourceRead()->getReadID();
							if(find(visitedNodes.begin(),visitedNodes.end(),nextRevEdgeSrc) != visitedNodes.end())    //both destinations visited, forms a loop, so break
								break;
							if(is_mergeable(nextRevEdge, currRevEdge) && !(nextRevEdge->isLoop()) && !(currRevEdge->isLoop()) )
							{
								indx = readIdMap.find(nextRevRead)->second;
								allMarked[indx]=true;
								visitedNodes.push_back(nextRevRead);
								nextRevEdge->setInvalid();
								nextRevEdge->getReverseEdge()->setInvalid();
								//currRevEdge->getReverseEdge()->updateEdge(nextRevEdge->getReverseEdge());
								Edge *temp_edge = Add(nextRevEdge, currRevEdge);
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
					Edge *new_edge = Add(currRevEdge, currFwdEdge);
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
			for (map<UINT64, t_edge_vec* >::iterator it=m_graph->find(prevReadID); it!=m_graph->end() && startIndx < nodeCount;it++)
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
	double endPSec = omp_get_wtime();
	FILE_LOG(logINFO)<<"Composite edge parallel section end: Time "<<double(endPSec - beginPSec)<<'\n';
	//Delete all edges marked as invalid
	for (map<UINT64, t_edge_vec* >::iterator it=m_graph->begin(); it!=m_graph->end();it++)
	{
		UINT64 i = 0;
		while(i<it->second->size())
		{
			if(it->second->at(i)->isInvalid())
			{
				removeEdge(it->second->at(i));
				delEdges++;
			}
			else
				i++;
		}
	}
	//Add all the combined edges
	for(size_t i=0;i<addEdgeList.size();i++)
	{
		if(!existsEdge(addEdgeList[i]))
		{
			insertEdge(addEdgeList[i]);
			counter++;
		}
		else
		{
			delete addEdgeList[i]->getReverseEdge();
			delete addEdgeList[i];
		}
	}
	//Delete nodes in the graph that have no edges anymore
	map<UINT64, t_edge_vec* >::iterator it=m_graph->begin();
	while(it!=m_graph->end())
	{
		if(it->second->empty())
		{
			delete it->second;
			m_graph->erase(it);
		}
		else
			it++;
	}
	FILE_LOG(logINFO) << setw(10) << addEdgeList.size() << " composite extra edges created." << "\n";
	FILE_LOG(logINFO) << setw(10) << delEdges << " edges deleted." << "\n";
	if(counter > 0){
		FILE_LOG(logINFO) << setw(10) << counter << " composite Edges merged." << "\n";
	}
	CLOCKSTOP;
	FILE_LOG(logINFO) << "numberOfEdges = " << m_numberOfEdges << "\n";
	return counter;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  contractCompositeEdges
 *  Description:  Contract composite paths in the overlap graph.
 *                u*-------->v>---------*w  => u*--------------------*w
 *                u*--------<v<---------*w  => u*--------------------*w
 * =====================================================================================
 */
UINT64 OverlapGraph::contractCompositeEdges(void)
{
	CLOCKSTART;
	UINT64 counter(0);
	for(map<UINT64, t_edge_vec* >::iterator it=m_graph->begin(); it!=m_graph->end();it++)
	{
		if(it->second->size() == 2) // Check if the node has only two edges.
		{
			// First edge, going into Read
			Edge *edge1 =it->second->at(0)->getReverseEdge();
			// Second edge, going out from Read
			Edge *edge2 = it->second->at(1);
			// One incoming edge and one outgoing edge.
			// And do not merge if either of the edges is a loop
			if( is_mergeable(edge1, edge2) && !(edge1->isLoop()) && !(edge2->isLoop()) )
			{
				Edge *new_edge = Add(edge1, edge2);
				insertEdge(new_edge);
				if(edge2 != edge1->getReverseEdge()){
					removeEdge(edge2);
				}
				removeEdge(edge1);
				++counter;	// Counter how many edges merged.
			}
		}

	}
	if(counter > 0){
		FILE_LOG(logINFO) << setw(10) << counter << " composite Edges merged." << "\n";
	}
	CLOCKSTOP;
	FILE_LOG(logINFO) << "numberOfEdges = " << m_numberOfEdges << "\n";
	return counter;
}

//Checks if an edge already exists
bool OverlapGraph::existsEdge(Edge *checkEdge)
{
	UINT64 sourceID = checkEdge->getSourceRead()->getReadID();
	auto it = m_graph->find(sourceID);
	if(it != m_graph->end())
	{
		t_edge_vec *edgeList = it->second;
		for(size_t i=0;i<edgeList->size();i++)
		{
			Edge *e = edgeList->at(i);
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
UINT64 OverlapGraph::removeDeadEndNodes(void)
{
	CLOCKSTART;
	vector< vector<UINT64>* > nodes_to_remove;

	#pragma omp parallel num_threads(p_ThreadPoolSize)
	{
		vector<UINT64> *nodes_to_remove_local = new vector<UINT64>;
		#pragma omp for schedule(dynamic)
		for(UINT64 i = 1; i <= m_dataset->size() ; i++) // For each read.
		{
			auto it = m_graph->find(i);
			if(it != m_graph->end() && !it->second->empty())	// If the read has some edges.
			{
				bool isDeadEnd = true;	// flag for dead end edge
				UINT64 inEdge = 0; 	// number of incoming edges to this node
				UINT64 outEdge = 0; 	// number of outgoing edges from this node

				// Find number of in- ane out- edges
				for(UINT64 j=0; j < it->second->size(); j++)
				{
					Edge * edge = it->second->at(j);
					/* Break case:
					 * 0. edge already marked as not dead end
					 * 1. composite edge with more than minReadsCountInEdgeToBeNotDeadEnd (default: 10)
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
						nodes_to_remove_local->push_back(it->first);
				}
			}
		}
		#pragma omp critical
		{
			nodes_to_remove.push_back(nodes_to_remove_local);
		}
	}
	UINT64 deleted_edges(0), dead_nodes(0);
	// Now delete the edges incident to these nodes
	for(auto itOuter = nodes_to_remove.cbegin(); itOuter != nodes_to_remove.cend(); ++itOuter){
		for(auto it = (*itOuter)->cbegin(); it != (*itOuter)->cend(); ++it){
			UINT64 nodeID = *it;
			t_edge_vec* eList = m_graph->at(nodeID);
			dead_nodes++;
			while(!(eList->empty()))
			{
				removeEdge(eList->front());
				++deleted_edges;
			}
		}
	}
	FILE_LOG(logDEBUG) << "number of dead end nodes found: " << dead_nodes << "\n";
	//Delete heal allocated vector for each thread
	for(auto itOuter = nodes_to_remove.cbegin(); itOuter != nodes_to_remove.cend(); ++itOuter){
		delete (*itOuter);
	}

	//Delete nodes in the graph that have no edges anymore
	map<UINT64, t_edge_vec* >::iterator it=m_graph->begin();
	while(it!=m_graph->end())
	{
		if(it->second->empty())
		{
			delete it->second;
			m_graph->erase(it);
		}
		else
			it++;
	}

	FILE_LOG(logDEBUG) << "number of edges deleted: " << deleted_edges << "\n";

	CLOCKSTOP;
	FILE_LOG(logINFO) << "numberOfEdges = " << m_numberOfEdges << "\n";
	return deleted_edges;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  removeShortBranches
 *  Description:  After flow analysis, at late stage of graph simplification.
 *  Now only look at nodes with 1 incident edge, if this edge is much shorter comparing
 *  to other edges incident to its neighbor, then remove this relatively short branch.
 * =====================================================================================
 */
UINT64 OverlapGraph::removeShortBranches(void)
{
	if(!(this->m_flowComputed)) {
		return 0;
	}
	CLOCKSTART;
	UINT64 num_nodes_rm(0);

	typedef vector<UINT64> LongBrLens;
	map<UINT64, LongBrLens> long_brlens_map; /* map from read number to length vector of length 2 */
	for(map<UINT64, t_edge_vec* >::iterator it=m_graph->begin(); it!=m_graph->end();it++)  // For each read.
	{
		/* If this read has exactly 1 incident edge, and this edge is not a loop */
		if(it->second->size()==1 && !(it->second->front()->isLoop())){

			/* The incident edge, going out from its only neighbor */
			Edge* one_edge =it->second->at(0)->getReverseEdge();

			UINT64 neighbor = one_edge->getSourceRead()->getReadID();
			UINT64 neighbor_degree = m_graph->at(neighbor)->size();

			/* If neighbor has more edges besides the "short branch" */
			if(neighbor_degree > 1){ 
				UINT64 one_length = one_edge->getOverlapOffset();
				// edge is going in (0) or out (1) from the read
				UINT8 in_out = ((one_edge->getOrientation() >> 1) & 1 );	
/* 				FILE_LOG(logDEBUG1) << "Short branch from " << neighbor << " to " << i 
 * 					<< " with orientation " << in_out << " and length " << one_length << "\n";
 */
				// If the longest branch lengths at this neighor hasn't been found yet,
				// look for it
				if(long_brlens_map.count(neighbor)==0){
/* 					FILE_LOG(logDEBUG1) << "Now look for the longest branches at node " 
 * 						<< neighbor << "\n";
 */
					LongBrLens long_brlens;
					// Initialize both longest length in and out to be 0
					long_brlens.push_back(0);
					long_brlens.push_back(0);
					long_brlens.at(in_out) = one_length;
					/* Find the longest branches at neighbor in both directions */
					for(UINT64 j = 0; j < neighbor_degree; ++j){
						Edge* e = m_graph->at(neighbor)->at(j);
						UINT8 direction = (e->getOrientation() >> 1) & 1; 
						if(e->getOverlapOffset() > long_brlens.at(direction))
							long_brlens.at(direction) = e->getOverlapOffset();

					}
					long_brlens_map[neighbor] = long_brlens;
				}
				// overlap is small 
				// edge is not already in the deleting list
				// edge is smaller than its reverse edge
				if(one_length * minFoldToBeShortBranch < long_brlens_map.at(neighbor).at(in_out) && one_length < minSizeToBeShortBranch){
					removeEdge(one_edge);
					++num_nodes_rm;
					FILE_LOG(logDEBUG1) << "Delete this edge, length: " << one_length << " and " << long_brlens_map.at(neighbor).at(in_out) << "\n";
				}
			}
		}
	}

	if(num_nodes_rm > 0){
		FILE_LOG(logINFO) << "short-branch nodes removed: " << num_nodes_rm << "\n";
	}
	CLOCKSTOP;
	FILE_LOG(logINFO) << "numberOfEdges = " << m_numberOfEdges << "\n";
	return num_nodes_rm;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  removeSimilarEdges
 *  Description:  Remove similar length edges between same pair of nodes
 * =====================================================================================
 */
UINT64 OverlapGraph::removeSimilarEdges(void)
{
	CLOCKSTART;
	UINT64 counter = 0;
	// pairs of edges that are found, the two lists have the same length, a pair are in the same index

	for(map<UINT64, t_edge_vec* >::iterator it=m_graph->begin(); it!=m_graph->end();it++){
		if(it->second->size() > 1){
			UINT64 num_edges = it->second->size();
			for(UINT64 j = 0; j < num_edges; j++)	// For all the edges.
			{
				Edge * e1 = it->second->at(j);
				UINT64 destination1 = e1->getDestinationRead()->getReadID();

				// Only check if edge's source is smaller than destination, 
				// and not already in the to remove list
				if( !(e1->isLoop()) && !e1->isInvalid()) {
					for(UINT64 k = j + 1; k < num_edges; k++) {

						Edge * e2 = it->second->at(k);
						// Source read is the same, check if 
						// 1. destination read is the same
						// 2. orientation is the same
						if(destination1 == e2->getDestinationRead()->getReadID() ){
							if(e1->getOrientation() == e2->getOrientation()) {
								// The lengths are more than 95% similar
								if(abs((int)(e1->getOverlapOffset() - e2->getOverlapOffset())) < (int)(e2->getOverlapOffset()/20)) 
								{
									//FILE_LOG(logDEBUG1) << *e1 << " and\n " << *e2 << " are similar\n\n";
									e1->updateBaseByBaseCoverageStat(m_dataset);
									e2->updateBaseByBaseCoverageStat(m_dataset);
									UINT64 e1_reads = (e1->isListofReads() ? e1->getListofReadsSize() : 0);
									UINT64 e2_reads = (e2->isListofReads() ? e2->getListofReadsSize() : 0);
									// Check coverage depth and number of reads 
									// to decide which one to keep
									if(e1->getCovDepth() < e2->getCovDepth() ||	
											(e1->getCovDepth() == e2->getCovDepth() && 
											 e1_reads < e2_reads))
									{
										// remove e1
										e1->setInvalid();
										break;
									}
									else{
										// remove e2
										e2->setInvalid();
									}
									++counter;
								}
							}
						}
						// destination is not the same any more
						else{
							break;
						}
					}
				}
			}
			//Delete all edges marked as invalid
			size_t j=0;
			while(j<it->second->size())
			{
				Edge * e1 = it->second->at(j);
				if(e1->isInvalid())
				{
					removeEdge(e1);
				}
				else{
					j++;
				}
			}
		}
	}
	if(counter > 0){
		FILE_LOG(logINFO) << counter << " edges removed during bubble popping." << "\n";
	}
	CLOCKSTOP;
	FILE_LOG(logINFO) << "numberOfEdges = " << m_numberOfEdges << "\n";
	return counter;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  clipBranches
 *  Description:  When a noqde has multiple in- or out- edges, clip the edges with small
 *  		  overlap lengths
 * =====================================================================================
 */
UINT64 OverlapGraph::clipBranches(void)
{
	CLOCKSTART;
	UINT64 num_clip_branches = 0;

	UINT64 max_in_ovl, max_out_ovl, ovl;

	for(map<UINT64, t_edge_vec* >::iterator it=m_graph->begin(); it!=m_graph->end();it++){
		if(it->second->size() > 1){
			max_in_ovl = 0; max_out_ovl = 0;
			t_edge_vec inEdges, outEdges;
			vector<UINT64> inOvls, outOvls;
			for(UINT64 j = 0; j < it->second->size(); j++){
				Edge* e = it->second->at(j);
				// Find the first overlap length
				ovl = e->getOverlapLen();
				// Do not consider loop for now TODO might later
				if(!e->isLoop()){
					// In edges
					if(!((e->getOrientation() >> 1) & 1)){
						inEdges.push_back(e);
						inOvls.push_back(ovl);
						if(ovl > max_in_ovl){
							max_in_ovl = ovl;
						}
					}
					// Out edges
					else{
						outEdges.push_back(e);
						outOvls.push_back(ovl);
						if(ovl > max_out_ovl){
							max_out_ovl = ovl;
						}
					}
				}
			}
			if(inEdges.size() > 1){
				for(UINT64 k = 0; k < inEdges.size(); k++){
					if((inOvls.at(k) + minOvlDiffToClip) < max_in_ovl){
//						FILE_LOG(logDEBUG1) << "Break edge " << *(inEdges.at(k)) << " with overlap length " << inOvls.at(k) << " comparing to " << max_in_ovl << "\n";
						t_edge_vec sub_edges = inEdges.at(k)->breakEdge(0,m_dataset);
						removeEdge(inEdges.at(k));
						for(auto it = sub_edges.begin(); it != sub_edges.end(); ++it)
							insertEdge(*it);
						++num_clip_branches;
					}
				}
			}
			if(outEdges.size() > 1){
				for(UINT64 k = 0; k < outEdges.size(); k++){
					if((outOvls.at(k) + minOvlDiffToClip) < max_out_ovl){
//						FILE_LOG(logDEBUG1) << "Break edge " << *(outEdges.at(k)) << " with overlap length " << outOvls.at(k) << " comparing to " << max_out_ovl << "\n";
						t_edge_vec sub_edges = outEdges.at(k)->breakEdge(0,m_dataset);
						removeEdge(outEdges.at(k));
						for(auto it = sub_edges.begin(); it != sub_edges.end(); ++it)
							insertEdge(*it);
						++num_clip_branches;
					}
				}
			}
//			FILE_LOG(logDEBUG1) << "Node " << i << " done with clipBranches" << "\n" << "\n";
		}
	}
	FILE_LOG(logINFO) << "Short overlap branches clipped: " << num_clip_branches << "\n";
	CLOCKSTOP;
	FILE_LOG(logINFO) << "numberOfEdges = " << m_numberOfEdges << "\n";
	return num_clip_branches;
}



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  reduceLoops
 *  Description:  This function remove loops
 *                a>--->b>--->b>--->c
 *                a<---<b<--->b>--->c
 * =====================================================================================
 */
UINT64 OverlapGraph::reduceLoops(void)
{
	if(this->m_flowComputed == false) {
		return 0;
	}
	CLOCKSTART;
	UINT64 counter = 0, remove_counter = 0;
	Edge *ab,*bb,*bc;
	for(map<UINT64, t_edge_vec* >::iterator it=m_graph->begin(); it!=m_graph->end();it++){
		if(it->second->size() == 4) // only four edges. The loop is counted twice.
		{
			UINT64 loopCount = 0, incomingEdgeCount = 0, outgoingEdgeCount = 0;
			for(UINT64 j = 0; j< it->second->size(); j++)
			{
				if(it->second->at(j)->isLoop()) // This is a loop
				{
					loopCount++;
					bb = it->second->at(j);
				}
				else if(((it->second->at(j)->getOrientation() >> 1) & 1 ) == 0) // incoming edge
				{
					incomingEdgeCount++;
					ab = it->second->at(j)->getReverseEdge();
				}
				else if(((it->second->at(j)->getOrientation() >> 1) & 1) == 1) // outgoing edge
				{
					outgoingEdgeCount++;
					bc = it->second->at(j);
				}
			}
			/* Be aware that the loops without flow are not removed in function removeAllEdgesWithoutFlow.
			 * Therefore if merge an edge with flow and a loop with flow 0, the resulted edge will have flow 0,
			 * which is wrong. Here we will fix it, reassign the correct number of flow to the loop.
			 * If later reduce tree with flow 0 and two other edges, 
			 * then one edge with flow 0 would have been deleted before 
			 * it's contracted with the other branch, wrong!
			 * *** Now do not consider flow any more since it's not reliable with the removal of edges with flow
			 * *** Whenever an edge with flow is removed, it creates an imbalance of flow at its incident nodes.
			 */
			
			// two in the loop and one incoming and one outgoing
			if(loopCount==2 && incomingEdgeCount == 1 && outgoingEdgeCount == 1)
			{  
				if(bb->getOrientation() == 0)
				{
					++counter;
					Edge *new_edge = Add(ab,  bb->getReverseEdge());
					insertEdge(new_edge);
					removeEdge(ab);
					removeEdge(bb);
				}
				else if(bb->getOrientation() == 3)
				{
					++counter;
					Edge *new_edge = Add(ab, bb);
					insertEdge(new_edge);
					removeEdge(ab);
					removeEdge(bb);
				}
				else{
					FILE_LOG(logDEBUG1) << "Loop has wrong orientation, cannot be used to reduce the graph" << "\n";
					++remove_counter;
					removeEdge(bb);
				}
			}
			/* two in the loop and two incoming
			 * qstat  in the case: *--->b>---<b<---* */
			else if (loopCount==2 && incomingEdgeCount == 2 && outgoingEdgeCount == 0 && bb->getOrientation() == 2){ // 
				counter++;
				Edge *new_edge = Add(ab, bb);
				insertEdge(new_edge);
				removeEdge(ab);
				removeEdge(bb);
			}
			/* two in the loop and two incoming
			 * reduce in the case: *---<b<--->b>---* */
			else if (loopCount==2 && incomingEdgeCount == 0 && outgoingEdgeCount == 2 && bb->getOrientation() == 1){ // 
				counter++; 
				Edge *new_edge = Add(bb, bc);
				insertEdge(new_edge);
				removeEdge(bc);
				removeEdge(bb);
			}
			/* Cannot reduce graph */
			else if (loopCount == 2){
				remove_counter++;
				removeEdge(bb);
			}
		}
	}
	if(counter > 0 || remove_counter > 0){
		FILE_LOG(logINFO) <<" Loops reduced: " << counter << "\n"; // Number of loop we were able to reduce
		FILE_LOG(logINFO) <<" Loops removed: " << remove_counter << "\n"; // Number of loop we were able to reduce
	}
	CLOCKSTOP;
	FILE_LOG(logINFO) << "numberOfEdges = " << m_numberOfEdges << "\n";
	return (counter + remove_counter);
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  calculateBoundAndCost
 *  Description:  This function calculates the cost and bounds for an edge in the overlap graph.
 *                This function is very sensitive to the assembled contigs and to the cost function parameters
 * =====================================================================================
 */
void OverlapGraph::calculateBoundAndCost(const Edge *edge, INT64* FLOWLB, INT64* FLOWUB, INT64* COST)
{
	for(UINT64 i = 0; i < 3; i++)		// For the simple edges we put very high cost
	{
		FLOWLB[i] = 0; FLOWUB[i] = 10; COST[i] = 500000;
	}

	if(edge->isListofReads() && !(edge->getListofReadsSize()==0)) // Composite Edge
	{
		// Case1: Composite Edge of at least minFlowReadCountThreshold (default: 20) reads. Must have at least one unit of flow.
		// Case2: Composite Edge length is longer than minFlowEdgeLengthThreshold (default: 1000)
		if(edge->getListofReadsSize() >= minReadsCountInEdgeToBe1MinFlow || edge->getEdgeLength() > minEdgeLengthToBe1MinFlow)
		{
			s_nGoodEdges++;
			s_nReads_in_goodEdges += edge->getListofReadsSize();
			// the first 1 flow must be used and has a cost of 1
			// each additional flow up to 8 flows has a cost of 100000

			// this edge carry the first 1 flow
			FLOWLB[0] = 1; FLOWUB[0] = 1; COST[0] = 1;
			// this edge carries the second 1 flow with high cost
			FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000;
			// this edge provides additional flow after the second flow
			FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000;
		}
		else // Short composite edge containing less than 20 reads. May have zero flow.
		{
			// the first 1 flow may not be required, but has a low cost of 1
			// each additional flow up to 8 flows has a cost of 100000

			// this edge carries the first 1 flow
			FLOWLB[0] = 0; FLOWUB[0] = 1; COST[0] = 1;
			// this edge carries the second unit of flow with high cost
			FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000;
			// this edge provides additional flow after the two units of flow.
			FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000;
		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  findEdge
 *  Description:  Return edge with most reads on it, between source and destination
 * =====================================================================================
 */
Edge* OverlapGraph::findEdge(const UINT64 & source, const UINT64 & destination)
{
	t_edge_vec edges = findEdges(source, destination);
	if(edges.empty())
		return nullptr;
	else
		return edges.front();
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  findEdges
 *  Description:  Return all edges between source and destination, sorted by number of reads
 *  		  in decreasing order.
 * =====================================================================================
 */
t_edge_vec OverlapGraph::findEdges(const UINT64 & source, const UINT64 & destination)
{
	t_edge_vec edges;
	if(m_graph->find(source) != m_graph->end())
	{
		t_edge_vec *eList = m_graph->at(source);
		for(UINT64 i = 0; i < eList->size(); i++) // For the list of edges of the source node.
		{
			// check if there is an edge to destination
			if(eList->at(i)->getDestinationRead()->getReadID() == destination){
				edges.push_back(eList->at(i));	// return the edge.
			}
		}
	}
	if(edges.empty()){
		FILE_LOG(logDEBUG1) << "Cannot find edge from " << source << " to " << destination << "\n";
	}
	// Sort the edges by number of reads contained in decreasing order
	else
		sort(edges.begin(),edges.end(),compareEdgesByReads);
	return edges;
}
/***
 * Constructor to load simple edges and perform partial simplification
 */
OverlapGraph::OverlapGraph(const vector<std::string> &edge_files, string simplifyPartialPath, DataSet *dataSet,
		UINT64 minOvl = 40, UINT64 parallelThreadPoolSize = 1)
	: m_minOvl(minOvl), p_ThreadPoolSize(parallelThreadPoolSize)
{
	CLOCKSTART;
	m_numberOfNodes	= 0;
	m_numberOfEdges	= 0;
	m_flowComputed 	= false;

	m_dataset=dataSet;

	m_graph = new map<UINT64, t_edge_vec* >;
	// loop edgeFilenameList and load partially simplified graphs
	if(edge_files.size()>0)
	{
		//check if all partial graph files exist...
		bool reCalcParSimpEdges=true;
		for(size_t i = 0;i< edge_files.size();i++)
		{
			string composite_out_edge_file = outputFilenamePrefix + "_" + SSTR(i) +"_ParSimpleEdges.txt";
			if(!ifstream(composite_out_edge_file))
			{
				reCalcParSimpEdges=false;
				break;
			}
		}
		if(reCalcParSimpEdges)
		{
			FILE_LOG(logINFO) << "Partial graphs already exist. Loading those..."<<endl;
			string prev_composite_out_edge_file = outputFilenamePrefix + "_" + SSTR(0) +"_ParSimpleEdges.txt";
			readParEdges(prev_composite_out_edge_file);
			for(size_t i = 1;i < edge_files.size();i++)
			{
				string composite_out_edge_file = outputFilenamePrefix + "_" + SSTR(i) +"_ParSimpleEdges.txt";
				readParEdges(composite_out_edge_file);
			}
		}
		else
		{
			string prev_composite_out_edge_file = outputFilenamePrefix + "_" + SSTR(0) +"_ParSimpleEdges.txt";
			string runSimplifyExeStr = simplifyPartialPath + "/parsimplify " + edge_files[0] + " " + prev_composite_out_edge_file
					+ " " + SSTR(minOvl) + " " + SSTR(parallelThreadPoolSize) ;
			//Perform the first partial graph simplification
			int retStatus = system(runSimplifyExeStr.c_str());
			if(retStatus!=0)
				MYEXIT("Error in executing partial graph simplification.Could not invoke partial simplification process.Please check file paths and -simPth option.");

			for(size_t i = 1;i< edge_files.size();i++)
			{
				//Setup arguments for the next partial graph simplification
				string composite_out_edge_file = outputFilenamePrefix + "_" + SSTR(i) +"_ParSimpleEdges.txt";
				string simplifyExeStr = simplifyPartialPath + "/parsimplify";

				char *path = new char[simplifyExeStr.length()+1];
				strcpy(path, simplifyExeStr.c_str());

				char* tokens[6];

				string exeName = "parsimplify";

				tokens[0] = new char[exeName.length() + 1];
				strcpy(tokens[0], exeName.c_str());

				tokens[1]=new char[edge_files[i].length() + 1];
				strcpy(tokens[1], edge_files[i].c_str());

				tokens[2]=new char[composite_out_edge_file.length() + 1];
				strcpy(tokens[2], composite_out_edge_file.c_str());

				tokens[3]=new char[SSTR(minOvl).length() + 1];
				strcpy(tokens[3], SSTR(minOvl).c_str());

				tokens[4]=new char[SSTR(parallelThreadPoolSize).length() + 1];
				strcpy(tokens[4], SSTR(parallelThreadPoolSize).c_str());

				tokens[5]={NULL};

				pid_t child_pid;

				//Spawn next partial graph simplification
				int status = posix_spawn(&child_pid, path, NULL, NULL, tokens, environ);

				if(status == 0)   // Successfully launched next partial graph simplification
				{
					//Read the previous partial graph before waiting
					readParEdges(prev_composite_out_edge_file);
					// If I want to wait for the child before continuing:
					int   child_status;
					pid_t wait_result = waitpid( child_pid, &child_status, 0 );
					if(wait_result != child_pid)
						MYEXIT("Error in executing partial graph simplification.");
				}
				else
				{
					MYEXIT("Could not invoke partial simplification process.Please check file paths and -simPth option.");
					break;
				}
				delete tokens[1];
				delete tokens[2];
				delete tokens[3];
				delete tokens[4];
				delete path;
				prev_composite_out_edge_file=composite_out_edge_file;
			}
			//Read the last one...
			readParEdges(prev_composite_out_edge_file);
		}
	}
	else
		MYEXIT("No edge files provided. Therefore there is no graph to simplify. So we are done!!!");

	sortEdgesByDestID();
	FILE_LOG(logINFO) << "graph has " << m_graph->size() << " unique number of nodes.\n";
	FILE_LOG(logINFO) << "numberOfEdges loaded from edge file(s): " << m_numberOfEdges << "\n";

	UINT64 counter(0);
	do {
		counter = contractCompositeEdgesPar();
	} while (counter > 0);
	FILE_LOG(logINFO) << "numberOfEdges = " << m_numberOfEdges << "\n";

	CLOCKSTOP;
}

/*
 * Constructor to load composite edge graph that has been simplified
 */
OverlapGraph::OverlapGraph(const string parGlobalGraph, string simplifyPartialPath, DataSet *dataSet,
		UINT64 minOvl = 40, UINT64 parallelThreadPoolSize = 1)
	: m_minOvl(minOvl), p_ThreadPoolSize(parallelThreadPoolSize)
{
	CLOCKSTART;
	m_numberOfNodes	= 0;
	m_numberOfEdges	= 0;
	m_flowComputed 	= false;
	m_dataset=dataSet;
	m_graph = new map<UINT64, t_edge_vec* >;
	// loop partial global graph and load it
	FILE_LOG(logINFO) << "Loading partially simplified graph...\n";
	readParEdges(parGlobalGraph);
	FILE_LOG(logINFO) << "graph has " << m_graph->size() << " unique number of nodes.\n";
	FILE_LOG(logINFO) << "numberOfEdges loaded from edge file(s): " << m_numberOfEdges << "\n";
	UINT64 counter(0);
	do {
		counter = contractCompositeEdgesPar();
	} while (counter > 0);
	FILE_LOG(logINFO) << "numberOfEdges = " << m_numberOfEdges << "\n";
	CLOCKSTOP;
}


/*
 * ===  FUNCTION  ======================================================================
 *         Name:  simplifyGraph1
 *  Description:  simplify graph right after graph construction
 * =====================================================================================
 */
void OverlapGraph::graphPathFindInitial()
{
	FILE_LOG(logINFO) << "Initial simplification: contract composite edges, remove dead end nodes,"
		<< " and clip branches with very short overlap length.\n";
	// Composite edge contraction with remove dead end nodes
	CLOCKSTART;
	double prev = omp_get_wtime();
	UINT64 counter(0);
	do {
		removeDeadEndNodes();
		counter = contractCompositeEdgesPar();
		double curr = omp_get_wtime();
		double lapsed = curr - prev;
		if(lapsed>=DISK_GRAPH_UPDATE)
		{
			//Update checkpoint graph
			printAllEdges(outputFilenamePrefix+"_CurrGraph_.txt");
			prev = omp_get_wtime();
		}
	} while (counter > 1000);
	FILE_LOG(logERROR) << "numberOfEdges = " << m_numberOfEdges << "\n";
	/* disconnect the edges incident to nodes and have small overlap lengths */
	removeSimilarEdges();
	clipBranches();
	//calculateMeanAndSdOfInnerDistance();
	//findSupportByMatepairsAndMerge();
	FILE_LOG(logINFO) << "Initial simplification done.\n";
	CLOCKSTOP;
}

//=============================================================================
// Default destructor
//=============================================================================
OverlapGraph::~OverlapGraph()
{
	// Free the memory used by the overlap graph.
//	if (m_deletedEdges != nullptr){
//		for(UINT64 i = 0 ; i < m_deletedEdges->size(); ++i){
//			if (m_deletedEdges->at(i) != nullptr){
//				delete m_deletedEdges->at(i);
//				m_deletedEdges->at(i) = nullptr;
//			}
//		}
//		delete m_deletedEdges;
//		m_deletedEdges = nullptr;
//	}
	if (m_graph!= nullptr){
		for(map<UINT64, t_edge_vec* >::iterator it=m_graph->begin(); it!=m_graph->end();it++){
			for(UINT64 j = 0; j< it->second->size(); j++) {
					if (it->second->at(j) != nullptr){
						delete it->second->at(j);
						it->second->at(j) = nullptr;
					}
				}
				delete it->second;
				it->second = nullptr;
		}
		delete m_graph;
		m_graph = nullptr;
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  sortEdgesByLength
 *  Description:  For each node, sort its incident edges by their lengths in increasing order
 * =====================================================================================
 */
void OverlapGraph::sortEdgesByDestID()
{
	for(map<UINT64, t_edge_vec* >::iterator it=m_graph->begin(); it!=m_graph->end();it++){
			sort(it->second->begin(), it->second->end(), compareEdgesByDestID);
	}
}
void OverlapGraph::sortEdgesByLength()
{
	for(map<UINT64, t_edge_vec* >::iterator it=m_graph->begin(); it!=m_graph->end();it++){
			sort(it->second->begin(), it->second->end(), compareEdgesByLength);
	}
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  simplifyScaffoldGraph
 *  Description:  simplify scaffolded graph iteratively
 * =====================================================================================
 */
void OverlapGraph::simplifyScaffoldGraph(void)
{
	UINT64 counter = 0;
	do
	{
		counter = contractCompositeEdgesPar();
		counter += removeSimilarEdges();
		counter += removeDeadEndNodes();
		// the three steps below requires flow to be computed
		// if simplifyGraph is called in the unitig graph, these two functions will just return.
		counter += reduceLoops();	// Reduce loops

	} while (counter > 0);
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  simplifyGraph
 *  Description:  simplify graph iteratively
 * =====================================================================================
 */
void OverlapGraph::simplifyGraph(void)
{
	CLOCKSTART;
	double prev = omp_get_wtime();
	UINT64 counter = 0;
	do
	{
		counter = contractCompositeEdgesPar();
		counter += removeSimilarEdges();
		counter += removeDeadEndNodes();
		// the three steps below requires flow to be computed
		// if simplifyGraph is called in the unitig graph, these two functions will just return.
		counter += removeShortBranches();	// Remove short branches
		counter += reduceLoops();	// Reduce loops

		double curr = omp_get_wtime();
		double lapsed = curr - prev;
		if(lapsed>=DISK_GRAPH_UPDATE)
		{
			//Update checkpoint graph
			printAllEdges(outputFilenamePrefix+"_CurrGraph_.txt");
			prev = omp_get_wtime();
		}
	} while (counter > 0);
	CLOCKSTOP;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  calculateFlowStream
 *  Description:  Calculate min cost flow
 *  An sample input to the CS2 algorithm
 * 
 *  p min       3840      13449	// p min numberOfNodes numberOfEdges
 *  n          1         0	// initial flow of 0 in node 1, node 1 is the supersource
 *  n       3840         0	// initial flow of 0 in node 3840, which is the supersink (equal to the number of nodes)
 *  // edge from supersink to supersource, LowBound(1), UpperBound(1000000), Cost per unit of flow(1000000)
 *  a       3840          1          1    1000000    1000000	
 *  // edge from supersource to node 2, with the defined cost function
 *  a          1          2          0    1000000          0	
 *
 *  connect each node to supersource and supersink
 *  connect every edge in the original graph
 *
 * =====================================================================================
 */
void OverlapGraph::calculateFlowStream(void) 
{
	CLOCKSTART;
	// Add super source and super sink nodes, add edge from super sink to super source with very big cost
	// Add edge from super source to every node in the graph, also from every node in the graph to the super sink
	// Every edge will be assigned a lower bound and an upper bound of the flow (capacity), and the cost associated with the edge

	//Delete nodes in the graph that have no edges anymore
	map<UINT64, t_edge_vec* >::iterator it=m_graph->begin();
	while(it!=m_graph->end())
	{
		if(it->second->empty())
		{
			delete it->second;
			m_graph->erase(it);
		}
		else
			it++;
	}
	FILE_LOG(logDEBUG) <<"Real graph size:"<<m_graph->size()<<endl;
	FILE_LOG(logDEBUG) <<"Recorded graph size:"<<m_numberOfNodes<<endl;

	UINT64 V = m_numberOfNodes*2 + 2, E = m_numberOfEdges * 3 + m_numberOfNodes * 4 + 1 , SUPERSOURCE = 1, SUPERSINK = V;

	// Flow bounds and cost of the edges, 
	// cost function originally is a piecewise function with 3 segments
	INT64 FLOWLB[3], FLOWUB[3], COST[3];            
	stringstream ss;
	// Problem description: Number of nodes and edges in the graph.
	ss << "p min " << setw(10) << V << " " << setw(10) << E << "\n";    

	// Flow in the super source
	ss << "n " << setw(10) << SUPERSOURCE << setw(10) << " 0" << "\n";  

	// Flow in the super sink.
	ss << "n " << setw(10) << SUPERSINK << setw(10) << " 0" << "\n";    

	// Add an edge from super sink to super source with very high cost (almost infinity), also at most can be used once
	FLOWLB[0] = 1; 
	FLOWUB[0] = std::numeric_limits<UINT64>::max(); 
	COST[0]   = 1000000;
	ss << "a " << setw(10) << SUPERSINK << " " << setw(10) << SUPERSOURCE 
		<< " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] 
		<< " " << setw(10) << COST[0] << "\n"; 


	// For n nodes in the graph, CS2 requires that the nodes are numbered from 1 to n. 
	// In the overlap graph, the nodes does not have sequential ID. We need to convert them to 1 - n
	// Mapping between original node ID and cs2 node ID
	// Mapping between original node ID and cs2 node ID
	//Note here that the index is starts from zero. But, CS2 index starts from 1. So we add 1
	UINT64 indx=0;
	map<UINT64, size_t> readIdMap;
	map<UINT64, size_t> idReadMap;
	for (map<UINT64, t_edge_vec* >::iterator it=m_graph->begin(); it!=m_graph->end();it++){
		idReadMap.insert(make_pair(indx,it->first));
		readIdMap.insert(make_pair(it->first,indx++));
	}

	// This loop set lower bound and upper bound from super source to each node, and from each node to super sink. All costs are 0.
	UINT64 currentIndex = 1;
	for(map<UINT64, t_edge_vec* >::iterator it=m_graph->begin(); it!=m_graph->end();it++){
		// edges to and from the super source and super sink

		FLOWLB[0] = 0; FLOWUB[0] = 1000000; COST[0] = 0;
		ss << "a " << setw(10) << SUPERSOURCE << " " << setw(10) << 2 * currentIndex << " "
			<< setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";

		ss << "a " << setw(10) << SUPERSOURCE << " " << setw(10) << 2 * currentIndex + 1 << " "
			<< setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";

		ss << "a " << setw(10) << 2 * currentIndex << " " << setw(10) << SUPERSINK << " "
			<< setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";

		ss << "a " << setw(10) << 2 * currentIndex + 1 << " " << setw(10) << SUPERSINK << " "
			<< setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";
		currentIndex++;
	}
	FILE_LOG(logDEBUG) <<"Set source sink edges for each node."<<endl;
	// This loop converts the original bi-directed edges to directed edges (1 becomes 6).
	for(map<UINT64, t_edge_vec* >::iterator it=m_graph->begin(); it!=m_graph->end();it++){
		// edges to and from the super source and super sink
		for(UINT64 j = 0; j < it->second->size(); j++)
		{
			Edge *edge = it->second->at(j);
			UINT64 u = 0, v=0;
			//Get source index from rid
			auto nodeIndxIt = readIdMap.find(edge->getSourceRead()->getReadID());
			if(nodeIndxIt != readIdMap.end())
				u = nodeIndxIt->second + 1;
			else
				MYEXIT("ERROR: Source read ID:"+ SSTR(edge->getSourceRead()->getReadID()) +" not in map. This is an error!!!");

			//Get destination index from rid
			nodeIndxIt = readIdMap.find(edge->getDestinationRead()->getReadID());
			if(nodeIndxIt != readIdMap.end())
				v = nodeIndxIt->second + 1;
			else
				MYEXIT("ERROR: Destination read ID:"+ SSTR(edge->getDestinationRead()->getReadID()) +" not in map. This is an error!!!");

			if(u < v || (u == v && edge < edge->getReverseEdge())) {
				calculateBoundAndCost(edge, FLOWLB, FLOWUB, COST);
				// Here for each edge we add three edges with different values of cost and bounds.
				// Total 6 edges considering the reverse edge too.
				// For details on how to convert the edges off different types please see my thesis.
				UINT64 u1 = 2 * u, u2 = 2 * u + 1, v1 =  2 * v, v2 = 2 * v + 1;
				if(edge->getOrientation() == 0)
				{
					ss << "a " << setw(10) << v1 << " " << setw(10) << u1 << " " << setw(10) <<
						FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";
					ss << "a " << setw(10) << u2 << " " << setw(10) << v2 << " " << setw(10) <<
						FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";

					ss << "a " << setw(10) << v1 << " " << setw(10) << u1 << " " << setw(10) <<
						FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << "\n";
					ss << "a " << setw(10) << u2 << " " << setw(10) << v2 << " " << setw(10) <<
						FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << "\n";

					ss << "a " << setw(10) << v1 << " " << setw(10) << u1 << " " << setw(10) <<
						FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << "\n";
					ss << "a " << setw(10) << u2 << " " << setw(10) << v2 << " " << setw(10) <<
						FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << "\n";
				}
				else if(edge->getOrientation() == 1)
				{
					ss << "a " << setw(10) << v2 << " " << setw(10) << u1 << " " << setw(10) <<
						FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";
					ss << "a " << setw(10) << u2 << " " << setw(10) << v1 << " " << setw(10) <<
						FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";

					ss << "a " << setw(10) << v2 << " " << setw(10) << u1 << " " << setw(10) <<
						FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << "\n";
					ss << "a " << setw(10) << u2 << " " << setw(10) << v1 << " " << setw(10) <<
						FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << "\n";

					ss << "a " << setw(10) << v2 << " " << setw(10) << u1 << " " << setw(10) <<
						FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << "\n";
					ss << "a " << setw(10) << u2 << " " << setw(10) << v1 << " " << setw(10) <<
						FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << "\n";

				}
				else if(edge->getOrientation() == 2)
				{
					ss << "a " << setw(10) << u1 << " " << setw(10) << v2 << " " << setw(10) <<
						FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";
					ss << "a " << setw(10) << v1 << " " << setw(10) << u2 << " " << setw(10) <<
						FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";

					ss << "a " << setw(10) << u1 << " " << setw(10) << v2 << " " << setw(10) <<
						FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << "\n";
					ss << "a " << setw(10) << v1 << " " << setw(10) << u2 << " " << setw(10) <<
						FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << "\n";

					ss << "a " << setw(10) << u1 << " " << setw(10) << v2 << " " << setw(10) <<
						FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << "\n";
					ss << "a " << setw(10) << v1 << " " << setw(10) << u2 << " " << setw(10) <<
						FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << "\n";

				}
				else if(edge->getOrientation() == 3)
				{
					ss << "a " << setw(10) << u1 << " " << setw(10) << v1 << " " << setw(10) <<
						FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";
					ss << "a " << setw(10) << v2 << " " << setw(10) << u2 << " " << setw(10) <<
						FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";

					ss << "a " << setw(10) << u1 << " " << setw(10) << v1 << " " << setw(10) <<
						FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << "\n";
					ss << "a " << setw(10) << v2 << " " << setw(10) << u2 << " " << setw(10) <<
						FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << "\n";

					ss << "a " << setw(10) << u1 << " " << setw(10) << v1 << " " << setw(10) <<
						FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << "\n";
					ss << "a " << setw(10) << v2 << " " << setw(10) << u2 << " " << setw(10) <<
						FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << "\n";
				}
			}
		}
	}
	FILE_LOG(logINFO) <<"Finished initializing flow to edges.\n";

	FILE_LOG(logDEBUG) << "Number of edges with flow 1 set is " << OverlapGraph::s_nGoodEdges << "\n";
	FILE_LOG(logDEBUG) << "Number of reads contained in these edges is " << OverlapGraph::s_nReads_in_goodEdges << "\n";
	FILE_LOG(logINFO) << "Calling CS2 for flow analysis\n";
	stringstream oss;
	main_cs2(&ss, oss);
	FILE_LOG(logINFO) << "Flow analysis finished\n";

	if(loglevel > 2){
		std::string flowfile = outputFilenamePrefix + "_init.flow";
		FILE_LOG(logDEBUG) << "Print result flow in graph to " << flowfile << "\n";
		ofstream flowin(flowfile.c_str());
		flowin << ss.str();
		flowin.close();

		flowfile = outputFilenamePrefix + "_calc.flow";
		FILE_LOG(logDEBUG) << "Print result flow in graph to " << flowfile << "\n";
		ofstream flowout(flowfile.c_str());
		flowout << oss.str();
		flowout.close();
	}

	FILE_LOG(logINFO) <<"Start assigning calculated flow to edges.\n";
	std::string s, d, f;
	UINT64 lineNum = 0;

	while(!oss.eof()) {
		lineNum++;
		UINT64 source, destination, flow;
		// get the flow from CS2
		oss >> source >> destination >> flow;	

		if(source != SUPERSINK && source != SUPERSOURCE && destination != SUPERSOURCE && destination != SUPERSINK && flow!=0)
		{
			//FILE_LOG(logINFO) << "Processing: "<< lineNum<<" Src:"<<source<<" Dest:"<<destination<<" Flow:"<<flow<<endl;
			double sourceIndx = ceil(source/2);
			//Note here that the index in the map is starts from zero. But, CS2 index starts from 1. So we subtract 1
			sourceIndx--;
			// Map the source to the original graph
			UINT64 mySource = idReadMap[sourceIndx];

			double destIndx = ceil(destination/2);
			//Note here that the index in the map is starts from zero. But, CS2 index starts from 1. So we subtract 1
			destIndx--;
			// Map the destination in the original graph
			UINT64 myDestination = idReadMap[destIndx];

			Edge *edge = findEdge(mySource, myDestination);	// Find the edge in the original graph.
			if(edge){
				edge->m_flow += flow;	// Add the flow in the original graph.
				edge->getReverseEdge()->m_flow += flow;	// Also add the same flow to the twin edge to make sure that the flow is the same for the twin edges (flow is doubled if both directions have same flow)
			}
			else{
				FILE_LOG(logDEBUG1) << "edge does not exist in graph!\n";
			}
		}
	}
	this->m_flowComputed = true;

	CLOCKSTOP;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  removeAllEdgesWithoutFlow
 *  Description:  After flow analysis, remove all edges with flow 0
 * =====================================================================================
 */
UINT64 OverlapGraph::removeAllEdgesWithoutFlow()
{
	if(!m_flowComputed)
		return 0;
	CLOCKSTART;
	UINT64 num_edge_rm(0);
	for (map<UINT64, t_edge_vec* >::iterator it=m_graph->begin(); it!=m_graph->end();it++) // For each read.
	{
		if(!it->second->empty())	// If the read has some edges.
		{
			for(UINT64 j=0; j < it->second->size(); j++) // For each edge
			{
				Edge * edge = it->second->at(j);

				//Also remove loops without flow. This means, by default, 
				//the edges formed by loops that do not contain many reads will not be used.
				if(edge->m_flow == 0 && !edge->isLoop() && edge->getListofReadsSize() <= minReadsCountToHave0Flow &&
						edge->getEdgeLength() < minEdgeLengthToHave0Flow) {
					removeEdge(edge);
					++num_edge_rm;
				}
			}
		}
	}
	//Remove nodes that have become empty...
	//Delete nodes in the graph that have no edges anymore
	map<UINT64, t_edge_vec* >::iterator it=m_graph->begin();
	while(it!=m_graph->end())
	{
		if(it->second->empty())
		{
			delete it->second;
			m_graph->erase(it);
		}
		else
			it++;
	}
	FILE_LOG(logINFO) <<"No flow edges removed: " << num_edge_rm << "\n";
	CLOCKSTOP;
	return num_edge_rm;
}

//=============================================================================
// Function to compare two edges by the destination read number. 
// The read number was determined by the read std::string lexicographically.
// Used for sorting.
//=============================================================================
bool compareEdgesByDestID (const Edge *edge1, const Edge* edge2)
{
	if (edge1->getDestinationRead()->getReadID() < edge2->getDestinationRead()->getReadID())
	{
		return true;
	}
	else if (edge1->getDestinationRead()->getReadID() == edge2->getDestinationRead()->getReadID()){
		return (compareEdgesByLength(edge1, edge2));
	}
	else
		return false;
}


//=============================================================================
// Function to compare two edges by overlap offset. 
// Used for transitive edge removal (which is not included any more)
//=============================================================================
bool compareEdgesByLength (const Edge *edge1, const Edge* edge2)
{
	return (edge1->getEdgeLength() < edge2->getEdgeLength());
}



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  compareEdgesByReads
 *  Description:  Compare edges by number of reads contained in the them
 * =====================================================================================
 */
bool compareEdgesByReads (const Edge *edge1, const Edge* edge2)
{
	UINT64 num_reads1 = (edge1->isListofReads() ? edge1->getListofReadsSize() : 0);
	UINT64 num_reads2 = (edge2->isListofReads() ? edge2->getListofReadsSize() : 0);
	return (num_reads1 > num_reads2);
}


//=============================================================================
// Function to compare two std::string length
//=============================================================================
bool compareStringByLength (const std::string & seq1, const std::string & seq2)
{
	return (seq1.length() < seq2.length());
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getEdges
 *  Description:  save all the edges in the graph in a vector of pointers to these edges
 * =====================================================================================
 */
void OverlapGraph::getEdges(t_edge_vec & contigEdges) const
{
	CLOCKSTART;
	contigEdges.clear();
	for(UINT64 i = 1; i<= m_dataset->size(); i++)
	{
		auto it = m_graph->find(i);
		if(it != m_graph->end() && !it->second->empty()) // if this read has some edge(s) going out of it (since now the graph is directed)
		{
			t_edge_vec *eList = it->second;
			for(UINT64 j = 0; j < eList->size(); j++)
			{
				Edge * e = eList->at(j);
				if(e->isSmallerEdge()){  /* Only consider the edges between non-contained reads */
					contigEdges.push_back(e); // List of contigs.
				}
			}
		}
	}
	FILE_LOG(logINFO) << "Number of edges(contigs) in the graph to print: " << contigEdges.size() << "\n";
	// Sort the contigs by lengths in decreasing order
	std::sort(contigEdges.begin(), contigEdges.end(), compareEdgesByLength);
	std::reverse(contigEdges.begin(), contigEdges.end());
	FILE_LOG(logINFO) << "Number of edges(contigs) to print: " << contigEdges.size() << "\n";
	CLOCKSTOP;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  getEdges
 *  Description:  Print all the edges in the graph in a vector of pointers to these edges
 * =====================================================================================
 */

void OverlapGraph::printEdge(Edge *contigEdge, ostream & filePointer) const
{
	UINT64 source = contigEdge->getSourceRead()->getReadID();
	UINT64 destination = contigEdge->getDestinationRead()->getReadID();
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
		filePointer<<0<<",";					//no edits
		filePointer<<contigEdge->m_flow<<"\t";					//flow value
		if(contigEdge->isListofReads())
		{
			for(size_t j=0;j<contigEdge->getListofReadsSize();j++)
			{
				orientation=contigEdge->getInnerOrientation(j);
				filePointer<<"("<<contigEdge->getInnerReadID(j)<<","<<orientation<<","
						<<contigEdge->getInnerOverlapOffset(j)<<")";

			}
		}
		filePointer<<'\n';
	}
}


/*
 * ===  FUNCTION  ======================================================================
 *         Name:  printEdge
 *  Description:  Print all the edge sequence in the graph in a vector of pointers to these edges used by contig writing function
 * =====================================================================================
 */

void OverlapGraph::printEdge(Edge *contigEdge, ostream & edgeFilePointer,ostream & fileUsedReadPointer, UINT64 edgeNameID) const
{
	UINT64 source = contigEdge->getSourceRead()->getReadID();
	UINT64 destination = contigEdge->getDestinationRead()->getReadID();
	UINT64 orientation=contigEdge->getOrientation();
	if(source < destination || (source == destination && contigEdge < contigEdge->getReverseEdge()))
	{
		UINT64 offsetSum = 0, lastOffset=contigEdge->getOverlapOffset();
		if(contigEdge->isListofReads())
		{
			for(size_t j=0;j<contigEdge->getListofReadsSize();j++)
				offsetSum+=contigEdge-> getInnerOverlapOffset(j);
			lastOffset=contigEdge->getInnerOverlapOffset(contigEdge->getListofReadsSize()-1);
		}
		edgeFilePointer<< "contig_" << setfill('0') << setw(10) << edgeNameID<< "\t";
		edgeFilePointer<<source<<"\t";		// store the source read information
		edgeFilePointer<<destination<<"\t";	// store the destination read information
		fileUsedReadPointer<<source<<'\n';		//Write the source read as used up
		m_dataset->at(source)->setUsedRead(true);	//set as used
		fileUsedReadPointer<<destination<<'\n';		//Write the destination read as used up
		m_dataset->at(destination)->setUsedRead(true);	//set as used
		edgeFilePointer<<orientation<<",";
		edgeFilePointer<<contigEdge->getOverlapOffset()-offsetSum<<",";  //first overlap offset
		edgeFilePointer<<offsetSum+(contigEdge->getDestinationRead()->getReadLength()-lastOffset)<<","; //overlap length
		edgeFilePointer<<0<<",";					//no substitutions
		edgeFilePointer<<0<<"\t";					//no edits
		if(contigEdge->isListofReads())
		{
			for(size_t j=0;j<contigEdge->getListofReadsSize();j++)
			{
				orientation=contigEdge->getInnerOrientation(j);
				edgeFilePointer<<"("<<contigEdge->getInnerReadID(j)<<","<<orientation<<","
						<<contigEdge->getInnerOverlapOffset(j)<<")";
				fileUsedReadPointer<<contigEdge->getInnerReadID(j)<<'\n';		//Write the read as used up
				m_dataset->at(contigEdge->getInnerReadID(j))->setUsedRead(true);
			}
		}
		edgeFilePointer<<'\n';
	}
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  printEdge
 *  Description:  Print all the edge sequence in the graph in a vector of pointers to these edges used by contig writing function
 * =====================================================================================
 */

void OverlapGraph::printEdgeCoverage(Edge *contigEdge, ostream & filePointer, UINT64 edgeNameID) const
{
	UINT64 source = contigEdge->getSourceRead()->getReadID();
	UINT64 destination = contigEdge->getDestinationRead()->getReadID();
	if(source < destination || (source == destination && contigEdge < contigEdge->getReverseEdge()))
	{
		vector<UINT64> coverageVals = contigEdge->getBaseByBaseCoverageValues(m_dataset);
		filePointer<< "contig_" << setfill('0') << setw(10) << edgeNameID<< ",";
		for(size_t j=0;j<coverageVals.size();j++)
		{
			filePointer<<coverageVals[j]<<",";

		}
		filePointer<<'\n';
	}
}

/*===== FUNCTION: operator<< ==================================================
 *  This function prints the overlap graph in overlap_graph->cytoscape file. The graph can be viewed by
 *  Cytoscape (free software available at http://cytoscape.org/)
 *  It also stores the contigs in a file.
 *
 *graph: {
 *layoutalgorithm :forcedir
 *fdmax:704
 *tempmax:254
 *tempmin:0
 *temptreshold:3
 *tempscheme:3
 *tempfactor:1.08
 *randomfactor:100
 *gravity:0.0
 *repulsion:161
 *attraction:43
 *ignore_singles:yes
 *node.fontname:"helvB10"
 *edge.fontname:"helvB10"
 *node.shape:box
 *node.width:80
 *node.height:20
 *node.borderwidth:1
 *node.bordercolor:31
 *node: { title:"43" label: "43" }	// node, Title and label are both node ID 43 (and read number)
 *node: { title:"65" label: "65" }
 *............................................
 *............................................
 * edges from source node 43 to destination node 32217, thickness of 3 means composite edge, thickness of 1 for simple edge
 * edge type of backarrowstyle:solid arrowstyle:solid color: green is >----------------<
 * edge type of arrowstyle:solid color: red is <----------------<
 * edge type of arrowstyle: none color: blue  is <------------------->
 * (1,0x,206,30) means (Flow, coverageDepth, OverlapOffset, numberOfReads)
 *
 *edge: { source:"43" target:"32217" thickness: 3 backarrowstyle:solid arrowstyle:solid color: green label: "(1,0x,206,30)" }
 *edge: { source:"65" target:"38076" thickness: 3 arrowstyle:solid color: red label: "(0,0x,75,11)" }
 *edge: { source:"280" target:"47580" thickness: 3 arrowstyle: none color: blue label: "(0,0x,123,11)" }
 *}
 *=============================================================================
 */
ostream& operator<< (ostream &out, const OverlapGraph & graph)
{
	CLOCKSTART;
	t_edge_vec contig_edges;
	graph.getEdges(contig_edges);
	std::string format = "cytoscape";

	if (format == "gdl"){
		UINT64 thickness;
		// Graph specification before the nodes and edges
		out << "graph: {\n\
			layoutalgorithm :forcedir\n\
			fdmax:704\n\
			tempmax:254\n\
			tempmin:0\n\
			temptreshold:3\n\
			tempscheme:3\n\
			tempfactor:1.08\n\
			randomfactor:100\n\
			gravity:0.0\n\
			repulsion:161\n\
			attraction:43\n\
			ignore_singles:yes\n\
			node.fontname:\"helvB10\"\n\
			edge.fontname:\"helvB10\"\n\
			node.shape:box\n\
			node.width:80\n\
			node.height:20\n\
			node.borderwidth:1\n\
			node.bordercolor:31\n";

		// All the nodes, title and label
		for(UINT64 i = 1; i< graph.m_dataset->size(); i++)
		{
			if(!graph.m_graph->at(i)->empty())
				// Print nodes even if there is some edge incident to it
				out << "node: { title:\""<< i <<"\" label: \"" << i << "\" }\n";
		}

		// All the edges
		for (UINT64 i = 0; i < contig_edges.size(); i++)
		{
			Edge *e = contig_edges.at(i);
			int num_reads = 0;
			if(e->isListofReads())
				num_reads = e->getListofReadsSize();
			thickness = (num_reads > 0) ? 1: 3;
			UINT64 source = e->getSourceRead()->getReadID(), destination = e->getDestinationRead()->getReadID();
			// Edge label: (first overlap length, edge length, number of reads, overlap offset, last overlap length)
			if(e->isSmallerEdge()) {
				if(e->getOrientation() == 0)
					out << "edge: { source:\"" << source << "\" target:\"" << destination
						<< "\" thickness: " << thickness << " arrowstyle: none backarrowstyle: solid color: red label: \"("
						<< e->getOverlapLen() << "," << e->getEdgeLength() << "," << num_reads << "," << e->getOverlapOffset()
						<< "," << e->getReverseEdge()->getOverlapLen()
						<< ")\" }" << "\n";
				else if(e->getOrientation() == 1)
					out << "edge: { source:\"" << source << "\" target:\"" << destination << "\" thickness: "
						<< thickness << " backarrowstyle:solid arrowstyle:solid color: green label: \"("
						<< e->getOverlapLen() << "," << e->getEdgeLength() << "," << num_reads << "," << e->getOverlapOffset()
						<< "," << e->getReverseEdge()->getOverlapLen()
						<< ")\" }" << "\n";
				else if(e->getOrientation() == 2)
					out << "edge: { source:\"" << source << "\" target:\"" << destination << "\" thickness: "
						<< thickness << " arrowstyle: none color: blue label: \"("
						<< e->getOverlapLen() << "," << e->getEdgeLength() << "," << num_reads << "," << e->getOverlapOffset()
						<< "," << e->getReverseEdge()->getOverlapLen()
						<< ")\" }" << "\n";
				else if(e->getOrientation() == 3)
					out << "edge: { source:\"" << source << "\" target:\"" << destination << "\" thickness: "
						<< thickness << " arrowstyle:solid color: red label: \"("
						<< e->getOverlapLen() << "," << e->getEdgeLength() << "," << num_reads << "," << e->getOverlapOffset()
						<< "," << e->getReverseEdge()->getOverlapLen()
						<< ")\" }" << "\n";
			}
		}
		out << "}";
	}
	else if (format == "cytoscape"){
		out << "source\ttarget\tfirtOvl\tcontigLen\tnumReads\toffset\tlastOvl\tedgeType\n";
		for (UINT64 i = 0; i < contig_edges.size(); i++)
		{
			Edge *e = contig_edges.at(i);
			int num_reads = 0;
			if(e->isListofReads())
				num_reads = e->getListofReadsSize();
			UINT64 source = e->getSourceRead()->getReadID(), destination = e->getDestinationRead()->getReadID();
			// Edge label: (first overlap length, edge length, number of reads, overlap offset, last overlap length)
			if(e->isSmallerEdge()) {
				out << source << "\t" << destination << "\t" << e->getOverlapLen() << "\t" << e->getEdgeLength() << "\t"
					<< num_reads << "\t" << e->getOverlapOffset() << "\t" << e->getReverseEdge()->getOverlapLen() << "\t"
					<< static_cast<int>(e->getOrientation()) << "\n";
			}
		}

	}
	CLOCKSTOP;
	return out;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  printContigs
 *  Description:  Print all the edges in the graph
 * =====================================================================================
 */
void OverlapGraph::printAllEdges(string edge_file) const
{
	CLOCKSTART;
	ofstream filePointer;
	filePointer.open(edge_file.c_str());

	for(UINT64 i = 1; i<= m_dataset->size(); i++)
	{
		auto it = m_graph->find(i);
		if(it != m_graph->end() && !it->second->empty()) // if this read has some edge(s) going out of it (since now the graph is directed)
		{
			t_edge_vec *eList = it->second;
			for(UINT64 j = 0; j < eList->size(); j++)
			{
				Edge * e = eList->at(j);
				if(e->isSmallerEdge()){
					printEdge(e,filePointer);
				}
			}
		}
	}
	filePointer.close();
	CLOCKSTOP;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  printContigs
 *  Description:  Print all the edges in the graph
 * =====================================================================================
 */
void OverlapGraph::readParEdges(string edge_file)
{
	CLOCKSTART;
	// Open file
	ifstream filePointer;
	filePointer.open(edge_file.c_str());
	if(!filePointer.is_open() )
		MYEXIT("Unable to open file: "+edge_file);

	FILE_LOG(logINFO) << "Reading partial edge file:"<<edge_file<<endl;

	// Read file
	std::string edge_text;
	while(getline(filePointer,edge_text)) {

		vector<string> tok = Utils::split(edge_text,'\t');
		Read *source = m_dataset->at(std::stoull(tok[0],nullptr,0));
		Read *destination = m_dataset->at(std::stoull(tok[1],nullptr,0));

		vector<string> infoTok = Utils::split(tok[2],',');
		UINT64 orientationForward = std::stoull(infoTok[0],nullptr,0);
		UINT64 overlapOffsetForward = std::stoull(infoTok[1],nullptr,0);
		UINT64 flowValue=0;
		//If there is flow information; load it
		if(infoTok.size()>5)
			flowValue=std::stoull(infoTok[5],nullptr,0);

		// Make the forward edge list
		UINT64 *listReadsForward = nullptr;
		UINT64 innerReadCount=0;
		UINT64 usedReadCtr=0;
		UINT64 unUsedMate=0;
		if(tok.size()>3)			//If composite edge load the inner reads
			usedReadCtr = createFwdList(tok[3], &listReadsForward, innerReadCount, unUsedMate, m_dataset);		//Returns the count of inner reads that have been used before

		//Do not load edge into graph if the edge was used in previous assemblies...
		if(isUsedEdge(innerReadCount,usedReadCtr,unUsedMate,source,destination))
		{
			delete[] listReadsForward;
			listReadsForward = nullptr;
		}
		else
		{
			Edge *edgeForward = new Edge(source,destination,orientationForward,
					overlapOffsetForward, listReadsForward, innerReadCount, flowValue);

			// Make the reverse edge
			UINT64 *listReadsReverse = nullptr;
			UINT64 lRSize=0;
			if(tok.size()>3)
				createRevList(edgeForward, &listReadsReverse, lRSize, m_dataset);
			UINT64 overlapOffsetReverse = edgeForward->getEdgeLength()-source->getReadLength();
			Edge *edgeReverse = new Edge(destination,source,twinEdgeOrientation(orientationForward),
					overlapOffsetReverse, listReadsReverse, lRSize, flowValue);

			edgeForward->setReverseEdge(edgeReverse);
			edgeReverse->setReverseEdge(edgeForward);

			//Insert edge in graph
			if(!existsEdge(edgeForward))
				insertEdge(edgeForward);
		}

	}
	filePointer.close();
	CLOCKSTOP;
}
/*
 * Check if an edge can be considered as used.
 */
bool OverlapGraph::isUsedEdge(UINT64 lFSize, UINT64 usedReadCtr,UINT64 unUsedMate, Read *source, Read *destination)
{
	//Removed used edges; If all the inner reads, source, destination have been used in previous assembly; do not load this edge
	if(lFSize > 0 && usedReadCtr>0 && usedReadCtr > (lFSize*0.5) && unUsedMate < (usedReadCtr*0.5))	//Composite edge
	{
		return true;
	}
	else if(lFSize == 0)		//Single edge
	{
		//Removed used edges; If both source and destination and corresponding
		//mates have been used in previous assembly; do not load this edge
		UINT64 sourceMate = m_dataset->getMatePair(source->getReadID());
		UINT64 destinationMate = m_dataset->getMatePair(destination->getReadID());
		if(sourceMate==0 || destinationMate==0)
		{
			if(source->isUsedRead() && destination->isUsedRead())
				return true;
		}
		else if(sourceMate > 0 && destinationMate > 0)
		{
			if(source->isUsedRead() && destination->isUsedRead() &&
					m_dataset->at(sourceMate)->isUsedRead() && m_dataset->at(destinationMate)->isUsedRead())
				return true;
		}
		else if(sourceMate > 0)
		{
			if(source->isUsedRead() && destination->isUsedRead() &&
					m_dataset->at(sourceMate)->isUsedRead())
				return true;
		}
		else
		{
			if(source->isUsedRead() && destination->isUsedRead() &&
					m_dataset->at(destinationMate)->isUsedRead())
				return true;
		}
	}
	return false;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  loadStringFromReadsFile
 *  Description:  Fill the bases for all the edges in the graph, by streaming the reads in the
 *  		  reads file. Whenever a read is read from the reads file, its nucleotide sequence
 *  		  is streamed in, and the edges where this read resides will have the corresponding
 *  		  bases populated(determined). In the meantime, the positions with mismatches will
 *  		  have a count of each base option.
 * =====================================================================================
 */
void OverlapGraph::loadStringFromReadsFile(const std::string &read_file, UINT64 & readID)
{
	CLOCKSTART;
	FILE_LOG(logINFO) << "load read strings and fill the bases for the edges: " << read_file << "\n";
	// To count of reads in this file
	UINT64 readCount = 0;

	if(read_file.substr(read_file.length() - 3 )==".gz")
	{
#ifdef INCLUDE_READGZ
			gzFile fp;
			kseq_t *seq;
			int l;
			fp = gzopen(read_file.c_str(), "r");
			seq = kseq_init(fp);
			while ((l = kseq_read(seq)) >= 0) {
				string line1=seq->seq.s;
				transform(line1.begin(), line1.end(), line1.begin(), ::toupper);
				if(!testRead(line1)) // Test the read is of good quality. If not replace all N's
					std::replace( line1.begin(), line1.end(), 'N', 'A');

				populate_read(readID, line1);
				++readID;
				++readCount;
				if(readCount % 1000000 == 0){
					FILE_LOG(logDEBUG) << setw(10) << (readCount/1000000)  << ",000,000"
						<< " reads streamed, "
						<< setw(7) << checkMemoryUsage() << " MB\n";
				}
			}
			kseq_destroy(seq);
			gzclose(fp);
#else
			MYEXIT("Unknown input file format. Looks like the file is in gzip compressed format."
					"The Omega3 code was not built with ZLIB using READGZ=1. To assemble either uncompress"
					"the file or build Omega3 with ZLIB library using make \"make READGZ=1\".");
#endif
	}
	else
	{
		// Open file
		ifstream filePointer;
		filePointer.open(read_file.c_str());
		if(!filePointer.is_open()){
			FILE_LOG(logERROR) << "Unable to open file: " << read_file << "\n";
			return;
		}
		vector<std::string> line;
		std::string line0,line1, text;
		enum FileType {FASTA, FASTQ, UNDEFINED};
		FileType fileType = UNDEFINED;

		while(!filePointer.eof())
		{
			// Check FASTA and FASTQ
			if(fileType == UNDEFINED) {
				getline (filePointer,text);
				if (text.length() > 0){
					if(text[0] == '>')
						fileType = FASTA;
					else if(text[0] == '@')
						fileType = FASTQ;
					else{
						FILE_LOG(logERROR) << "Unknown input file format."<<"\n";
						break;
					}
					filePointer.seekg(0, ios::beg);
				}
			}

			line.clear();

			// FASTA file read
			if(fileType == FASTA) {
				getline (filePointer,line0);	// get ID line
				getline (filePointer,line1,'>');	// get string line

				line1.erase(std::remove(line1.begin(), line1.end(), '\n'),
						line1.end());
			}
			// FASTQ file read
			else if(fileType == FASTQ) {
				getline(filePointer, line0);	// ID
				getline(filePointer, line1);	// String
				// Ignore the next two lines
				filePointer.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				filePointer.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			}

			// Get ReadID after removing the > or @ identifier and convert string to UINT64
			populate_read(readID, line1);
			++readID;
			++readCount;
			if(readID % 1000000 == 0 && readID > 0){
				FILE_LOG(logINFO) << setw(10) << (readID / 1000000) << ",000,000"
					<< " reads streamed, "
					<< setw(7) << checkMemoryUsage() << " MB\n";
			}
		}
		filePointer.close();
	}
	FILE_LOG(logINFO) << setw(10) << readCount << " reads streamed from this read file\n";
	CLOCKSTOP;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  populate_read
 *  Description:  Given a read in the dataset, and its sequence, populate the strings for
 *   		  the edges.
 * =====================================================================================
 */
void OverlapGraph::populate_read(const UINT64 &readID, const std::string & read_str)
{
	Read *read = m_dataset->at(readID);
	std::string read_str_reverse = Utils::reverseComplement(read_str);
	// Edges with read as source or destination
	auto srcIt = m_graph->find(readID);
	if(srcIt != m_graph->end() && !srcIt->second->empty()) // if this read has some edge(s).
	{
		for (auto it = srcIt->second->begin(); it != srcIt->second->end(); ++it){
			if((*it)->isSmallerEdge()){
				Edge *edge = *it;
				if (((edge->getOrientation() >> 1) & 1))
					edge->loadReadString(read_str, -1);
				else
					edge->loadReadString(read_str_reverse, -1);
			}
			else{
				Edge *edge = (*it)->getReverseEdge();
				if ((edge->getOrientation() & 1))
					edge->loadReadString(read_str, -2);
				else
					edge->loadReadString(read_str_reverse, -2);
			}
		}
	}
	// Edges with read on it
	vector< t_edge_loc_pair > fwd_edges = read->getFwdEdges();
	vector< t_edge_loc_pair > bwd_edges = read->getBwdEdges();

	for(auto it = fwd_edges.cbegin(); it != fwd_edges.cend(); ++it){
		it->first->loadReadString(read_str, it->second);
	}
	for(auto it = bwd_edges.cbegin(); it != bwd_edges.cend(); ++it){
		it->first->loadReadString(read_str_reverse, it->second);
	}
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  printContigs
 *  Description:  Print contigs/scaffolds for all the edges in the graph, by streaming all the reads files.
 * =====================================================================================
 */
void OverlapGraph::streamContigs(const vector<std::string> &read_SingleFiles,const vector<std::string> &read_PairFiles,
		vector<std::string> &read_PairInterFiles, string contig_file, string edge_file,string edge_cov_file,
		string usedReadFileName, string namePrefix, UINT64 &printed_contigs)
{
	CLOCKSTART;
	//streaming reads
	UINT64 readID = 1;
	for(auto it = read_PairFiles.cbegin(); it != read_PairFiles.cend(); ++it){
		loadStringFromReadsFile(*it, readID);
	}
	for(auto it = read_PairInterFiles.cbegin(); it != read_PairInterFiles.cend(); ++it){
		loadStringFromReadsFile(*it, readID);
	}
	for(auto it = read_SingleFiles.cbegin(); it != read_SingleFiles.cend(); ++it){
		loadStringFromReadsFile(*it, readID);
	}

	t_edge_vec contigEdges; 
	getEdges(contigEdges);

	//Open contig file
	ofstream contigFilePointer;
	contigFilePointer.open(contig_file.c_str());
	if(!contigFilePointer)
		MYEXIT("Unable to open file: "+contig_file);

	//Open edge file
	ofstream edgeFilePointer;
	edgeFilePointer.open(edge_file.c_str());
	if(!edgeFilePointer)
		MYEXIT("Unable to open file: "+edge_file);

	//Open coverage file
	ofstream fileCoveragePointer;
	fileCoveragePointer.open(edge_cov_file.c_str());
	if(!fileCoveragePointer)
		MYEXIT("Unable to open file: "+edge_cov_file);

	//Open used read file
	ofstream fileUsedReadPointer;
	fileUsedReadPointer.open(usedReadFileName.c_str(), std::ofstream::trunc);
	if(!fileUsedReadPointer)
		MYEXIT("Unable to open file: "+usedReadFileName);

	#pragma omp parallel for schedule(dynamic) num_threads(p_ThreadPoolSize)
	for(auto it = contigEdges.begin(); it < contigEdges.end(); ++it)
	{
		if((*it)->getEdgeLength() >= minContigLengthTobeReported && (*it)->getListofReadsSize() >= minNumberofReadsTobePrinted ){
			string contigString = (*it)->getEdgeString();
			#pragma omp critical
			{
				++printed_contigs;
				printEdge(*it,edgeFilePointer,fileUsedReadPointer,printed_contigs);
				printEdgeCoverage(*it,fileCoveragePointer,printed_contigs);
				(*it)->updateBaseByBaseCoverageStat(m_dataset);
				contigFilePointer << ">"<<namePrefix<<"_" << setfill('0') << setw(10) << printed_contigs
				<<" Coverage: " << (*it)->getCovDepth()
				<<" Length: " << contigString.length() << "\n";

				UINT32 start=0;
				do
				{
					contigFilePointer << contigString.substr(start, 100) << "\n";  // save 100 BP in each line.
					start+=100;
				} while (start < contigString.length());
			}
		}
	}
	edgeFilePointer.close();
	fileCoveragePointer.close();
	contigFilePointer.close();
	fileUsedReadPointer.close();
	FILE_LOG(logINFO) << "Total number of contigs printed: " << printed_contigs << "\n";
	CLOCKSTOP;
}


/*
 * ===  FUNCTION  ======================================================================
 *         Name:  printContigs
 *  Description:  Print contigs/scaffolds for all the edges in the graph, from already loaded read files.
 * =====================================================================================
 */
/*void OverlapGraph::printContigs(string contig_file, string edge_file,string edge_cov_file,
		string usedReadFileName, string namePrefix, UINT64 &printed_contigs)
{
	CLOCKSTART;

	t_edge_vec contigEdges;
	getEdges(contigEdges);

	//Open contig file
	ofstream contigFilePointer;
	contigFilePointer.open(contig_file.c_str());
	if(!contigFilePointer)
		MYEXIT("Unable to open file: "+contig_file);

	//Open edge file
	ofstream edgeFilePointer;
	edgeFilePointer.open(edge_file.c_str());
	if(!edgeFilePointer)
		MYEXIT("Unable to open file: "+edge_file);

	//Open coverage file
	ofstream fileCoveragePointer;
	fileCoveragePointer.open(edge_cov_file.c_str());
	if(!fileCoveragePointer)
			MYEXIT("Unable to open file: "+edge_cov_file);

	//Open used read file
	ofstream fileUsedReadPointer;
	fileUsedReadPointer.open(usedReadFileName.c_str(), std::ofstream::trunc);
	if(!fileUsedReadPointer)
			MYEXIT("Unable to open file: "+usedReadFileName);

	#pragma omp parallel for schedule(dynamic) num_threads(p_ThreadPoolSize)
	for(auto it = contigEdges.begin(); it < contigEdges.end(); ++it)
	{
		if((*it)->getEdgeLength() >= minContigLengthTobeReported && (*it)->getListofReadsSize() >= minNumberofReadsTobePrinted ){
			populate_edge(*it);
			string contigString = (*it)->getEdgeString();
			#pragma omp critical
			{
				++printed_contigs;
				printEdge(*it,edgeFilePointer,fileUsedReadPointer,printed_contigs);
				printEdgeCoverage(*it,fileCoveragePointer,printed_contigs);
				//contigFilePointer << ">"<<namePrefix<<"_" << setfill('0') << setw(10) << printed_contigs
				//<< " Edge ("  << (*it)->getSourceRead()->getReadID() << ", "
				//<< (*it)->getDestinationRead()->getReadID()
				//<< ") String Length: " << contigString.length() << "\n";
				(*it)->updateBaseByBaseCoverageStat(m_dataset);
				contigFilePointer << ">"<<namePrefix<<"_" << setfill('0') << setw(10) << printed_contigs
				<<" Coverage: " << (*it)->getCovDepth()
				<<" Length: " << contigString.length() << "\n";

				UINT32 start=0;
				do
				{
					contigFilePointer << contigString.substr(start, 100) << "\n";  // save 100 BP in each line.
					start+=100;
				} while (start < contigString.length());
			}
		}
	}
	edgeFilePointer.close();
	fileCoveragePointer.close();
	contigFilePointer.close();
	fileUsedReadPointer.close();
	FILE_LOG(logINFO) << "Total number of contigs printed: " << printed_contigs << "\n";
	CLOCKSTOP;
}
*/
/**********************************************************************************************************************
	Returns true if the read contains only {A,C,G,T} and does not contain more than 80% of the same nucleotide
**********************************************************************************************************************/
bool OverlapGraph::testRead(const string & read)
{
	UINT64 readLength = read.length();
	for(UINT64 i = 0; i < readLength; i++) // Count the number of A's, C's , G's and T's in the string.
	{
		if(read[i]!= 'A' && read[i] != 'C' && read[i] != 'G' && read[i] != 'T')
			return false;
	}
	return true;
}


/*
 * ===  FUNCTION  ======================================================================
 *         Name:  populate_edge
 *  Description:  Given a edge in the graph, populate the string for
 *   		  the edge.
 * =====================================================================================
 */
/*void OverlapGraph::populate_edge(Edge *edge)
{
	string srcReadStr = edge->getSourceRead()->getStringForward();
	if (((edge->getOrientation() >> 1) & 1))
		edge->loadReadString(srcReadStr, -1);
	else
	{
		std::string srcReadStr_rev = Utils::reverseComplement(srcReadStr);
		edge->loadReadString(srcReadStr_rev, -1);
	}

	string destReadStr = edge->getDestinationRead()->getStringForward();
	if ((edge->getOrientation() & 1))
		edge->loadReadString(destReadStr, -2);
	else
	{
		std::string destReadStr_rev = Utils::reverseComplement(destReadStr);
		edge->loadReadString(destReadStr_rev, -2);
	}

	for(size_t i=0;i<edge->getListofReadsSize();i++)
	{
		string read_str = m_dataset->at(edge->getInnerReadID(i))->getStringForward();
		UINT64 orient = edge->getInnerOrientation(i);
		if(orient==0)
			read_str = Utils::reverseComplement(read_str);
		edge->loadReadString(read_str, i);
	}
}
*/

/**********************************************************************************************************************
	Merge two edges in the overlap graph.
	CP: the flow of the new edge is the minimum flow of the two old edges and the flow of the new edge is deducted from those of the old edges
	CP: Remove the old edges that doesn't have any flow left, but keep the old edges that still have flow left.
**********************************************************************************************************************/
void OverlapGraph::merge2Edges(Edge *edge1, Edge *edge2)
{
	Edge *new_edge = Add(edge1, edge2);
	UINT16 flow = min(edge1->m_flow,edge2->m_flow);		// We take the minimum of the flow in the new edge.

	new_edge->m_flow = flow;						// Modify the flow in the new forward edge.
	new_edge->getReverseEdge()->m_flow = flow;		// Modify the flow in the new reverse edge.

	insertEdge(new_edge);						    // Insert the new new edge in the graph.

	edge1->m_flow = edge1->m_flow - flow;				// Remove the used flow from edge1.
	edge1->getReverseEdge()->m_flow = edge1->m_flow;	// Remove the used flow from the reverse of edge1.

	edge2->m_flow = edge2->m_flow - flow;				// Remove the used flow from edge2.
	edge2->getReverseEdge()->m_flow = edge2->m_flow;	// Remove the used flow from the reverse of edge2.

	if(edge2 != edge1->getReverseEdge() && (edge2->m_flow == 0 || flow == 0))
		removeEdge(edge2);
	if(edge1->m_flow == 0 || flow == 0)				// If no flow left in edge1
		removeEdge(edge1);							// edge1 is deleted from the graph.
}


/**********************************************************************************************************************
	Check for paths in the overlap graph between each mate-pair and calculate the support by using the paths.
**********************************************************************************************************************/
UINT64 OverlapGraph::findSupportByMatepairsAndMerge(void)
{
	CLOCKSTART;

	/***TO DO: Change logic to check if paired-end data exists***/
	// if the file set is not mate-pair, then just skip
	//if(meanOfInsertSizes.size() == 0) // no mate-pair
	//  return 0;
	UINT64 noPathsFound = 0, pathsFound = 0, mpOnSameEdge=0;

	vector <pairedEdges> listOfPairedEdges;

	vector <vector <pairedEdges>* > listOfPairedEdges_localList;
	for(UINT64 i = 0; i < p_ThreadPoolSize; i++)	// for each thread create a vector of pair edges to populate
	{
		vector <pairedEdges> *listOfPairedEdges_local = new vector <pairedEdges>;
		listOfPairedEdges_localList.push_back(listOfPairedEdges_local);
	}
	#pragma omp parallel for schedule(guided) reduction(+:noPathsFound), reduction(+:pathsFound), reduction(+:mpOnSameEdge) num_threads(p_ThreadPoolSize)
	for(UINT64 i = 1; i <= m_dataset->size(); i++)	// for each read in the dataset
	{
		UINT64 threadID = omp_get_thread_num();
		vector <pairedEdges> *listOfPairedEdges_local = listOfPairedEdges_localList[threadID];
		Read * read1 = m_dataset->at(i);		// Get a read read1 in the graph.
		vector<UINT64> r1MPList = m_dataset->getMatePairList(read1);
		for(UINT64 j = 0; j < r1MPList.size(); j++)		// For each mate-pair read1
		{
			Read * read2 = m_dataset->at(r1MPList[j]);
			// CP: read1 and read2 are paired-ends
			if(read1->getReadID() > read2->getReadID()) 		// To avoid finding path twice
				continue;
			// CP: ignore this pair if their dataset has an average insert size of 0
			//if(meanOfInsertSizes.at(read1->getMatePairList()->at(j).datasetNumber) == 0)
			//	continue;
			//Orientation is F/R so its always 2
			// a path is defined by a vector of edges,
			vector <Edge *> copyOfPath;
			// This list of flags is used to mark the edges that are common in all paths found.
			vector <UINT64> copyOfFlags;
			bool findPaths=false;
			findPaths = findPathBetweenMatepairs(read1, read2, 2,
				m_dataset->getDataSetNumber(r1MPList[j]), copyOfPath, copyOfFlags);
			if(findPaths == true)
			{
				// Matepair on different edge
				if(copyOfPath.size() == 0)
					noPathsFound++;
				else
					pathsFound++;
			}
			else // Mate pair on the same edge
				mpOnSameEdge++;

			if(copyOfPath.size() > 1)	// Path found
			{
				for(UINT64 k = 0; k < copyOfFlags.size(); k++)
				{
					// edge at k and k+1 is supported. We need to add it in the list if not present. If already the pair of edges present then increase the support
					if(copyOfFlags.at(k) == 1)
					{
						UINT64 l;
						for(l = 0; l < listOfPairedEdges_local->size(); l++)
						{
							// already present in the list
							if(listOfPairedEdges_local->at(l).edge1 == copyOfPath.at(k) &&
									listOfPairedEdges_local->at(l).edge2 == copyOfPath.at(k+1))
							{
								listOfPairedEdges_local->at(l).uniqSupport++;	// only increase the support
								break;
							}
							// already present in the list
							else if(listOfPairedEdges_local->at(l).edge2->getReverseEdge() == copyOfPath.at(k)
									&& listOfPairedEdges_local->at(l).edge1->getReverseEdge() == copyOfPath.at(k+1))
							{
								listOfPairedEdges_local->at(l).uniqSupport++;	// only increase the support
								break;
							}
						}
						if(l == listOfPairedEdges_local->size()) // not present in the list
						{
								// add in the list with support 1
							if(copyOfPath.at(k)->getSourceRead()->getReadID() != copyOfPath.at(k)->getDestinationRead()->getReadID()
									|| copyOfPath.at(k+1)->getSourceRead()->getReadID() != copyOfPath.at(k+1)->getDestinationRead()->getReadID())
								// do not want to add support between edge (a,a) and (a,a)
							{
								pairedEdges newPair;
								newPair.edge1 = copyOfPath.at(k);
								newPair.edge2 = copyOfPath.at(k+1);
								newPair.uniqSupport = 1;
								newPair.isFreed = false;
								listOfPairedEdges_local->push_back(newPair);
							}
						}
					}
				}
			}
		}
	}//End of For loop

	#pragma omp parallel num_threads(p_ThreadPoolSize)
	{
		UINT64 threadID = omp_get_thread_num();
		FILE_LOG(logINFO) << "Thread "<< threadID << ": Edge Pairs = " << listOfPairedEdges_localList[threadID]->size()<<endl;
		UINT64 numReadsPerThread = (m_dataset->size()/p_ThreadPoolSize);
		UINT64 startIndex = (numReadsPerThread*threadID) +1;
		UINT64 endIndex = m_dataset->size();
		if((threadID+1)<p_ThreadPoolSize)
			endIndex = (numReadsPerThread*(threadID+1));
		vector <pairedEdges> listOfPairedEdges_lFinal;
		for(UINT64 i = 0; i < listOfPairedEdges_localList.size(); i++)	//Check each partial edge pair list
		{
			for(UINT64 k = 0; k < listOfPairedEdges_localList[i]->size(); k++) // for each edge pair
			{
				//get common read of pair
				UINT64 comReadID = listOfPairedEdges_localList[i]->at(k).edge1->getDestinationRead()->getReadID();
				if(comReadID>=startIndex && comReadID<endIndex) //If true add to local list
				{
					UINT64 l;
					for(l = 0; l < listOfPairedEdges_lFinal.size(); l++)
					{
						// already present in the list
						if(listOfPairedEdges_lFinal.at(l).edge1 == listOfPairedEdges_localList[i]->at(k).edge1 &&
								listOfPairedEdges_lFinal.at(l).edge2 == listOfPairedEdges_localList[i]->at(k).edge2)
						{
							listOfPairedEdges_lFinal.at(l).uniqSupport+=listOfPairedEdges_localList[i]->at(k).uniqSupport;	// only increase the support
							break;
						}
						// already present in the list as reverse
						else if(listOfPairedEdges_lFinal.at(l).edge2->getReverseEdge() == listOfPairedEdges_localList[i]->at(k).edge1
								&& listOfPairedEdges_lFinal.at(l).edge1->getReverseEdge() == listOfPairedEdges_localList[i]->at(k).edge2)
						{
							listOfPairedEdges_lFinal.at(l).uniqSupport+=listOfPairedEdges_localList[i]->at(k).uniqSupport;	// only increase the support
							break;
						}
					}
					if(l == listOfPairedEdges_lFinal.size()) // not present in the list; add to global list
						listOfPairedEdges_lFinal.push_back(listOfPairedEdges_localList[i]->at(k));
				}
			}
		}
		#pragma omp critical
		{
			listOfPairedEdges.reserve(listOfPairedEdges.size() + listOfPairedEdges_lFinal.size() );
			listOfPairedEdges.insert(listOfPairedEdges.end(), listOfPairedEdges_lFinal.begin(),listOfPairedEdges_lFinal.end());
		}
	}
	for(UINT64 i = 0; i < p_ThreadPoolSize; i++)	// for each read in the thread delete the local lists
	{
		delete listOfPairedEdges_localList[i];
	}

	FILE_LOG(logINFO) << "No paths found between " << noPathsFound << " matepairs that are on different edge.\n";
	FILE_LOG(logINFO) << "Paths found between " << pathsFound << " matepairs that are on different edge.\n";
	FILE_LOG(logINFO) << "Total matepairs on different edges " << pathsFound+ noPathsFound << '\n';
	FILE_LOG(logINFO) << "Total matepairs on same edge " << mpOnSameEdge << '\n';
	FILE_LOG(logINFO) << "Total matepairs " << pathsFound+noPathsFound+mpOnSameEdge << '\n';

	FILE_LOG(logINFO)<<"List of pair edges:"<<listOfPairedEdges.size()<<endl;

	UINT64 pairsOfEdgesMerged = 0;

	//Sorting ensures we are merging edge pairs with higher support first
	sort(listOfPairedEdges.begin(), listOfPairedEdges.end());

	for(UINT64 i = 0; i<listOfPairedEdges.size(); i++)
	{
		if(listOfPairedEdges.at(i).isFreed == false && listOfPairedEdges.at(i).uniqSupport >= minUinqSupport &&
				listOfPairedEdges.at(i).edge1->getEdgeLength()>=minSizeToBeShortBranch &&
				listOfPairedEdges.at(i).edge2->getEdgeLength()>=minSizeToBeShortBranch)
		{
			pairsOfEdgesMerged++;
			FILE_LOG(logINFO) << setw(4) << i + 1 << " Merging (" << setw(10) << listOfPairedEdges.at(i).edge1->getSourceRead()->getReadID()
					<< "," << setw(10) <<listOfPairedEdges.at(i).edge1->getDestinationRead()->getReadID() << ") Length: "
					<< setw(8) << listOfPairedEdges.at(i).edge1->getOverlapOffset() + listOfPairedEdges.at(i).edge1->getDestinationRead()->getReadLength() << " Flow: " << setw(3)
					<< listOfPairedEdges.at(i).edge1->m_flow << " and (" << setw(10) << listOfPairedEdges.at(i).edge2->getSourceRead()->getReadID()
					<< "," << setw(10) << listOfPairedEdges.at(i).edge2->getDestinationRead()->getReadID() << ") Length: "
					<< setw(8) << listOfPairedEdges.at(i).edge2->getOverlapOffset() + listOfPairedEdges.at(i).edge2->getDestinationRead()->getReadLength() << " Flow: " << setw(3) << listOfPairedEdges.at(i).edge2->m_flow
					<< " are supported " << setw(4) << listOfPairedEdges.at(i).uniqSupport<<" times"<< '\n';

			Edge * e1f = listOfPairedEdges.at(i).edge1, *e1r = listOfPairedEdges.at(i).edge1->getReverseEdge();
			Edge * e2f = listOfPairedEdges.at(i).edge2, *e2r = listOfPairedEdges.at(i).edge2->getReverseEdge();

			merge2Edges(listOfPairedEdges.at(i).edge1, listOfPairedEdges.at(i).edge2);

			// BH: once we merge edge1 and edge2. We make sure that we do not try to merge these edges again.
			// We mark all the pair of edges that contains edge1 and edge2 or their reverse edges as freed.
			#pragma omp parallel for schedule(dynamic) num_threads(p_ThreadPoolSize)
			for(UINT64 j = i + 1; j<listOfPairedEdges.size(); j++)
			{
				if( listOfPairedEdges.at(j).edge1 == e1f || listOfPairedEdges.at(j).edge1 == e1r
						|| listOfPairedEdges.at(j).edge1 == e2f || listOfPairedEdges.at(j).edge1 == e2r )
					listOfPairedEdges.at(j).isFreed = true;
				if( listOfPairedEdges.at(j).edge2 == e1f || listOfPairedEdges.at(j).edge2 == e1r
						|| listOfPairedEdges.at(j).edge2 == e2f || listOfPairedEdges.at(j).edge2 == e2r )
					listOfPairedEdges.at(j).isFreed = true;
			}
		}
	}
	FILE_LOG(logINFO) << pairsOfEdgesMerged <<" Pairs of Edges merged out of " << listOfPairedEdges.size() << " supported pairs of edges" <<endl;
	CLOCKSTOP;
	return pairsOfEdgesMerged;

}

/**********************************************************************************************************************
	This function finds all the paths between two reads in the overlap graph.
**********************************************************************************************************************/

// CP: inputs: read1 and read2 are the pair, orient of the pair is defined in the MPlist struct, datasetNumber retrievs the mean and SD of insert size of the dataset
// CP: outputs: copyOfPath is a path between the two reads and copyOfFlags indicates whether the connection between this edge
// to the next edge is supported by all possible paths or not
// CP: return true if valid paths are found between them.
bool OverlapGraph::findPathBetweenMatepairs(const Read * read1, const Read * read2,
		UINT8 orient, UINT8 datasetNumber, vector <Edge *> &copyOfPath, vector <UINT64> &copyOfFlags)
{
	UINT64 pathFound = 0;			// CP: the total number of paths found between the two reads
	// CP: two variables passed to the exploreGraph function for return values
	vector <Edge *> firstPath;
	vector <UINT64> flags;

	vector<t_edge_loc_pair> listRead1, listRead2;

	// mate pair orientation
	// 0 = 00 means the reverse of this read and the reverse of the second read are matepairs.
	// 1 = 01 means the reverse of this read and the forward of the second read are matepairs.
	// 2 = 10 means the forward of this read and the reverse of the second read are matepairs.
	// 3 = 11 means the forward of this read and the forward of the second read are matepairs.

	listRead1 = (orient == 2 || orient == 3) ? read1->getFwdEdges() : read1->getBwdEdges();

	// CP: Explain how this forward-forward/forward-reverse is related to mate-pair orientation listed above?
	// BH: if the first read is a forward read and the second read is a reverse read then we have to look for the forward read
	// in the first edge and the reverse read in the second edge.

	// CP: use this, if the matepairs are forward - forward
	// CP: Ted, we need to add this to the config file, instead of hard-coding it here
	//listRead2 = (orient == 1 || orient == 3) ? read2->getListOfEdgesForward() : read2->getListOfEdgesReverse(); // if the matepairs are forward - forward
	//locationOnEdgeRead2 = (orient == 1 || orient == 3) ? read2->getLocationOnEdgeForward() : read2->getLocationOnEdgeReverse(); // if the matepairs are forward - forward

	// CP: use this, if the matepairs are forward - reverse
	listRead2 = (orient == 0 || orient == 2) ? read2->getBwdEdges() : read2->getFwdEdges(); // if the matepairs are forward - reverse

	// CP: return false, if the two reads are not part of any edge or they are on the same edge
	if(listRead1.size()==0 || listRead2.size()==0)	// Not  present in any edges.
	{
		return false;
	}
	else
	{
		for(UINT64 i = 0; i < listRead1.size(); i++)
		{
			for(UINT64 j = 0; j < listRead2.size(); j++)
			{
				Edge * firstEdge = listRead1[i].first;
				Edge * lastEdge = listRead2[j].first;
				if(firstEdge ==lastEdge || firstEdge == lastEdge->getReverseEdge()) // Both are on the same edge
				{
					return false;
				}
			}
		}
	}

	for(UINT64 i = 0; i < listRead1.size(); i++)	// If not on the same edge.
	{
		for(UINT64 j = 0; j < listRead2.size(); j++)
		{
			Edge * firstEdge = listRead1[i].first;
			Edge * lastEdge = listRead2[j].first;
			if(firstEdge !=lastEdge && firstEdge != lastEdge->getReverseEdge()) //not on the same edge
			{

				int r1Index=listRead1[i].second;
				int r2Index=listRead2[j].second;
				UINT64 r1Offset=0,r2Offset=0;

				if(r1Index < static_cast<int>(listRead1[i].first->getListofReadsSize() - 1))  // Not the last read on edge
					r1Offset=listRead1[i].first->getInnerOverlapSum(0,r1Index+1);
				else
					r1Offset=listRead1[i].first->getInnerOverlapSum(0,listRead1[i].first->getListofReadsSize());

				if(r2Index < static_cast<int>(listRead2[j].first->getListofReadsSize() - 1))  // Not the last read on edge
					r2Offset = listRead2[j].first->getInnerOverlapSum(0,r2Index+1);
				else
					r2Offset = listRead2[j].first->getInnerOverlapSum(0,listRead2[j].first->getListofReadsSize());



				// distance from the end of the read1 to the end of the last read of the edge,
				int distanceOnFirstEdge = firstEdge->getOverlapOffset() - r1Offset - read1->getReadLength();
				// distance from the beginning of the first read of the edge to the beginning of the read2
				int distanceOnLastEdge = r2Offset;
				// the two reads can't be too far apart within their own edges
				if((distanceOnFirstEdge + distanceOnLastEdge) <
						(m_dataset->getDataSetInfo()->at(datasetNumber).avgInnerDistance + insertSizeRangeSD * m_dataset->getDataSetInfo()->at(datasetNumber).avgInnerDistanceSD))
				{
					// from first edge  try to find a path to the last edge.
					UINT64 newPaths = 0;
					vector <Edge *> listOfEdges;
					vector <UINT64> pathLengths;
					exploreGraph(firstEdge, lastEdge, distanceOnFirstEdge, distanceOnLastEdge, datasetNumber,
							0, firstPath, flags,newPaths,listOfEdges,pathLengths);

					if(newPaths > 0)	// If some path found.
					{
						pathFound+=newPaths;	// How many paths found.
						if(copyOfPath.empty())	// Path found for the first time.
						{
							for(UINT64 k = 0; k < firstPath.size(); k++)	// Add the paths in the list.
								copyOfPath.push_back(firstPath.at(k));
							for(UINT64 k = 0; k < firstPath.size() - 1; k++)	// Also add the flag.
								copyOfFlags.push_back(flags.at(k));
						}
						else		// Not the first path
						{
							// CP: compare each supported pair of edge in copyOfPath with all pairs in firstpath
							// CP: if this supported pair of edge in copyOfPath is still supported by a pair in firstPath, it remains supported.
							// CP: if this supported pair of edge in copyOfPath is not supported by any pair in firstPath, it is changed to not supported.
							UINT64 k , l;
							// Check if the previous path contains the same pair of edges adjacent to the new path and has flag 1.
							for( k = 0; k< copyOfPath.size() - 1; k++)
							{
								for( l = 0; l < firstPath.size() - 1; l++)
								{
									if(copyOfPath.at(k) == firstPath.at(l) &&  copyOfPath.at(k+1) == firstPath.at(l+1) && flags.at(l) == 1)
										break;
								}
								if(l == firstPath.size() - 1)	// This pair is not uniquely supported
									copyOfFlags.at(k) = 0;
							}
						}
					}
				}
			}
		}
	}
	return true;
}

/**********************************************************************************************************************
	Explore all paths starting from firstEdge and tries to reach lastEdge using depth first search.
	Depth is limited to 100 to reduce running time
**********************************************************************************************************************/
// CP: this a recursive function
// CP: inputs: firstEdge and lastEdge in the path, distanceOnFirstEdge and distanceOnLastEdge are the relative position of the paired reads on the two edges
// CP: inputs: datasetNumber retrieves the mean and SD of insert size of the dataset
// CP: outputs: level is the level of the depth first search, firstPath and flags are ???
// BH: when we find the first path we save it and flag is used to indicate that all the edges in the first path are supported.
//     When we find another path, we unmark the flag for the pair of edges in the first path that are not present in the second path and so on.
// CP: return the number of paths found
void OverlapGraph::exploreGraph(Edge* firstEdge, Edge * lastEdge, int distanceOnFirstEdge,
		int distanceOnLastEdge, UINT64 datasetNumber, UINT64 level, vector <Edge *> &firstPath, vector <UINT64> &flags,
		UINT64 &pathFound, vector <Edge *> &listOfEdges, vector <UINT64> &pathLengths)
{
	// length of the path from the beginning of the source read of the first edge to the current edge's source read's beginning
	INT64 meanDist = m_dataset->getDataSetInfo()->at(datasetNumber).avgInnerDistance;
	INT64 meanSD = m_dataset->getDataSetInfo()->at(datasetNumber).avgInnerDistanceSD;
	if(level == 0)
	{
		// clear the variables and resize their capacity to free memory
		pathFound = 0;
		firstPath.resize(0);
		flags.resize(0);
		listOfEdges.resize(0);
		pathLengths.resize(0);
	}
	else
	{
		listOfEdges.resize(level);
		pathLengths.resize(level);
	}
	// BH: we do not go deeper than EXPLORE_DEPTH levels. We can put this in the config file.
	// CP: when reaching the maximum depth, return 0 path found and exit the recursive call loop
	if(level > EXPLORE_DEPTH) return; // Do not go very deep.


	if(level == 0)
	{
		// CP: at level 0, start from the end of the first edge
		listOfEdges.push_back(firstEdge);
		pathLengths.push_back(distanceOnFirstEdge);
	}
	else
	{
		if(firstEdge == lastEdge) // Destination found.
		{
			// If we read our destination read, we check if the distance is within 3 sd of the mean.
			// mean - 3*sd can be negative.
			if((INT64)(distanceOnLastEdge + pathLengths.at(level - 1)) >=  (meanDist - insertSizeRangeSD * meanSD)
					&& (INT64)(distanceOnLastEdge + pathLengths.at(level - 1)) <= meanDist + insertSizeRangeSD * meanSD)
			{
				// CP: the path length is within the insert size range
				listOfEdges.push_back(firstEdge);
				pathLengths.push_back(distanceOnLastEdge + pathLengths.at(level - 1));
				pathFound++;
				if(pathFound == 1)	// First path
				{

					for(UINT64 i = 0; i <listOfEdges.size() ; i++) // Copy the path.
						firstPath.push_back(listOfEdges.at(i));

					for(UINT64 i = 0; i <listOfEdges.size() - 1 ; i++) // All adjacent edges in the path are supported.
						flags.push_back(1);
				}
				else		// Not the first path.
				{
					UINT64 i , j;
					for( i = 0; i< firstPath.size() - 1; i++)		// Compare to the first path.
					{
						for( j = 0; j < listOfEdges.size() - 1; j++)
						{
							// A pair of edges adjacent in the first path is also adjacent in this path. We keep the support.
							if(firstPath.at(i) == listOfEdges.at(j) &&  firstPath.at(i+1) == listOfEdges.at(j+1))
								break;
						}
						// A pair of adjacent edge in the first path is not adjacent in this path. So this pair is not supported anymore. So we clear the flag.
						if(j == listOfEdges.size() - 1)
							flags.at(i) = 0;
					}
				}
				// CP: when found a valid path, return 1 path found and exit the recursive call
				return;
			}
			else		// add the new edge in the path
			{
				// CP: this length of this path is not within the valid range of the insert size and keep going deeper
				listOfEdges.push_back(firstEdge);
				pathLengths.push_back( distanceOnFirstEdge + pathLengths.at(level - 1) );
			}
		}
		else	// add the new edge in the path.
		{
			// CP: this path has not reached the destination edge yet and keep going deeper
			listOfEdges.push_back(firstEdge);
			pathLengths.push_back( distanceOnFirstEdge + pathLengths.at(level - 1) );
		}
	}
	auto it = m_graph->find(firstEdge->getDestinationRead()->getReadID());
	if(it != m_graph->end())
	{
		for(UINT64 i = 0 ; i < it->second->size(); i++ )	// go deeper in the graph to reach destination edge
		{
			// CP: going to each edge that the last read of the firstEdge is connected to
			Edge * nextEdge =  it->second->at(i);
			// CP: if the two edges are compatible and the current path is not already too long, then call self again
			if(matchEdgeType(firstEdge, nextEdge) && (INT64)pathLengths.at(level) < meanDist + insertSizeRangeSD * meanSD)
				exploreGraph(nextEdge, lastEdge, nextEdge->getOverlapOffset(),
						distanceOnLastEdge, datasetNumber, level+1, firstPath, flags,pathFound,listOfEdges,pathLengths);
		}
	}
}
/**********************************************************************************************************************
	Original Scaffolder function.
	*******************************************************************************************************************/
UINT64 OverlapGraph::scaffolder(void)
{
	CLOCKSTART;
	UINT64 pairsOfEdgesMerged = 0;
	vector<pairedEdges> listOfPairedEdges;

	// listOfCompositeEdges contains all the composite edges
	vector<Edge *> *listOfCompositeEdges = new vector<Edge *>;

	for(UINT64 i = 1; i<= m_dataset->size(); i++)
	{
		auto it = m_graph->find(i);
		if(it != m_graph->end() && !it->second->empty()) // if this read has some edge(s) going out of it (since now the graph is directed)
		{
			t_edge_vec *edgeList = it->second;
			for(UINT64 j = 0; j < edgeList->size(); j++)
			{
				Edge * edge = edgeList->at(j);
				if(edge->isListofReads() && !(edge->getLastOverlapOffset()==0)) // composite edge
				{
					listOfCompositeEdges->push_back(edge);	// List of all the composite edges.
				}
			}
		}
	}
	for(UINT64 i = 0 ; i < listOfCompositeEdges->size(); i++) // For each composite edge in the graph
	{
		Edge *edge1 = listOfCompositeEdges->at(i);
		// Find the list of other edges that share unique matepairs with the current edge. Only check one of the endpoints of the current edge.
		vector<Edge *> *listOfFeasibleEdges = getListOfFeasibleEdges(edge1);
		for(UINT64 j = 0; j < listOfFeasibleEdges->size(); j++ ) // Check the current edge vs the list of edges for support for scaffolder
		{
			Edge *edge2 =listOfFeasibleEdges->at(j);
			INT64 gapDistance;
			INT64 support = checkForScaffold(edge1,edge2,gapDistance); // check the support and distance
			if(support>0)	// If there are support then add the current pair in the list.
			{
				pairedEdges newPair;
				newPair.edge1 = edge1;
				newPair.edge2 = edge2;
				newPair.uniqSupport = support;
				newPair.distance = gapDistance;
				newPair.isFreed = false;
				listOfPairedEdges.push_back(newPair);

			}
		}
		delete listOfFeasibleEdges;	// Free the memory
	}
	delete listOfCompositeEdges;

	sort(listOfPairedEdges.begin(), listOfPairedEdges.end());		// Sort the list according to support.

	for(UINT64 i = 0; i < listOfPairedEdges.size(); i++)
	{
		if(listOfPairedEdges.at(i).isFreed == false && listOfPairedEdges.at(i).uniqSupport >= minUinqSupport &&
				listOfPairedEdges.at(i).edge1->getEdgeLength()>=minSizeToBeShortBranch &&
				listOfPairedEdges.at(i).edge2->getEdgeLength()>=minSizeToBeShortBranch)
		{
			pairsOfEdgesMerged++;
			FILE_LOG(logINFO) << setw(4) << i + 1 << " (" << setw(10) << listOfPairedEdges.at(i).edge1->getSourceRead()->getReadID()
					<< "," << setw(10) <<listOfPairedEdges.at(i).edge1->getDestinationRead()->getReadID() << ") Length: "
					<< setw(8) << listOfPairedEdges.at(i).edge1->getOverlapOffset() << " Flow: " << setw(3) << listOfPairedEdges.at(i).edge1->m_flow
					<< " and (" << setw(10) << listOfPairedEdges.at(i).edge2->getSourceRead()->getReadID() << "," << setw(10)
					<< listOfPairedEdges.at(i).edge2->getDestinationRead()->getReadID() << ") Length: " << setw(8)
					<< listOfPairedEdges.at(i).edge2->getOverlapOffset() << " Flow: " << setw(3) << listOfPairedEdges.at(i).edge2->m_flow << " are supported "
					<< setw(4) << listOfPairedEdges.at(i).uniqSupport << " times. Average distance: "<< setw(4) << listOfPairedEdges.at(i).distance << '\n';
			Edge * e1f = listOfPairedEdges.at(i).edge1, *e1r = listOfPairedEdges.at(i).edge1->getReverseEdge();
			Edge * e2f = listOfPairedEdges.at(i).edge2, *e2r = listOfPairedEdges.at(i).edge2->getReverseEdge();
			mergeEdgesDisconnected(listOfPairedEdges.at(i).edge1, listOfPairedEdges.at(i).edge2,listOfPairedEdges.at(i).distance);		// Merge the edges.
			// BH: if an edge is merged already, I make sure that I will not try to merge it again with other edges.
			for(UINT64 j = i + 1; j<listOfPairedEdges.size(); j++)
			{
				if( listOfPairedEdges.at(j).edge1 == e1f || listOfPairedEdges.at(j).edge1 == e1r
						|| listOfPairedEdges.at(j).edge1 == e2f || listOfPairedEdges.at(j).edge1 == e2r )
					listOfPairedEdges.at(j).isFreed = true;
				if( listOfPairedEdges.at(j).edge2 == e1f || listOfPairedEdges.at(j).edge2 == e1r
						|| listOfPairedEdges.at(j).edge2 == e2f || listOfPairedEdges.at(j).edge2 == e2r )
					listOfPairedEdges.at(j).isFreed = true;
			}
		}
	}
	CLOCKSTOP;
	return pairsOfEdgesMerged;
}

/**********************************************************************************************************************
 	 This functions returns a list of edges that might be joined to "edge"
 	 CP: For the input edge, find all the feasible edges that are linked with the input edge by a pair of edge-unique reads with appropriate distance
 	 CP: edge-unique reads are reads that are only present on one edge
***********************************************************************************************************************/
vector<Edge *> * OverlapGraph::getListOfFeasibleEdges(const Edge *edge)
{
	// We want to find if there are other edges that share matepairs. current edge (u,v) we check the matepairs near the node v. That's why we took the reverse edge.
	// CP: why do you only check near v, not u?
	// BH: Here we are checkin if we can merge (u,v) followed by another edge. That why we only take the reads near v.
	// BH: At some point we will call this function with the reverse edge (v,u) as well. In that case we will look at the reads near vertex u.
	Edge * rEdge=edge->getReverseEdge();
	vector<Edge *> * feasibleListOfEdges = new vector<Edge *>;
	UINT64 dist = 0;
	for(UINT64 i = 0; i <rEdge->getListofReadsSize(); i++) // for each read in the edge
	{
		// CP: dist is the distance from the beginning of rEdge to the beginning of this read.
		dist+=rEdge->getInnerOverlapOffset(i);	// offset. We do not have to go much deeper. we need to make sure that we atleast go upto the longest insert size.
		// CP: if this read is too far inside rEdge, then no need to go further inside. Stop and return the list
		if(dist > 2*longestMeanOfInsertDistance)	// longest insert size mean
			 break;
		// CP: retrieve the current read: r1
		UINT64 mp1=rEdge->getInnerReadID(i); // mate pair 1
		Read *r1 = m_dataset->at(mp1); // read1

		if(r1->getFwdEdges().size() == 1) // only present in this current edge
		{
			// CP: r1 is a unique read that is present in this current edge only. Ignore non-unique reads
			// CP: this for loop considers r1 and forward edges of r2: r1->.......r2->
			vector<UINT64> r1MPList = m_dataset->getMatePairList(r1);
			for(UINT64 j = 0; j < r1MPList.size(); j++) // for each matepair of current read1
			{
				// CP: r2 is the paired read of r1
				UINT64 mp2 = r1MPList[j]; // matepair 2
				Read* r2 = m_dataset->at(mp2); // read2
				UINT8 orient = 2;
				vector<t_edge_loc_pair > list;

				if(orient == 0  || orient == 2)
					list = r2->getBwdEdges(); // edge contain read backward
				else
					list = r2->getFwdEdges(); // edge contain read forward


				// CP: use read2 if it's on one and only one edge and it's not on the input forward/reverse edge
				if(list.empty() || list.size() > 1 || list[0].first == edge ||
						list[0].first == edge->getReverseEdge())
					// Not present uniquely on the edge
					continue;
				//AB: We do a summation to get the offset location of the read on the edge.
				//AB: See detailed justification of this in calculateMeanAndSdOfInnerDistance()
				int r2Index=list[0].second;				//There is only one index 0 as the read is only on one edge
				UINT32 r2Offset=0;

				if(r2Index < static_cast<int>(list[0].first->getListofReadsSize() - 1))  // Not the last read on edge
					r2Offset=list[0].first->getInnerOverlapSum(0,r2Index+1);
				else
					r2Offset=list[0].first->getInnerOverlapSum(0,list[0].first->getListofReadsSize());

				// AB: Check distance is adequate
				if(r2Offset > 2*longestMeanOfInsertDistance)
					// not within the distance of longest inner distance size.
					continue;

				UINT64 k;
				for(k = 0; k<feasibleListOfEdges->size();k++)		// add in the list of feasible edges. This list is expected to be small.
				{
					if(feasibleListOfEdges->at(k) == list[0].first)	// already in the list.
						break;
				}
				if(k == feasibleListOfEdges->size())	// Not present in the list.
				{
					feasibleListOfEdges->push_back(list[0].first);	// insert the edge in the list.
				}
			}
		}
	}

	return feasibleListOfEdges;	// list of edges that might be joined with the current edge for scaffolding
}

/**********************************************************************************************************************
	Calculate Mean and Standard Deviation of insert size of each dataset.
**********************************************************************************************************************/
bool OverlapGraph::calculateMeanAndSdOfInnerDistance(void)
{
	CLOCKSTART;
	longestMeanOfInsertDistance = 0;
	vector<INT64> *innerDistSizes = new vector<INT64>;

	for(UINT64 d = 0; d < m_dataset->getDataSetInfo()->size(); d++)	// For each dataset.
	{
		if(m_dataset->getDataSetInfo()->at(d).isPaired)
		{
			FILE_LOG(logINFO) << "Calculating mean and SD of dataset: " << d << endl;
			innerDistSizes->clear();
			#pragma omp parallel for schedule(dynamic) num_threads(p_ThreadPoolSize)
			for(UINT64 i = m_dataset->getDataSetInfo()->at(d).r1Start; i <= m_dataset->getDataSetInfo()->at(d).r1End; i++)	// for each read in the dataset
			{
				Read *read1 = m_dataset->at(i), *read2;	// Get a read read1 in the graph.
				vector<UINT64> r1MPList = m_dataset->getMatePairList(read1);
				for(UINT64 j = 0; j < r1MPList.size(); j++) 	// For each mate-pair read2
				{
					if(m_dataset->getDataSetNumber(r1MPList[j]) == d)		// if it is in the current dataset.
					{
						vector<t_edge_loc_pair> listOfEdgesRead1, listOfEdgesRead2;
						// All the edges that contain forward string of read1
						listOfEdgesRead1 = read1->getFwdEdges();

						read2 = m_dataset->at(r1MPList[j]);	// Get the read2 object.
						// All the edges that contain reverse string of read2
						listOfEdgesRead2 = read2->getBwdEdges();

						for(UINT64 k = 0; k < listOfEdgesRead1.size(); k++)				// For each edge containing read1
						{
							for(UINT64 l = 0; l < listOfEdgesRead2.size(); l++)			// For each edge containing read2
							{
								/*
								 * AB: the listOfEdgesRead1 has index of the read in the edge not offset location so, we need to do a little dance to get the offset values
								 * AB: First get the indexes. fairly easy
								 * AB: Next get the summation of the offset till index+1 as the read and offset list does not include the source and destination
								 * and the first read in the list is the second read in the alignment.
								 */
								int r1Index=listOfEdgesRead1.at(k).second;
								int r2Index=listOfEdgesRead2.at(l).second;
								UINT32 r1Offset=0,r2Offset=0;

								if(r1Index < static_cast<int>(listOfEdgesRead1.at(k).first->getListofReadsSize() - 1))  // Not the last read on edge
									r1Offset=listOfEdgesRead1.at(k).first->getInnerOverlapSum(0,r1Index+1);
								else
									r1Offset=listOfEdgesRead1.at(k).first->getInnerOverlapSum(0,listOfEdgesRead1.at(k).first->getListofReadsSize());

								if(r2Index < static_cast<int>(listOfEdgesRead2.at(l).first->getListofReadsSize() - 1))  // Not the last read on edge
									r2Offset = listOfEdgesRead2.at(l).first->getInnerOverlapSum(0,r2Index+1);
								else
									r2Offset = listOfEdgesRead2.at(l).first->getInnerOverlapSum(0,listOfEdgesRead2.at(l).first->getListofReadsSize());

								INT64 mpDist = (INT64)(r2Offset - (r1Offset + read1->getReadLength()));
								// Both reads are on the same edge
								if(listOfEdgesRead1.at(k).first == listOfEdgesRead2.at(l).first && mpDist>0 && mpDist < MAX_INNER_DIST_TRESH)
								{
									// Distance between the two edges is less than 1000. Some times some mate pairs are far off the actual value.
									// This upper bound is used to get only good mate pairs.
									// We may need to change the threshold for datasets with longer insert size.
									// Insert the distance between the  reads in the list.
										#pragma omp critical(insertDist)
										{
											innerDistSizes->push_back(mpDist);
										}
								}
							}
						}
					}
				}
			}
		}

		double sum = 0, variance=0, mean=0, sdiv=0, smpSize=0;

		if(innerDistSizes->size() == 0) // If no insert size found
		{
			FILE_LOG(logINFO) << "No insert-size found for dataset: " << d << endl;
			m_dataset->getDataSetInfo()->at(d).avgInnerDistance=0;
			m_dataset->getDataSetInfo()->at(d).avgInnerDistanceSD=0;
			continue;
		}
		smpSize=innerDistSizes->size();

		for(UINT64 i = 0; i < innerDistSizes->size(); i++)
			sum += innerDistSizes->at(i);			// Calculate the sum of all the insert sizes of the current dataset.

		mean = sum/smpSize;
		m_dataset->getDataSetInfo()->at(d).avgInnerDistance=mean;	// Calculate the mean. In push it in the variable of the OverlapGraph object.


		for(UINT64 i = 0; i < innerDistSizes->size(); i++)	// Calculate the variance of the current dataset.
			variance += (mean - innerDistSizes->at(i)) * (mean - innerDistSizes->at(i));

		sdiv=sqrt(variance/innerDistSizes->size());
		m_dataset->getDataSetInfo()->at(d).avgInnerDistanceSD=sdiv;  // Calculate and insert the standard deviation.

		// Print the values of the current dataset.
		FILE_LOG(logINFO) << "Mean set to: " << mean << '\n';
		FILE_LOG(logINFO) << "SD set to: " << sdiv << '\n';
		FILE_LOG(logINFO) << "Reads on same edge: " << innerDistSizes->size() << endl;

		if((INT64)longestMeanOfInsertDistance < mean)
		{
			longestMeanOfInsertDistance = mean;
		}

	}
	FILE_LOG(logINFO) << "Mean of longest inner distance : " << longestMeanOfInsertDistance << endl;

	delete innerDistSizes;
	CLOCKSTOP;
	return true;
}

INT64 OverlapGraph::checkForScaffold(const Edge *edge1, const Edge *edge2, INT64 &averageGapDistance)
{
	UINT64 dist = 0;
	INT64 support = 0, oppose=0;
	averageGapDistance = 0;
	Edge *rEdge1 = edge1->getReverseEdge();
	vector<t_edge_loc_pair> listRead1, listRead2;		//  This is the lists of edges that contain read1 and read2
	// CP: listOfReads contains all the reads in the end section of edge1
	vector<UINT64> listOfReads;
	for(UINT64 i = 0; i <rEdge1->getListofReadsSize(); i++)
	{
		dist+=rEdge1->getInnerOverlapOffset(i);
		if(dist>2*longestMeanOfInsertDistance)
			 break;
		listOfReads.push_back(rEdge1->getInnerReadID(i));
	}
	for(UINT64 i = 0; i < listOfReads.size(); i++)
	{
		Read *read1 = m_dataset->at(listOfReads.at(i));
		vector<UINT64> r1MPList = m_dataset->getMatePairList(read1);
		for(UINT64 j = 0; j < r1MPList.size(); j++)// For each matepair
		{
			Read *read2 = m_dataset->at(r1MPList[j]);	// Get the read object of the matepair.
			//if(read1->getReadNumber() > read2->getReadNumber()) // To avoid duplicate computation
				//	continue;
			UINT64 orient = 2;						//Orientation for illumina is always 2

			UINT64 datasetNumber = m_dataset->getDataSetNumber(r1MPList[j]);		// Get the dataset number

			// 0 = 00 means the reverse of r1 and the reverse of r2 are matepairs.
			// 1 = 01 means the reverse of r1 and the forward of r2 are matepairs.
			// 2 = 10 means the forward of r1 and the reverse of r2 are matepairs.
			// 3 = 11 means the forward of r1 and the forward of r2 are matepairs.
			// To calculate distance of forward read, flip the read and get the location of the offset.
			listRead1 = (orient == 0 || orient == 1) ? read1->getFwdEdges() : read1->getBwdEdges();
			// To calculate distance of reverse read, flip the read and get the location of the offset.
			listRead2 = (orient == 0 || orient == 2) ? read2->getBwdEdges() : read2->getFwdEdges();

			//if either of the lists are empty or not uniquely mapped DONOT proceed
			if(listRead1.size() != 1 || listRead2.size() != 1)
				continue;
			/*
			 * AB: the listOfRead1 has index of the read in the edge not offset location so, we need to do a little dance to get the offset values
			 * AB: First get the indexes. fairly easy
			 * AB: Next get the summation of the offset till index+1 as the read and offset list does not include the source and destination
			 * and the first read in the list is the second read in the alignment.
			 */
			int r1Index=listRead1[0].second;
			int r2Index=listRead2[0].second;
			UINT32 r1Offset=0,r2Offset=0;

			if(r1Index < static_cast<int>(listRead1[0].first->getListofReadsSize() - 1))  // Not the last read on edge
				r1Offset=listRead1[0].first->getInnerOverlapSum(0,r1Index+1);
			else
				r1Offset=listRead1[0].first->getInnerOverlapSum(0,listRead1[0].first->getListofReadsSize());

			if(r2Index < static_cast<int>(listRead2[0].first->getListofReadsSize() - 1))  // Not the last read on edge
				r2Offset = listRead2[0].first->getInnerOverlapSum(0,r2Index+1);
			else
				r2Offset = listRead2[0].first->getInnerOverlapSum(0,listRead2[0].first->getListofReadsSize());

			// Only consider if the distance is less than mean+3*SD
			if(listRead1[0].first == edge1->getReverseEdge()
					&& listRead2[0].first == edge2
					&& (r1Offset + r2Offset) < (m_dataset->getDataSetInfo()->at(datasetNumber).avgInnerDistance + insertSizeRangeSD * m_dataset->getDataSetInfo()->at(datasetNumber).avgInnerDistanceSD) )
				// Both the reads are present on only on edge and the distance is less that mean+3*sd
			{
				dist = r1Offset + r2Offset;
				// if there are already in the same edge, don't do anything
				if(listRead1[0].first == listRead2[0].first ||  listRead1[0].first == listRead2[0].first->getReverseEdge()) // Not on the same edge
					continue;
				averageGapDistance += static_cast<int>(m_dataset->getDataSetInfo()->at(datasetNumber).avgInnerDistance - dist);
				support++;
			}
			else
			{
				oppose++;
			}
		}
	}
	if(support)
		averageGapDistance = (INT64)(averageGapDistance/(INT64)(support));
	return (support - oppose);
}

/**********************************************************************************************************************
	Check if two strings overlap. At least 10 bp must overlap.
	CP: return the length of the overlap. 0 if no overlap
**********************************************************************************************************************/
UINT64 OverlapGraph::findOverlap(const string & string1, const string & string2)
{
	UINT64 minimum = min(string1.length(), string2.length());
	for(UINT64 i = minimum - 1; i >= 10; i--)
	{
		if(string1.substr(string1.length() - i, i) == string2.substr(0,i))
		{
			return i;
		}
	}
	return 0;
}

/**********************************************************************************************************************
	Merge two edges that do not share any node.
**********************************************************************************************************************/
bool OverlapGraph::mergeEdgesDisconnected(Edge *edge1, Edge *edge2, INT64 gapLength)
{
	if(edge1->getDestinationRead()->getReadID() == edge2->getSourceRead()->getReadID() && matchEdgeType (edge1, edge2))
		// If the two edges are already connected. A --->B and B---->C. They share a common node
	{
		merge2Edges(edge1,edge2); // Merge the edges.
		return true;
	}

	// A------>B and C------>D. We need to check if the nodes B and C overlaps or not
	//string string1, string2;
	//string1 = ( edge1->getOrientation() == 1 || edge1->getOrientation() == 3 )
	//		? edge1->getDestinationRead()->getStringForward() : edge1->getDestinationRead()->getStringReverse(); // We will check if the two nodes overlap or not
	//string2 = ( edge2->getOrientation() == 2 || edge2->getOrientation() == 3 )
	//		? edge2->getSourceRead()->getStringForward() : edge2->getSourceRead()->getStringReverse();

	// Find the overlap between B and C. If they do not overlap then the return will be zero. We check for at least 10 bp overlap
	//UINT64 overlapLength = findOverlap(string1,string2);
	UINT64 overlapLength = 0;

	UINT64 overlapOffset1, overlapOffset2;
	Read *read1 = edge1->getSourceRead(), *read2 = edge2->getDestinationRead(); // Get the reads.
	UINT8 orientationForward = mergedEdgeOrientationDisconnected(edge1,edge2); // Orientation of the forward edge based on the orientations of edge1 and edge2
	UINT8 orientationReverse = twinEdgeOrientation(orientationForward);			// Orientation of the reverse edge.

	if(overlapLength == 0) // Strings in the read B and C do not overlap
	{
		// CP: do you insert Ns if they don't overlap? It's important to insert Ns
		// BH: We do not add N here. But when we generate the string from the edges, we insert N's there.
		if(gapLength>0)
		{
			overlapOffset1 = edge1->getDestinationRead()->getReadLength() + gapLength;	// In this case we concatenate the strings in the edges. So the offset is the length of the read B
			overlapOffset2 = edge2->getSourceRead()->getReadLength() + gapLength;		// This is the overlap offset of the reverse edge.
		}
		else //If gapLength<0 we should have string overlap. If not, then the gap length estimation is incorrect and we add a gap of 10
		{
			overlapOffset1 = edge1->getDestinationRead()->getReadLength() + 10;	// In this case we concatenate the strings in the edges. So the offset is the length of the read B
			overlapOffset2 = edge2->getSourceRead()->getReadLength() + 10;		// This is the overlap offset of the reverse edge.
		}
	}
	else	// Strings in the read B and C do overlap
	{
		overlapOffset1 = edge1->getDestinationRead()->getReadLength() - overlapLength; // Overlap offset of the forward edge is taken according to the length of the string in B
		overlapOffset2 = edge2->getSourceRead()->getReadLength() - overlapLength; // overlap offset of the reverse edge is taken according to the length of the string in C
	}


	// CP: merge the forward edge
	UINT64 * listReadsForward = nullptr;		// List of reads in the forward edge.
	UINT32 lFSize=0;

	mergeListDisconnected(edge1, edge2, overlapOffset1, &listReadsForward, lFSize);
	Edge *edgeForward = new Edge(read1,read2,orientationForward, edge1->getOverlapOffset() + edge2->getOverlapOffset() + overlapOffset1,
			listReadsForward, lFSize);

	// CP: merge the reverse edge
	UINT64 *listReadsReverse = nullptr;			// List of reads in the reverse edge.
	UINT32 lRSize=0;
	mergeListDisconnected(edge2->getReverseEdge(),edge1->getReverseEdge(), overlapOffset2, &listReadsReverse, lRSize);
	// BH: lengthReverseEdge is the overlap offset of the new edge.
	UINT64 lengthReverseEdge = edge1->getReverseEdge()->getOverlapOffset() + edge2->getReverseEdge()->getOverlapOffset() + overlapOffset2;
	Edge *edgeReverse = new Edge(read2, read1, orientationReverse, lengthReverseEdge,
			listReadsReverse, lRSize);

//	edgeReverse->makeEdge(read2, read1, orientationReverse, edge1->getReverseEdge()->getOverlapOffset() + edge2->getReverseEdge()->getOverlapOffset() + overlapOffset2,
//	listReadsReverse, listOverlapOffsetsReverse, listOrientationsReverse);

	edgeForward->setReverseEdge(edgeReverse);	// set the pointer of reverse edge
	edgeReverse->setReverseEdge(edgeForward);	// set the pointer of reverse edge


	UINT16 flow = min(edge1->m_flow,edge2->m_flow);	// Take the minimum of the flow from the two original edges.
//	UINT64 coverage = min(edge1->coverageDepth, edge2->coverageDepth);	// not used
	edgeForward->m_flow = flow;	// set the flow in the forward edge.
//	edgeForward->coverageDepth = coverage;

	edgeReverse->m_flow = flow;	// set the flow in the reverse edge.
//	edgeReverse->coverageDepth = coverage;

	//if(flowComputed && flow == 0 && edgeForward->getOverlapOffset() > 1000)
	//{
	//	cout << "Check for error inserting edge between " << edgeForward->getSourceRead()->getReadNumber() << " and
	// " << edgeForward->getDestinationRead()->getReadNumber() << " Length: " << edgeForward->getOverlapOffset() << endl;
	//}
	insertEdge(edgeForward); // insert forward and reverse the edge in the graph.


	edge1->m_flow = edge1->m_flow - flow;		// decrease the flow in the original edge.
	edge1->getReverseEdge()->m_flow = edge1->getReverseEdge()->m_flow - flow; // decrease the flow in the original edge.
//	edge1->coverageDepth = edge1->coverageDepth - coverage;
//	edge1->getReverseEdge()->coverageDepth = edge1->getReverseEdge()->coverageDepth - coverage;

	edge2->m_flow = edge2->m_flow - flow; // decrease the flow in the original edge.
	edge2->getReverseEdge()->m_flow = edge2->getReverseEdge()->m_flow - flow; // decrease the flow in the original edge.
//	edge2->coverageDepth = edge2->coverageDepth - coverage;
//	edge2->getReverseEdge()->coverageDepth = edge2->getReverseEdge()->coverageDepth - coverage;

	if(edge2 != edge1->getReverseEdge() && (edge2->m_flow == 0 || flow == 0))
			removeEdge(edge2);
	if(edge1->m_flow == 0 || flow == 0)				// If no flow left in edge1
		removeEdge(edge1);							// edge1 is deleted from the graph.

	return true;
}



/**********************************************************************************************************************
	Merge the list of reads, list of overlap offsets and list of orientations of two edges.
	CP: input: edge1 and edge2, NOT modified. All these pointers need to be changed to const Edge *edge1
	CP: the overlapOffset is the overlapOffset of the new edge, which is from the beginning of source read of the old edge1
	to the beginning of the destination read of the old edge2.
	CP: return-by-pointer: *listReads,  *listOverlapOffsets, *listOrientations
**********************************************************************************************************************/

bool OverlapGraph::mergeListDisconnected(Edge *edge1, Edge *edge2, UINT64 overlapOffset,
		UINT64 **returnListReads, UINT32 &lSize)
{
	UINT64 sum = 0;
	lSize=edge1->getListofReadsSize()+edge2->getListofReadsSize()+2;
	size_t lCtr=0;
	UINT64 *listReads = new UINT64[lSize];
	// CP: Add all the reads in the first edge, NOT including its source read
	for(UINT64 i = 0; i < edge1->getListofReadsSize(); i++)	// Get the list from the first edge.
	{
		listReads[lCtr++]=edge1->getInnerReadInfo(i);
		sum += edge1->getInnerOverlapOffset(i);
	}

	// CP: Add the destination read of the first edge and set its offset and orientation
	// CP: Add the destination read of the first edge.
	// CP: calculate the overlapOffset of the destination read

	UINT64 overlapOff = ((edge1->getOverlapOffset()-sum) << 32);
	UINT64 orient=0;
	if(edge1->getOrientation() == 1 || edge1->getOrientation() == 3)		// CP: calculate the orientation of the destination read
		orient=1;
	orient = orient << 63;
	UINT64 destRID =  edge1->getDestinationRead()->getReadID() | overlapOff | orient;
	listReads[lCtr++] = destRID;

	// CP: Add the source read of the second edge.
	// Add the source read of the second edge.
	// the overlapOff input from the function parameter
	overlapOff = (overlapOffset << 32);
	orient=0;
	if(edge2->getOrientation() == 2 || edge2->getOrientation() == 3)		// CP: calculate the orientation of the destination read
		orient = 1;
	orient = orient << 63;
	UINT64 srcRID = edge2->getSourceRead()->getReadID() | overlapOff | orient;
	listReads[lCtr++] = srcRID;

	// CP: Add all the reads in the second edge, NOT including its destination read
	for(UINT64 i = 0; i < edge2->getListofReadsSize(); i++)	// Get the list from the second edge.
	{
		listReads[lCtr++]=edge2->getInnerReadInfo(i);
	}
	*returnListReads = listReads;
	return true;
}




/**********************************************************************************************************************
	Orientation of the merged edge
**********************************************************************************************************************/
UINT8 OverlapGraph::mergedEdgeOrientationDisconnected(const Edge *edge1, const Edge *edge2)
{
	UINT8 or1 = edge1->getOrientation(), or2 = edge2->getOrientation(),returnValue;
	     if( (or1 == 1 || or1 == 0) && (or2 == 0 || or2 == 2) )		// <---------*  and *-----------<
		returnValue = 0;
	else if( (or1 == 1 || or1 == 0) && (or2 == 1 || or2 == 3) )		// <---------*  and *----------->
		returnValue = 1;
	else if( (or1 == 2 || or1 == 3) && (or2 == 0 || or2 == 2) )		// >---------*  and *-----------<
		returnValue = 2;
	else if( (or1 == 2 || or1 == 3) && (or2 == 1 || or2 == 3) )		// >---------*  and *----------->
		returnValue = 3;
	else
	{
		FILE_LOG(logINFO)<<(int)or1<<" "<<(int)or2<<endl;
		MYEXIT("Unable to merge.")
	}
	return returnValue;
}



/**********************************************************************************************************************
	Orientation of a reverse edge;
	Twin edge of Orientation 0 = <-------< is Orientation 3 = >------->
	Twin edge of Orientation 1 = <-------> is Orientation 1 = <------->
	Twin edge of Orientation 2 = >-------< is Orientation 2 = >-------<
	Twin edge of Orientation 3 = >-------> is Orientation 0 = <-------<
**********************************************************************************************************************/
UINT8 OverlapGraph::twinEdgeOrientation(UINT8 orientation)
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

void OverlapGraph::updateReadsLocations(Edge *edge, EdgeOP operation, DataSet *d)
{
	if(edge->isListofReads()){
		for(UINT64 i = 0; i < edge->getListofReadsSize(); ++i){
			Read *updateRead = m_dataset->at(edge->getInnerReadID(i));
			updateEdgeInfo(updateRead, edge, i, operation);
		}
	}
}



/*
 * ===  FUNCTION  ======================================================================
 *         Name:  updateEdgeInfo
 *  Description:  Update a read's residing edge information according to the
 *  		  insertion of new edge of deletion of existing edge
 *  		  TODO: in the case of loop reduction, it's possible that one read
 *  		  appears more than once in an edge. Should this function be called more
 *  		  than once in such cases?
 * =====================================================================================
 */
void OverlapGraph::updateEdgeInfo(Read * updateRead, Edge *edge, UINT32 read_index, EdgeOP operation)
{
	//CLOCKSTART;
	// Insert an edge with this read included on it
	if (operation == INSERTION){
		if((edge->getInnerOrientation(read_index) & 1) == 1)
			updateRead->setEdge(edge,read_index,0);
		else
			updateRead->setEdge(edge,read_index,1);
	}
	// Delete an edge with this read included on it
	else{
		if((edge->getInnerOrientation(read_index) & 1) == 1)
			updateRead->delEdge(edge, read_index,0);
		else
			updateRead->delEdge(edge, read_index,1);
	}
	//CLOCKSTOP;
}


void OverlapGraph::generateGFAOutput(ostream & gfaFilePointer)
{
	gfaFilePointer << "H\tVN:Z:2.0\n";
	UINT64 path_id=0;
	for(UINT64 i = 1; i<= m_dataset->size(); i++) //For each read
	{
		//Write segments
		//gfaFilePointer << "S\t" << i <<"\t"<<m_dataset->at(i)->getStringForward().length()<< "\t"<<m_dataset->at(i)->getStringForward() << std::endl;
		//No read sequence in GFA file
		gfaFilePointer << "S\t" << i <<"\t"<<m_dataset->at(i)->getReadLength()<< "\t*\n";
		auto it = m_graph->find(i);
		if(it != m_graph->end() && !it->second->empty()) // if this read has some edge(s) going out of it (since now the graph is directed)
		{
			t_edge_vec *eList = it->second;
			for(UINT64 j = 0; j < eList->size(); j++)	//For each edge
			{
				Edge *e = eList->at(j);
				//Write links
				UINT64 source = e->getSourceRead()->getReadID();
				UINT64 destination = e->getDestinationRead()->getReadID();
				string fOrientation=(e->getOrientation()==2||e->getOrientation()==3)?"+":"-";		//Orientation of source read in the first link of the edge
				string lOrientation=(e->getOrientation()==1||e->getOrientation()==3)?"+":"-";		//Orientation of destination read in the last link of the edge
				if(source < destination || (source == destination && e < e->getReverseEdge()))
				{
					if(e->getListofReadsSize()>0)
					{
						//Add first link with source
						gfaFilePointer << "L\t" << source << "\t"<<fOrientation<<"\t";
						string orientation=(e->getInnerOrientation(0)==0)?"-":"+";
						gfaFilePointer << e->getInnerReadID(0) << "\t"<<orientation<<"\t"<<
							e->getSourceRead()->getReadLength() - e->getInnerOverlapOffset(0)<<"M\n";
						string pathStr=SSTR(source)+fOrientation+",";
						string pathStrOverlap=SSTR(e->getSourceRead()->getReadLength()-e->getInnerOverlapOffset(0))+"M"+",";
						//Add middle links
						for(size_t j=1;j<e->getListofReadsSize();j++)
						{
							string orientation=(e->getInnerOrientation(j-1)==0)?"-":"+";
							gfaFilePointer << "L\t" << e->getInnerReadID(j-1) << "\t"<<orientation<<"\t";
							pathStr=pathStr+SSTR(e->getInnerReadID(j-1))+orientation+",";

							orientation=(e->getInnerOrientation(j)==0)?"-":"+";
							gfaFilePointer << e->getInnerReadID(j) << "\t"<<orientation<<"\t"<<
									m_dataset->at(e->getInnerReadID(j-1))->getReadLength() - e->getInnerOverlapOffset(j)<<"M\n";
							pathStrOverlap=pathStrOverlap+SSTR(m_dataset->at(e->getInnerReadID(j-1))->getReadLength() - e->getInnerOverlapOffset(j))+"M"+",";
						}
						//Add last link
						size_t lastInnerReadID=e->getInnerReadID(e->getListofReadsSize()-1);		//Get ID of last inner read
						orientation=(e->getInnerOrientation(e->getListofReadsSize()-1)==0)?"-":"+";			//Get the orientation
						gfaFilePointer << "L\t" << lastInnerReadID << "\t"<<orientation<<"\t";
						gfaFilePointer << destination << "\t"<<lOrientation<<"\t"<<
								m_dataset->at(lastInnerReadID)->getReadLength() - (e->getOverlapOffset()-e->getInnerOverlapSum(0,e->getListofReadsSize()))<<"M\n";
						pathStr=pathStr+SSTR(destination)+lOrientation;
						pathStrOverlap=pathStrOverlap.substr(0,pathStrOverlap.length()-1);
						path_id++;
						gfaFilePointer << "P\t" << path_id << "\t"<< pathStr <<"\t" << pathStrOverlap <<'\n';
					}
					else
					{
						//Add first link with source-destination (simple edge)
						gfaFilePointer << "L\t" << source << "\t"<<fOrientation<<"\t";
						gfaFilePointer << destination << "\t"<<lOrientation<<"\t"<<e->getOverlapOffset() << '\n';
					}

				}
			}
		}
	}
}

void OverlapGraph::generateGFA2Edge(ostream & gfaFilePointer, UINT64 edge_id, UINT64 source, string sOri,
		UINT64 destination,string dOri, UINT64 offset)
{
	if(sOri=="+" && dOri=="+")
	{
		gfaFilePointer << "E\t" << edge_id <<"\t"<< source << "\t"<<sOri<<"\t"<<destination<<"\t";

		UINT64 ovlLength = m_dataset->at(source)->getReadLength()-offset;
		gfaFilePointer << offset <<"\t"<< m_dataset->at(source)->getReadLength() <<"$\t0\t"<<
				ovlLength<<"\t"<<ovlLength<<"M\n";
	}
	if(sOri=="+" && dOri=="-")
	{
		gfaFilePointer << "E\t" << edge_id <<"\t"<< source << "\t"<<dOri<<"\t"<<destination<<"\t";

		UINT64 ovlLength = m_dataset->at(source)->getReadLength()-offset;
		gfaFilePointer << offset <<"\t"<< m_dataset->at(source)->getReadLength() <<"$\t"
				<<m_dataset->at(source)->getReadLength()-ovlLength<<"\t"
				<< m_dataset->at(source)->getReadLength()<<"$\t"<<ovlLength<<"M\n";
	}
	else if(sOri=="-" && dOri=="+")
	{
		gfaFilePointer << "E\t" << edge_id <<"\t"<< destination << "\t"<<sOri<<"\t"<<source<<"\t";

		UINT64 ovlLength = m_dataset->at(source)->getReadLength()-offset;
		gfaFilePointer << "0\t"<< ovlLength <<"$\t"
						<<m_dataset->at(source)->getReadLength()-ovlLength<<"\t"
						<< m_dataset->at(source)->getReadLength()<<"$\t"<<ovlLength<<"M\n";
	}
	else if(sOri=="-" && dOri=="-")
	{
		gfaFilePointer << "E\t" << edge_id <<"\t"<< source << "\t"<<dOri<<"\t"<<destination<<"\t";

		UINT64 ovlLength = m_dataset->at(source)->getReadLength()-offset;
		gfaFilePointer << "0\t"<< ovlLength <<"$\t"
						<<m_dataset->at(destination)->getReadLength()<<"$\t"
						<< m_dataset->at(destination)->getReadLength()-ovlLength<<"$\t"<<ovlLength<<"M\n";
	}
}

void OverlapGraph::generateGFA2Output(ostream & gfaFilePointer)
{
	gfaFilePointer << "H\tVN:Z:2.0\n";
	UINT64 path_id=0, edge_id=0;

	for(UINT64 i = 1; i<= m_dataset->size(); i++) //For each read
	{
		//Write segments
		gfaFilePointer << "S\t" << i <<"\t"<<m_dataset->at(i)->getReadLength()<< "\t*\n";
		auto it = m_graph->find(i);
		if(it != m_graph->end() && !it->second->empty()) // if this read has some edge(s) going out of it (since now the graph is directed)
		{
			t_edge_vec *eList = it->second;
			for(UINT64 j = 0; j < eList->size(); j++)	//For each edge
			{
				Edge *e = eList->at(j);
				//Write links
				UINT64 source = e->getSourceRead()->getReadID();
				UINT64 destination = e->getDestinationRead()->getReadID();
				string fOrientation=(e->getOrientation()==2||e->getOrientation()==3)?"+":"-";		//Orientation of source read in the first link of the edge
				string lOrientation=(e->getOrientation()==1||e->getOrientation()==3)?"+":"-";		//Orientation of destination read in the last link of the edge
				if(source < destination || (source == destination && e < e->getReverseEdge()))
				{
					if(e->getListofReadsSize()>0)
					{
						//Add first link with source
						edge_id++;
						string orientation=(e->getInnerOrientation(0)==0)?"-":"+";
						generateGFA2Edge(gfaFilePointer,edge_id,source,fOrientation,e->getInnerReadID(0),orientation,e->getInnerOverlapOffset(0));
						string pathStr=SSTR(edge_id)+"\t";
						//Add middle links
						for(size_t j=1;j<e->getListofReadsSize();j++)
						{
							string sorientation=(e->getInnerOrientation(j-1)==0)?"-":"+";
							string dorientation=(e->getInnerOrientation(j)==0)?"-":"+";
							edge_id++;
							generateGFA2Edge(gfaFilePointer,edge_id,e->getInnerReadID(j-1),sorientation,e->getInnerReadID(j),dorientation,e->getInnerOverlapOffset(j));
							pathStr=pathStr+SSTR(edge_id)+"\t";
						}
						//Add last link
						size_t lastInnerReadID=e->getInnerReadID(e->getListofReadsSize()-1);		//Get ID of last inner read
						orientation=(e->getInnerOrientation(e->getListofReadsSize()-1)==0)?"-":"+";			//Get the orientation
						edge_id++;
						generateGFA2Edge(gfaFilePointer,edge_id,lastInnerReadID,orientation,destination,lOrientation,
								(e->getOverlapOffset()-e->getInnerOverlapSum(0,e->getListofReadsSize())));
						pathStr=pathStr+SSTR(edge_id);
						path_id++;
						gfaFilePointer << "PO\t" << path_id << "\t"<< pathStr <<'\n';
					}
					else
					{
						//Add first link with source-destination (simple edge)
						edge_id++;
						generateGFA2Edge(gfaFilePointer,edge_id,source,fOrientation,destination,lOrientation,e->getOverlapOffset());
					}
				}
			}
		}
	}
}


/*
 * OverlapGraph.cpp
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 * Author: Abhishek Biswas
 */

#include "OverlapGraph.h"

#include "Common.h"

/**********************************************************************************************************************
	Function to compare two edges. Used for sorting.
**********************************************************************************************************************/
bool compareEdges (Edge *edge1, Edge* edge2)
{
	return (edge1->getOverlapOffset() < edge2->getOverlapOffset());
}

/**********************************************************************************************************************
	Constructor. Build the overlap graph using the hash table.
**********************************************************************************************************************/
OverlapGraph::OverlapGraph(HashTable *ht, UINT64 maxThreads,UINT64 maxParGraph,
		UINT64 maxMemSizeGB, string fnamePrefixGraph,string fnamePrefixSimplify,
		string simplifyPartialPath, int myid, int numprocs)
{
	// Initialize the variables.
	numberOfNodes = 0;
	numberOfEdges = 0;
	parallelThreadPoolSize=maxThreads;
	writeParGraphSize=maxParGraph;
	myProcID=myid;
	hashTable = ht;
	dataSet = ht->getDataset();
	UINT64 numElements=dataSet->getNumberOfUniqueReads()+1;
	myMarked=new int[numElements];								//Contained reads are considered already marked to remove them form further consideration...
	std::memset(myMarked, 0, numElements*sizeof(int));				//0 not marked; >0 already marked

	UINT64 mem_used = checkMemoryUsage();
	UINT64 maxMemSizeMB = maxMemSizeGB*1024;
	INT64 memPerThdMB = (maxMemSizeMB-mem_used)/parallelThreadPoolSize;

	if(memPerThdMB<=0)
		MYEXIT("Error!!! User did not provide enough memory for graph construction.");

	if(memPerThdMB>=MEGA_THD_MEMORY_SIZE)						//[20-INF)GB per thread available
		writeParGraphSize=MEGA_PAR_GRAPH_SIZE;
	else if(memPerThdMB>=MAX_THD_MEMORY_SIZE && memPerThdMB<MEGA_THD_MEMORY_SIZE)		//[10-20)GB per thread available
		writeParGraphSize=MAX_PAR_GRAPH_SIZE;
	else if(memPerThdMB>=MID_THD_MEMORY_SIZE && memPerThdMB<MAX_THD_MEMORY_SIZE)     	//[5,10)GB per thread available
		writeParGraphSize=MID_PAR_GRAPH_SIZE;
	else
		writeParGraphSize=MIN_PAR_GRAPH_SIZE;		// (0,5)GB per thread available


	buildOverlapGraphFromHashTable(ht,fnamePrefixGraph,fnamePrefixSimplify,simplifyPartialPath,numprocs);
}

/**********************************************************************************************************************
	Default destructor.
**********************************************************************************************************************/
OverlapGraph::~OverlapGraph()
{
	delete[] myMarked;
	// Free the memory used by the overlap graph.
}


/**********************************************************************************************************************
	Build the overlap graph from hash table
**********************************************************************************************************************/
bool OverlapGraph::buildOverlapGraphFromHashTable(HashTable *ht, string fnamePrefixGraph,string fnamePrefixSimplify,
		string simplifyPartialPath, int numprocs)
{
	CLOCKSTART;
	//Contained reads are considered already marked to remove them form further consideration...
	markContainedReads(fnamePrefixGraph, dataSet->getFRMap(), numprocs);
	UINT64 numNodes=dataSet->getNumberOfUniqueReads()+1;
	//Check if partial previous run data exists... Load partial graph data and mark reads.
	bool prevResultExists=false;
	#pragma omp parallel num_threads(parallelThreadPoolSize)
	{
		int threadID = omp_get_thread_num();
		string parFileName = fnamePrefixGraph + "_" + SSTR(myProcID) + "_" + SSTR(threadID) + "_parGraph.txt";
		if(ifstream(parFileName.c_str()))
		{
			prevResultExists=true;
			cout << "Thread:" << threadID << " Partial graph file exists. Loading marked reads." <<endl;
			ifstream filePointer;
			filePointer.open(parFileName.c_str());
			string text;
			UINT64 procCtr=0;
			while(getline(filePointer,text))
			{
				procCtr++;
				vector<string> toks = splitTok(text,'\t');
				//Get source destination IDs
				UINT64 sourceReadFindex = std::stoull(toks[0],nullptr,0);
				UINT64 destReadFindex = std::stoull(toks[1],nullptr,0);
				auto sourceIt = dataSet->getFRMap()->find(sourceReadFindex);
				auto destIt = dataSet->getFRMap()->find(destReadFindex);
				if(sourceIt == dataSet->getFRMap()->end() || destIt == dataSet->getFRMap()->end() ) // If the reads do not exist go to the next record
					continue;

				//Check if both marked or not
				vector<string> toks2 = splitTok(toks[2],',');
				UINT64 markFlag = std::stoull(toks2[toks2.size()-1],nullptr,0);

				//Check if destination is also marked by this thread.
				// 0: Only source is marked
				// 1: Only destination is marked
				// 2: Both source and destination are marked
				if(markFlag==0)
					myMarked[sourceIt->second]=1;
				else if (markFlag==1)
					myMarked[destIt->second]=1;
				else
				{
					myMarked[sourceIt->second]=1;
					myMarked[destIt->second]=1;
				}
				if(procCtr%1000000==0)
					cout<<"Proc"<<myProcID<< " Thread:" << threadID << " " <<procCtr<<" marked reads loaded ..."<<endl;
			}
		}
	}
	if(numprocs>1 && prevResultExists)			//Previous result exists and not only one process as it will deadlock
	{
		MPI_Request request[numprocs];
		for(int i=0;i<numprocs;i++)
		{
			if(i==myProcID)
			{
				for(int j=0;j<numprocs;j++)
					if(i!=j)
						MPI_Isend(myMarked, numNodes, MPI_INT, j, 0, MPI_COMM_WORLD, &request[j]);
			}
			else
			{
				int* readIDBuf = new int[numNodes];
				// Now receive the message with the allocated buffer
				MPI_Recv(readIDBuf, numNodes, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				//Set contained read values
				for(UINT64 i = 1; i < numNodes; i++) // For each readid
				{
					if(readIDBuf[i]>0 && myMarked[i]==0)
						myMarked[i]=1;
				}
				delete[] readIDBuf;
			}
		}
		for(int i=0;i<numprocs;i++)
			if(i!=myProcID)
				MPI_Wait(&request[i],MPI_STATUS_IGNORE);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	//Restart operations complete. Delete file index to read ID map
	dataSet->freeFindexReadIDMAP();
	//Starting graph construction

	hashTable->setLockAll();
	#pragma omp parallel num_threads(parallelThreadPoolSize)
	{
		int threadID = omp_get_thread_num();
		if(threadID==0 && numprocs>1)
		{
			bool allCompleteFlag=false;
			bool allRemoteFinish=false;
			// Allocate a buffer to hold the incoming numbers
			int *readIDBuf = new int[numNodes];
			size_t usedCounter=1;
			while(!allCompleteFlag)
			{
				//Put main thread to sleep before next round of MPI...
				std::this_thread::sleep_for(std::chrono::seconds(10));
				MPI_Request sendRequest[numprocs-1];
				size_t reqCtr=0;
				//Send all my data
				for(int j=0;j<numprocs;j++)
				{
					if(myProcID!=j)
					{
						#pragma omp critical(getRemoteData)
						{
							MPI_Isend(myMarked, numNodes, MPI_INT, j, 0, MPI_COMM_WORLD, &sendRequest[reqCtr++]);
						}
					}
				}
				//Receive all my data
				for(int i=0;i<numprocs;i++)
				{
					if(myProcID!=i)
					{
						std::memset(readIDBuf, 0, numNodes*sizeof(int));
						#pragma omp critical(getRemoteData)
						{
							MPI_Recv(readIDBuf, numNodes, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						}
						for(UINT64 i = usedCounter; i < numNodes; i++) // For each readid
						{
							if(readIDBuf[i]>0 && myMarked[i]==0)
							{
								#pragma omp atomic write
									myMarked[i]=1;
							}
						}
					}
				}
				//Wait till all the sending is done...
				#pragma omp critical(getRemoteData)
				{
					MPI_Waitall(numprocs-1, sendRequest,MPI_STATUS_IGNORE);
				}
				//Test completion
				int myFinFlag=1;
				for(;usedCounter<numNodes;usedCounter++)
				{
					int usedFlag=0;
					#pragma omp atomic read
						usedFlag=myMarked[usedCounter];
					if(usedFlag==0)
					{
						myFinFlag=0;
						break;
					}
				}
				//cout<<"Proc:"<<myProcID<<" Main communication fin flag: "<<myFinFlag << " stuck at: "<<loopCtr<<endl;
				//Inform others of status
				reqCtr=0;
				allRemoteFinish=true;
				for(int j=0;j<numprocs;j++)
				{
					if(myProcID!=j)
					{
						#pragma omp critical(getRemoteData)
						{
							MPI_Isend(&myFinFlag, 1, MPI_INT, j, 0, MPI_COMM_WORLD, &sendRequest[reqCtr++]);
						}
					}

				}
				//Get status from others
				for(int i=0;i<numprocs;i++)
				{
					if(myProcID!=i)
					{
						int remoteFin=0;
						#pragma omp critical(getRemoteData)
						{
							MPI_Recv(&remoteFin, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						}
						allRemoteFinish = (allRemoteFinish && remoteFin);
					}
				}
				//Wait till all the sending is done...
				#pragma omp critical(getRemoteData)
				{
					MPI_Waitall(numprocs-1, sendRequest,MPI_STATUS_IGNORE);
				}
				//If this process is finished and all remote processes are finished end while loop
				if(allRemoteFinish && myFinFlag)
					allCompleteFlag=true;
			}//end of while
			delete[] readIDBuf;
			cout<<"Proc:"<<myProcID<<" Main communication thread complete!!!"<<endl;
		}
		else
		{
			UINT64 startReadID=0,prevReadID=0;
			UINT64 workBlock = numNodes/(numprocs*(parallelThreadPoolSize-1));
			prevReadID=startReadID=(myProcID*(parallelThreadPoolSize-1)*workBlock)+threadID;
			while(startReadID!=0 && startReadID<numNodes) // Loop till all nodes marked
			{
				map<UINT64,nodeType> *exploredReads = new map<UINT64,nodeType>;							//Record of nodes processed
				queue<UINT64> *nodeQ = new queue<UINT64>;												//Queue
				map<UINT64, vector<Edge*> * > *parGraph = new map<UINT64, vector<Edge*> * >;			//Partial graph

				vector<Edge *> *newList = new vector<Edge *>;
				parGraph->insert( std::pair<UINT64, vector<Edge*> * >(startReadID, newList)); // Insert start node
				UINT64 writtenMakedNodes=0;
				nodeQ->push(startReadID);  											// // Initialize queue start and end.
				while(!nodeQ->empty() && writtenMakedNodes<writeParGraphSize) 													// This loop will explore all connected component starting from read startReadID.
				{
					UINT64 read1 = nodeQ->front();										//Pop from queue...
					nodeQ->pop();
					int usedFlag=0;
					#pragma omp atomic read
						usedFlag=myMarked[read1];
					if(usedFlag==0)
					{
						#pragma omp atomic write
							myMarked[read1]=1;
						if(exploredReads->find(read1) ==  exploredReads->end()) //if node is UNEXPLORED
						{
							insertAllEdgesOfRead(read1, exploredReads, parGraph);					// Explore current node.
							exploredReads->insert( std::pair<UINT64,nodeType>(read1,EXPLORED) );
						}
						if(parGraph->at(read1)->size() != 0) 								// Read has some edges (required only for the first read when a new queue starts.
						{
							if(exploredReads->at(read1) == EXPLORED) 					// Explore unexplored neighbors first.
							{
								for(UINT64 index1 = 0; index1 < parGraph->at(read1)->size(); index1++ )
								{
									UINT64 read2 = parGraph->at(read1)->at(index1)->getDestinationRead()->getReadNumber();
									if(exploredReads->find(read2) ==  exploredReads->end()) 			// Not explored.
									{
										nodeQ->push(read2);  						// Put in the queue.
										insertAllEdgesOfRead(read2, exploredReads, parGraph);
										exploredReads->insert( std::pair<UINT64,nodeType>(read2,EXPLORED) );
									}
								}
								markTransitiveEdges(read1, parGraph); // Mark transitive edges
								exploredReads->at(read1) = EXPLORED_AND_TRANSITIVE_EDGES_MARKED;
							}
							if(exploredReads->at(read1) == EXPLORED_AND_TRANSITIVE_EDGES_MARKED)
							{
								for(UINT64 index1 = 0; index1 < parGraph->at(read1)->size(); index1++) 				// Then explore all neighbour's neighbors
								{
									UINT64 read2 = parGraph->at(read1)->at(index1)->getDestinationRead()->getReadNumber();
									if(exploredReads->at(read2) == EXPLORED)
									{
										for(UINT64 index2 = 0; index2 < parGraph->at(read2)->size(); index2++) 		// Explore all neighbors neighbors
										{
											UINT64 read3 = parGraph->at(read2)->at(index2)->getDestinationRead()->getReadNumber();
											if(exploredReads->find(read3) ==  exploredReads->end()) 				// Not explored
											{
												nodeQ->push(read3);  					// Put in the queue
												insertAllEdgesOfRead(read3, exploredReads, parGraph);
												exploredReads->insert( std::pair<UINT64,nodeType>(read3,EXPLORED) );
											}
										}
										markTransitiveEdges(read2, parGraph); // Mark transitive edge
										exploredReads->at(read2) = EXPLORED_AND_TRANSITIVE_EDGES_MARKED;
									}
								}
								removeTransitiveEdges(read1, parGraph); // Remove the transitive edges
								exploredReads->at(read1) = EXPLORED_AND_TRANSITIVE_EDGES_REMOVED;
								writtenMakedNodes++;
							}
						}
					}
				}
				INT64 mem_used = checkMemoryUsage();
				if(writtenMakedNodes>MIN_LOG_SIZE)
					cout<<"Proc:"<<myProcID<<" Thread:"<<threadID<<" Start Read ID:"<<startReadID<<" Reads Processed:"<<writtenMakedNodes<<" Memory Used:" << mem_used << '\n';
				saveParGraphToFile(fnamePrefixGraph + "_" + SSTR(myProcID) + "_" + SSTR(threadID) + "_parGraph.txt" , exploredReads, parGraph);
				for (map<UINT64, vector<Edge*> * >::iterator it=parGraph->begin(); it!=parGraph->end();it++)
				{
					UINT64 readID = it->first;
					for(UINT64 j = 0; j< parGraph->at(readID)->size(); j++)
					{
						delete parGraph->at(readID)->at(j);
					}
					delete parGraph->at(readID);
				}
				delete parGraph;
				delete exploredReads;
				delete nodeQ;
				startReadID=0;
				for(UINT64 i=prevReadID;i<numNodes;i++)
				{
					int usedFlag=0;
					#pragma omp atomic read
						usedFlag=myMarked[i];
					if(usedFlag==0){
						startReadID=prevReadID=i;
						break;
					}
				}
			}
			string edge_file=fnamePrefixGraph + "_" + SSTR(myProcID) + "_" + SSTR(threadID) + "_parGraph.txt" ;
			string prev_composite_out_edge_file = fnamePrefixSimplify + "_" + SSTR((myProcID*(parallelThreadPoolSize-1))+(threadID-1)) +"_ParSimpleEdges.txt";
			string runSimplifyExeStr = simplifyPartialPath + "/parsimplify " + edge_file + " " + prev_composite_out_edge_file
									+ " " + SSTR(dataSet->getMinimumOverlapLength()) + " " + SSTR(parallelThreadPoolSize) ;
			//Perform the first partial graph simplification
			int retStatus = system(runSimplifyExeStr.c_str());
		}
	}
	hashTable->unLockAll();
	hashTable->endEpoch();
	cout<<endl<<"Graph Construction Complete"<<endl;
	cout<<"Process:"<<myProcID<<" RMA OPS:"<<hashTable->getRMACount()<<endl;
	CLOCKSTOP;
	return true;
}


/**********************************************************************************************************************
	This function check if a read contains other small reads. If a read is contained in more than one super read
	then it is assigned to the longest such super read. Also duplicate reads are marked
**********************************************************************************************************************/
void OverlapGraph::markContainedReads(string fnamePrefix, map<UINT64, UINT64> *fIndxReadIDMap, int numprocs)
{
	CLOCKSTART;
	ofstream filePointer;
	UINT64 nonContainedReads = 0;
	string testContainedReadFile = fnamePrefix+"_"+SSTR(myProcID)+"_"+SSTR(0)+"_containedReads.txt";
	if(ifstream(testContainedReadFile.c_str()))
	{
		#pragma omp parallel num_threads(parallelThreadPoolSize)
		{
			int threadID = omp_get_thread_num();
			cout << "Contained read file already exists. Using this file." <<endl;
			string containedReadFile = fnamePrefix+ "_" + SSTR(myProcID) + "_" + SSTR(threadID) +"_containedReads.txt";
			ifstream filePointer;
			filePointer.open(containedReadFile.c_str());
			if(!filePointer)
				MYEXIT("Unable to open contained reads file: +"+containedReadFile);
			string text;
			UINT64 procCtr=0;
			while(getline(filePointer,text))
			{
				procCtr++;
				vector<string> toks = splitTok(text,'\t');
				UINT64 containedReadFindex = std::stoull(toks[0],nullptr,0);
				UINT64 containingReadFindex = std::stoull(toks[1],nullptr,0);
				auto it = fIndxReadIDMap->find(containedReadFindex);
				Read *r = dataSet->getReadFromID(it->second); // Get the read
				r->superReadID=containingReadFindex;
				if(procCtr%1000000==0)
					cout<<procCtr<<" contained reads processed..."<<endl;
			}
			filePointer.close();
		}
	}
	else
	{
		vector< shared_ptr<ofstream> > filePointerList;
		for(UINT64 i = 0; i < parallelThreadPoolSize; i++) // For each thread
		{
			string containedReadFile = fnamePrefix+ "_" + SSTR(myProcID) + "_" + SSTR(i) +"_containedReads.txt";
			shared_ptr<ofstream> filePointer = make_shared<ofstream>(containedReadFile);
			//filePointer.open(containedReadFile.c_str());
			if(!*(filePointer))
				MYEXIT("Unable to open contained read file: +"+containedReadFile);
			filePointerList.push_back(filePointer);
		}
		UINT64 nextReadOffset=0;
		hashTable->setLockAll();
		#pragma omp parallel num_threads(parallelThreadPoolSize)
		{
			int threadID = omp_get_thread_num();
			UINT64 thisReadOffset=0;
			//assign initial start points...
			#pragma omp critical(assignNextOffset)
			{
				thisReadOffset=nextReadOffset;
				if(nextReadOffset<hashTable->getMemoryMaxLocalOffset(myProcID))			//Read offset exceeds local limit...
					nextReadOffset = hashTable->getLocalNextOffset(nextReadOffset,myProcID); // Get the next read
			}
			while(thisReadOffset<hashTable->getMemoryMaxLocalOffset(myProcID))
			{
				UINT64 read1ID = hashTable->getLocalReadID(thisReadOffset,myProcID); // Get the read
				Read *read1 = dataSet->getReadFromID(read1ID);
				string read1String = hashTable->getLocalStringForward(thisReadOffset,myProcID); // Get the forward of the read
				UINT64 read1Len=read1String.length();
				//assign next start point...
				#pragma omp critical(assignNextOffset)
				{
					thisReadOffset=nextReadOffset;
					if(nextReadOffset<hashTable->getMemoryMaxLocalOffset(myProcID))			//Read offset exceeds local limit...
						nextReadOffset = hashTable->getLocalNextOffset(nextReadOffset,myProcID); // Get the next read
				}
				if(read1->superReadID!=0 || !hashTable->needsProcessing(read1ID,read1String,myProcID) || read1->containmentChecked)		//If read is already marked as contained, there is no need to look for contained reads within it
																											//If read will be processed by another node then don't proceed
					continue;

				vector<UINT64*> *localReadHits;
				/*Start global communication epoch*/
				localReadHits = hashTable->setLocalHitList_nocache(read1String,myProcID); // Search the substring in the hash table
				/*End global communication epoch*/
				string subString;
				read1->containmentChecked=true;
				ostringstream strstream;
				for(UINT64 j = 0; j < read1String.length() - hashTable->getHashStringLength(); j++) // fGr each substring of read1 of length getHashStringLength
				{
					subString = read1String.substr(j,hashTable->getHashStringLength()); // Get the substring from read1
					map<UINT64,string> listOfReads = hashTable->getLocalHitList_nocache(localReadHits, subString, j);
					if(!listOfReads.empty()) // If other reads contain the substring as prefix or suffix
					{
						for(map<UINT64,string>::iterator it=listOfReads.begin() ; it!=listOfReads.end(); ++it) // For each read in the list.
						{
							UINT64 data = it->first; // We used bit operation in the hash table to store read ID and orientation
							UINT64 read2ID = data & 0X3FFFFFFFFFFFFFFF;
							string read2String = it->second;
							Read *read2 = dataSet->getReadFromID(read2ID); 	// Least significant 62 bits store the read number.
																								// Most significant 2 bits store the orientation.
																								// Orientation 0 means prefix of forward of the read
																								// Orientation 1 means suffix of forward of the read
																								// Orientation 2 means prefix of reverse of the read
																								// Orientation 3 means prefix of reverse of the read
							if(read2->superReadID!=0)		//If read is already marked as contained, there is no need to check if its contained again
								continue;
							UINT64 read2Len = read2String.length();
							if(read1->getReadNumber() != read2->getReadNumber() && checkOverlapForContainedRead(read1String,read2String,(data >> 62),j)) // read1 need to be longer than read2 in order to contain read2
																																					 // Check if the remaining of the strings also match
							{
								if(read1String.length() > read2Len)
								{
									UINT64 overlapLen=0;
									UINT64 orientation=1;
									switch (data >> 62) // Most significant 2 bit represents  00 - prefix forward, 01 - suffix forward, 10 -  prefix reverse, 11 -  suffix reverse.
									{
										case 0: orientation = 3; overlapLen = read1Len - j; break; 				// 3 = r1>------->r2
										case 1: orientation = 0; overlapLen = hashTable->getHashStringLength() + j; break; 		// 0 = r1<-------<r2
										case 2: orientation = 2; overlapLen = read1Len - j; break; 				// 2 = r1>-------<r2
										case 3: orientation = 1; overlapLen = hashTable->getHashStringLength() + j; break; 		// 1 = r2<------->r2
									}
									if(read2->superReadID == 0) // This is the first super read found. we store the ID of the super read.
											read2->superReadID = read1ID;
									//Write contained read information regardless as it is a super read has been identified
									strstream<<read2->getFileIndex()<<"\t"<<read1->getFileIndex()<<"\t"<<orientation<<","
											<<read2Len<<","
											<<"0"<<","<<"0"<<","								//No substitutions or edits
											<<read2Len<<","					//Cointained Read (len,start,stop)
											<<"0"<<","
											<<read2Len<<","
											<<read1Len<<","					//Super Read (len,start,stop)
											<<read1Len-overlapLen<<","
											<<read1Len-overlapLen+read2Len
											<<'\n';
								}
								else if(read1String.length() == read2Len && read1->getReadNumber() < read2->getReadNumber())
								{
									UINT64 overlapLen=0;
									UINT64 orientation=1;
									switch (data >> 62) // Most significant 2 bit represents  00 - prefix forward, 01 - suffix forward, 10 -  prefix reverse, 11 -  suffix reverse.
									{
										case 0: orientation = 3; overlapLen = read1Len - j; break; 				// 3 = r1>------->r2
										case 1: orientation = 0; overlapLen = hashTable->getHashStringLength() + j; break; 		// 0 = r1<-------<r2
										case 2: orientation = 2; overlapLen = read1Len - j; break; 				// 2 = r1>-------<r2
										case 3: orientation = 1; overlapLen = hashTable->getHashStringLength() + j; break; 		// 1 = r2<------->r2
									}

									if(read2->superReadID==0)
										read2->superReadID = read1ID;
									if(read1->getReadNumber() < read2->superReadID)
										read2->superReadID = read1ID;

									//Write duplicate read information regardless as it is a super read has been identified
									strstream<<read2->getFileIndex()<<"\t"<<read1->getFileIndex()<<"\t"<<orientation<<","
											<<read2Len<<","
											<<"0"<<","<<"0"<<","								//No substitutions or edits
											<<read2Len<<","					//Duplicate Read (len,start,stop)
											<<"0"<<","
											<<read2Len<<","
											<<read1Len<<","					//Super Read (len,start,stop)
											<<read1Len-overlapLen<<","
											<<read1Len-overlapLen+read2Len
											<<'\n';
								}
							}
						}
					}
				}//End of inner for
				hashTable->deleteLocalHitList(localReadHits);
				*(filePointerList[threadID])<<strstream.str();
				filePointerList[threadID]->flush();
			}//end of while
		}
		for(UINT64 i = 0; i < parallelThreadPoolSize; i++) // For each thread close file pointer
		{
			filePointerList[i]->close();
		}
		cout<<"Proc:"<<myProcID<<" Completed contained read computation."<<endl;
		hashTable->unLockAll();
		hashTable->endEpoch();
	}

	if(numprocs>1)			//Exchange computed contained reads; Only one process will deadlock
	{

		int ctdReads=0;
		//Get contained read count
		#pragma omp parallel for schedule(guided) reduction(+:ctdReads) num_threads(parallelThreadPoolSize)
		for(UINT64 i = 1; i <= dataSet->getNumberOfUniqueReads(); i++) // For each read
		{
			Read *read1 = dataSet->getReadFromID(i); // Get the read
			if(read1->superReadID!=0)		//Check if read is marked as contained
				ctdReads++;
		}
		UINT64 *buf = new UINT64[ctdReads];
		size_t j=0;
		//Populate buffer of contained reads
		for(UINT64 i = 1; i <= dataSet->getNumberOfUniqueReads(); i++) // For each read
		{
			Read *read1 = dataSet->getReadFromID(i); // Get the read
			if(read1->superReadID!=0)		//Check if read is marked as contained
			{
				buf[j]=i;
				j++;
			}
		}

		MPI_Request request[numprocs];
		for(int i=0;i<numprocs;i++)
		{
			if(i==myProcID)
			{
				for(int j=0;j<numprocs;j++)
					if(i!=j)
						MPI_Isend(buf, ctdReads, MPI_UINT64_T, j, 0, MPI_COMM_WORLD, &request[j]);
			}
			else
			{
				int number_amount;
				MPI_Status status;
				// Probe for an incoming message from process i
				MPI_Probe(i, 0, MPI_COMM_WORLD, &status);
				// When probe returns, the status object has the size and other
				// attributes of the incoming message. Get the message size
				MPI_Get_count(&status, MPI_UINT64_T, &number_amount);
				// Allocate a buffer to hold the incoming numbers
				UINT64* readIDBuf = new UINT64[number_amount];
				// Now receive the message with the allocated buffer
				MPI_Recv(readIDBuf, number_amount, MPI_UINT64_T, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				//Set contained read values
				for(int i = 0; i < number_amount; i++) // For each readid
				{
					Read *read1 = dataSet->getReadFromID(readIDBuf[i]); // Get the read
					read1->superReadID=1;
				}
				delete[] readIDBuf;
			}
		}
		for(int i=0;i<numprocs;i++)
			if(i!=myProcID)
				MPI_Wait(&request[i],MPI_STATUS_IGNORE);
		delete[] buf;
	}
	//Delete buffer of contained reads
	#pragma omp parallel for schedule(guided) reduction(+:nonContainedReads) num_threads(parallelThreadPoolSize)
	for(UINT64 i = 1; i <= dataSet->getNumberOfUniqueReads(); i++) // For each read
	{
		Read *read1 = dataSet->getReadFromID(i); // Get the read
		if(read1->superReadID==0)		//If read is already marked as contained, there is no need to look for contained reads within it
			nonContainedReads++;
		else
			myMarked[i-1]=1;			//Once marked this read will not be marked by the graph construction...
	}
	cout<< '\n' << setw(10) << nonContainedReads << " Non-contained reads. (Keep as is)\n";
	cout<< setw(10) << dataSet->getNumberOfUniqueReads()-nonContainedReads << " contained reads. (Need to change their mate-pair information)\n";
	cout<<"Process:"<<myProcID<<" RMA OPS:"<<hashTable->getRMACount()<<endl;
	CLOCKSTOP;
}


/**********************************************************************************************************************
	Hash table search found that a proper substring of read1 is a prefix or suffix of read2 or reverse complement of
	read2 (dependents on the orient).
	orient 0 means prefix of forward of the read2
	orient 1 means suffix of forward of the read2
	orient 2 means prefix of reverse of the read2
	orient 3 means prefix of reverse of the read2
	We need to check if the remaining of the stings match to see if read2 is contained in read1.
**********************************************************************************************************************/
bool OverlapGraph::checkOverlapForContainedRead(string read1, string read2, UINT64 orient, UINT64 start)
{
	UINT64 hashStringLength = hashTable->getHashStringLength(), lengthRemaining1, lengthRemaining2;
	string string2 = (orient == 0 || orient== 1) ? read2 : reverseComplement(read2); // Get the string in read2 based on the orientation.
	if(orient == 0 || orient == 2)
									// orient 0
									//   >--------MMMMMMMMMMMMMMM*******------> read1      M means match found by hash table
									//            MMMMMMMMMMMMMMM*******>       read2      * means we need to check these characters for match
									//				OR
									// orient 2
									//	 >---*****MMMMMMMMMMMMMMM*******------> read1
									//		      MMMMMMMMMMMMMMM*******<	    Reverese complement of read2
	{
		lengthRemaining1 = read1.length() - start - hashStringLength; 	// This is the remaining of read1
		lengthRemaining2 = string2.length() - hashStringLength; 	// This is the remaining of read2
		if(lengthRemaining1 >= lengthRemaining2)
		{
			return read1.substr(start + hashStringLength, lengthRemaining2) == string2.substr(hashStringLength, lengthRemaining2); // If the remaining of the string match, then read2 is contained in read1
		}
	}
	else							// orient 1
									//   >---*****MMMMMMMMMMMMMMM-------------> read1      M means match found by hash table
									//      >*****MMMMMMMMMMMMMMM       		read2      * means we need to check these characters for match
									//				OR
									// orient 3
									//	 >---*****MMMMMMMMMMMMMMM-------------> read1
									//		<*****MMMMMMMMMMMMMMM				Reverse Complement of Read2
	{
		lengthRemaining1 = start;
		lengthRemaining2 = string2.length() - hashStringLength;
		if(lengthRemaining1 >= lengthRemaining2)
		{
			return read1.substr(start - lengthRemaining2, lengthRemaining2) == string2.substr(0, lengthRemaining2); // If the remaining of the string match, then read2 is contained in read1
		}
	}
	return false;

}




/**********************************************************************************************************************
	Checks if two read overlaps.
	Hash table search found that a proper substring of read1 is a prefix or suffix of read2 or reverse complement of
	read2 (depents on the orient).
	Orientation 0 means prefix of forward of the read2
	Orientation 1 means suffix of forward of the read2
	Orientation 2 means prefix of reverse of the read2
	Orientation 3 means prefix of reverse of the read2
	We need to check if the remaining of the stings match to see if read1 and read2 overlap.
**********************************************************************************************************************/
bool OverlapGraph::checkOverlap(string read1, string read2, UINT64 orient, UINT64 start)
{
	UINT64 hashStringLength = hashTable->getHashStringLength();
	string string2 = (orient == 0 || orient== 1) ? read2 : reverseComplement(read2); // Get the string in read2 based on the orientation.
	if(orient == 0 || orient == 2)		// orient 0
										//   >--------MMMMMMMMMMMMMMM*************> 			read1      M means match found by hash table
										//            MMMMMMMMMMMMMMM*************------->      read2      * means we need to check these characters for match
										//				OR
										// orient 2
										//	 >---*****MMMMMMMMMMMMMMM*************> 			read1
										//		      MMMMMMMMMMMMMMM*************-------<	    Reverese complement of read2
	{
		if(read1.length()- start - hashStringLength >= string2.length() - hashStringLength) // The overlap must continue till the end.
			return false;
		return read1.substr(start + hashStringLength, read1.length()-(start + hashStringLength)) == string2.substr(hashStringLength,  read1.length()-(start + hashStringLength)); // If the remaining strings match.
	}
	else								// orient 1
										//   	>********MMMMMMMMMMMMMMM-------------> 			read1      M means match found by hash table
										//  >----********MMMMMMMMMMMMMMM       		    		read2      * means we need to check these characters for match
										//				OR
										// orient 3
										//	 	>********MMMMMMMMMMMMMMM-------------> 			read1
										//	<----********MMMMMMMMMMMMMMM						Reverse Complement of Read2
	{
		if(string2.length()-hashStringLength < start)
			return false;
		return read1.substr(0, start) == string2.substr(string2.length()-hashStringLength-start, start); // If the remaining strings match.
	}
}

/**********************************************************************************************************************
	Insert an edge in the overlap graph.
**********************************************************************************************************************/
bool OverlapGraph::insertEdge(Edge * edge, map<UINT64, vector<Edge*> * > *parGraph)
{

	UINT64 ID = edge->getSourceRead()->getReadNumber(); // This is the source read.
	if(parGraph->find(ID) == parGraph->end()){ 			// If there is no edge incident to the node
		vector<Edge *> *newList = new vector<Edge *>;
		parGraph->insert( std::pair<UINT64, vector<Edge*> * >(ID, newList));
	}
	parGraph->at(ID)->push_back(edge);						// Insert the edge in the list of edges of ID
	return true;
}

/**********************************************************************************************************************
	Insert an edge in the graph.
**********************************************************************************************************************/
bool OverlapGraph::insertEdge(Read *read1, Read *read2, UINT64 r1Len, UINT64 r2Len, UINT8 orient, UINT16 overlapOffset, map<UINT64, vector<Edge*> * > *parGraph)
{
	Edge * edge1 = new Edge(read1,r1Len,read2,r2Len,orient,overlapOffset);								// Create a new edge in the graph to insert.
	UINT16 overlapOffsetReverse = r2Len + overlapOffset - r1Len;	// Set the overlap offset accordingly for the reverse edge. Note that read lengths are different.
																						// If read lengths are the same. Then the reverse edge has the same overlap offset.
	Edge * edge2 = new Edge(read2,r2Len,read1,r1Len,twinEdgeOrientation(orient),overlapOffsetReverse);		// Create a new edge for the reverses string.

	edge1->setReverseEdge(edge2);		// Set the reverse edge pointer.
	edge2->setReverseEdge(edge1);		// Set the reverse edge pinter.
	insertEdge(edge1, parGraph);					// Insert the edge in the overlap graph.
	insertEdge(edge2, parGraph);					// Insert the edge in the overlap graph.
	return true;
}

/**********************************************************************************************************************
	Insert all edges of a read in the overlap graph
**********************************************************************************************************************/
bool OverlapGraph::insertAllEdgesOfRead(UINT64 readNumber, map<UINT64,nodeType> * exploredReads, map<UINT64, vector<Edge*> * > *parGraph)
{
	Read *read1 = dataSet->getReadFromID(readNumber); 	// Get the current read read1.
	string readString="";
	/*Start global communication epoch*/
	readString = hashTable->getStringForward(read1,myProcID); 		// Get the forward string of read1.
	/*End global communication epoch*/
	string subString;
	vector<UINT64> insertedEdgeList;
	UINT64 read1Len = readString.length();

	/*Start global communication epoch*/
	//vector<UINT64*> *localReadHits;
	//localReadHits = hashTable->setLocalHitList_nocache(readString,myProcID); // Search the substring in the hash table
	vector<shared_ptr<UINT64> > localReadHits = hashTable->setLocalHitList(readString,myProcID); // Search the substring in the hash table
	/*End global communication epoch*/
	for(UINT64 j = 1; j < read1Len-hashTable->getHashStringLength(); j++) // For each proper substring of length getHashStringLength of read1
	{
		subString = readString.substr(j,hashTable->getHashStringLength());  // Get the proper substring s of read1.
		map<UINT64,string> listOfReads = hashTable->getLocalHitList(localReadHits, subString, j); // Search the string in the hash table.
		//map<UINT64,string> listOfReads = hashTable->getLocalHitList_nocache(localReadHits, subString, j);
		if(!listOfReads.empty()) // If other reads contain the substring as prefix or suffix
		{
			int insertCtr=0;
			for(map<UINT64,string>::iterator it=listOfReads.begin() ; it!=listOfReads.end() && insertCtr < MAX_EDGE_PER_KMER; ++it) // For each read in the list.
			{
				UINT64 data = it->first; // We used bit operations in the hash table. Most significant 2 bits store orientation and least significant 62 bits store read ID.
				UINT64 read2ID = data & 0X3FFFFFFFFFFFFFFF;
				string read2String = it->second;
				UINT16 overlapLen=0;
				UINT8 orientation=1;
				Read *read2 = dataSet->getReadFromID(read2ID); 	// Least significant 62 bits store the read number.
				if(exploredReads->find(read2ID) !=  exploredReads->end())
					continue;                                                       // No need to discover the same edge again. All edges of read2 is already inserted in the graph.

				if(readNumber != read2ID 											//Must not be a loop
						&& find(insertedEdgeList.begin(), insertedEdgeList.end(), read2ID)==insertedEdgeList.end()     //Must not have already added an edge with greater overlap
						&& read1->superReadID == 0 && read2->superReadID == 0		// Both read need to be non contained.
						&& checkOverlap(readString,read2String,(data >> 62),j)) 				// Must overlap
				{
					switch (data >> 62) // Most significant 2 bit represents  00 - prefix forward, 01 - suffix forward, 10 -  prefix reverse, 11 -  suffix reverse.
					{
						case 0: orientation = 3; overlapLen = read1Len - j; break; 				// 3 = r1>------->r2
						case 1: orientation = 0; overlapLen = hashTable->getHashStringLength() + j; break; 		// 0 = r1<-------<r2
						case 2: orientation = 2; overlapLen = read1Len - j; break; 				// 2 = r1>-------<r2
						case 3: orientation = 1; overlapLen = hashTable->getHashStringLength() + j; break; 		// 1 = r2<------->r2
					}
					insertEdge(read1,read2,read1Len, read2String.length(),orientation,read1Len-overlapLen, parGraph); 			// Insert the edge in the graph.
					insertedEdgeList.push_back(read2ID);
					insertCtr++;
				}
			}
		}
	}
	//hashTable->deleteLocalHitList(localReadHits);
	for(UINT64 j = 0; j < localReadHits.size(); j++)
		localReadHits[j].reset();

	if(parGraph->at(readNumber)->size() != 0)
		sort(parGraph->at(readNumber)->begin(),parGraph->at(readNumber)->end(), compareEdges); // Sort the list of edges of the current node according to the overlap offset (ascending).
	return true;
}

/**********************************************************************************************************************
	Mark all the transitive edges of a read.
	For Details: E.W. Myers. The fragment assembly string graph. Bioinformatics, 21(suppl 2):ii79-ii85, 2005.
**********************************************************************************************************************/
bool OverlapGraph::markTransitiveEdges(UINT64 readNumber, map<UINT64, vector<Edge*> * > *parGraph)
{
	map<UINT64,markType> *markedNodes = new map<UINT64,markType>();
	for(UINT64 i = 0; i < parGraph->at(readNumber)->size(); i++){ // Mark all the neighbors of the current read as INPLAY
		markedNodes->insert(std::pair<UINT64,markType>(parGraph->at(readNumber)->at(i)->getDestinationRead()->getReadNumber(), INPLAY));
	}
	for(UINT64 i = 0; i < parGraph->at(readNumber)->size(); i++) // Traverse through the list of edges according to their overlap offset.
	{
		UINT64 read2 = parGraph->at(readNumber)->at(i)->getDestinationRead()->getReadNumber(); // For each neighbor
		if(markedNodes->at(read2) == INPLAY) 										// If the neighbor is marked as INPLAY
		{
			for(UINT64 j = 0; j < parGraph->at(read2)->size(); j++)
			{
				UINT64 read3 = parGraph->at(read2)->at(j)->getDestinationRead()->getReadNumber(); // Get the neighbors neighbors
				if(markedNodes->find(read3) != markedNodes->end() && markedNodes->at(read3) == INPLAY)
				{
					UINT8 type1 = parGraph->at(readNumber)->at(i)->getOrientation();
					UINT8 type2 = parGraph->at(read2)->at(j)->getOrientation();
					if((type1 == 0 ||  type1 == 2) && (type2==0 || type2==1)) 	// Check edge orientation
						markedNodes->at(read3) = ELIMINATED; 					// Mark as ELIMINATED
					else if((type1==1||type1==3) && (type2==2 || type2==3)) 	// Check edge orientation
						markedNodes->at(read3) = ELIMINATED; 					// Mark as ELIMINATED
				}
			}
		}
	}
	for(UINT64 i = 0;i < parGraph->at(readNumber)->size(); i++)
	{
		if(markedNodes->at(parGraph->at(readNumber)->at(i)->getDestinationRead()->getReadNumber()) == ELIMINATED) // Current read to a node marked as ELIMINATED
		{
			parGraph->at(readNumber)->at(i)->transitiveRemovalFlag = true; 					// Mark this edge as transitive edge. Will remove this edge later.
			parGraph->at(readNumber)->at(i)->getReverseEdge()->transitiveRemovalFlag = true;	// Mark also the reverse edge. Will remove this edge later.
		}
	}
	delete markedNodes;
	return true;
}



/**********************************************************************************************************************
	Remove all transitive edges of a given read.
	For Details: E.W. Myers. The fragment assembly string graph. Bioinformatics, 21(suppl 2):ii79-ii85, 2005.
**********************************************************************************************************************/
bool OverlapGraph::removeTransitiveEdges(UINT64 readNumber, map<UINT64, vector<Edge*> * > *parGraph)
{
	for(UINT64 index = 0; index < parGraph->at(readNumber)->size(); index++)  		// Go through the list of edges of the current read.
	{
		if(parGraph->at(readNumber)->at(index)->transitiveRemovalFlag == true)		// This edge is marked as transitive. We will first remove the reverese edge.
		{
			Edge *twinEdge = parGraph->at(readNumber)->at(index)->getReverseEdge();
			UINT64 ID = twinEdge->getSourceRead()->getReadNumber();
			for(UINT64 index1 = 0; index1 < parGraph->at(ID)->size(); index1++) 	// Get the reverse edge first
			{
				if(parGraph->at(ID)->at(index1) == twinEdge)
				{
					delete twinEdge;
					parGraph->at(ID)->at(index1) = parGraph->at(ID)->at(parGraph->at(ID)->size()-1); // Move the transitive edges at the back of the list and remove.
					parGraph->at(ID)->pop_back();
					break;
				}
			}
		}
	}
	UINT64 j=0;
	for(UINT64 index=0; index < parGraph->at(readNumber)->size(); index++) // Then we will remove all the transitive edges of the current read.
	{
		if(parGraph->at(readNumber)->at(index)->transitiveRemovalFlag == false)		// We move all the non-transitive edges at the beginning of the list
			parGraph->at(readNumber)->at(j++) = parGraph->at(readNumber)->at(index);
		else		// Free the transitive edge
			delete parGraph->at(readNumber)->at(index);
	}
	parGraph->at(readNumber)->resize(j);
	return true;
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


/**********************************************************************************************************************
	Save the partial  graph in a text file
**********************************************************************************************************************/
bool OverlapGraph::saveParGraphToFile(string fileName, map<UINT64,nodeType> * exploredReads,map<UINT64, vector<Edge*> * > *parGraph)
{
	//CLOCKSTART;
	ostringstream strstream;
	for (map<UINT64, vector<Edge*> * >::iterator it=parGraph->begin(); it!=parGraph->end();)
	{
		UINT64 readID = it->first;
		if(!it->second->empty() && exploredReads->find(readID) !=  exploredReads->end())
		{
			if(exploredReads->at(readID) == EXPLORED_AND_TRANSITIVE_EDGES_REMOVED)
			{
				for(UINT64 j = 0; j < it->second->size(); j++)	// for each edge of the node
				{
					vector<UINT64> list;
					Edge * e = it->second->at(j);
					Edge *twinEdge = it->second->at(j)->getReverseEdge();
					UINT64 source = e->getSourceRead()->getReadNumber();
					UINT64 destination = e->getDestinationRead()->getReadNumber();
					if(source < destination || (source == destination && e < e->getReverseEdge()))
					{
						UINT64 srcLen =e->getSrcLen();
						list.push_back(e->getSourceRead()->getFileIndex());	// store the edge information first
						list.push_back(e->getDestinationRead()->getFileIndex());
						list.push_back(e->getOrientation());
						list.push_back(srcLen - e->getOverlapOffset());  //overlap length
						list.push_back(0);					//no substitutions
						list.push_back(0);					//no edits
						//Source Read (len,start,stop)
						list.push_back(srcLen);
						list.push_back(e->getOverlapOffset());
						list.push_back(srcLen-1);
						//Destination Read (len,start,stop)
						list.push_back(e->getDestLen());
						list.push_back(0);
						list.push_back(srcLen - e->getOverlapOffset()-1);

						//Check if destination is also marked by this thread.
						// 0: Only source is marked
						// 1: Only destination is marked
						// 2: Both source and destination are marked
						if(exploredReads->at(destination) == EXPLORED_AND_TRANSITIVE_EDGES_REMOVED)
							list.push_back(2);
						else
							list.push_back(0);
					}
					else
					{
						UINT64 srcLen = twinEdge->getSrcLen();
						list.push_back(e->getDestinationRead()->getFileIndex());	// store the edge information first
						list.push_back(e->getSourceRead()->getFileIndex());
						list.push_back(twinEdge->getOrientation());
						list.push_back(srcLen - twinEdge->getOverlapOffset());  //overlap length
						list.push_back(0);					//no substitutions
						list.push_back(0);					//no edits
						//Source Read (len,start,stop)
						list.push_back(srcLen);
						list.push_back(twinEdge->getOverlapOffset());
						list.push_back(srcLen-1);
						//Destination Read (len,start,stop)
						list.push_back(twinEdge->getDestLen());
						list.push_back(0);
						list.push_back(srcLen - twinEdge->getOverlapOffset()-1);
						//Check if destination is also marked by this thread.
						// 0: Only source is marked
						// 1: Only destination is marked
						// 2: Both source and destination are marked
						if(exploredReads->at(destination) == EXPLORED_AND_TRANSITIVE_EDGES_REMOVED)
							list.push_back(2);
						else
							list.push_back(1);
					}
					//write to file
					if(list.size()>0)
					{
						strstream<<list.at(0)<<"\t";
						strstream<<list.at(1)<<"\t";
						for(UINT64 i = 2; i < list.size()-1; i++)	// store in a file for future use.
							strstream<<list.at(i)<<",";
						strstream<<"NA,"<<list.at(list.size()-1)<<'\n';
					}
					//remove twin edges
					UINT64 twinID = twinEdge->getSourceRead()->getReadNumber();
					for(UINT64 index1 = 0; index1 < parGraph->at(twinID)->size(); index1++) 	// Get the reverse edge first
					{
						if(parGraph->at(twinID)->at(index1) == twinEdge)
						{
							delete twinEdge;
							parGraph->at(twinID)->at(index1) = parGraph->at(twinID)->at(parGraph->at(twinID)->size()-1); // Move the transitive edges at the back of the list and remove.
							parGraph->at(twinID)->pop_back();
							break;
						}
					}

				}
				//remove edges
				for(UINT64 j = 0; j< parGraph->at(readID)->size(); j++)
				{
					delete parGraph->at(readID)->at(j);
				}
				delete parGraph->at(readID);
				parGraph->erase(it++);
				exploredReads->at(readID) = EXPLORED_AND_TRANSITIVE_EDGES_WRITTEN;
			}
			else
				++it;
		}
		else
			++it;
	}
	ofstream filePointer;
	filePointer.open(fileName.c_str(), std::ios_base::app);
	if(!filePointer)
		MYEXIT("Unable to open file: "+fileName);
	filePointer<<strstream.str();
	//CLOCKSTOP;
	return true;
}

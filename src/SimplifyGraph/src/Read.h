#ifndef READ_H
#define READ_H

//============================================================================
// Name        : Read.h
// Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey, Abhishek Biswas
// Version     : v1.2
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : Read header file
//============================================================================

#include "Config.h"
#include "dna.h"

class Edge;

enum EdgeOP {INSERTION, DELETION};
typedef std::pair<Edge* , UINT32> t_edge_loc_pair;
class Read
{
	private:
		/* ====================  DATA MEMBERS  ======================================= */

		//dna_bitset *m_seq;                   // String representation of the read.


		/* array of pairs, with each pair storing
		 * first: pointer to an edge that includes this read
		 * second: index(location) of this read on the edge and orientation
		*/
		Edge ** edgeP;						//Array of pointers to the edges containing the read
		UINT32 * edgeOriIndex;				//MSB: Orientation of the read sequence in the Edge,
											//31 LSB : Index of the read on the edge starting from source
		UINT64 *containedReads;				//Array of contained reads(2MSB:Read Orientation, next 30 MSB: overlap start, 32LSB: contained read ID, )

		UINT64 m_readID;         // Unique Identification of the read.

		UINT32 readLen;          //We store only length of the read in the read class to save memory.
		                         //Actual read sequence can be stored using the dna_bitset above that is commented out
								 //Right now all read sequence is streamed

		UINT16 noOfEdges;					//Number of edges this read belongs to
		UINT16 noOfAllocEdgeMemAvail;		//Number of allocated edge locations available that have been deleted
											//This is used to avoid allocation and deallocation when the edges are to be updated...


		UINT16 noOfConReads;				//Number of contained reads

		bool containedReadFlag;				// True if read is contained...

		bool usedRead;						// True if read is used by previous assembly...

		/* ====================  METHODS      ======================================= */
		void initEdgeInfo();

	public:
		/* ====================  LIFECYCLE     ======================================= */
		Read();

		Read(const UINT64 length);

		//Expose only is read sequence is being stored
		//Read(const std::string & seq);

		~Read(void);

		/* ====================  OPERATORS     ======================================= */
		Read& operator=(const Read &s_read);

		/* ====================  DATA MEMBERS  ======================================= */
		friend class Edge;

		/* ====================  MUTATORS      ======================================= */
		void setReadID(const UINT64 & ID){m_readID = ID;}

		void setSeq(const std::string & seq);

		void setIsContained(bool isCon) { containedReadFlag=isCon; }

		void setConRead(UINT64 conReadID, UINT32 conOvlStart, UINT64 orient);

		void setEdge(Edge *edge, UINT32 readIndx, UINT32 orient);

		void delEdge(Edge *edge, UINT64 readIndx, UINT64 orient);

		void ClearEdgeInfo();

		void setUsedRead(bool val){
			usedRead=val;
		}

		/* ====================  ACCESSORS     ======================================= */ 

		//If you are going to store the read sequence then use these 3 functions
		//std::string getStringForward(void) const {return m_seq->toString();}
		//std::string getStringReverse(void) const {return m_seq->toRevComplement();}
		//UINT32 getReadLength(void) const {return m_seq->getLength();}

		UINT32 getReadLength(void) const {return readLen;}

		UINT64 getReadID(void) const {return m_readID;}

		bool isContainedRead() const { return containedReadFlag; }

		UINT16 getContainedReadCount() const { return noOfConReads; }

		UINT64 getContainedReadID(UINT64 indx) const
		{
			return (containedReads[indx] & 0X00000000FFFFFFFF);
		}
		UINT32 getContainedReadOvlStart(UINT64 indx) const
		{
			return ((containedReads[indx] >> 32) & 0X000000003FFFFFFF);
		}
		UINT8 getContainedReadOrientation(UINT64 indx) const
		{
			return (containedReads[indx] >> 62);
		}

		UINT32 getEdgeReadIndx(UINT64 indx) const
		{
			return (edgeOriIndex[indx] & 0X7FFFFFFF);
		}
		UINT8 getEdgeOriIndx(UINT64 indx) const
		{
			return (edgeOriIndex[indx] >> 31);
		}

		vector<t_edge_loc_pair> getFwdEdges() const;

		vector<t_edge_loc_pair> getBwdEdges() const;

		bool isUsedRead() {
			return usedRead;
		}

};/* End of class Read */

std::ostream &operator<<(std::ostream & out, const Read & read);

#endif /* READS_H */

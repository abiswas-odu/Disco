/*
 * HashTable.h
 *
 *  Created on: Apr 22, 2013
 *      Author: b72
 */

#ifndef HASHTABLE_H_
#define HASHTABLE_H_
#include "../src/Common.h"
#include "../src/Dataset.h"

#define BIG_CONSTANT(x) (x##LLU)

/* useful constants */
enum
{
	BASE_A = 0x0,	/* binary: 00 */
	BASE_C = 0x1,	/*'binary: 01 */
	BASE_G = 0x2,	/* binary: 10 */
	BASE_T = 0x3,	/* binary: 11 */
};

#define BASE_MASK 0x00000003	/* binary: 11 */


/**********************************************************************************************************************
	Class to store hashtable.
**********************************************************************************************************************/
class HashTable{
	private:
		Dataset *dataSet;							// Pointer of the dataset. this is NOT modified here
		UINT64 hashTableSize; 						// Ted: Size of the hash table. This is the prime number of mod operation.
		UINT64 hashDataTableSize; 						// Ted: Size of the hash data table. This is based on the number of reads.
		UINT64 *hashTable; 							// AB: List of hash offset of each hash key location
		UINT64 *hashData; 							// AB: List of hash keys of read number (ID) and orientation.
		UINT16 hashStringLength;					// Ted: Length of prefix and suffix of the reads to hash. This is equal to the minumum overlap length.
		mutable UINT64 numberOfHashCollision;		// Counter to count the number of hash collisions. For debugging only.
													// It's mutable such that it can be modified in the const member function, getListOfReads
		bool insertIntoTable(Read *read, string forwardRead, UINT64 *hashDataLengths);	// Insert a string in the hash table.
		bool hashReadLengths(string forwardRead); 					// Ted: Hash prefix and suffix of the read and its reverse complement in the hash table. Turn over to the constant
		void setHashTableSize(UINT64 size); 		// Set the size of the hash table.
		void setHashTableDataSize(UINT64 size);		// Set the size of the hash data table.
		string reverseComplement(const std::string & seq) const;
		UINT64 getHashIndex(const string & subString) const;	//Convert a set of data bytes from hash data table to DNA string
		string toString(UINT64 hashDataIndex,UINT64 stringLen) const;
	public:
		HashTable(void);							// Default constructor.
		~HashTable();								// Destructor.
		void createOffsetTable();
		bool insertDataset(Dataset *d, UINT64 minOverlapLength,UINT64 parallelThreadPoolSize);	// Insert the dataset in the hash table.
		vector<UINT64> * getListOfReads(const string & subString) const; 			// Get the list of reads that contain subString as prefix or suffix.
		UINT64 hashFunction(const string & subString) const; 						// Hash function.
		UINT64 getHashTableSize(void) const {return hashTableSize;}		// Get the size of the hash table.
		UINT64 getHashStringLength() const {return hashStringLength;}		// Get the hash string length.
		Dataset * getDataset(void) const {return dataSet;}					// Get the pointer to the dataset.

		string getStringForward(UINT64 offset) const; 									// Get the forward string of the read at offset.
		string getStringReverse(UINT64 offset) const;  								// Get the reverse string of the read at offset.
		UINT64 getReadLength(UINT64 offset) const; 								// Get the length of the string in the read at offset.

		void readReadLengthsFromFile(string fileName, UINT64 minOverlap);
		void populateReadLengths();												//Populate the read lengths in the hash table for future offset calculation
		void populateReadData();												//Populate the read sequence in the hash data
		void readReadSequenceFromFile(string fileName, UINT64 minOverlap, UINT64 *hashDataLengths, UINT64 &readID);

};


#endif /* HASHTABLE_H_ */

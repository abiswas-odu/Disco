/*
 * Reads.h
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 */


#ifndef READ_H_
#define READ_H_

#include "../src/Common.h"

/**********************************************************************************************************************
	Class to store a read.
**********************************************************************************************************************/

class Read
{
	private:
		UINT64 readNumber; 						// Unique Identification of the read.
		UINT64 fileIndex; 					// The sequence number of the read in the files...
		UINT64 readHashOffset; 					// The offset number of the read in the hash data...
	public:
		UINT64 superReadID;						// 0 = not a contained read
												// otherwise superReadID contains the ID of the uniqe super read.
		Read(void);								// Default constructor.
		Read(UINT64 fIndx);					// Another constructor.
		~Read(void);							// Destructor.

		bool setReadNumber(UINT64 id); 			// Set the read number.
		void setFileIndex(UINT64 id); 			// Set the file index number.
		void setReadHashOffset(UINT64 offset); 			// Set the hash data table offset.
		UINT64 getReadNumber(void) const {return readNumber;} 								// Get the read number of the current read.
		UINT64 getFileIndex(void) const {return fileIndex;} 								// Get the fileIndex number of the current read.
		UINT64 getReadHashOffset(void)const {return readHashOffset;}  			// Get the hash data table offset.
};

#endif /* READS_H_ */

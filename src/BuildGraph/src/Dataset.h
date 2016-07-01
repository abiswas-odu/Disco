/*
 * Dataset.h
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 */


#ifndef DATASET_H_
#define DATASET_H_
#include "../src/Common.h"
#include "../src/Read.h"


/**********************************************************************************************************************
	Class to store the dataset
**********************************************************************************************************************/
class Dataset
{
	private:
		UINT64 numberOfReads;								// Number of total reads present in the dataset.
		UINT64 numberOfUniqueReads; 						// number of unique reads in the dataset.
		UINT64 minimumOverlapLength;						// Length of the shortest read in the dataset.
		vector<Read *> *reads; 								// List of reads in the dataset.
		string reverseComplement(const string & read); 		// Get the reverse complement of a string.
															// The dataset contains only good quality reads.

	public:
		vector<string> pairedEndDatasetFileNames;
		vector<string> singleEndDatasetFileNames;
		UINT64 shortestReadLength;
		UINT64 longestReadLength;

		bool testRead(const string & read); 				// Test if the read is good. Contains only {A,C,G,T} and does not contain more than 80% of same base.

		Dataset(void); 										// Default constructor.
		Dataset(vector<string> pairedEndFileNames, vector<string> singleEndFileNames,string fileNamePrefix, UINT64 minOverlap);// anotherconstructor.
		~Dataset(void);										// Default destructor.
		bool readDataset(string fileName, UINT64 minOverlap, UINT64 datasetNumber, UINT64 &fIndx); // Read the database from fasta/fastq file. Matepairs should be one after another in the file.
		UINT64 getNumberOfReads(void); 						// Get the number of total reads in the database.
		UINT64 getNumberOfUniqueReads(void); 				// Get the number of unique reads in the database.
		bool printDataset(void); 							// Print few the reads in the dataset. For debuggin only.
        Read * getReadFromID(UINT64 ID); 					// Find a read in the database given the ID in constant time.
        Read * getReadFromFileIndex(UINT64 fID);
};


#endif /* DATASET_H_ */

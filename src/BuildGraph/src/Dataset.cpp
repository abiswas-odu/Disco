/*
 * Dataset.cpp
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider, Abhishek Biswas
 */

#include "Dataset.h"
#include "Common.h"
#ifdef INCLUDE_READGZ
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)
#endif

/**********************************************************************************************************************
	Default constructor
**********************************************************************************************************************/
Dataset::Dataset(void)
{
	// Initialize the variables.
	numberOfUniqueReads = 0;
	numberOfReads = 0;
	minimumOverlapLength = 0;
	shortestReadLength = 0XFFFFFFFFFFFFFFFF;
	longestReadLength = 0X0000000000000000;
	reads = new vector<Read *>;

}


/**********************************************************************************************************************
	Another constructor
**********************************************************************************************************************/
Dataset::Dataset(vector<string> pairedEndFileNames, vector<string> singleEndFileNames,
		string fileNamePrefix, UINT64 minOverlap, UINT64 maxThreads)
{
	// Initialize the variables.
	numberOfUniqueReads = 0;
	numberOfReads = 0;
	shortestReadLength = 0XFFFFFFFFFFFFFFFF;
	longestReadLength = 0X0000000000000000;
	reads = new vector<Read *>;
	pairedEndDatasetFileNames = pairedEndFileNames;
	singleEndDatasetFileNames = singleEndFileNames;
	minimumOverlapLength = minOverlap;
	UINT64 counter = 0, fileIndex=0;

	filterStrings = {"ACACACACACACACACACACACACACACA",
			"AGAGAGAGAGAGAGAGAGAGAGAGAGAGA",
			"ATATATATATATATATATATATATATATA",
			"CGCGCGCGCGCGCGCGCGCGCGCGCGCGC",
			"CTCTCTCTCTCTCTCTCTCTCTCTCTCTC",
			"AAGAAGAAGAAGAAGAAGAAGAAGAAGAA",
			"ATAATAATAATAATAATAATAATAATAAT",
			"TAATAATAATAATAATAATAATAATAATA",
			"AACAACAACAACAACAACAACAACAACAA",
			"ACAACAACAACAACAACAACAACAACAAC",
			"CAACAACAACAACAACAACAACAACAACA",
			"AAGAAGAAGAAGAAGAAGAAGAAGAAGAA",
			"AGAAGAAGAAGAAGAAGAAGAAGAAGAAG",
			"GAAGAAGAAGAAGAAGAAGAAGAAGAAGA",
			"TTCTTCTTCTTCTTCTTCTTCTTCTTCTT",
			"AAATAAATAAATAAATAAATAAATAAATA",
			"TAAATAAATAAATAAATAAATAAATAAAT",
			"ATAAATAAATAAATAAATAAATAAATAAA",
			"AATAAATAAATAAATAAATAAATAAATAA",
			"AATTAATTAATTAATTAATTAATTAATTA",
			"ATTAATTAATTAATTAATTAATTAATTAA",
			"TTAATTAATTAATTAATTAATTAATTAAT",
			"TAATTAATTAATTAATTAATTAATTAATT",
			"AAAGAAAGAAAGAAAGAAAGAAAGAAAGA",
			"AAAGAAAGAAAGAAAGAAAGAAAGAAAGA",
			"AGAAAGAAAGAAAGAAAGAAAGAAAGAAA",
			"GAAAGAAAGAAAGAAAGAAAGAAAGAAAG",
			"TACATACATACATACATACATACATACAT",
			"ACATACATACATACATACATACATACATA",
			"CATACATACATACATACATACATACATAC",
			"ATACATACATACATACATACATACATACA",
			"GTTTGTTTGTTTGTTTGTTTGTTTGTTTG",
			"TGTTTGTTTGTTTGTTTGTTTGTTTGTTT",
			"TTTGTTTGTTTGTTTGTTTGTTTGTTTGT",
			"AGGGAGGGAGGGAGGGAGGGAGGGAGGGA",
			"GAGGGAGGGAGGGAGGGAGGGAGGGAGGG",
			"GGAGGGAGGGAGGGAGGGAGGGAGGGAGG",
			"GGGAGGGAGGGAGGGAGGGAGGGAGGGAG"};

	merCheckStrings={"AC","AG","AT","CG","CT","GT","AAT","ATA","TAA","AAC","ACA","CAA","AAG","AGA","GAA","GGGGCC"};

//					"TGAGCCCGGCC","TTGGCCCGGCC","TTAGCCCGGCC","TAAGGCCGGGC","TCAGCCCGGCC","TGAGGCCGGGC","TAGCCCGGCCC",
//	                "TAGGGCCGGGC","TCAGGCCGGGC","TGGCCCGGCCC","TACGGCCGGGGC","TAGCCCCGGCCG","AGCCCGGCCCG","TCGGGCCGGGC",
//	                "AAGGCCGGGCC","TAAGGCCGCGCCCC","TAGGCCCGCGGC","TAGGGCCGAGGC","TTAGGCCGCGGGCCGCCGCCC","TAAGGCCCGGGC",
//	                "TAAGGCCGGGGC","TTAGCCCCGGGCC","TAGCCCGGCCG","TCAGGCCGGGGC","TACGGCCGGGC","TTAGGGGCGCGGCC","TAGCCCGGCGGCGG",
//	                "TGCGGCCGGGC","AGCCCGGCCGC","TGAGGCCCGGGC","TAGCCTCGGGCGGGCC","TAAGGCCGCGCCCCA","TTATGGGGCGCGGCC","TCGGCCCCAGCC",
//	                "TAAGGCCGCCCCGCCGGGGCC","TTAGGCCCCGGCGGGGCGGCC","TTACGGGGCCGGGCC","AGGGCCGGGCG","AGGGCCGGGCC","TCGGGGGCGCGGCC",
//					"TGAGGCCCCGGGC","AAGGCCGGGGC","TAACGGCCGGGC","TGAGGGCCGGGC","AAGGCCGGGCG","TTCGCCCGGCCT","TGAAGGCCGGGC","TTACGCCCGGCC",
//	                "TAGCCGCC","TTAGCCCGGGCC","TAAGGGCCGCGCCCC","TAAGGGCCGCGGGC","TGAGGCCGCCGGGC","TCAGCCCGGCGGCC","TAAGGGCCGGGC",
//	                "TAGCCCGGGCCG","TTAGGGGCGCGGCGCCGGCC","TTGGGGCCCGGCC","TAAGGCTAGGCC","TTAGGCCTAGCC","TTAGGGCCGGGGCGCGC","TAAGGCCACGGGC",
//	                "TAAGGCCGCGGCGCTGCGGGC","TCGAGGCTAGGC","TTAGGGGCCGGCC","TAAGGCCGCCCC","TCAGCGCCCGGCC","TTAGCCCGGGGCC","TACGGCCGGCCGGGC",
//	                "TTGGGCCCGGCC","TGAGGCCGCGC","TCAGCGCGGCC","TGGGCCCGGCCG","TCAGGGGCGCGGCC","TAGCCCCGGGCCC","TTAGGGGCGGCCGGCC"};

	parallelThreadPoolSize=maxThreads;

	ofstream filePointer;
	string fileName = fileNamePrefix+"_ReadIDMap.txt";
	filePointer.open(fileName.c_str());
	if(!filePointer)
		MYEXIT("Unable to open file: " + fileName);


	UINT64 startReadID=rQStart=fileIndex;
	readDataset(singleEndDatasetFileNames.at(0), minimumOverlapLength, counter++, fileIndex);

	if(fileIndex <= startReadID)
		MYEXIT("File empty. No reads loaded from "+ singleEndDatasetFileNames.at(0));
	filePointer<<singleEndDatasetFileNames.at(0)<<": Singleton file "<<1<<"\nReadID Range: ("<<startReadID+1<<",";
	filePointer<<fileIndex<<")\n";
	rQEnd=fileIndex-1;


	startReadID=rSStart=fileIndex;
	readDataset(singleEndDatasetFileNames.at(1), minimumOverlapLength, counter++, fileIndex);

	if(fileIndex <= startReadID)
		MYEXIT("File empty. No reads loaded from "+ singleEndDatasetFileNames.at(1));
	filePointer<<singleEndDatasetFileNames.at(1)<<": Singleton file "<<2<<"\nReadID Range: ("<<startReadID+1<<",";
	filePointer<<fileIndex<<")\n";
	rSEnd=fileIndex-1;

	filePointer.close();
	cout << endl << "Shortest read length in all datasets: " << setw(5) << shortestReadLength<<endl;
	cout << " Longest read length in all datasets: " << setw(5) << longestReadLength <<endl;

	for(UINT64 i = 0 ; i < reads->size(); i++) 		// Assing ID's to the reads.
		reads->at(i)->setReadNumber(i + 1);
	numberOfUniqueReads=reads->size();
	reads->shrink_to_fit();

	if(numberOfUniqueReads<=0)						//Exit of no reads found...
		MYEXIT("No reads found in the read files provided! Please check if the filename(s) and path(s) are correct.");

	//Create a file index to readID lookup table. Used to load previous partial results in case of a restart...
	fIndxReadIDMap = new map<UINT64, UINT64>;
	for(UINT64 i = 1; i <= getNumberOfUniqueReads(); i++)
	{
		UINT64 fIndx = getReadFromID(i)->getFileIndex();
		auto it = fIndxReadIDMap->end();
		fIndxReadIDMap->insert(it, pair<UINT64,UINT64>(fIndx,i));
	}
}

/*
 * Deallocate the memory storing the map between file index and read id after its not needed
 */
void Dataset::freeFindexReadIDMAP()
{
	delete fIndxReadIDMap;
}
/**********************************************************************************************************************
	This function reads the dataset from FASTA/FASTQ files
**********************************************************************************************************************/
bool Dataset::readDataset(string fileName, UINT64 minOverlap, UINT64 datasetNumber, UINT64 &fIndx)
{
	CLOCKSTART;
	cout << "Reading dataset: " << datasetNumber << " from file: " << fileName << endl;
	UINT64 goodReads = 0, badReads = 0;
	map<UINT64,string> readList;
	if(fileName.substr( fileName.length() - 3 )==".gz")
	{
#ifdef INCLUDE_READGZ
		gzFile fp;
		kseq_t *seq;
		int l;
		fp = gzopen(fileName.c_str(), "r");
		seq = kseq_init(fp);
		#pragma omp parallel num_threads(parallelThreadPoolSize)
		{
			#pragma omp single
			{
				while ((l = kseq_read(seq)) >= 0) {
					if( (goodReads + badReads ) != 0 && (goodReads + badReads)%1000000 == 0)
						cout<< setw(10) << goodReads + badReads << " read processed added in hashtable. " << setw(10) << goodReads << " good reads." << setw(10) << badReads << " bad reads." << endl;
					string line1=seq->seq.s;
					fIndx++;							//Increment file index of the read
					if(readList.size()>=READ_TASK_BLOCK)
					{
						readList.insert(std::pair<UINT64, string>(fIndx, line1));
						#pragma omp task firstprivate(readList)
						{
							for(auto iterator = readList.begin(); iterator != readList.end(); iterator++) {
								UINT64 fileIndexOfRead=iterator->first;
								string seqLine=iterator->second;
								for (std::string::iterator p = seqLine.begin(); seqLine.end() != p; ++p) // Change the case
									*p = toupper(*p);
								if(seqLine.length() > minOverlap && testRead(seqLine) ) // Test the read is of good quality.
								{
									UINT64 len = seqLine.length();
									Read *r1=new Read(fileIndexOfRead,len);


									if(len > longestReadLength)
									{
										#pragma omp atomic write
											longestReadLength = len;
									}
									if(len < shortestReadLength)
									{
										#pragma omp atomic write
											shortestReadLength = len;
									}

									#pragma omp critical
									{
										reads->push_back(r1);						// Store the first string in the dataset.
									}
									#pragma omp atomic
										goodReads++;
								}
								else
								{
									#pragma omp atomic
										badReads++;
								}
							}
						}
						readList.clear();
					}
					else
						readList.insert(std::pair<UINT64, string>(fIndx, line1));
				}
				#pragma omp taskwait
			}
		}
		kseq_destroy(seq);
		gzclose(fp);
#else
		MYEXIT("Unknown input file format. Looks like the file is in gzip compressed format."
				"The Disco code was not built with ZLIB using READGZ=1. To assemble either uncompress"
				"the file or build Disco with ZLIB library using make \"make READGZ=1\".");
#endif
	}
	else
	{
		ifstream myFile;
		myFile.open (fileName.c_str());
		if(!myFile)
			MYEXIT("Unable to open file: "+fileName)

		vector<string> line;
		string text;
		enum FileType { FASTA, FASTQ, UNDEFINED};
		FileType fileType = UNDEFINED;
		#pragma omp parallel num_threads(parallelThreadPoolSize)
		{
			#pragma omp single
			{
				while(getline(myFile,text))
				{
					string line1="",line0="";
					if( (goodReads + badReads ) != 0 && (goodReads + badReads)%1000000 == 0)
						cout<< setw(10) << goodReads + badReads << " reads processed in dataset " << setw(2) << datasetNumber <<". " << setw(10) << goodReads << " good reads." << setw(10) << badReads << " bad reads." << endl;
					if(fileType == UNDEFINED)
					{
						if(text[0] == '>')
							fileType = FASTA;
						else if(text[0] == '@')
							fileType = FASTQ;
						else
							MYEXIT("Unknown input file format.");
					}
					line.clear();
					if(fileType == FASTA) 			// Fasta file
					{
						line.push_back(text);
						getline (myFile,text,'>');
						line.push_back(text);

						line.at(1).erase(std::remove(line.at(1).begin(), line.at(1).end(), '\n'), line.at(1).end());
						line.at(0).erase(std::remove(line.at(0).begin(),line.at(0).end(),'>'),line.at(0).end());			//Sequence name
						line0=line.at(0);
						line1 = line.at(1);								// The first string is in the 2nd line.

					}
					else if(fileType == FASTQ) 					// Fastq file.
					{
						line.push_back(text);
						for(UINT64 i = 0; i < 3; i++) 	// Read the remaining 3 lines. Total of 4 lines represent one sequence in a fastq file.
						{
							getline (myFile,text);
							line.push_back(text);
						}
						line.at(0).erase(std::remove(line.at(0).begin(),line.at(0).end(),'>'),line.at(0).end());			//Sequence name
						line0=line.at(0);
						line1 = line.at(1); 			// The first string is in the 2nd line.
					}
					fIndx++;							//Increment file index of the read
					if(readList.size()>=READ_TASK_BLOCK)
					{
						readList.insert(std::pair<UINT64, string>(fIndx, line1));
						#pragma omp task firstprivate(readList)
						{
							for(auto iterator = readList.begin(); iterator != readList.end(); iterator++) {
								UINT64 fileIndexOfRead=iterator->first;
								string seqLine=iterator->second;
								for (std::string::iterator p = seqLine.begin(); seqLine.end() != p; ++p) // Change the case
									*p = toupper(*p);
								if(seqLine.length() > minOverlap && testRead(seqLine) ) // Test the read is of good quality.
								{
									UINT64 len = seqLine.length();
									Read *r1=new Read(fileIndexOfRead,len);


									if(len > longestReadLength)
									{
										#pragma omp atomic write
											longestReadLength = len;
									}
									if(len < shortestReadLength)
									{
										#pragma omp atomic write
											shortestReadLength = len;
									}

									#pragma omp critical
									{
										reads->push_back(r1);						// Store the first string in the dataset.
									}
									#pragma omp atomic
										goodReads++;
								}
								else
								{
									#pragma omp atomic
										badReads++;
								}
							}
						}
						readList.clear();
					}
					else
						readList.insert(std::pair<UINT64, string>(fIndx, line1));
				}
				#pragma omp taskwait
			}
		}
		myFile.close();
	}
	//Insert remaining reads left in the map
	for(auto iterator = readList.begin(); iterator != readList.end(); iterator++) {
		UINT64 fileIndexOfRead=iterator->first;
		string seqLine=iterator->second;
		for (std::string::iterator p = seqLine.begin(); seqLine.end() != p; ++p) // Change the case
			*p = toupper(*p);
		if(seqLine.length() > minOverlap && testRead(seqLine) ) // Test the read is of good quality.
		{
			UINT64 len = seqLine.length();
			Read *r1=new Read(fileIndexOfRead,len);

			if(len > longestReadLength)
			{
				longestReadLength = len;
			}
			if(len < shortestReadLength)
			{
				shortestReadLength = len;
			}
			reads->push_back(r1);						// Store the first string in the dataset.
			goodReads++;
		}
		else
		{
				badReads++;
		}
	}
	numberOfReads+=goodReads;				// Counter of the total number of reads.
    cout << endl << "Dataset: " << setw(2) << datasetNumber << endl;
    cout << "File name: " << fileName << endl;
	cout << setw(10) << goodReads << " good reads in current dataset."  << endl;
	cout << setw(10) << badReads << " bad reads in current dataset." << endl;
	cout << setw(10) << goodReads + badReads << " total reads in current dataset." << endl;
	cout << setw(10) << numberOfReads << " good reads in all datasets." << endl << endl;;
	CLOCKSTOP;
	return true;
}


/**********************************************************************************************************************
	This function returns the number of reads
**********************************************************************************************************************/
UINT64 Dataset::getNumberOfReads(void)
{
	return numberOfReads;
}


/**********************************************************************************************************************
	This function returns the number of unique reads
**********************************************************************************************************************/
UINT64 Dataset::getNumberOfUniqueReads(void)
{
	return numberOfUniqueReads;
}

/**********************************************************************************************************************
	Returns true if the read contains only {A,C,G,T} and does not contain more than 80% of the same nucleotide
**********************************************************************************************************************/
bool Dataset::testRead(const string & read)
{
	UINT64 cnt[4] = {0,0,0,0};
	UINT64 readLength = read.length();
	if(readLength < MIN_READ_SIZE)
		return false;
	for(UINT64 i = 0; i < readLength; i++) // Count the number of A's, C's , G's and T's in the string and check no other characters exist
	{
		if(read[i]!= 'A' && read[i] != 'C' && read[i] != 'G' && read[i] != 'T')
			return false;
		cnt[(read[i] >> 1) & 0X03]++; // Least significant 2nd and 3rd bits of ASCII value used here
	}
	UINT64 threshold = read.length()*.7;	// 70% of the length.
	if(cnt[0] >= threshold || cnt[1] >= threshold || cnt[2] >= threshold || cnt[3] >= threshold)
		return false;	// If 80% bases are the same base.

	//This loop filters out reads with unusually micro-repeating kmers at start or end
	for(size_t i=0;i<filterStrings.size();i++)
	{
		UINT64 len=filterStrings[i].size();
		if(read.size()<len)				//Read must be at-least as long as the filter strings
			return false;
		if(filterStrings[i] == read.substr(0,len))
			return false;
		if(filterStrings[i] == read.substr(readLength-len))
			return false;
	}
	//Check if dimers or trimers exist in very large numbers
	threshold = readLength*.5;	// 50% of the length.
	for(size_t i=0;i<merCheckStrings.size();i++)
	{
		size_t repeatCtr = countSubstring(read, merCheckStrings[i]);
		repeatCtr = repeatCtr * merCheckStrings[i].length();
		if(repeatCtr >= threshold)		//Check if repeat appears over more than 50% of the read's sequence
			return false;				//One of the dimers/trimers exist in large micro repeats
	}
//	//Check if some longer k-mers exist in very large numbers
//	for(size_t k=4;k < minimumOverlapLength;k++) //for substring upto overlap len
//	{
//		for(UINT64 j = 0; j < readLength - k; j++) //for each substring
//		{
//			string subString = read.substr(j,k); // Get the substring
//			size_t repeatLenCtr = countSubstring(read, subString);
//			repeatLenCtr = repeatLenCtr * k;
//			if(repeatLenCtr >= threshold)		//Check if repeat appears over more than 50% of the read's sequence
//				return false;				//One of the k-mers exist in large micro repeats
//		}
//	}
	return true;
}

/**********************************************************************************************************************
	Returns the reverse complement of a read
**********************************************************************************************************************/
string Dataset::reverseComplement(const string & read)
{
	UINT64 stringLength = read.length();
	string reverse(stringLength,'0');
	for(UINT64 i = 0;i < stringLength; i++)						// Then complement the string. Change A to T, C to G, G to C and T to A.
	{
		if( read[i] & 0X02 ) // C or G
			reverse.at(stringLength -  i - 1 ) = read[i] ^ 0X04;
		else // A <==> T
			reverse.at(stringLength -  i - 1 ) = read[i] ^ 0X15;
	}
	return reverse; // return the reverse complement as a string
}



/**********************************************************************************************************************
	This function returns a read from its ID.
**********************************************************************************************************************/
Read * Dataset::getReadFromID(UINT64 ID)
{
	if(ID < 1 || ID > numberOfUniqueReads)	// ID outside the range.
	{
		stringstream ss;
		ss << "ID " << ID << " out of bound.";
		string s = ss.str();
		MYEXIT(s);
	}
	return reads->at(ID - 1);
}

/**********************************************************************************************************************
	This function returns a read from its FileIndex.
**********************************************************************************************************************/
Read * Dataset::getReadFromFileIndex(UINT64 fID)
{
	UINT64 readID=0;
	size_t i = fID<reads->size()?fID:reads->size();
	for(;i>0;i--)
	{
		Read *r = getReadFromID(i);
		if(r->getFileIndex()==fID)
		{
			readID = r->getReadNumber();
			break;
		}
	}
	return getReadFromID(readID);
}

/**********************************************************************************************************************
	Default destructor
**********************************************************************************************************************/
Dataset::~Dataset(void)
{
	// Free the memory used by the dataset.
	for(UINT64 i = 0; i < reads->size(); i++)
	{
		delete reads->at(i);
	}
	delete reads;
}


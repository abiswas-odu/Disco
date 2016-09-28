/*
 * Dataset.cpp
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider, Abhishek Biswas
 */

#include "Dataset.h"
#include "Common.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

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
Dataset::Dataset(vector<string> pairedEndFileNames, vector<string> singleEndFileNames, string fileNamePrefix, UINT64 minOverlap)
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

	ofstream filePointer;
	string fileName = fileNamePrefix+"_ReadIDMap.txt";
	filePointer.open(fileName.c_str());
	if(!filePointer)
		MYEXIT("Unable to open file: " + fileName);

	for(UINT64 i = 0; i < pairedEndDatasetFileNames.size(); i++)						// Read the paired-end datasets.
	{
		UINT64 startReadID=fileIndex;
		readDataset(pairedEndDatasetFileNames.at(i), minimumOverlapLength, counter++, fileIndex);
		if(fileIndex <= startReadID)
			MYEXIT("File empty. No reads loaded from "+ pairedEndDatasetFileNames.at(i));
		filePointer<<pairedEndDatasetFileNames.at(i)<<": Paired-end file "<<i+1<<"\nReadID Range: ("<<startReadID+1<<",";
		filePointer<<fileIndex<<")\n";
	}

	for(UINT64 i = 0; i < singleEndDatasetFileNames.size(); i++)						// Read the single-end datasets.
	{
		UINT64 startReadID=fileIndex;
		readDataset(singleEndDatasetFileNames.at(i), minimumOverlapLength, counter++, fileIndex);

		if(fileIndex <= startReadID)
			MYEXIT("File empty. No reads loaded from "+ singleEndDatasetFileNames.at(i));
		filePointer<<singleEndDatasetFileNames.at(i)<<": Singleton file "<<i+1<<"\nReadID Range: ("<<startReadID+1<<",";
		filePointer<<fileIndex<<")\n";
	}
	filePointer.close();
	cout << endl << "Shortest read length in all datasets: " << setw(5) << shortestReadLength<<endl;
	cout << " Longest read length in all datasets: " << setw(5) << longestReadLength <<endl;

	for(UINT64 i = 0 ; i < reads->size(); i++) 		// Assing ID's to the reads.
		reads->at(i)->setReadNumber(i + 1);
	numberOfUniqueReads=reads->size();
	reads->shrink_to_fit();
}
/**********************************************************************************************************************
	This function reads the dataset from FASTA/FASTQ files
**********************************************************************************************************************/
bool Dataset::readDataset(string fileName, UINT64 minOverlap, UINT64 datasetNumber, UINT64 &fIndx)
{
	CLOCKSTART;
	cout << "Reading dataset: " << datasetNumber << " from file: " << fileName << endl;
	UINT64 goodReads = 0, badReads = 0;
	if(fileName.substr( fileName.length() - 3 )==".gz")
	{
		gzFile fp;
		kseq_t *seq;
		int l;
		fp = gzopen(fileName.c_str(), "r");
		seq = kseq_init(fp);
		while ((l = kseq_read(seq)) >= 0) {
			if( (goodReads + badReads ) != 0 && (goodReads + badReads)%1000000 == 0)
				cout<< setw(10) << goodReads + badReads << " read processed added in hashtable. " << setw(10) << goodReads << " good reads." << setw(10) << badReads << " bad reads." << endl;
			string line1=seq->seq.s;
			fIndx++;							//Increment file index of the read
			for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
				*p = toupper(*p);
			if(line1.length() > minOverlap && testRead(line1) ) // Test the read is of good quality.
			{
				Read *r1=new Read(fIndx);
				UINT64 len = line1.length();
				if(len > longestReadLength)
					longestReadLength = len;
				if(len < shortestReadLength)
					shortestReadLength = len;
				reads->push_back(r1);						// Store the first string in the dataset.
				numberOfReads++;							// Counter of the total number of reads.
				goodReads++;
			}
			else
				badReads++;
		}
		kseq_destroy(seq);
		gzclose(fp);
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
			for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
				*p = toupper(*p);
			if(line1.length() > minOverlap && testRead(line1) ) // Test the read is of good quality.
			{
				Read *r1=new Read(fIndx);
				UINT64 len = line1.length();
				if(len > longestReadLength)
					longestReadLength = len;
				if(len < shortestReadLength)
					shortestReadLength = len;

				reads->push_back(r1);						// Store the first string in the dataset.
				numberOfReads++;							// Counter of the total number of reads.
				goodReads++;
			}
			else
				badReads++;
		}

		myFile.close();
	}
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
	for(UINT64 i = 0; i < readLength; i++) // Count the number of A's, C's , G's and T's in the string.
	{
		if(read[i]!= 'A' && read[i] != 'C' && read[i] != 'G' && read[i] != 'T')
			return false;
		cnt[(read[i] >> 1) & 0X03]++; // Least significant 2nd and 3rd bits of ASCII value used here
	}
	UINT64 threshold = read.length()*.8;	// 80% of the length.
	if(cnt[0] >= threshold || cnt[1] >= threshold || cnt[2] >= threshold || cnt[3] >= threshold)
		return false;	// If 80% bases are the same base.
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

//============================================================================
// Name        : DataSet.cpp
// Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey
// Version     : v1.2
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : DataSet cpp file
//============================================================================

#include "DataSet.h"

void DataSet::loadReadsFromReadFile(const std::string &read_file)
{
	CLOCKSTART;
	FILE_LOG(logINFO) << "load reads from read file: " << read_file << "\n";
	UINT64 readID(m_vec_reads->size() + 1);

	// To count of reads in this file
	UINT64 readCount = 0;

	// Open file
	ifstream filePointer;
	filePointer.open(read_file.c_str());
	if(!filePointer.is_open()){
		FILE_LOG(logWARNING) << "Unable to open file: " << read_file << "\n";
		return;
	}

	// Variables
	string text;
	enum FileType {FASTA, FASTQ, UNDEFINED};
	FileType fileType = UNDEFINED;

	while(getline(filePointer,text))
	{
		string line0="",line1="";
		// Check FASTA and FASTQ
		if(fileType == UNDEFINED) {
			if (text.length() > 0){
				if(text[0] == '>'){
					FILE_LOG(logINFO) << "Input reads file format: FASTA\n";
					fileType = FASTA;
				}
				else if(text[0] == '@'){
					FILE_LOG(logINFO) << "Input reads file format: FASTA\n";
					fileType = FASTQ;
				}
				else{
					FILE_LOG(logERROR) << "Unknown input file format."<<endl;
					break;
				}
			}
		}
		line0=text;	// get ID line
		// FASTA file read
		if(fileType == FASTA) {
			getline (filePointer,line1,'>');	// get string line
			line1.erase(std::remove(line1.begin(), line1.end(), '\n'), 
					line1.end());
			//std::transform(line1.begin(), line1.end(), line1.begin(), ::toupper);
		}
		// FASTQ file read
		else if(fileType == FASTQ) {
			getline(filePointer, line1);	// String
			// Ignore the next two lines
			filePointer.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			filePointer.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}
		transform(line1.begin(), line1.end(), line1.begin(), ::toupper);
		if(!testRead(line1)) // Test the read is of good quality. If not replace all N's
			std::replace( line1.begin(), line1.end(), 'N', 'A');
		Read *r = new Read(line1);
		r->setReadID(readID);
		// add read to the dataset
		addRead(r);
		++readID;
		++readCount;
		if(readCount % 1000000 == 0){
			FILE_LOG(logDEBUG) << setw(10) << (readCount/1000000)  << ",000,000"  
				<< " reads loaded to memory, "
				<< setw(7) << checkMemoryUsage() << " MB\n";
		}
	}
	filePointer.close();
	FILE_LOG(logDEBUG) << setw(10) << readCount << " reads loaded from this read file\n";
	CLOCKSTOP;
}
DataSet::DataSet(const vector<std::string> &read_SingleFiles,const vector<std::string> &read_PairFiles,
		vector<std::string> &read_PairInterFiles)
{
	CLOCKSTART;
	m_vec_reads = new vector<Read *>;
	dataSetInfo = new vector<DataSetInfo>;
	isSequenceLoaded=false;
	UINT64 dsCount=0;
	//Read seperated paired files
	for(auto it = read_PairFiles.cbegin(); it != read_PairFiles.cend(); ++it){
		DataSetInfo newDataSet;
		newDataSet.datasetNumber=dsCount++;
		newDataSet.isPaired=1;
		newDataSet.isReadInterleaved=0;
		newDataSet.matePairOrientation=0;
		newDataSet.r1Start = m_vec_reads->size()+1;
		//Load read information
		loadReadsFromReadFile(*it);
		newDataSet.r1End = m_vec_reads->size();
		newDataSet.r2Start = m_vec_reads->size() + 1;
		newDataSet.r1FileName = *it;
		it++;
		//Load read information
		loadReadsFromReadFile(*it);
		newDataSet.r2End = m_vec_reads->size();
		newDataSet.r2FileName = *it;
		dataSetInfo->push_back(newDataSet);
	}
	//Read interleaved paired files
	for(auto it = read_PairInterFiles.cbegin(); it != read_PairInterFiles.cend(); ++it){
			DataSetInfo newDataSet;
			newDataSet.datasetNumber=dsCount++;
			newDataSet.isPaired=1;
			newDataSet.isReadInterleaved=1;
			newDataSet.matePairOrientation=0;
			newDataSet.r1Start = m_vec_reads->size()+1;
			//Load read information
			loadReadsFromReadFile(*it);
			newDataSet.r1End = m_vec_reads->size();
			newDataSet.r1FileName = *it;
			newDataSet.r2FileName = "";
			newDataSet.r2Start = 0;
			newDataSet.r2End = 0;
			dataSetInfo->push_back(newDataSet);
	}
	//Read single read files
	for(auto it = read_SingleFiles.cbegin(); it != read_SingleFiles.cend(); ++it){
		DataSetInfo newDataSet;
		newDataSet.datasetNumber=dsCount++;
		newDataSet.isPaired=0;
		newDataSet.matePairOrientation=0;
		newDataSet.r1Start = m_vec_reads->size()+1;
		loadReadsFromReadFile(*it);
		newDataSet.r1End = m_vec_reads->size();
		dataSetInfo->push_back(newDataSet);
		newDataSet.r2Start = 0;
		newDataSet.r2End = 0;
	}
	m_vec_reads->shrink_to_fit();
	CLOCKSTOP;
}
/*
 * Load used read from previous rounds of assembly.
 */
UINT64 DataSet::LoadUsedReads(string usedReadFileName)
{
	// Open file and read used reads is available...
	ifstream usedReadFilePointer;
	UINT64 usedReadCtr=0;
	usedReadFilePointer.open(usedReadFileName.c_str());
	if(!usedReadFilePointer.is_open()){
		FILE_LOG(logWARNING) << "No used read file present: " << usedReadFileName << "\n";
	}
	else
	{
		string text;
		while(getline(usedReadFilePointer,text))
		{
			UINT64 readID = stoi(text);
			if(!at(readID)->isUsedRead())
			{
				at(readID)->setUsedRead(true);
				usedReadCtr++;
			}
		}
		FILE_LOG(logINFO)<< SSTR(usedReadCtr) << " used reads loaded."<<endl;
	}
	return usedReadCtr;
}
DataSet::~DataSet()
{
	if (m_vec_reads != nullptr){
		// The DataSet class should not delete reads
		for(UINT64 i = 0; i < m_vec_reads->size(); i++)
		{
			if (m_vec_reads->at(i) != nullptr){
				delete m_vec_reads->at(i);
				m_vec_reads->at(i) = nullptr;
			}
		}
		delete m_vec_reads;
		m_vec_reads = nullptr;
	}
}

std::ostream& operator<< (std::ostream &out, DataSet & a_data_set)
{
	out << "Dataset with size " << a_data_set.size() << endl;
	return out;
}

DataSet& DataSet::operator= (const DataSet &s_dataset)
{
	if(this == &s_dataset)
		return *this;
	delete m_vec_reads;
	m_vec_reads = new vector<Read*>;
	m_vec_reads->reserve((s_dataset.m_vec_reads)->size());
	for(auto it = s_dataset.m_vec_reads->cbegin(); it != s_dataset.m_vec_reads->cend(); ++it){
		m_vec_reads->push_back(*it);
	}
	return *this;
}

void DataSet::rmRead(Read *r)
{
	m_vec_reads->erase(std::remove(m_vec_reads->begin(), m_vec_reads->end(), r),m_vec_reads->end());
}


Read * DataSet::at(UINT64 ID) const
{
	assert(ID >= 1 && ID <= m_vec_reads->size());
	return m_vec_reads->at(ID - 1);
}

/**********************************************************************************************************************
	Returns true if the read contains only {A,C,G,T} and does not contain more than 80% of the same nucleotide
**********************************************************************************************************************/
bool DataSet::testRead(const string & read)
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
	This function reads the contained read file generated during graph construction and store matepair information.
**********************************************************************************************************************/
void DataSet::storeContainedReadInformation(vector<string> containedReadFile)
{
	CLOCKSTART;
	FILE_LOG(logINFO)<< "Store contained read information..."<< endl;
	for(UINT64 i=0;i<containedReadFile.size();i++)
	{
		ifstream myFile;
		myFile.open(containedReadFile[i].c_str());
		if(!myFile)
			MYEXIT("Unable to open file: "+containedReadFile[i]);

		string text;
		while(getline (myFile,text))
		{
			vector<string> toks = Utils::split(text,'\t');
			UINT64 containedReadID = atoi(toks[0].c_str());
			UINT64 containingReadID = atoi(toks[1].c_str());

			vector<string> info = Utils::split(toks[2],',');
			UINT64 containedReadOri = atoi(info[0].c_str());
			UINT64 containedReadOverlapStart = atoi(info[8].c_str());

			at(containedReadID)->setIsContained(true);

			at(containingReadID)->setConRead(containedReadID, containedReadOverlapStart, containedReadOri);
		}
		myFile.close();
	}
	CLOCKSTOP;
}

/*
 * Returns dataset number for a given read ID
 *
 */
UINT64 DataSet::getDataSetNumber(UINT64 readID)
{
	for(size_t i=0; i < dataSetInfo->size();i++)
	{
		if((readID>=dataSetInfo->at(i).r1Start && readID<=dataSetInfo->at(i).r1End)
				|| (readID>=dataSetInfo->at(i).r2Start && readID<=dataSetInfo->at(i).r2End))
			return dataSetInfo->at(i).datasetNumber;
	}
	MYEXIT(readID + "not found in any dataset.");
	return 0;
}


/*
 * Get coverage of a given base index of a given read
 */
UINT32 DataSet::getReadCoverage(UINT64 readID, UINT64 indx) const
{
	Read *read = at(readID);
	UINT32 cov=1;
	if(!read->isContainedRead() && !read->getContainedReadCount()==0)
	{
		for(size_t i=0;i<read->getContainedReadCount();i++)
		{
			UINT64 conRID = read->getContainedReadID(i);
			UINT64 overlapStart = read->getContainedReadOvlStart(i);
			UINT32 conReadLen = at(conRID)->getReadLength();
			if(overlapStart>indx && indx < (overlapStart+conReadLen))
				cov++;
		}
	}
	return cov;
}
/****
 * Returns the mate id of a read, 0 otherwise
 */
UINT64 DataSet::getMatePair(UINT64 r1ID)
{
	UINT64 r2ID=0;
	for(size_t i=0; i < dataSetInfo->size();i++)
	{
		if(dataSetInfo->at(i).isPaired)
		{
			if(r1ID>=dataSetInfo->at(i).r1Start && r1ID<=dataSetInfo->at(i).r1End)
			{
				if(dataSetInfo->at(i).isReadInterleaved)			//Check if the read is from an interleaved read file.
				{
					int testForR = (r1ID-dataSetInfo->at(i).r1Start);
					if (testForR % 2)  //If odd reverse read
						r2ID = r1ID-1;
					else				//If odd forward read
						r2ID = r1ID+1;
					if(r2ID!=0 && !at(r2ID)->isContainedRead())
						return r2ID;
					break;
				}
				else
				{
					UINT64 r2ID = (r1ID-dataSetInfo->at(i).r1Start) + dataSetInfo->at(i).r2Start;
					if(!at(r2ID)->isContainedRead())
						return r2ID;
					break;
				}
			}
			else if(r1ID>=dataSetInfo->at(i).r2Start && r1ID<=dataSetInfo->at(i).r2End)  //Interleaved file check not required as r2Start and r2End are 0s in interleaved files
			{
				r2ID = (r1ID-dataSetInfo->at(i).r2Start) + dataSetInfo->at(i).r1Start;
				if(!at(r2ID)->isContainedRead())
					return r2ID;
				break;
			}
		}
	}
	return r2ID;
}

vector<UINT64> DataSet::getMatePairList(Read *read)
{
	vector<UINT64> mpList;
	UINT64 r1ID = read->getReadID();
	UINT64 r2ID = getMatePair(r1ID);
	if(r2ID!=0)
		mpList.push_back(r2ID);
	for(UINT64 i=0;i<read->getContainedReadCount();i++)
	{
		UINT64 c1ID = read->getContainedReadID(i);
		UINT64 c2ID = getMatePair(c1ID);
		if(c2ID!=0)
			mpList.push_back(c2ID);
	}
	return mpList;
}


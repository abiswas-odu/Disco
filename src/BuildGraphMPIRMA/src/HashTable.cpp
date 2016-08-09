/*
 * HashTable.cpp
 *
 *  Created on: Apr 22, 2013
 *      Author: b72
 */


#include "../../BuildGraphMPIRMA/src/HashTable.h"

#include "../../BuildGraphMPIRMA/src/Common.h"



/**********************************************************************************************************************
	Function to get a prime number larger than a given number.
	For efficiency of the hash function the hash table size is set to a prime number.
	This will reduce the number of hash collision.
**********************************************************************************************************************/

UINT64 getPrimeLargerThanNumber(UINT64 number)
{
	// Pre-computed list of prime number.
	UINT64 array[450] = {1114523, 1180043, 1245227, 1310759, 1376447, 1442087, 1507379, 1573667, 1638899, 1704023, 1769627, 1835027, 1900667, 1966127, 2031839, 2228483, 2359559, 2490707, 2621447, 2752679, 2883767, 3015527, 3145739, 3277283, 3408323, 3539267, 3670259, 3801143, 3932483, 4063559, 4456643, 4718699, 4980827, 5243003, 5505239, 5767187, 6029603, 6291563, 6553979, 6816527, 7079159, 7340639, 7602359, 7864799, 8126747, 8913119, 9437399, 9962207, 10485767, 11010383, 11534819, 12059123, 12583007, 13107923, 13631819, 14156543, 14680067, 15204467, 15729647, 16253423, 17825999, 18874379, 19923227, 20971799, 22020227, 23069447, 24117683, 25166423, 26214743, 27264047, 28312007, 29360147, 30410483, 31457627, 32505983, 35651783, 37749983, 39845987, 41943347, 44040383, 46137887, 48234623, 50331707, 52429067, 54526019, 56623367, 58720307, 60817763, 62915459, 65012279, 71303567, 75497999, 79691867, 83886983, 88080527, 92275307, 96470447, 100663439, 104858387, 109052183, 113246699, 117440699, 121635467, 125829239, 130023683, 142606379, 150994979, 159383759, 167772239, 176160779, 184549559, 192938003, 201327359, 209715719, 218104427, 226493747, 234882239, 243269639, 251659139, 260047367, 285215507, 301989959, 318767927, 335544323, 352321643, 369100463, 385876703, 402654059, 419432243, 436208447, 452986103, 469762067, 486539519, 503316623, 520094747, 570425399, 603979919, 637534763, 671089283, 704643287, 738198347, 771752363, 805307963, 838861103, 872415239, 905971007, 939525143, 973079279, 1006633283, 1040187419, 1140852767, 1207960679, 1275069143, 1342177379, 1409288183, 1476395699, 1543504343, 1610613119, 1677721667, 1744830587, 1811940419, 1879049087, 1946157419, 2013265967, 2080375127, 2281701827, 2415920939, 2550137039, 2684355383, 2818572539, 2952791147, 3087008663, 3221226167, 3355444187, 3489661079, 3623878823, 3758096939, 3892314659, 4026532187, 4160749883, 4563403379, 4831838783, 5100273923, 5368709219, 5637144743, 5905580687, 6174015503, 6442452119, 6710886467, 6979322123, 7247758307, 7516193123, 7784629079, 8053065599, 8321499203, 9126806147, 9663676523, 10200548819, 10737418883, 11274289319, 11811160139, 12348031523, 12884902223, 13421772839, 13958645543, 14495515943, 15032386163, 15569257247, 16106127887, 16642998803, 18253612127, 19327353083, 20401094843, 21474837719, 22548578579, 23622320927, 24696062387, 25769803799, 26843546243, 27917287907, 28991030759, 30064772327, 31138513067, 32212254947, 33285996803, 36507222923, 38654706323, 40802189423, 42949673423, 45097157927, 47244640319, 49392124247, 51539607599, 53687092307, 55834576979, 57982058579, 60129542339, 62277026327, 64424509847, 66571993199, 73014444299, 77309412407, 81604379243, 85899346727, 90194314103, 94489281203, 98784255863, 103079215439, 107374183703, 111669150239, 115964117999, 120259085183, 124554051983, 128849019059, 133143986399, 146028888179, 154618823603, 163208757527, 171798693719, 180388628579, 188978561207, 197568495647, 206158430447, 214748365067, 223338303719, 231928234787, 240518168603, 249108103547, 257698038539, 266287975727, 292057776239, 309237645803, 326417515547, 343597385507, 360777253763, 377957124803, 395136991499, 412316861267, 429496730879, 446676599987, 463856468987, 481036337207, 498216206387, 515396078039, 532575944723, 584115552323, 618475290887, 652835029643, 687194768879, 721554506879, 755914244627, 790273985219, 824633721383, 858993459587, 893353198763, 927712936643, 962072674643, 996432414899, 1030792152539, 1065151889507, 1168231105859, 1236950582039, 1305670059983, 1374389535587, 1443109012607, 1511828491883, 1580547965639, 1649267441747, 1717986918839, 1786706397767, 1855425872459, 1924145348627, 1992864827099, 2061584304323, 2130303780503, 2336462210183, 2473901164367, 2611340118887, 2748779070239, 2886218024939, 3023656976507, 3161095931639, 3298534883999, 3435973836983, 3573412791647, 3710851743923, 3848290698467, 3985729653707, 4123168604483, 4260607557707, 4672924419707, 4947802331663, 5222680234139, 5497558138979, 5772436047947, 6047313952943, 6322191860339, 6597069767699, 6871947674003, 7146825580703, 7421703488567, 7696581395627, 7971459304163, 8246337210659, 8521215117407, 9345848837267, 9895604651243, 10445360463947, 10995116279639, 11544872100683, 12094627906847, 12644383722779, 13194139536659, 13743895350023, 14293651161443, 14843406975659, 15393162789503, 15942918604343, 16492674420863, 17042430234443, 18691697672867, 19791209300867, 20890720927823, 21990232555703, 23089744183799, 24189255814847, 25288767440099, 26388279068903, 27487790694887, 28587302323787, 29686813951463, 30786325577867, 31885837205567, 32985348833687, 34084860462083, 37383395344739, 39582418600883, 41781441856823, 43980465111383, 46179488367203, 48378511622303, 50577534878987, 52776558134423, 54975581392583, 57174604644503, 59373627900407, 61572651156383, 63771674412287, 65970697666967, 68169720924167, 74766790688867, 79164837200927, 83562883712027, 87960930223163, 92358976733483, 96757023247427, 101155069756823, 105553116266999, 109951162779203, 114349209290003, 118747255800179, 123145302311783, 127543348823027, 131941395333479, 136339441846019, 149533581378263, 158329674402959, 167125767424739, 175921860444599, 184717953466703, 193514046490343, 202310139514283, 211106232536699, 219902325558107, 228698418578879, 237494511600287, 246290604623279, 255086697645023, 263882790666959, 272678883689987, 299067162755363, 316659348799919, 334251534845303, 351843720890723, 369435906934019, 387028092977819, 404620279022447, 422212465067447, 439804651111103, 457396837157483, 474989023199423, 492581209246163, 510173395291199, 527765581341227, 545357767379483, 598134325510343, 633318697599023, 668503069688723, 703687441776707, 738871813866287, 774056185954967, 809240558043419, 844424930134187, 879609302222207, 914793674313899, 949978046398607, 985162418489267, 1020346790579903, 1055531162666507, 1090715534754863}; 	//list of 450 sorted prime numbers;
	vector<UINT64> primeNumbers(array, array + sizeof array / sizeof array[0]);
	for(UINT64 i = 0; i < primeNumbers.size(); i++)
		if(primeNumbers.at(i)> number)
			return primeNumbers.at(i);		// Return the smallest prime in the list larger than number.
	return number + 1;
}




/**********************************************************************************************************************
	Default Constructor
**********************************************************************************************************************/
HashTable::HashTable(UINT64 parallelProcessPoolSize)
{
	// Initialize the variables.
	hashStringLength = 0;
	numberOfHashCollision = 0;
	memoryDataPartitions = NULL;
	memoryReadCount=NULL;
	dataSet=NULL;
	hashTable = NULL;
	hashData = NULL;
	win = MPI_WIN_NULL;
	cachedHashTable=NULL;
	hashQ=NULL;
	cachedReadTable=NULL;
	readQ=NULL;
	rmaCtr=0;
}




/**********************************************************************************************************************
	Add Database in hashTable
**********************************************************************************************************************/
void HashTable::insertDataset(Dataset* d, UINT64 minOverlapLength, UINT64 parallelProcessPoolSize, int myid)
{
	CLOCKSTART;
	dataSet=d;
	hashStringLength = minOverlapLength - 1;
	numberOfHashCollision = 0;
	UINT64 size = getPrimeLargerThanNumber(d->getNumberOfUniqueReads() * 8 + 1);  // Size should be at least twice the number of entries in the hash table to reduce hash collision.
	setHashTableSize(size);

	//Init the caches
	cachedHashTable = new map<UINT64, std::shared_ptr<UINT64> >();
	hashQ = new queue<UINT64>;

	cachedReadTable = new map<UINT64, string>();
	readQ = new queue<UINT64>;

	populateReadLengths(); 				//Each hash entry is populated with the length of it's record

	//Change lengths of records to offset at each index...
	UINT64 nextOffset=hashTable[0];
	hashTable[0]=0;
	for(size_t i=1; i<hashTableSize; i++)
	{
		UINT64 tempOffset=nextOffset+hashTable[i-1];
		nextOffset=hashTable[i];
		hashTable[i]=tempOffset;
	}

	size_t hashDataTableSize = nextOffset + hashTable[hashTableSize-1];   //Calculate the size of the hash data table
	cout << "Hash Data Table size set to: " << hashDataTableSize << endl;
	// Defines the minimum number of hash records in each partition
	UINT64 minWindowSize =  hashDataTableSize/parallelProcessPoolSize;
	UINT64 currWinSize=0;
	// Populate partitioning vectors based on hash record and data boundaries.
	size_t rankIndx=0;
	memoryDataPartitions = new vector<UINT64>(parallelProcessPoolSize+1,0);
	memoryReadCount = new vector<UINT64>(parallelProcessPoolSize,0);
	memoryDataPartitions->at(rankIndx)=0;
	rankIndx++;
	for(size_t i=1; i<hashTableSize; i++)
	{
		currWinSize+=hashTable[i]-hashTable[i-1];
		if(currWinSize>=minWindowSize)
		{
			memoryDataPartitions->at(rankIndx)=hashTable[i];			//The data index upto which the MPI rank stores value (MPI rank == vector index)
			currWinSize=0;
			rankIndx++;
		}
	}
	//Update empty ranks with max index...
	for(;rankIndx<memoryDataPartitions->size();rankIndx++)
	{
		memoryDataPartitions->at(rankIndx)=hashDataTableSize;
	}
	cout<<"Hash table partition information:"<<endl;
	for(size_t i=0;i<memoryReadCount->size();i++)
	{
		cout<<"Data Part Offset:"<<memoryDataPartitions->at(i)<<" Read Count:"<<memoryReadCount->at(i)<<endl;
	}
	setHashTableDataSize(myid);
	populateReadData(myid);												//Populate the hash data table with reads
	CLOCKSTOP;
}

/**********************************************************************************************************************
	Insert a read lengths in the hashTable
**********************************************************************************************************************/
void HashTable::populateReadLengths()
{
	CLOCKSTART;
	for(UINT64 i = 0; i < dataSet->pairedEndDatasetFileNames.size(); i++)						// Read the paired-end datasets.
	{
		readReadLengthsFromFile(dataSet->pairedEndDatasetFileNames.at(i), hashStringLength+1);
	}
	for(UINT64 i = 0; i < dataSet->singleEndDatasetFileNames.size(); i++)						// Read the single-end datasets.
	{
		readReadLengthsFromFile(dataSet->singleEndDatasetFileNames.at(i), hashStringLength+1);
	}
	CLOCKSTOP;
}


/**********************************************************************************************************************
	Insert a read sequence in hashData
**********************************************************************************************************************/
void HashTable::populateReadData(int myid)
{
	CLOCKSTART;
	UINT64 *hashDataLengths = new UINT64[hashTableSize];
	std::memset(hashDataLengths, 0, hashTableSize*sizeof(UINT64));
	UINT64 readID=0;
	for(UINT64 i = 0; i < dataSet->pairedEndDatasetFileNames.size(); i++)						// Read the paired-end datasets.
	{
		readReadSequenceFromFile(dataSet->pairedEndDatasetFileNames.at(i), hashStringLength+1, hashDataLengths, readID ,myid);
	}

	for(UINT64 i = 0; i < dataSet->singleEndDatasetFileNames.size(); i++)						// Read the single-end datasets.
	{
		readReadSequenceFromFile(dataSet->singleEndDatasetFileNames.at(i), hashStringLength+1, hashDataLengths, readID, myid);
	}
	delete[] hashDataLengths;
	MPI_Win_fence(0, win);
	CLOCKSTOP;
}

/**********************************************************************************************************************
	Read sequence lengths into the hashTable
**********************************************************************************************************************/
void HashTable::readReadLengthsFromFile(string fileName, UINT64 minOverlap)
{
	CLOCKSTART;
	cout << "Reading read lengths from file: " << fileName << endl;
	ifstream myFile;
	myFile.open (fileName.c_str());
	if(!myFile)
		MYEXIT("Unable to open file: "+fileName)
	UINT64 goodReads = 0, badReads = 0;
	vector<string> line;
	string text;
	enum FileType { FASTA, FASTQ, UNDEFINED};
	FileType fileType = UNDEFINED;
	while(getline(myFile,text))
	{
		string line1="",line0="";
		if( (goodReads + badReads ) != 0 && (goodReads + badReads)%1000000 == 0)
			cout<< setw(10) << goodReads + badReads << " read lengths added in hashtable. " << setw(10) << goodReads << " good reads." << setw(10) << badReads << " bad reads." << endl;
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
		for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
		    *p = toupper(*p);
		if(line1.length() > minOverlap && dataSet->testRead(line1) ) // Test the read is of good quality.
		{
			hashReadLengths(line1); 								// Calculate the offset lengths of each hash table key in the hash table.
			goodReads++;
		}
		else
			badReads++;
	}

	myFile.close();
    cout << "File name: " << fileName << endl;
	cout << setw(10) << goodReads << " good reads in current file."  << endl;
	cout << setw(10) << badReads << " bad reads in current file." << endl;
	cout << setw(10) << goodReads + badReads << " total reads in current file." << endl;
	CLOCKSTOP;
}

/**********************************************************************************************************************
	Read sequence data into hashData
**********************************************************************************************************************/
void HashTable::readReadSequenceFromFile(string fileName, UINT64 minOverlap, UINT64 *hashDataLengths, UINT64 &readID, int myid)
{
	CLOCKSTART;
	cout << "Reading read data from file: " << fileName << endl;
	ifstream myFile;
	myFile.open (fileName.c_str());
	if(!myFile)
		MYEXIT("Unable to open file: "+fileName)
	UINT64 goodReads = 0, badReads = 0;
	vector<string> line;
	string text;
	enum FileType { FASTA, FASTQ, UNDEFINED};
	FileType fileType = UNDEFINED;
	while(getline(myFile,text))
	{
		string line1="",line0="";
		if( (goodReads + badReads ) != 0 && (goodReads + badReads)%1000000 == 0)
			cout<< setw(10) << goodReads + badReads << " read sequence added in hashtable. " << setw(10) << goodReads << " good reads." << setw(10) << badReads << " bad reads." << endl;
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
		for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
		    *p = toupper(*p);
		if(line1.length() > minOverlap && dataSet->testRead(line1) ) // Test the read is of good quality.
		{
			readID++;
			cout<<"Proc:"<<myid<<"Processing read:"<<readID<<endl;
			insertIntoTable(dataSet->getReadFromID(readID), line1 ,hashDataLengths, myid); 								// Calculate the offset lengths of each hash table key in the hash table.
			goodReads++;
		}
		else
			badReads++;
	}

	myFile.close();
    cout << "File name: " << fileName << endl;
	cout << setw(10) << goodReads << " good reads in current file."  << endl;
	cout << setw(10) << badReads << " bad reads in current file." << endl;
	cout << setw(10) << goodReads + badReads << " total reads in current file." << endl;
	CLOCKSTOP;
}
/**********************************************************************************************************************
	Insert a read length in the hashTable
**********************************************************************************************************************/
bool HashTable::hashReadLengths(string forwardRead)
{
	string prefixForward = forwardRead.substr(0,hashStringLength); 											// Prefix of the forward string.
	string suffixForward = forwardRead.substr(forwardRead.length() - hashStringLength,hashStringLength);	// Suffix of the forward string.

	/* number of bytes necessary to store dna_str as a bitset */
	size_t dna_word = (forwardRead.length() / 32) + (forwardRead.length() % 32 != 0);

	UINT64 index = getHashIndex(prefixForward);						// Get the index using the hash function.
	hashTable[index]+=(dna_word+1);

	index = getHashIndex(suffixForward);
	hashTable[index]+=(dna_word+1);

	return true;
}

/**********************************************************************************************************************
	Set the hashTableSize
**********************************************************************************************************************/
void HashTable::setHashTableSize(UINT64 size)
{
	cout << "Hash Table size set to: " << size << endl;
	hashTableSize=size;
    hashTable = new UINT64[hashTableSize];
	std::memset(hashTable, 0, hashTableSize*sizeof(UINT64));
}


/**********************************************************************************************************************
	Set the Hash Data Size
**********************************************************************************************************************/
void HashTable::setHashTableDataSize(int myid)
{
	int numElements=memoryDataPartitions->at(myid+1)-memoryDataPartitions->at(myid);
	//hashData = new UINT64[numElements];
	cout<<"MPI UINT64:"<<sizeof(MPI_LONG_LONG_INT)<<endl;
	hashData = NULL;
	MPI_Alloc_mem(numElements*sizeof(MPI_UINT64_T), MPI_INFO_NULL, (void **)&hashData);
	MPI_Win_create(hashData, numElements*sizeof(MPI_UINT64_T), sizeof(MPI_UINT64_T), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
 	//MPI_Win_allocate(numElements * sizeof(MPI_UINT64_T), sizeof(MPI_UINT64_T), MPI_INFO_NULL, MPI_COMM_WORLD, hashData, &win);
	std::memset(hashData, 0, numElements*sizeof(MPI_UINT64_T));
	MPI_Win_fence(0, win);
}

/**********************************************************************************************************************
	Returns the hash value of a subString
**********************************************************************************************************************/
UINT64 HashTable::getHashIndex(const string & subString) const
{
	UINT64 hashIndxFwd = hashFunction(subString);
	UINT64 hashIndxRev =  hashFunction(reverseComplement(subString));
	if(hashIndxFwd<=hashIndxRev)
		return hashIndxFwd;
	else
		return hashIndxRev;
}

/**********************************************************************************************************************
	Returns the hash value of a subString
**********************************************************************************************************************/
UINT64 HashTable::hashFunction(const string & subString) const
{
	UINT64 sum1 = 1, sum2 = 1, length = subString.length();
	for(UINT64 i = 0; i < length; i++)			// We take the bit representation of the string. A = 00, C = 01, G = 11 and T = 10
	{
		if(i < 32)
			// sum1 is for the first 32 bp. bit shifted to left 2 bits.
			// Change the character to integer. A=Ox41=01000001
			//                                  C=0x43=01000011
			//                                  G=0x47=01000111
			//                                  T=0x54=01010100
			// Then, shift to right way 1 bit.
			// Then, bit and operation with 00000011
			// Then, it just have A=00, C=01, G=11,T=10
			sum1 = (sum1 << 2) | ( ( (int)(subString[i] ) >> 1 ) & 0X03 );
		else
			sum2 = (sum2 << 2) | ( ( (int)(subString[i] ) >> 1 ) & 0X03 );
	}
	// multiply two mod results, and mod again.
	return ((sum1 % hashTableSize) * (sum2  % hashTableSize)) % hashTableSize; 	// Modulus operation to get the index in the hash table.
}

/**********************************************************************************************************************
	Insert a subString in a hashTable
**********************************************************************************************************************/
bool HashTable::insertIntoTable(Read *read, string forwardRead, UINT64 *hashDataLengths, int myid)
{
	string prefixForward = forwardRead.substr(0,hashStringLength); 											// Prefix of the forward string.
	string suffixForward = forwardRead.substr(forwardRead.length() - hashStringLength,hashStringLength);	// Suffix of the forward string.

	/* number of 64 bit words necessary to store dna_str as a bitset */
	UINT64 dna_word = (forwardRead.length() / 32) + (forwardRead.length() % 32 != 0);

	UINT64 prefixForwardID =  read->getReadNumber() | (forwardRead.length() << 48); 	// Most significant bit is used to store the orientation of the string in the read.
																							// 0 = 0 means prefix of the forward string.
																							// 1 = 1 means suffix of the forward string.
																							//Read string length is stored is next 15 most significant bits.
																							//Read ID is stored in the 48 least significant bits.
	UINT64 orient=1;
	UINT64 suffixForwardID =  read->getReadNumber() | (forwardRead.length() << 48) | (orient << 63);
																	// Most significant bit is used to store the orientation of the string in the read.
																	// 0 = 0 means prefix of the forward string.
																	// 1 = 1 means suffix of the forward string.
																	//Read string length is stored is next 15 most significant bits.
																	//Read ID is stored in the 48 least significant bits.

	/*store prefix forward data*/

	UINT64 index = getHashIndex(prefixForward);						// Get the index using the hash function.
	UINT64 baseOffset = hashTable[index];
	memoryReadCount->at(getOffsetRank(baseOffset))++;
	// Get the start offset of the hash value in the hashData table
	if(isGlobalOffsetInRange(baseOffset, myid))
	{
		UINT64 localBaseOffset = getLocalOffset(baseOffset,myid);
		UINT64 currentOffset = hashDataLengths[index];					// Get the current offset of the hash value in the hashData table
		hashData[localBaseOffset+currentOffset] = prefixForwardID;
		/* for each base of the DNA sequence */
		for (size_t i = 0; i < forwardRead.length(); i++)
		{
			UINT64 shift = (32*2-2) - 2*(i % 32);

			switch (forwardRead[i])
			{
			case 'A':
				hashData[localBaseOffset+currentOffset + 1 + i/32] |= (UINT64)BASE_A << shift;
				break;
			case 'C':
				hashData[localBaseOffset+currentOffset + 1 + i/32] |= (UINT64)BASE_C << shift;
				break;
			case 'G':
				hashData[localBaseOffset+currentOffset + 1 + i/32] |= (UINT64)BASE_G << shift;
				break;
			case 'T':
				hashData[localBaseOffset+currentOffset + 1 + i/32] |= (UINT64)BASE_T << shift;
				break;
			default:
				throw std::invalid_argument("invalid DNA base");
			}
		}
	}
	//store the offset of the read data in the read. This will be used to obtain the sequence instead of storing the read sequence in the read object.
	UINT64	readHashOffset = baseOffset+hashDataLengths[index];
	read->setReadHashOffset(readHashOffset);

	//Update data lengths to maintain global consistency of offset values
	hashDataLengths[index] +=  (dna_word + 1);

	/*store prefix reverse data*/
	index = getHashIndex(suffixForward);						// Get the index using the hash function.
	baseOffset = hashTable[index];							// Get the start offset of the hash value in the hashData table
	memoryReadCount->at(getOffsetRank(baseOffset))++;
	if(isGlobalOffsetInRange(baseOffset, myid))
	{
		UINT64 localBaseOffset = getLocalOffset(baseOffset,myid);
		UINT64 currentOffset = hashDataLengths[index];					// Get the current offset of the hash value in the hashData table
		hashData[localBaseOffset+currentOffset] = suffixForwardID;
		/* for each base of the DNA sequence */
		for (size_t i = 0; i < forwardRead.length(); i++)
		{
			UINT64 shift = (32*2-2) - 2*(i % 32);
			switch (forwardRead[i])
			{
			case 'A':
				hashData[localBaseOffset+currentOffset + 1 + i/32] |= (UINT64)BASE_A << shift;
				break;
			case 'C':
				hashData[localBaseOffset+currentOffset + 1 + i/32] |= (UINT64)BASE_C << shift;
				break;
			case 'G':
				hashData[localBaseOffset+currentOffset + 1 + i/32] |= (UINT64)BASE_G << shift;
				break;
			case 'T':
				hashData[localBaseOffset+currentOffset + 1 + i/32] |= (UINT64)BASE_T << shift;
				break;
			default:
				throw std::invalid_argument("invalid DNA base");
			}
		}
	}
	//Update data lengths to maintain global consistency of offset values
	hashDataLengths[index] += (dna_word + 1);
	return true;
}
/**********************************************************************************************************************
	Returns true id the global offset of the process id
**********************************************************************************************************************/
bool HashTable::isGlobalOffsetInRange(UINT64 globalOffset, int myid) const
{
	if(globalOffset>=memoryDataPartitions->at(myid) && globalOffset<memoryDataPartitions->at(myid+1))
		return true;
	return false;

}
/**********************************************************************************************************************
	Returns the local offset of the read in the data table.
**********************************************************************************************************************/
UINT64 HashTable::getLocalOffset(UINT64 globalOffset, int myid) const
{
	UINT64 thisProcBaseHash = memoryDataPartitions->at(myid);
	return globalOffset-thisProcBaseHash;
}

/**********************************************************************************************************************
	Returns rank for a given offset value
**********************************************************************************************************************/
int HashTable::getOffsetRank(UINT64 globalOffset) const
{
	int rank=-1;
	for(size_t i=0;i<memoryDataPartitions->size()-1;i++)
	{
		if(isGlobalOffsetInRange(globalOffset,i))
			rank = i;
	}
	return rank;
}
/**********************************************************************************************************************
	Populate a list of reads from the hash table containing the subString as prefix or suffix. Do not cache result...
**********************************************************************************************************************/
vector<UINT64*> * HashTable::setLocalHitList_nocache(const string readString, int myid)
{
	vector<UINT64*> *localReadHits;				//AB: store data locally
	#pragma omp critical(getRemoteData)
	{
		try{
			UINT64 hashQueryCount = readString.length() - getHashStringLength();
			localReadHits = new vector<UINT64*>(hashQueryCount);
			for(UINT64 j = 0; j < readString.length() - getHashStringLength(); j++) // for each substring of read1 of length getHashStringLength
			{
				string subString = readString.substr(j,getHashStringLength()); // Get the substring from read string
				UINT64 index = getHashIndex(subString);	// Get the index using the hash function.
				//Compute various indices
				UINT64 globalStartOffset=hashTable[index];
				UINT64 globalEndOffset= (index==hashTableSize-1)?hashDataTableSize:hashTable[index+1];
				UINT64 hash_block_len = globalEndOffset-globalStartOffset;
				localReadHits->at(j)=NULL;
				if(hash_block_len)
				{
					int t_rank = getOffsetRank(globalStartOffset);
					MPI_Aint localOffset = getLocalOffset(globalStartOffset,t_rank);
					//cout<<"oRank:"<<myid<<", tRank:"<<t_rank<<", len:"<<hash_block_len<<", Offset:"<<localOffset<<endl;
					//Allocate memory
					localReadHits->at(j) = new UINT64[hash_block_len];
					std::memset(localReadHits->at(j), 0, hash_block_len*sizeof(MPI_UINT64_T));
					// Get data from RMA
					void *buf = localReadHits->at(j);
					if((localOffset+hash_block_len)>getMemoryMaxLocalOffset(t_rank))
					{
						cout<<"Req. Start: oProcess:"<<myid<<" tRank:"<<t_rank<<" GlobalOffset:"<<globalStartOffset<<" LocalOffset"<<localOffset<<" blockLen:"<<hash_block_len<<endl;
					}
					MPI_Get(buf, hash_block_len, MPI_UINT64_T, t_rank, localOffset, hash_block_len, MPI_UINT64_T, win);
					//MPI_Win_flush(t_rank, win);
					rmaCtr++;
				}
			}
			MPI_Win_flush_all(win);
		}
		catch (const std::bad_alloc&) {
			cout<<"Proc"<<myid<<" Exception Caught getting read data!!!";
			exit(0);
		}
	}
	return localReadHits;
}
/**********************************************************************************************************************
	Populate a list of reads from the hash table containing the subString as prefix or suffix.
**********************************************************************************************************************/
vector<shared_ptr<UINT64> > HashTable::setLocalHitList(const string readString, int myid)
{
	UINT64 hashQueryCount = readString.length() - getHashStringLength();
	vector<shared_ptr<UINT64> > localReadHits(hashQueryCount);
	#pragma omp critical(getRemoteData)
	{
		for(UINT64 j = 0; j < hashQueryCount; j++) // for each substring of read1 of length getHashStringLength
		{
			string subString = readString.substr(j,getHashStringLength()); // Get the substring from read string
			UINT64 index = getHashIndex(subString);	// Get the index using the hash function.
			//Compute various indices
			UINT64 globalStartOffset=hashTable[index];
			UINT64 globalEndOffset= (index==hashTableSize-1)?hashDataTableSize:hashTable[index+1];
			UINT64 hash_block_len = globalEndOffset-globalStartOffset;
			if(hash_block_len)
			{
				//check local cache and entry NOT found
				if(!getCachedHashTable(index, hash_block_len, localReadHits[j]))
				{
					int t_rank = getOffsetRank(globalStartOffset);
					MPI_Aint localOffset = getLocalOffset(globalStartOffset,t_rank);
					//Allocate memory
					void *dataBlock = new UINT64[hash_block_len];
					std::memset(dataBlock, 0, hash_block_len*sizeof(MPI_UINT64_T));
					// Get data from RMA
					MPI_Get(dataBlock, hash_block_len, MPI_UINT64_T, t_rank, localOffset, hash_block_len, MPI_UINT64_T, win);
					std::shared_ptr<UINT64> sp( static_cast<UINT64 *>(dataBlock), std::default_delete<UINT64[]>());
					localReadHits[j] = sp;
					insertCachedHashTable(index, sp);
					sp.reset();
					rmaCtr++;
				}
			}
		}
		MPI_Win_flush_all(win);
	}
	return localReadHits;
}
/**********************************************************************************************************************
	Get  a list of read containing the subString as prefix or suffix. Release cache on completion.
**********************************************************************************************************************/
map<UINT64,string> HashTable::getLocalHitList(vector<shared_ptr<UINT64> > localReadHits, string subString, UINT64 subStringIndx)
{
	map<UINT64,string> listOfReads;
	UINT64 *dataBlock=(localReadHits[subStringIndx]).get();
	if(dataBlock)
	{
		UINT64 index = getHashIndex(subString);	// Get the index using the hash function.
		//Compute various indices
		UINT64 globalStartOffset=hashTable[index];
		UINT64 globalEndOffset= (index==hashTableSize-1)?hashDataTableSize:hashTable[index+1];
		UINT64 hash_block_len = globalEndOffset-globalStartOffset;
		UINT64 startOffset=0;
		while(startOffset<hash_block_len)
		{
			UINT64 readID = dataBlock[startOffset] & 0X0000FFFFFFFFFFFF;
			UINT64 orient = dataBlock[startOffset] >> 63;
			UINT64 stringLen = (dataBlock[startOffset] >> 48) & 0X0000000000007FFF;
			UINT64 dataLen = (stringLen / 32) + (stringLen % 32 != 0);
			string forwardRead = toStringMPI(dataBlock, stringLen, startOffset+1);
			//cout<<"rid:"<<readID<<", ori:"<<orient<<", len:"<<stringLen<<", Read:"<<forwardRead<<endl;
			if(orient==0){
				string prefixForward = forwardRead.substr(0,hashStringLength);
				string suffixReverse = reverseComplement(prefixForward);
				if(subString==prefixForward)
				{
					UINT64 orientation=0;
					UINT64 data = readID | (orientation << 62);
					listOfReads.insert(pair<UINT64,string>(data,forwardRead));
				}
				else if(subString==suffixReverse){
					UINT64 orientation=3;
					UINT64 data = readID | (orientation << 62);
					listOfReads.insert(pair<UINT64,string>(data,forwardRead));
				}
			}
			else {
				string suffixForward = forwardRead.substr(forwardRead.length() - hashStringLength,hashStringLength);
				string prefixReverse = reverseComplement(suffixForward);
				if(subString==suffixForward)
				{
					UINT64 orientation=1;
					UINT64 data = readID | (orientation << 62);
					listOfReads.insert(pair<UINT64,string>(data,forwardRead));
				}
				else if(subString==prefixReverse){
					UINT64 orientation=2;
					UINT64 data = readID | (orientation << 62);
					listOfReads.insert(pair<UINT64,string>(data,forwardRead));
				}
			}
			startOffset+=dataLen+1;
		}
	}
	return listOfReads;
}

/**********************************************************************************************************************
	Get  a list of read containing the subString as prefix or suffix. Release cache on completion.
**********************************************************************************************************************/
map<UINT64,string> HashTable::getLocalHitList_nocache(vector<UINT64*> *localReadHits, string subString, UINT64 subStringIndx)
{
	map<UINT64,string> listOfReads;
	UINT64 *dataBlock=localReadHits->at(subStringIndx);
	if(dataBlock)
	{
		UINT64 index = getHashIndex(subString);	// Get the index using the hash function.
		//Compute various indices
		UINT64 globalStartOffset=hashTable[index];
		UINT64 globalEndOffset= (index==hashTableSize-1)?hashDataTableSize:hashTable[index+1];
		UINT64 hash_block_len = globalEndOffset-globalStartOffset;
		UINT64 startOffset=0;
		while(startOffset<hash_block_len)
		{
			UINT64 readID = dataBlock[startOffset] & 0X0000FFFFFFFFFFFF;
			UINT64 orient = dataBlock[startOffset] >> 63;
			UINT64 stringLen = (dataBlock[startOffset] >> 48) & 0X0000000000007FFF;
			UINT64 dataLen = (stringLen / 32) + (stringLen % 32 != 0);
			if((startOffset+dataLen)<hash_block_len)
			{
				string forwardRead = toStringMPI(dataBlock, stringLen, startOffset+1);
				//cout<<"rid:"<<readID<<", ori:"<<orient<<", len:"<<stringLen<<", Read:"<<forwardRead<<endl;
				if(orient==0){
					string prefixForward = forwardRead.substr(0,hashStringLength);
					string suffixReverse = reverseComplement(prefixForward);
					if(subString==prefixForward)
					{
						UINT64 orientation=0;
						UINT64 data = readID | (orientation << 62);
						listOfReads.insert(pair<UINT64,string>(data,forwardRead));
					}
					else if(subString==suffixReverse){
						UINT64 orientation=3;
						UINT64 data = readID | (orientation << 62);
						listOfReads.insert(pair<UINT64,string>(data,forwardRead));
					}
				}
				else {
					string suffixForward = forwardRead.substr(forwardRead.length() - hashStringLength,hashStringLength);
					string prefixReverse = reverseComplement(suffixForward);
					if(subString==suffixForward)
					{
						UINT64 orientation=1;
						UINT64 data = readID | (orientation << 62);
						listOfReads.insert(pair<UINT64,string>(data,forwardRead));
					}
					else if(subString==prefixReverse){
						UINT64 orientation=2;
						UINT64 data = readID | (orientation << 62);
						listOfReads.insert(pair<UINT64,string>(data,forwardRead));
					}
				}
				startOffset+=dataLen+1;
			}
			else
			{
				cout<<"Data corruption!!! RID:"<<readID<<" String Len"<<stringLen<<endl;
				cout<<"Hash Index:"<<index<<" Block Offset:"<<globalStartOffset<<" Block End Offset"<<globalEndOffset<<endl;
				exit(0);
			}
		}
	}
	return listOfReads;
}
/**********************************************************************************************************************
	Get the read counts for a given process...
**********************************************************************************************************************/
UINT64 HashTable::getMemoryReadCount(int myid)
{
	return memoryReadCount->at(myid);
}

/**********************************************************************************************************************
	Get the maximum read counts...
**********************************************************************************************************************/
UINT64 HashTable::getMaxMemoryReadCount()
{
	return *max_element(memoryReadCount->begin(),memoryReadCount->end());

}
/**********************************************************************************************************************
	Get the max offset for this rank
**********************************************************************************************************************/
UINT64 HashTable::getMemoryMaxLocalOffset(int rank)
{
	return (memoryDataPartitions->at(rank+1)-memoryDataPartitions->at(rank));
}
/**********************************************************************************************************************
	Fence and end a epoch...
**********************************************************************************************************************/
void HashTable::endEpoch()
{
	MPI_Win_fence(0, win);
}
/**********************************************************************************************************************
	Destructor
**********************************************************************************************************************/
HashTable::~HashTable(void)
{
	// Free the memory used by the hash table.
	delete[] hashTable;
	delete memoryDataPartitions;
	delete memoryReadCount;
	MPI_Win_free(&win);
	MPI_Free_mem(hashData);

	//delete cache table
	clearHashTableCache();
	delete cachedHashTable;
	delete hashQ;
	delete cachedReadTable;
	delete readQ;
}

/*
 * Convert a set of data bytes from hash data table to string
 *
 */
string HashTable::toString(UINT64 hashDataIndex,UINT64 stringLen) const
{
	string dna_str="";
	string strArr[] = {"A","C","G","T"};
	/* for each base of the DNA sequence */
	for (size_t i = 0; i < stringLen; i++)
	{
		uint8_t shift = (32*2-2) - 2*(i % 32);
		/* get the i-th DNA base */
		UINT64 base = (hashData[hashDataIndex + i/32] & ((UINT64)BASE_MASK << shift)) >> shift;
		dna_str += strArr[base];
	}
	return dna_str;
}

UINT64 HashTable::getLocalReadID(UINT64 localOffset, int myid) const
{
	UINT64 readID;
	if(localOffset<(memoryDataPartitions->at(myid+1)-memoryDataPartitions->at(myid)))
		readID = hashData[localOffset] & 0X0000FFFFFFFFFFFF;
	else
		MYEXIT("Read ID out of range: ID"+SSTR(localOffset)+", Rank:"+SSTR(myid));
	return readID;
}
UINT64 HashTable::getLocalReadOrient(UINT64 localOffset, int myid) const
{
	UINT64 orient;
	if(localOffset<(memoryDataPartitions->at(myid+1)-memoryDataPartitions->at(myid)))
		orient = hashData[localOffset] >> 63;
	else
		MYEXIT("Read ID out of range: ID"+SSTR(localOffset)+", Rank:"+SSTR(myid));
	return orient;
}
UINT64 HashTable::getLocalReadLength(UINT64 localOffset, int myid) const
{
	UINT64 readLen;
	if(localOffset<(memoryDataPartitions->at(myid+1)-memoryDataPartitions->at(myid)))
		readLen = ((hashData[localOffset] >> 48) & 0X0000000000007FFF);
	else
		MYEXIT("Read ID out of range: ID"+SSTR(localOffset)+", Rank:"+SSTR(myid));
	return readLen;
}
string HashTable::getLocalStringForward(UINT64 localOffset, int myid) const
{
	string readStr;
	if(localOffset<(memoryDataPartitions->at(myid+1)-memoryDataPartitions->at(myid)))
	{
		UINT64 stringLen=getLocalReadLength(localOffset,myid);
		readStr = toString(localOffset+1,stringLen);
	}
	else
		MYEXIT("Read ID out of range: ID"+SSTR(localOffset)+", Rank:"+SSTR(myid));
	return readStr;
}
UINT64 HashTable::getLocalNextOffset(UINT64 localOffset, int myid) const
{
	UINT64 nextOffset=0;
	if(localOffset<(memoryDataPartitions->at(myid+1)-memoryDataPartitions->at(myid)))
	{
		UINT64 stringLen=getLocalReadLength(localOffset,myid);
		UINT64 dnaWordLen = (stringLen / 32) + (stringLen % 32 != 0);
		nextOffset = localOffset+dnaWordLen+1;
	}
	else
		MYEXIT("Read ID out of range: ID"+SSTR(localOffset)+", Rank:"+SSTR(myid));
	return nextOffset;
}

/*
 * Convert a set of data bytes from a data block to string
 *
 */
string HashTable::toStringMPI(UINT64 *hashDataBlock,UINT64 stringLen, UINT64 startOffset) const
{
	string dna_str="";
	string strArr[] = {"A","C","G","T"};
	/* for each base of the DNA sequence */
	for (size_t i = 0; i < stringLen; i++)
	{
		uint8_t shift = (32*2-2) - 2*(i % 32);
		/* get the i-th DNA base */
		UINT64 base = (hashDataBlock[startOffset + i/32] & ((UINT64)BASE_MASK << shift)) >> shift;
		dna_str += strArr[base];
	}
	return dna_str;
}
UINT64 HashTable::getReadLength(UINT64 globalOffset, UINT64 readID , int myid) const
{
	UINT64 dataRec;
	dataRec=0;
	int t_rank = getOffsetRank(globalOffset);
	MPI_Aint localOffset = getLocalOffset(globalOffset,t_rank);
	MPI_Get(&dataRec, 1, MPI_UINT64_T, t_rank, localOffset, 1, MPI_UINT64_T, win);
	MPI_Win_flush(t_rank,win);
	rmaCtr++;
	UINT64 rmaReadID = dataRec & 0X0000FFFFFFFFFFFF;
	if(readID!=rmaReadID)
	{
		cout<<"Proc:"<<myid<<" Error!!! Read ID do not match:"<<readID<<" "<<rmaReadID<<" Goff:"<<globalOffset<<" lOff:"<<localOffset<<" Len:"<<((dataRec >> 48) & 0X0000000000007FFF)<<endl;
		exit(0);
	}
	return ((dataRec >> 48) & 0X0000000000007FFF); 	//2nd MSB to 16th MSB are read length
	return 0;
}

string HashTable::getStringForward(Read *r, int myid)
{
	UINT64 globalOffset = r->getReadHashOffset();
	UINT64 rid = r->getReadNumber();
	UINT64 fid = r->getFileIndex();
	UINT64 *dataBlock = NULL;
	string seq="";
	#pragma omp critical(getRemoteData)
	{
		//check local cache
		if(!getCachedRead(rid,seq))
		{

			UINT64 stringLen=getReadLength(globalOffset, rid, myid);
			int t_rank = getOffsetRank(globalOffset);
			MPI_Aint localOffset = getLocalOffset(globalOffset,t_rank)+1;
			UINT64 dna_word_len = (stringLen / 32) + (stringLen % 32 != 0);
			dataBlock = new UINT64[dna_word_len];
			std::memset(dataBlock, 0, dna_word_len*sizeof(MPI_UINT64_T));
			try
			{
				//cout<<"Proc"<<myid<<" ReadID:"<<fid<<", trank:"<<t_rank<<", Goffset:"<<globalOffset<<", LOffset:"<<localOffset<<", Len:"<<stringLen<<endl;
				MPI_Get(dataBlock, dna_word_len, MPI_UINT64_T, t_rank, localOffset, dna_word_len, MPI_UINT64_T, win);
				MPI_Win_flush(t_rank,win);
				rmaCtr++;
				seq = toStringMPI(dataBlock,stringLen,0);
				insertCachedRead(rid, seq);
			}
			catch (const std::bad_alloc&) {
				cout<<"Proc"<<myid<<" Exception Caught getting read data!!!";
				//cout<<"ReadID:"<<rid<<", trank:"<<t_rank<<", Goffset:"<<globalOffset<<", LOffset:"<<localOffset<<", Len:"<<stringLen<<endl;
				exit(0);
			}
		}
	}
	if(dataBlock)
		delete[] dataBlock;
	return seq;
}

string HashTable::getStringReverse(Read *r, int myid)
{
	return reverseComplement(getStringForward(r,myid));
}
/*
 * Check if the given read belongs to this MPI node and and returns true if contained read has to be processed by this
 * MPI node
  */
bool HashTable::needsProcessing(UINT64 read1ID, string readString, int myid)
{
	string prefixForward = readString.substr(0,hashStringLength); 											// Prefix of the forward string.
	string suffixForward = readString.substr(readString.length() - hashStringLength,hashStringLength);	// Suffix of the forward string.

	UINT64 index1 = getHashIndex(prefixForward);						// Get the index using the hash function.
	UINT64 baseOffset1 = hashTable[index1];							// Get the start offset of the hash value in the hashData table
	int rank1 = getOffsetRank(baseOffset1);

	UINT64 index2 = getHashIndex(suffixForward);						// Get the index using the hash function.
	UINT64 baseOffset2 = hashTable[index2];							// Get the start offset of the hash value in the hashData table
	int rank2 = getOffsetRank(baseOffset2);
	int computeRank=0;

	if(read1ID%2==0)
		computeRank = rank1<rank2?rank1:rank2;
	else
	{
		computeRank = rank1<rank2?rank2:rank1;
	}
	return (myid==computeRank?true:false);
}

bool HashTable::getCachedHashTable(UINT64 index, UINT64 hash_block_len, shared_ptr<UINT64> &cacheBlock)
{
	map<UINT64, shared_ptr<UINT64> >::iterator it;
	it = cachedHashTable->find(index);
	if (it != cachedHashTable->end())
	{
		cacheBlock = it->second;
		return true;
	}
	else
		return false;
}
void HashTable::insertCachedHashTable(UINT64 index, shared_ptr<UINT64> &cacheBlock)
{
	if(cachedHashTable->size()<HASH_CACHE_SIZE)
	{
		std::shared_ptr<UINT64> loc(cacheBlock);
		cachedHashTable->insert(pair<UINT64, shared_ptr<UINT64>>(index,loc));
		hashQ->push(index);
	}
	else if(cachedHashTable->size()>0)
	{
		UINT64 oindex = hashQ->front();
		hashQ->pop();
		map<UINT64, shared_ptr<UINT64> >::iterator it = cachedHashTable->find(oindex);
		it->second.reset();
		cachedHashTable->erase(it);
		hashQ->push(index);
		std::shared_ptr<UINT64> loc(cacheBlock);
		cachedHashTable->insert(pair<UINT64, shared_ptr<UINT64>>(index,loc));
	}
}

void HashTable::clearHashTableCache()
{
	//delete cache entries
	for(map<UINT64,shared_ptr<UINT64> >::iterator it=cachedHashTable->begin() ; it!=cachedHashTable->end(); ++it) // For each read in the list.
	{
		hashQ->pop();
		it->second.reset();
	}
	cachedHashTable->clear();
}

bool HashTable::getCachedRead(UINT64 rid, string &readString)
{
	map<UINT64, string>::iterator it;
	it = cachedReadTable->find(rid);
	if (it != cachedReadTable->end())
	{
		readString=it->second;
		return true;
	}
	else
		return false;
}
void HashTable::insertCachedRead(UINT64 rid, string readString)
{
	if(cachedReadTable->size()<READ_CACHE_SIZE)
	{
		readQ->push(rid);
		cachedReadTable->insert(pair<UINT64, string>(rid,readString));
	}
	else if(cachedReadTable->size()>0)
	{
		UINT64 oRid = readQ->front();
		readQ->pop();
		map<UINT64, string>::iterator it = cachedReadTable->find(oRid);
		cachedReadTable->erase(it);
		readQ->push(rid);
		cachedReadTable->insert(pair<UINT64, string>(rid,readString));
	}
}

void HashTable::deleteLocalHitList(vector<UINT64*> *localReadHits)
{
	for(size_t i=0;i<localReadHits->size();i++)
	{
		if(localReadHits->at(i))
			delete[] localReadHits->at(i);
	}
	delete localReadHits;
	localReadHits=NULL;
}



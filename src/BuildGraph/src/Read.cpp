/*
 * Read.cpp
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 */

#include "Read.h"

#include "Common.h"



/**********************************************************************************************************************
	Default constructor
**********************************************************************************************************************/
Read::Read(void)
{
	// Initialize the variables.
	readNumber = 0;
	superReadID = 0;
	fileIndex=0;
}

/**********************************************************************************************************************
	Another constructor with file index
**********************************************************************************************************************/
Read::Read(UINT64 fIndx,int len)
{
	readNumber = 0;
	superReadID = 0;
	fileIndex=fIndx;
	readCoverage.insert(readCoverage.end(), len, 0);
}

/**********************************************************************************************************************
	Default destructor
**********************************************************************************************************************/
Read::~Read(void)
{
	// delete all the pointers.
}
/**********************************************************************************************************************
	This function assigns an ID to the read.
**********************************************************************************************************************/
bool Read::setReadNumber(UINT64 id)
{
	if(id <= 0) MYEXIT("ID less than 1.");
	readNumber = id;												// Set the read number.
	return true;
}
/**********************************************************************************************************************
	This function assigns an ID to the read.
**********************************************************************************************************************/
void Read::setFileIndex(UINT64 id)
{
	if(id <= 0) MYEXIT("ID less than 1.");
	fileIndex = id;												// Set the read number.
}
/**********************************************************************************************************************
	Returns the reverse complement of a read.
**********************************************************************************************************************/

void Read::setReadHashOffset(UINT64 offset)
{
	readHashOffset=offset;
}
/**********************************************************************************************************************
	This function assigns the super read ID to the read.
**********************************************************************************************************************/
void Read::setSuperReadID(UINT64 id)
{
	superReadID=id;
}

void Read::incrementReadCoverage(int indx)
{
	readCoverage.at(indx)=readCoverage.at(indx)+1;
}

UINT64 Read::getReadCoverage(int indx)
{
	return readCoverage.at(indx);
}



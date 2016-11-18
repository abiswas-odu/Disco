#ifndef __dna_h__
#define __dna_h__


#include <stdexcept>	/* std::invalid_argument */
#include <inttypes.h>
#include "Config.h"


#define BASE_MASK 0x0000000000000003	/* binary: 11 */

#define FRAG_SIZE 64
#define HALF_FRAG_SIZE 32	/* Half of the number of bits used as the bit storage block. Here it is 32 so it is set to 16.*/

/* useful constants */
enum
{
	BASE_A = 0x0,	/* binary: 00 */
	BASE_C = 0x1,	/*'binary: 01 */
	BASE_G = 0x2,	/* binary: 10 */
	BASE_T = 0x3,	/* binary: 11 */
};

class dna_bitset
{
private:
	uint64_t* m_data;
	size_t   m_len;
public:

	/**
		 * @brief returns 1 is the substrings are equal
		 * @param this pointer to first sequence to be compared (subject)
		 * @param start index in subject to compare
		 * @param length of substring of subject to compare
		 * @param pointer to second sequence to be compared (query)
		 * @param start index in query to compare
		 * @param length of substring of query to compare
	*/
	bool compareSubString (UINT64 seq1Start, UINT64 seq1Len, const dna_bitset *seq2, UINT64 seq2Start, UINT64 seq2Len, UINT64 orient) const;
	/**
		 * @brief default constructor
	 */
	dna_bitset ()
	{
		m_len=0;
		m_data=NULL;
	}
	/**
	 * @brief constructor for initializing only length
	 */
	dna_bitset (UINT64 seqLen)
	{
		m_len=seqLen;
		m_data=NULL;
	}
	/**
	 * @brief constructor
	 * @param dna_str a string containing a DNA sequence (e.g. "ATGCA...")
	 * @param dna_len length of the DNA sequence
	 */
	dna_bitset (const string &dna_str)
	{
		m_len = dna_str.length();

		/* number of bytes necessary to store dna_str as a bitset */
		size_t dna_word = (m_len / HALF_FRAG_SIZE) + (m_len % HALF_FRAG_SIZE != 0);

		m_data = new uint64_t[dna_word];

		std::memset(m_data, 0, dna_word*(FRAG_SIZE/8));

		/* for each base of the DNA sequence */
		for (size_t i = 0; i < m_len; i++)
		{
			uint64_t shift = (HALF_FRAG_SIZE*2-2) - 2*(i % HALF_FRAG_SIZE);

			switch (dna_str[i])
			{
			case 'A':
				m_data[i/HALF_FRAG_SIZE] |= (UINT64)BASE_A << shift;
				break;
			case 'C':
				m_data[i/HALF_FRAG_SIZE] |= (UINT64)BASE_C << shift;
				break;
			case 'G':
				m_data[i/HALF_FRAG_SIZE] |= (UINT64)BASE_G << shift;
				break;
			case 'T':
				m_data[i/HALF_FRAG_SIZE] |= (UINT64)BASE_T << shift;
				break;
			default:
				throw std::invalid_argument("invalid DNA base");
			}
		}
	}
	/**
	 * @brief destructor
	 */
	~dna_bitset ()
	{
		delete[] m_data;
	}
	/**
		 * @brief returns reverse complement of the stored DNA sequence as a C++ string
		 */
	string toRevComplement () const
	{
		string dna_str="";
		string strRevCArr[] = {"T","G","C","A"};
		/* for each base of the DNA sequence */
		for (int i = m_len-1; i >= 0; i--)
		{
			uint64_t shift = (HALF_FRAG_SIZE*2-2) - 2*(i % HALF_FRAG_SIZE);
			/* get the i-th DNA base */
			uint64_t base = (m_data[i/HALF_FRAG_SIZE] & ((UINT64)BASE_MASK << shift)) >> shift;
			dna_str += strRevCArr[base];
		}
		return dna_str;
	}
	/**
		* @brief returns the stored DNA sequence as a C++ string
	 */
	string toString () const
	{
		string dna_str="";
		string strArr[] = {"A","C","G","T"};
		/* for each base of the DNA sequence */
		for (size_t i = 0; i < m_len; i++)
		{
			uint64_t shift = (HALF_FRAG_SIZE*2-2) - 2*(i % HALF_FRAG_SIZE);
			/* get the i-th DNA base */
			uint64_t base = (m_data[i/HALF_FRAG_SIZE] & ((UINT64)BASE_MASK << shift)) >> shift;
			dna_str += strArr[base];
		}
		return dna_str;
	}
	int getLength(void) const {return m_len;}

};

#endif /* __dna_h__ */



#include "dna.h"

/**
	 * @brief returns 1 is the substrings are equal
	 * @param pointer to first sequence to be compared (subject)
	 * @param start index in subject to compare
	 * @param length of substring of subject to compare
	 * @param pointer to second sequence to be compared (query)
	 * @param start index in query to compare
	 * @param length of substring of query to compare
	 * @param overlap orientation
	 * // orient 0
	//   >--------MMMMMMMMMMMMMMM*******------> read1      M means match found by hash table
	//            MMMMMMMMMMMMMMM*******>       read2      * means we need to check these characters for match
	//				OR
	// orient 2
	//	 >---*****MMMMMMMMMMMMMMM*******------> read1
	//		      MMMMMMMMMMMMMMM*******<	    Reverese complement of read2
	 *
	 * // orient 1
	//   >---*****MMMMMMMMMMMMMMM-------------> read1      M means match found by hash table
	//      >*****MMMMMMMMMMMMMMM       		read2      * means we need to check these characters for match
	//				OR
	// orient 3
	//	 >---*****MMMMMMMMMMMMMMM-------------> read1
	//		<*****MMMMMMMMMMMMMMM				Reverse Complement of Read2
*/
bool dna_bitset::compareSubString (UINT64 seq1Start, UINT64 seq1Len, const dna_bitset *seq2, UINT64 seq2Start, UINT64 seq2Len, UINT64 orient) const
{
	if(seq1Len!=seq2Len)
		return false;
	if(orient == 0 || orient== 1)
	{
		for (UINT64 i = seq1Start, j = seq2Start; i < seq1Start+seq1Len; i++,j++)
		{
			uint8_t shift1 = (HALF_FRAG_SIZE*2-2) - 2*(i % HALF_FRAG_SIZE);
			///* get the i-th DNA base from subject
			uint8_t base1 = (m_data[i/HALF_FRAG_SIZE] & (BASE_MASK << shift1)) >> shift1;

			uint8_t shift2 = (HALF_FRAG_SIZE*2-2) - 2*(j % HALF_FRAG_SIZE);
			///* get the j-th DNA base from query
			uint8_t base2 = (seq2->m_data[j/HALF_FRAG_SIZE] & (BASE_MASK << shift2)) >> shift2;

			if(base1!=base2)
				return false;
		}
	}
	else
	{
		for (UINT64 i = seq1Start, j = seq2->m_len-seq2Start-1; i < seq1Start+seq1Len; i++,j--)
		{
			uint8_t shift1 = (HALF_FRAG_SIZE*2-2) - 2*(i % HALF_FRAG_SIZE);
			///* get the i-th DNA base from subject
			uint8_t base1 = (m_data[i/HALF_FRAG_SIZE] & (BASE_MASK << shift1)) >> shift1;

			uint8_t shift2 = (HALF_FRAG_SIZE*2-2) - 2*(j % HALF_FRAG_SIZE);
			///* get the j-th DNA base from query
			uint8_t base2 = ~((seq2->m_data[j/HALF_FRAG_SIZE] & (BASE_MASK << shift2)) >> shift2) & BASE_MASK;

			if(base1!=base2)
				return false;
		}
	}
	return true;
}

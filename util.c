#include <util.h>
#include <pairwise_align.h>
#include <string.h>
#include <unistd.h>

/*
  Fair note, this routine will cause memory utilities like Valgrind to whine viciously about
  reads in uninitalized memory. This is because it is correct, this utility will just read
  bytes in page size increments to fault in all pages and get real backing memory for them.
  This causes no real problem, since the uninitalized bytes read don't do anything.
*/
void touch_memory(void *mem, size_t size) {
  index_t page_size = sysconf(_SC_PAGESIZE);
  index_t *this_memory = (index_t *)mem;
  index_t size_increment = page_size / sizeof(index_t);
  index_t size_max = size / sizeof(index_t);
  for(index_t idx=0; idx < size_max; idx+=size_increment) {
    this_memory[idx] = 0;
  }
}

/* Fix up a sequence to remove the gaps
   Input-
        good_match_t *A    - the match struct for the sequence.  Used for the hyphen member
        int *dest          - destination array.  Needs to be allocated before it is passed in
        int *source        - sequence to be scrubbed
        int length         - length of the source.  Dest should be the same size
   Output-
        int *dest          - will be full of a gapless sequence
        int (return value) - size of dest
*/

index_t scrub_hyphens(good_match_t *A, codon_t *dest, codon_t *source, index_t length)
{
  index_t source_index=0, dest_index=0;
  while(source_index < length)
  {
    while(source_index < length && source[source_index] == A->simMatrix->hyphen) source_index++;
    dest[dest_index] = source[source_index];
    source_index++; dest_index++;
  }
  return dest_index-1;
}

/* Helper function for displaying acid chains.  This function will take a chain and make an
   ascii representation of the acids as defined in A
   Input-
       good_match_t *A    - the match struct for the sequence.
       char *result       - the resulting chain as a sequnce of acids
       int *chain         - the sequence to convert into ascii
       int length         - size of the chain
*/

void assemble_acid_chain(good_match_t *A, char *result, codon_t *chain, index_t length)
{
  memset(result, '\0', length);
  for(int idx=0; idx < length; idx++)
  {
    result[idx] = (char)A->simMatrix->aminoAcid[chain[idx]];
  }
  result[length] = '\0';  
}

/* Helper function for displaying codon chains.  This function will take a chain and make an
   ascii representation of the codon sequence out of it
   Input-
       good_match_t *A    - the match struct for the sequence.
       char *result       - the resulting chain as a sequnce of codons
       int *chain         - the sequence to convert into ascii
       int length         - size of the chain
*/

void assemble_codon_chain(good_match_t *A, char *result, codon_t *chain, index_t length)
{
  memset(result, '\0', length);
  for(int idx=0; idx < length; idx++)
  {
    if(chain[idx] == HYPHEN)
    {
      result[idx*3] = '-';
      result[idx*3+1] = '-';
      result[idx*3+2] = '-';
    }
    else
    {
      result[idx*3] = (char)A->simMatrix->codon[chain[idx]][0];
      result[idx*3+1] = (char)A->simMatrix->codon[chain[idx]][1];
      result[idx*3+2] = (char)A->simMatrix->codon[chain[idx]][2];
    }
  }
  result[length*3] = '\0';
}

/* Very simple scoring routine.  Compares main to match.  Both must be
   at least length long.
   Input-
       goot_match_t *A - the match structure used to generate main and match
       int main[]      - the main sequences
       int match[]     - the sequence to compare to main
       int length      - the size of smallest of main or match
  Output-
       int             - the score of the match
*/

score_t simple_score(good_match_t *A, codon_t main[], codon_t match[], index_t length)
{
  score_t score = 0;
  int mainMatch = 1;
  int matchMatch = 1;
  
  // recompute score by a brain dead simple method
  for(int i=0; i < length; i++)
  {
    if(main[i] == A->simMatrix->hyphen)
    {
      if(mainMatch == 1)
      {
        mainMatch = 0;
        score = score - A->simMatrix->gapStart;
      }
      score = score - A->simMatrix->gapExtend;
      continue;
    }
    if(match[i] == A->simMatrix->hyphen)
    {
      if(matchMatch == 1)
      {
        matchMatch = 0;
        score = score - A->simMatrix->gapStart;
      }
      score = score - A->simMatrix->gapExtend;
      continue;
    }
    mainMatch = 1;
    matchMatch = 1;
    score = score + A->simMatrix->similarity[main[i]][match[i]];
  }
  return score;
}

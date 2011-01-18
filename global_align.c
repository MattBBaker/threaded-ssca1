#include <util.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <global_align.h>

int inverse_identity_matrix[4][4] = {{0,1,1,1},{1,0,1,1},{1,1,0,1},{1,1,1,0}};

void print_square(short **square, int length)
{
  for(int idx=0; idx < length; idx++)
  {
    for(int jdx=0; jdx < length; jdx++)
    {
      printf("%2i ", square[idx][jdx]);
    }
    printf("\n");
  }
}

int verifyGlobal(ga_t G[], int size_g, int misPenalty, int gapPenalty, int maxDisplay)
{
  int retval=0;
  int D = size_g < maxDisplay ? size_g : maxDisplay;

  printf("\nDisplaying %i of %i sets of global alignment pair scores.\n\n", D, size_g);
  printf("       Penalty for mismatching bases: %i\n", misPenalty);
  printf("     Penalty for each space in a gap: %i\n", gapPenalty);

  if(D>0)
    printf("\nThe score for sequence x matched against sequence y.\n\n");

  for(int b=0; b < D; b++)
  {
    if(G[b].length > 0)
    {
      printf("Set %i:\n", b);
      print_square(G[b].globalScores, G[b].length);
    }
  }

  return retval;
}

/*
 * Function scorePair() - Compute the global alignment score for two
 * sequences using a variant of the Smith-Waterman algorithm.
 * Input-
 *     good_match_t *A        - [structure pointer] structure containing the results of kernel 1 and kernel 2
 *     seq_t *first_seq       - [structure pointer] structure containing the first sequence to be aligned, X
 *     seq_t *second_seq      - [structure pointer] structure containing the second sequence to be aligned, Y
 *     int codon2bases[3][64] - codon translation data
 *     int gapPenalty         - penelty for gaps, as defined in the parameters
 *     int disMatrix[4][4]    - used to score codons
 * Output-
 *     short                  - global alignment score
 */

short scorePair(good_match_t *A, seq_t *first_seq, seq_t *second_seq, int codon2bases[3][64], int gapPenalty, int disMatrix[4][4])
{
  int X[first_seq->length];
  int Xb[(first_seq->length) * 3];
  int length_x = scrub_hyphens(A, X, first_seq->main, first_seq->length );
  int n = length_x * 3;

  int Y[second_seq->length];
  int Yb[(second_seq->length) * 3];
  int length_y = scrub_hyphens(A, Y, second_seq->main, second_seq->length );
  int m = length_y * 3;

  int adder;

  short V[m];
  short G[m];

  for(int idx=0; idx < m; idx++)
  {
    V[idx] = (idx+1) * gapPenalty;
  }

  for(int idx=0; idx < length_x; idx++)
  {
    Xb[idx*3] = codon2bases[0][X[idx]];
    Xb[idx*3+1] = codon2bases[1][X[idx]];
    Xb[idx*3+2] = codon2bases[2][X[idx]];
  }

  for(int idx=0; idx < length_y; idx++)
  {
    Yb[idx*3] = codon2bases[0][Y[idx]];
    Yb[idx*3+1] = codon2bases[1][Y[idx]];
    Yb[idx*3+2] = codon2bases[2][Y[idx]];
  }

  for(int idx=0; idx < n; idx++)
  {
    for(int jdx=0; jdx < m; jdx++)
    {
      adder = (jdx == 0) ? idx * gapPenalty : V[jdx-1];
      G[jdx] = disMatrix[Xb[idx]][Yb[jdx]] + adder;
    }

    for(int jdx=0; jdx < m; jdx++)
    {
      V[jdx] = (V[jdx] + gapPenalty < G[jdx]) ? V[jdx] + gapPenalty : G[jdx];
    }

    for(int jdx=1; jdx < m; jdx++)
    {
      V[jdx] = (V[jdx-1] + gapPenalty < V[jdx]) ? V[jdx-1] + gapPenalty : V[jdx];
    }
  }

  return V[m-1];
}

/* Frees the data produced by this kernel
 * Input-
 *    ga_t *doomed - [ga_t[size]] Array of ga_structures
 *    int size     - [int] size of the ga_t array
 */

void release_ga(ga_t *doomed, int size)
{
  for(int idx=0; idx < size; idx++)
  {
    for(int jdx=0; jdx < doomed[idx].length; jdx++)
    {
      free(doomed[idx].globalScores[jdx]);
    }
    free(doomed[idx].globalScores);
  }
  free(doomed);
}

/*
 * This function uses a variant of the Smith-Waterman dynamic programming
 * algorithm to match each pair of subsequences in each of the sets of similar
 * subsequences reported by Kernel 3.  It converts the codon sequences to base
 * sequences and does global alignment, i.e., matches the two entire
 * subsequences.  It minimizes a distance function instead of maximizing a
 * similarity function as is done in Kernels 1-3.  There is no penalty for
 * each gap, the gapPenalty is for each space in the gap.
 *
 * This function returns a symmetric matrix of scores for each set of
 * similar sequences, where entry (x,y) holds the alignment score of
 * subsequence x with subsequence y.
 *
 * For a detailed description of the SSCA #1 Optimal Pattern Matching problem,
 * please see the SSCA #1 Written Specification.
 *
 * INPUT
 * A                - [good_match_t] results from pairwiseAlign
 *   seqData        - [seq_data_t] data sequences created by genScalData()
 *     main         - [char pointer] first codon sequence
 *     match        - [char pointer] second codon sequence
 *     maxValidation- [int] longest matching validation string.
 *   simMatrix      - [sim_matrix_t] codon similarity created by genSimMatrix()
 *     similarity   - [int [64][64]] 1-based codon/codon similarity table
 *     aminoAcid    - [char[65]] 1-based codon to aminoAcid table
 *     bases        - [char [4]] 1-based encoding to base letter table
 *     codon        - [char [64][3]] 1-based codon to base letters table
 *     encode       - [int [128]] aminoAcid character to last codon number
 *     hyphen       - [char] encoding representing a hyphen (gap or space)
 *     exact        - [int] value for exactly matching codons
 *     similar      - [int] value for similar codons (same amino acid)
 *     dissimilar   - [int] value for all other codons
 *     gapStart     - [int] penalty to start gap (>=0)
 *     gapExtend    - [int] penalty to for each codon in the gap (>0)
 *     matchLimit   - [int] longest match including hyphens
 *   goodEnds       - [int[M][2]] M matches; main/match endpoints
 *   goodScores     - [int[M]] the scores for the goodEnds
 * size_s           - [int] size of the S array.  Will also be the size of the returned array
 * S                - [good_match_t[size_s]] similar sequence sets for each bestSeq
 *   bestStarts     - [int[2][M]] main/match startpoints
 *   bestEnds       - [int[2][M]] main/match endpoints
 *   bestSeqs       - [seq_t[M]] main/match sequences
 *   bestScores     - [int[M]] the scores for the bestSeqs
 * misPenalty       - [int] penalty for each base/base mismatch
 * gapPenalty       - [int] penalty for each space in a gap
 *
 * OUTPUT
 * GA               - [ga_t[size_s]] array of alignments for set g
 *   length         - [int] size of the globalScores matrix
 *   globalScores   - [short[length][length]] alignment score (x,y)
 *
 */

ga_t *globalAlign(good_match_t *A, int size_s, good_match_t *S[size_s], int misPenalty, int gapPenalty)
{
  int dis_matrix[4][4];
  int codon2bases[3][64];
  int M;
  ga_t *GA = malloc(sizeof(ga_t) * size_s);
  memset(GA, 0, sizeof(ga_t *) * size_s);
  for(int idx=0; idx<4; idx++)
    for(int jdx=0; jdx<4; jdx++) 
      dis_matrix[idx][jdx] = inverse_identity_matrix[idx][jdx] * misPenalty;

  for(int idx=0; idx < 64; idx++)
  {
    codon2bases[0][idx] = idx / 16;
    codon2bases[1][idx] = (idx / 4) % 4;
    codon2bases[2][idx] = idx % 4;
  }

#ifdef _OPENMP
#pragma omp parallel for private(M)
#endif
  for(int idx=0; idx < size_s; idx++)
  {
    M = S[idx]->bestLength;
    if(0==M)
    {
      GA[idx].globalScores = NULL;
      continue;
    }

    GA[idx].length=M;

    GA[idx].globalScores = malloc(sizeof(short *) * M);
    for(int jdx=0; jdx < M; jdx++)
    {
      GA[idx].globalScores[jdx] = malloc(sizeof(short) * M);
      memset(GA[idx].globalScores[jdx], 0, M * sizeof(short));
    }
    for(int x=0; x < M; x++)
    {
      for(int y=x+1; y < M; y++)
      {
      	GA[idx].globalScores[x][y] = scorePair(A, &(S[idx]->bestSeqs[x]), &(S[idx]->bestSeqs[y]), codon2bases, gapPenalty, dis_matrix);
	      GA[idx].globalScores[y][x] = GA[idx].globalScores[x][y];
      }
    }
  }

  return GA;
}

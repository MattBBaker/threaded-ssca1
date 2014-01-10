#include <stdio.h>
#include <util.h>
#include <stdlib.h>
#include <string.h>
#include <gen_scal_data.h>
#include <pairwise_align.h>
#include <global_align.h>
#include <multiple_align.h>
#include <sort.h>
#include <gen_sim_matrix.h>

#define REALLOC_INCREMENT 1024

int verifyMultiple(ma_t MA[], index_t size_ma, int maxDisplay)
{
  int retval=0, t=0;
  int D = size_ma < maxDisplay ? size_ma : maxDisplay;

  printf("\nDisplaying %i of %lu sets of global alignments.\n", D, size_ma);

  if(D > 0)
  {
    printf("\nWithin each set the aligned center-sequence is shown first,");
    printf("\nthen the other aligned sequences ordered by score.\n");
  }

  for(int b=0; b < D; b++)
  {
    if(MA[b].length == 0)
    {
      printf("\nSet %i is empty.\n", b);
    }
    else
    {
      t=0;
      for(int idx=0; idx < MA[b].length; idx++)
      {
        t += MA[b].scores[idx];
      }
      printf("\nSet %i; size %lu, total score %i, average score %.2f\n",
             b, MA[b].length, t, (float)t/(float)MA[b].length);
      for(int i=0; i < MA[b].length; i++)
      {
        char print_buffer[MA[b].alignment[i].length+1];
        for(int idx=0; idx < MA[b].alignment[i].length; idx++)
        {
          if(MA[b].alignment[i].main[idx] == HYPHEN)
            print_buffer[idx] = '-';
          else
            print_buffer[idx] = MA[b].alignment[i].main[idx];
        }
        print_buffer[MA[b].alignment[i].length] = '\0';
        printf("  %4i  %s\n", (int)MA[b].scores[i], print_buffer);
      }
    }
  }
    
  return retval;
}

void print_seq(seq_t *sequence)
{
  char *main_buffer=malloc(sequence->length+1);
  char *match_buffer=malloc(sequence->length+1);
  main_buffer[0] = '\0';
  match_buffer[0] = '\0';

  if(sequence->main != NULL)
  {
    for(int idx=0; idx < sequence->length; idx++)
    {
      if(sequence->main[idx] == HYPHEN)
        main_buffer[idx] = '-';
      else
        main_buffer[idx] = sequence->main[idx];
    }
    main_buffer[sequence->length] = '\0';
  }

  if(sequence->match != NULL)
  {
    for(int idx=0; idx < sequence->length; idx++)
    {
      if(sequence->match[idx] == HYPHEN)
        match_buffer[idx] = '-';
      else
        match_buffer[idx] = sequence->match[idx];
    }
    main_buffer[sequence->length] = '\0';
  }


  printf(" main: %s;\nmatch: %s;\n", main_buffer, match_buffer);
}

/*
 * Generate a typical matching pair.
 *
 * With a small change this recursive approach can be used to generate all
 * possible matches.  We only need to identify a single match.  C implementation
 * is limited in recursion only by memory.
 * slightly different from the one found in scan_backwards.c and locate_similar.
 * Maybe consider trying to merge them?
 */

void tracepath_x(good_match_t *A, int n, int m, score_t T[n][m], codon_t Yb[m], codon_t Xb[n], int ci, int cj, int dir, seq_t *sequences, int rs[2], int depth)
{
  int Ci, Cj, size_data;

  if(ci==-1)
  {
    size_data = depth + cj + 1;
    sequences->main = (codon_t*)malloc(sizeof(codon_t)*(size_data));
    sequences->match = (codon_t*)malloc(sizeof(codon_t)*(size_data));
    for(int idx=0; idx <= cj; idx++)
    {
      sequences->main[idx] = A->simMatrix->hyphen;
      sequences->match[idx] = A->simMatrix->bases[Yb[idx]];
    }
    rs[0] = 0;
    rs[1] = 0;
    sequences->length = size_data;
    sequences->backing_memory = size_data;
    return;
  }

  if(cj==-1)
  {
    size_data = depth + ci + 1;
    sequences->main = (codon_t*)malloc(sizeof(codon_t)*(size_data));
    sequences->match = (codon_t*)malloc(sizeof(codon_t)*(size_data));
    for(int idx=0; idx <= ci; idx++)
    {
      sequences->main[idx] = A->simMatrix->bases[Xb[idx]];
      sequences->match[idx] = A->simMatrix->hyphen;
    }
    rs[0] = 0;
    rs[1] = 0;
    sequences->length = size_data;
    sequences->backing_memory = size_data;
    return;
  }

  if(dir == 0 || dir & T[ci][cj])
    dir = T[ci][cj] >> 2;

  if(dir == 0)
  {
    rs[0] = -1;
    rs[1] = -1;
    sequences->length = 0;
    sequences->backing_memory = 0;
    return;
  }

  Ci = A->simMatrix->bases[Xb[ci]];
  Cj = A->simMatrix->bases[Yb[cj]];

  if(dir & 4)
  {
    tracepath_x(A, n, m, T, Yb, Xb, ci-1, cj-1, 0, sequences, rs, depth+1);
    if(rs[0] != -1 && rs[1] != -1)
    {
      sequences->main[sequences->length-depth-1] = Ci;
      sequences->match[sequences->length-depth-1] = Cj;
      return;
    }
  }
  if(dir & 2)
  {
    tracepath_x(A, n, m, T, Yb, Xb, ci-1, cj, 2, sequences, rs, depth+1);
    if(rs[0] != -1 && rs[1] != -1)
    {
      sequences->main[sequences->length-depth-1] = Ci;
      sequences->match[sequences->length-depth-1] = A->simMatrix->hyphen;
      return;
    }
  }
  if(dir & 1)
  {
    tracepath_x(A, n, m, T, Yb, Xb, ci, cj-1, 1, sequences, rs, depth+1);
    if(rs[0] != -1 && rs[1] != -1)
    {
      sequences->main[sequences->length-depth-1] = A->simMatrix->hyphen;
      sequences->match[sequences->length-depth-1] = Cj;
      return;
    }
  }
}

/*
 * Function alignPair() - Compute the global alignment for two sequences
 * Input-
 *     good_match_t *A - A structure used to generate this match
 *     int codon2bases[3][64] - codon translation data
 *     int gapPenalty - penelty for gaps, as defined in the parameters
 *     int disMatrix[4][4] - used to score codons
 *     seq_t *first_seq - Sequence X
 *     seq_t *second_seq - Sequence Y
 *     seq_t *sequences - Storage area for output.  Needs to be allocated before calling
 * Output-
 *     seq_t *sequences - main has X aligned against Y and match has Y aligned against X
 */

int alignPair(good_match_t *A, int codon2bases[3][64], int gapPenalty, int disMatrix[4][4], seq_t *first_seq, seq_t *second_seq, seq_t *sequences)
{
  codon_t X[first_seq->length];
  codon_t Xb[(first_seq->length) * 3];
  index_t length_x = scrub_hyphens(A, X, first_seq->main, first_seq->length );
  index_t n = length_x * 3;

  codon_t *Y = second_seq->main;
  codon_t Yb[(second_seq->length) * 3];

  index_t m = second_seq->length * 3;

  score_t T[n][m];

  memset(T, 0, sizeof(score_t)*m*n);

  score_t V[m], Vg[m], F[m], E[m+1], G[m];
  E[m] = 0;
  int rs[2];
  int adder;

  for(int idx=0; idx < m; idx++)
  {
    V[idx] = (idx+1) * gapPenalty;
    F[idx] = (idx+1) * gapPenalty;
  }

  for(int idx=0; idx < length_x; idx++)
  {
    Xb[idx*3] = codon2bases[0][X[idx]];
    Xb[idx*3+1] = codon2bases[1][X[idx]];
    Xb[idx*3+2] = codon2bases[2][X[idx]];
  }

  for(int idx=0; idx < second_seq->length; idx++)
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
      V[jdx] = (F[jdx] < G[jdx]) ? F[jdx] : G[jdx];
    }

    E[0] = V[0] + gapPenalty;

    for(int jdx=1; jdx < m; jdx++)
    {
      E[jdx] = ((E[jdx-1] < V[jdx-1]) ? E[jdx-1] : V[jdx-1]) + gapPenalty;
      V[jdx] = (E[jdx] < V[jdx]) ? E[jdx] : V[jdx];
    }

    for(int jdx=0; jdx < m; jdx++)
    {
      T[idx][jdx] = 4*(V[jdx] == E[jdx]) + 8*(V[jdx] == F[jdx]) + 16*(V[jdx] == G[jdx]);
      Vg[jdx] = V[jdx] + gapPenalty;
      F[jdx] = ((F[jdx] + gapPenalty) < Vg[jdx]) ? (F[jdx] + gapPenalty) : Vg[jdx];
    }

    for(int jdx=0; jdx < m; jdx++)
    {
      // The last value in E should always be 0
      T[idx][jdx] = T[idx][jdx] + (E[jdx+1] == Vg[jdx]) + 2 * (F[jdx] == Vg[jdx]);
    }
  }

  // Debug, dumpt matrix T
  /*
    for(int tx=0; tx < n; tx++)
    {
    for(int ty=0; ty < m; ty++)
    printf("%2hi ", T[tx][ty]);
    printf("\n");
    }
  */

  tracepath_x(A, n, m, T, Yb, Xb, n-1, m-1, 0, sequences, rs, 0);

  if(rs[0] != 0 || rs[1] != 0)
  {
    printf("multipleAlign tracepath bug");
    abort();
  }

  int ret_val = (int)V[m-1];

  return ret_val;
}

/*
 * finds the index of the lowest value in 'scan'
 * Input-
 *    int *scan - array to scan
 *    int size  - size of scan
 * Output-
 *    int       - first index of the lowest value
 */

index_t lowest_index(score_t *scan, index_t size)
{
  index_t current_index = 0;

  for(index_t idx=1; idx < size; idx++)
  {
    if(scan[idx] < scan[current_index])
      current_index = idx;
  }
 
  return current_index;
}

/*
 * Check to see if the two sequences are identical
 * Input-
 *    seq_t *source - Sequence to check
 *    seq_t *dest   - Sequence to check
 * Output-
 *    int           - 1 if the sequences are identical, 0 otherwise
 */


//int is_aligned(seq_t *source, seq_t *dest) __attribute__((optimize(3)));


int is_aligned(seq_t *source, seq_t *dest)
{
  if(source->length != dest->length) return 0;
  for(int idx=0; idx < source->length; idx++)
  {
    if(source->main[idx] != dest->match[idx])
    {
      return 0;
    }
  }
  return 1;
}

// grab this from global_align.c
extern int inverse_identity_matrix[4][4];

/*
 * Releases the results of kernel 5
 * Input-
 *    ma_t *doomed - [ma_t[size_m]] array of ma_t structures
 *    int size_m   - [int] size of the array
 */

void release_ma(ma_t *doomed, index_t size_m)
{
  for(int idx=0; idx < size_m; idx++)
  {
    for(int jdx=0; jdx < doomed[idx].length; jdx++)
    {
      free(doomed[idx].alignment[jdx].main);
    }
    free(doomed[idx].alignment);
    free(doomed[idx].scores);
  }
  free(doomed);
}

/*
 * Function multipleAlign() - Kernel 5 - Multiple Sequence Alignment.
 *
 * This function uses the 'center star' method for 'sum of pairs' multiple
 * alignment, as described in Gusfield, Algorithms on Stings, Trees, and
 * Sequences, pg 348.  Using the alginment score data from Kernel 4, it finds
 * an ideal center sequence, one which minimizes the sum of its score versus
 * all of the other sequences.  (The sequences are sorted by those scores to
 * clarify the results.)  Then each sequence is aligned with the center
 * sequence and the alignments are merged by adding extra spaces as necessary.
 *
 * This function uses the same variant of the Smith-Waterman dynamic
 * programming algorithm as is used in Kernel 4, but includes the tracepath
 * function which generates an actual alignment.
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
 * maxReports       - [int] maximum number of endpoints reported
 * minSeparation    - [int] minimum startpoint separation in codons
 * size_s           - [int] size of the S array.  Will also be the size of the returned array
 * S                - [good_match_t[]] results from Kernel 3, similar sequences
 *   bestStarts       - [int[M][2]] main/match startpoints
 *   bestEnds         - [int[M][2]] main/match endpoints
 *   bestSeqs         - [seq_t[M]] main/match sequences
 *     main           - [int[M]] codon sequence from main sequence
 *     match          - [int[M]] codon sequence from match sequence
 *     length         - length of the codons, including gaps
 *   bestScores       - [int[M]] the scores for the bestSeqs
 *   bestLength       - [int] value of M.
 * GA               - [ga_t] alignments for set g
 *   length         - [int] size of the globalScores matrix
 *   globalScores   - [short[length][length]] alignment score (x,y)
 * OUTPUT
 * MA               - [1-D Data structure] alignments for set g
 *   scores         - [1-D integer vector] score for each alignment
 *   alignments     - [1-D cell array of strings] final alignments
 */

//#define DEBUGGER

ma_t *multipleAlign(good_match_t *A, index_t size_s, good_match_t *S[size_s], ga_t *GA, int misPenalty, int gapPenalty)
{
  int dis_matrix[4][4];
  int codon2bases[3][64];
  int M, center, global_score, i, flag;
  seq_t *C;
  seq_t *alignment;
  codon_t *realloc_temp;
  //int *sumScores;
  //int *order;

  ma_t *MA = (ma_t*)malloc(sizeof(ma_t)*size_s);
  memset(MA, 0, sizeof(ma_t)*size_s);

  //write_out_inputs(A, size_s, S, GA, misPenalty, gapPenalty);

  //return NULL;
  // create our distance matrix
  for(int idx=0; idx<4; idx++)
    for(int jdx=0; jdx<4; jdx++)
      dis_matrix[idx][jdx] = inverse_identity_matrix[idx][jdx] * misPenalty;

  // create the codon decoding matrix
  for(int idx=0; idx < 64; idx++)
  {
    codon2bases[0][idx] = idx / 16;
    codon2bases[1][idx] = (idx / 4) % 4;
    codon2bases[2][idx] = idx % 4;
  }

  // main loop
  for(int g=0; g < size_s; g++)
  {
    M=S[g]->bestLength;
    MA[g].scores = NULL;
    MA[g].alignment = NULL;
    MA[g].length = 0;

    // if this entry is empty, skip it
    if(M==0)
    {
      continue;
    }

    score_t sumScores[M];
    index_t order[M];
    memset(sumScores, 0, sizeof(score_t)*M);

    // sum over all sequence matches
    for(int idx=0; idx < GA[g].length; idx++)
      for(int jdx=0; jdx < GA[g].length; jdx++)
        sumScores[idx] += GA[g].globalScores[idx][jdx];

    // find the minimum
    center = lowest_index(sumScores, M);

    MA[g].scores = (score_t*)malloc(M*sizeof(score_t));
    MA[g].alignment = (seq_t*)malloc(M*sizeof(seq_t));
    MA[g].length = M;
    memset(MA[g].scores, 0, M*sizeof(score_t));
    memset(MA[g].alignment, 0, M*sizeof(seq_t));

    // make a copy of the scores
    for(int idx=0; idx < M; idx++) MA[g].scores[idx] = GA[g].globalScores[center][idx];
    index_sort(MA[g].scores, order, M);
    C = (seq_t*)malloc(sizeof(seq_t));
    C->main = (codon_t*)malloc(sizeof(codon_t)*S[g]->bestSeqs[center].length);
    C->match = NULL;
    C->length = scrub_hyphens(A, C->main, S[g]->bestSeqs[center].main, S[g]->bestSeqs[center].length);

    // Decode the codon sequence
    MA[g].alignment[0].main = (codon_t*)malloc(sizeof(codon_t) * C->length * 3);
    for(int idx=0; idx < C->length; idx++)
    {
      MA[g].alignment[0].main[idx*3] = A->simMatrix->codon[C->main[idx]][0];
      MA[g].alignment[0].main[idx*3+1] = A->simMatrix->codon[C->main[idx]][1];
      MA[g].alignment[0].main[idx*3+2] = A->simMatrix->codon[C->main[idx]][2];
    }
    MA[g].alignment[0].length = C->length * 3;
    MA[g].alignment[0].backing_memory = C->length * 3;

    // Align each other sequence with the center sequence.
    for(int x=1; x < M; x++)
    {
      alignment = (seq_t*)malloc(sizeof(seq_t));
      global_score = alignPair(A, codon2bases, gapPenalty, dis_matrix, &(S[g]->bestSeqs[order[x]]), C, alignment);
      alignment->backing_memory = alignment->length;

      if(global_score!=MA[g].scores[x]) // debug statement
        printf("Waring!  Newly computed score does not agree with previous score! global_score=%i MA[%i].scores[%i]=%i\n", 
               (int)global_score, g, x, (int)MA[g].scores[x]);

      // add spaces to make the centers match
      i=0;
#ifdef DEBUGGER
      printf("Entering is_aligned loop\n");
#endif
      while(!is_aligned(&(MA[g].alignment[0]), alignment))
      {
        //skip = 0;
        if(i >= alignment->length) // extend cAligned
        {
          flag = 1;
        }
        else if(i >= MA[g].alignment[0].length) // extend alignment table
        {
          flag = 0;
        }
        else if(alignment->match[i] == MA[g].alignment[0].main[i]) // this one matches
        {
          i++;
          continue;
          //skip = 1;
        }
        else  // extend the one that is not a gap
        {
          flag = (MA[g].alignment[0].main[i] == A->simMatrix->hyphen || MA[g].alignment[0].main[i] == 45);
        }


        //if(!skip)
        //{
        #ifdef DEBUGGER
        if(g==64)
        {
          printf("Not aligned: g=%i flag=%i alignment->length=%i i=%i\n", g, flag, alignment->length, i);
          print_seq(&(MA[g].alignment[0]));
          print_seq(alignment);
          getchar();
        }
        #endif

        if(flag)
        {
          if(alignment->backing_memory == alignment->length)
          {
            realloc_temp = (codon_t*)realloc(alignment->main, sizeof(codon_t) * (alignment->backing_memory + REALLOC_INCREMENT));
            if(realloc_temp == NULL)
            {
              printf("Realloc error\n");
              abort();
            }
            alignment->main = realloc_temp;

            realloc_temp = (codon_t*)realloc(alignment->match, sizeof(codon_t) * (alignment->backing_memory + REALLOC_INCREMENT));
            if(realloc_temp == NULL)
            {
              printf("Realloc error\n");
              abort();
            }
            alignment->match = realloc_temp;
            alignment->backing_memory = alignment->backing_memory + REALLOC_INCREMENT;
          }

          for(int idx=alignment->length; idx > i; idx--)
          {
            alignment->main[idx] = alignment->main[idx-1];
            alignment->match[idx] = alignment->match[idx-1];
          }

          alignment->main[i] = A->simMatrix->hyphen;
          alignment->match[i] = A->simMatrix->hyphen;

          alignment->length++;
        }
        else
        {
          for(int j=0; j < x; j++)
          {
            if(MA[g].alignment[j].backing_memory == MA[g].alignment[j].length)
            {
              realloc_temp = (codon_t*)realloc(MA[g].alignment[j].main, sizeof(codon_t) * (MA[g].alignment[j].backing_memory + REALLOC_INCREMENT));
              if(realloc_temp == NULL)
              {
                printf("Realloc error\n");
                abort();
              }
              MA[g].alignment[j].main = realloc_temp;
              MA[g].alignment[j].backing_memory += REALLOC_INCREMENT;
            }

            for(int idx=MA[g].alignment[j].length; idx > i; idx--)
            {
              MA[g].alignment[j].main[idx] = MA[g].alignment[j].main[idx-1];
            }

            MA[g].alignment[j].main[i] = A->simMatrix->hyphen;
            MA[g].alignment[j].length++;
          }
        }
        //}
      }
#ifdef DEBUGGER
      printf("done aligning\n");
      print_seq(&(MA[g].alignment[0]));
      print_seq(alignment);
#endif
      MA[g].alignment[x].main = alignment->main;
      MA[g].alignment[x].match = NULL;
      free(alignment->match);
      MA[g].alignment[x].length = alignment->length;
      MA[g].alignment[x].backing_memory = alignment->backing_memory;
      free(alignment);
    }
    free(C->main);
    free(C);
  }

  return MA;
}

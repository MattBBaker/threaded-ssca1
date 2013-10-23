/*
*********************************************

Copyright 2008, UT-Battelle, LLC.
All rights Reserved.
See file LICENSING for licensing information.

*********************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <sort.h>
#include <pairwise_align.h>
#include <string.h>

typedef struct
{
  int *goodEnds[2];
  int *goodScores;
  int report;
  int size;
} current_ends_t;

/* consider adding this endpoint, possibly eliminating near by endpoints */

void considerAdding(int score, int minSeparation, int main_index, int match_index, 
                    int maxReports, current_ends_t *score_ends)
{
  int elements_to_copy;
  //printf("Considering\n");

  for(int r=score_ends->report-1; r>=0; r--)
  {
    if((main_index - score_ends->goodEnds[0][r]) >= minSeparation) break; // retain point r
    if(abs(match_index - score_ends->goodEnds[1][r]) >= minSeparation) continue;  // if not near by
    if(score_ends->goodScores[r] > score) return; // discard new point, maybe others

    // discard point r
    elements_to_copy = score_ends->report;

    for(int i = r; i < elements_to_copy; i++)
    {
      score_ends->goodScores[i]=score_ends->goodScores[i+1];
      score_ends->goodEnds[0][i]=score_ends->goodEnds[0][i+1];
      score_ends->goodEnds[1][i]=score_ends->goodEnds[1][i+1];
    }
    score_ends->report--;
  }

  // debug code, V[j] typically should not be this high

  //if(V[j] > 100)
  //  printf("adding element: %i at index: %i\n", V[j], *report);

  //printf("adding element: %i at index: %i\n", V[j], report);
  
  // add a new point
  score_ends->goodScores[score_ends->report]=score;
  score_ends->goodEnds[0][score_ends->report]=main_index;
  score_ends->goodEnds[1][score_ends->report]=match_index;
  score_ends->report++;

  // When the table is full, sort and discard all but the best end points.
  // Keep the table in entry order, just compact-out the discarded entries.
  if(score_ends->report == score_ends->size)
  {
    int worst_keeper;
    int new_best_index=0;

    int *index_array = malloc(score_ends->size*sizeof(int));
    int *sorted_array = malloc(score_ends->size*sizeof(int));
    int *best_index = malloc(score_ends->size*sizeof(int));

    sorted_array[0]=0;

    memcpy(sorted_array, score_ends->goodScores, sizeof(int)*score_ends->size);
    index_sort(sorted_array, index_array, score_ends->report);

    worst_keeper = score_ends->size - maxReports;
    //*minScore = sorted_array[worst_keeper] + 1;

    for(int index_for_index=worst_keeper; index_for_index < score_ends->size; index_for_index++)
    {
      best_index[new_best_index] = index_array[index_for_index];
      new_best_index++;
    }

    sort(best_index, new_best_index);
    for(int idx=0; idx < new_best_index; idx++)
    {
      score_ends->goodScores[idx] =score_ends->goodScores[best_index[idx]];
      score_ends->goodEnds[0][idx]=score_ends->goodEnds[0][best_index[idx]];
      score_ends->goodEnds[1][idx]=score_ends->goodEnds[1][best_index[idx]];
    }
    score_ends->report = maxReports;

    free(index_array);
    free(sorted_array);
    free(best_index);
  } 
}

/* release_good_match:
 * Free the goot_match_t structure generated by Kernel 1 and Kernel 2
 * Note: You can still release a matrix that has not been through Kernel 2,
 *       there are not ill side effects.
 * Input:
 *     good_match_t *doomed - the structure to be freed
 * Output:
 *     None
 */

void release_good_match(good_match_t *doomed)
{
  if(doomed==NULL) return;
  free(doomed->goodScores);
  free(doomed->goodEnds[1]);
  free(doomed->goodEnds[0]);
  free(doomed->bestStarts[0]);
  free(doomed->bestStarts[1]);
  free(doomed->bestEnds[0]);
  free(doomed->bestEnds[1]);
  free(doomed->bestScores);
  for(int idx=0; idx<doomed->bestLength; idx++)
  {
    free(doomed->bestSeqs[idx].main);
    free(doomed->bestSeqs[idx].match);
  }
  free(doomed->bestSeqs);
  free(doomed);
}

/* pairwise_align 
 * real meat of the program, this function finds codon similarities in seq_data using the matrix sim_matrix
 * Input:
 *   seq_data_t *seq_data     - Sequence data generated by genScalData()
 *   sim_matrix_t *sim_matrix - Codon similarity matrix generated by genSimMatrix()
 *   int minScore             - Minimum end point score, from the init_parameters() function
 *   int maxReports           - Maximum number of reports to keep, from the init_parameters() function
 *   int minSeparation        - Minimum end point seperation in codons, from the init_parameters() function
 *
 *  Output:
 *    good_matrix_t * - a matrix of good matches
 *       ->simMatrix  - a pointer to the sim_matrix_t used
 *       ->seqData    - a pointer to the seq_data_t used
 *	 ->goodEnds   - a [2][maxReports] matrix with main/match endpoints
 *       ->goodScores - a [maxReports] good scores for upto maxReports endpoints
 *       ->numReports - an integer, the number of reports represented
 */

typedef struct
{
  seq_data_t *seq_data;
  sim_matrix_t *sim_matrix;
  current_ends_t *good_ends;
  int search_length;
  int max_match;
  int start_offset;
  int min_score;
  int min_separation;
  int max_reports;
} payload_t;

#define index2d(x,y,stride) ((y) + ((x) * (stride)))

//void *pairwise_worker(void *data) {
good_match_t *pairwise_align(seq_data_t *seq_data, sim_matrix_t *sim_matrix, int minScore, int maxReports, int minSeparation) {
  const int sortReports = maxReports * 3;
  const int *mainSeq = seq_data->main;
  const int *matchSeq = seq_data->match;
  const int gapExtend = sim_matrix->gapExtend;
  const int gapFirst = sim_matrix->gapStart + gapExtend;
  current_ends_t *good_ends = malloc(sizeof(current_ends_t));
  good_ends->size = sortReports;
  good_ends->report = 0;
  good_ends->goodScores = malloc(sizeof(int)*sortReports);
  good_ends->goodEnds[0] = malloc(sizeof(int)*sortReports);
  good_ends->goodEnds[1] = malloc(sizeof(int)*sortReports);

  int *score_matrix = malloc(sizeof(int)*3*seq_data->matchLen);
  int *match_gap_matrix = malloc(sizeof(int)*2*seq_data->matchLen);
  int *main_gap_matrix = malloc(sizeof(int)*2*seq_data->matchLen);
  long score_start, score_end;
  int G, W, E, F, cmp_a, cmp_b;
  long m, n;
  int max_values=0;
  int *index_array;
  int *sort_array;
  good_match_t *answer;

  W = sim_matrix->similarity[mainSeq[0]][matchSeq[0]];
  score_matrix[index2d(0,0,seq_data->matchLen)] = 0 > W ? 0 : W;
  main_gap_matrix[0] = -gapFirst + W;
  match_gap_matrix[0] = -gapFirst + W;

  W = sim_matrix->similarity[mainSeq[0]][matchSeq[1]];
  G = W;
  E = main_gap_matrix[index2d(0,0,seq_data->matchLen)];
  cmp_a = 0 > E ? 0 : E;
  cmp_a = cmp_a > G ? cmp_a : G;
  score_matrix[index2d(1,1,seq_data->matchLen)] = cmp_a;
  cmp_a = E - gapExtend;
  cmp_b = G - gapFirst;
  main_gap_matrix[index2d(1,0,seq_data->matchLen)] = cmp_a > cmp_b ? cmp_a : cmp_b;
  match_gap_matrix[index2d(1,0,seq_data->matchLen)] = -gapFirst > cmp_b ? -gapFirst : cmp_b;

  W = sim_matrix->similarity[mainSeq[1]][matchSeq[0]];
  G = W;
  F = match_gap_matrix[index2d(0,0,seq_data->matchLen)];
  cmp_a = 0 > F ? 0 : F;
  cmp_a = cmp_a > G ? cmp_a : G;
  score_matrix[index2d(1,0,seq_data->matchLen)] = cmp_a;
  cmp_a = F - gapExtend;
  cmp_b = G - gapFirst;
  main_gap_matrix[index2d(1,1,seq_data->matchLen)] = -gapFirst > cmp_b ? -gapFirst : cmp_b;
  match_gap_matrix[index2d(1,1,seq_data->matchLen)] = cmp_a > cmp_b ? cmp_a : cmp_b;
  
  for(int idx=2; idx < seq_data->matchLen * 2; idx++) {
    score_start = (idx-(seq_data->matchLen-1)) > 0 ? (idx-(seq_data->matchLen-1)) : 0;
    score_end = (idx) < (seq_data->matchLen-1) ? (idx) : (seq_data->matchLen-1);
    if(idx < seq_data->matchLen) {
      m = 0;
      n = idx;
      W = sim_matrix->similarity[mainSeq[m]][matchSeq[n]];
      G = W;
      F = match_gap_matrix[index2d((idx-1)%2,n-1,seq_data->matchLen)];
      cmp_a = F > 0 ? F : 0;
      cmp_a = cmp_a > G ? cmp_a : G;
      score_matrix[index2d((idx%3),m,seq_data->matchLen)] = cmp_a;
      if((W > 0 && cmp_a > minScore && cmp_a == G) &&
         ((m == seq_data->matchLen - 1) || (n == seq_data->matchLen - 1) ||
          (sim_matrix->similarity[mainSeq[m+1]][matchSeq[n+1]] <= 0))) {
        considerAdding(score_matrix[index2d((idx%3),m,seq_data->matchLen)], minSeparation, m, n, maxReports, good_ends);
      }
      cmp_a = F - gapExtend;
      cmp_b = G - gapFirst;
      match_gap_matrix[index2d(idx%2,n,seq_data->matchLen)] = cmp_a > cmp_b ? cmp_a : cmp_b;
      score_start = score_start+1;

      m = idx;
      n = 0;
      E = main_gap_matrix[index2d((idx-1)%2,m-1,seq_data->matchLen)];
      cmp_a = E > 0 ? E : 0;
      cmp_a = cmp_a > G ? cmp_a : G;
      score_matrix[index2d((idx%3),m,seq_data->matchLen)] = cmp_a;
      if((cmp_a > minScore && W > 0 && cmp_a == G) &&
         ((m == seq_data->matchLen - 1) || (n == seq_data->matchLen - 1) ||
          (sim_matrix->similarity[mainSeq[m+1]][matchSeq[n+1]] <= 0))) {
        considerAdding(score_matrix[index2d((idx%3),m,seq_data->matchLen)], minSeparation, m, n, maxReports, good_ends);
      }
      cmp_a = E - gapExtend;
      cmp_b = G - gapFirst;
      main_gap_matrix[index2d(idx%2,m,seq_data->matchLen)] = cmp_a > cmp_b ? cmp_a : cmp_b;
      score_end = score_end - 1;
    }

#pragma omp parallel for private(m,n,W,G,F,E,cmp_a,cmp_b) firstprivate(idx, score_matrix, sim_matrix, main_gap_matrix, match_gap_matrix, gapFirst, gapExtend, seq_data, minScore, minSeperation, maxReports, good_ends) schedule(static)
    for(int antidiagonal = score_start; antidiagonal <= score_end; antidiagonal++) {
      m = antidiagonal;
      n = idx - m;

      W = sim_matrix->similarity[mainSeq[m]][matchSeq[n]];
      G = score_matrix[index2d(((idx-2)%3),m-1,seq_data->matchLen)] + W;
      F = match_gap_matrix[index2d((idx-1)%2,n-1,seq_data->matchLen)];
      E = main_gap_matrix[index2d((idx-1)%2,m-1,seq_data->matchLen)];
      cmp_a = 0;
      cmp_a = cmp_a > E ? cmp_a : E;
      cmp_a = cmp_a > F ? cmp_a : F;
      cmp_a = cmp_a > G ? cmp_a : G;
      score_matrix[index2d((idx%3),m,seq_data->matchLen)] = cmp_a;
      if((cmp_a > minScore && W > 0 && cmp_a == G) &&
         ((m == seq_data->matchLen - 1) || (n == seq_data->matchLen - 1) ||
          (sim_matrix->similarity[mainSeq[m+1]][matchSeq[n+1]] <= 0))) {
#pragma omp critical
        considerAdding(score_matrix[index2d((idx%3),m,seq_data->matchLen)], minSeparation, m, n, maxReports, good_ends);
      }
      cmp_a = E - gapExtend;
      cmp_b = G - gapFirst;
      main_gap_matrix[index2d((idx)%2,m,seq_data->matchLen)] = cmp_a > cmp_b ? cmp_a : cmp_b;
      cmp_a = F - gapExtend;
      match_gap_matrix[index2d((idx)%2,n,seq_data->matchLen)] = cmp_a > cmp_b ? cmp_a : cmp_b;
    }
  }

  answer = malloc(sizeof(good_match_t));
  answer->simMatrix = sim_matrix;
  answer->seqData = seq_data;
  answer->goodEnds[0] = malloc(sizeof(int)*maxReports);
  answer->goodEnds[1] = malloc(sizeof(int)*maxReports);
  answer->goodScores = malloc(sizeof(int)*maxReports);

  memset(answer->goodEnds[0], 0, sizeof(int)*maxReports);
  memset(answer->goodEnds[1], 0, sizeof(int)*maxReports);
  memset(answer->goodScores, 0, sizeof(int)*maxReports);

  answer->bestEnds[0] = NULL;
  answer->bestStarts[0] = NULL;
  answer->bestEnds[1] = NULL;
  answer->bestStarts[1] = NULL;
  answer->bestSeqs = NULL;
  answer->bestScores = NULL;

  max_values=good_ends->report;

  if(max_values > maxReports) max_values = maxReports;

  index_array = malloc(sortReports * sizeof(int));
  sort_array = malloc(sortReports * sizeof(int));
  memset(sort_array, 0, sizeof(int)*sortReports);
  memcpy(sort_array, good_ends->goodScores, good_ends->report * sizeof(int));
  memset(index_array, 0, sizeof(int)*sortReports);
  index_sort(sort_array, index_array, good_ends->report);

  for(int idx=0; idx < max_values; idx++) {
    answer->goodScores[idx] = good_ends->goodScores[index_array[good_ends->report-(idx+1)]];
    answer->goodEnds[0][idx] = good_ends->goodEnds[0][index_array[good_ends->report-(idx+1)]];
    answer->goodEnds[1][idx] = good_ends->goodEnds[1][index_array[good_ends->report-(idx+1)]];
  }


  free(score_matrix);
  free(main_gap_matrix);
  free(match_gap_matrix);
  free(good_ends->goodScores);
  free(good_ends->goodEnds[0]);
  free(good_ends->goodEnds[1]);
  free(good_ends);

  free(sort_array);
  free(index_array);
  answer->numReports = max_values;
  return answer;
}

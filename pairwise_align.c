#include <stdio.h>
#include <stdlib.h>
#include <sort.h>
#include <pairwise_align.h>
#include <string.h>

int *indexes[2];
int index_size;
int target;

typedef int score_t;

score_t *E[3];
score_t *F[3];

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
  if(doomed == NULL) return;
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

int add_match(int **main_ends, int **match_ends, score_t **scores, int match_count, 
              score_t score, int main_end, int match_end)
{
  (*main_ends)[match_count] = main_end;
  (*match_ends)[match_count] = match_end;
  (*scores)[match_count] = score;
  match_count++;

  if(match_count == target)
  {
    target *= 2;
    int *new_main = realloc((*main_ends), sizeof(int) * target);
    int *new_match = realloc((*match_ends), sizeof(int) * target);
    score_t *new_scores = realloc((*scores), sizeof(score_t) * target);

    if(new_main == NULL || new_match == NULL || new_scores == NULL)
    {
      printf("Could not reallocate memory for matches!\n");
      abort();
    }

    (*main_ends) = new_main;
    (*match_ends) = new_match;
    (*scores) = new_scores;
  }

  return match_count;
}

#define max(x,y) ((x) > (y) ? (x) : (y) )
void score_diagnal(seq_data_t *seq_data, sim_matrix_t *sim_matrix, score_t *score_matrix[], int idx)
{
  int main_gap_start, match_gap_start, match_score, main_gap_extend, match_gap_extend;
  score_t cmp;
  int top, bottom;

  index_size=0;
  bottom = 1 > idx-seq_data->mainLen+1 ? 1 : idx-seq_data->mainLen+1;
  top = idx > seq_data->mainLen ? seq_data->mainLen : idx;

  for(int jdx=bottom; jdx < top; jdx++)
  {
    indexes[0][index_size] = jdx;
    indexes[1][index_size] = (idx - jdx);
    index_size++;
  }

#ifdef _OPENMP
#pragma omp parallel for private(main_gap_score,match_gap_score,match_score,cmp) shared(score_matrix,seq_data,sim_matrix,idx,indexes)
#endif
  for(int jdx=0; jdx < index_size; jdx++)
  {
    main_gap_start = score_matrix[(idx-1)%3][indexes[0][jdx]] -
      (sim_matrix->gapStart+sim_matrix->gapStart);
    main_gap_extend = E[(idx-1)%3][indexes[0][jdx]] - sim_matrix->gapExtend;
    E[idx%3][indexes[0][jdx]] = max(main_gap_start,main_gap_extend);

    match_gap_start = score_matrix[(idx-1)%3][indexes[0][jdx]-1] -
      (sim_matrix->gapStart+sim_matrix->gapStart);
    match_gap_extend = F[(idx-1)%3][indexes[0][jdx]-1] - sim_matrix->gapExtend;
    F[idx%3][indexes[0][jdx]] = max(match_gap_start,match_gap_extend);

    match_score = score_matrix[(idx-2)%3][indexes[0][jdx]-1] +
      sim_matrix->similarity[seq_data->main[indexes[0][jdx]]][seq_data->match[indexes[1][jdx]]];
    cmp = max(E[idx%3][indexes[0][jdx]], F[idx%3][indexes[0][jdx]]);
    cmp = max(cmp,match_score);
    cmp = max(cmp, 0);
    score_matrix[idx%3][indexes[0][jdx]] = cmp;
  }
}

int harvest_diagnal(seq_data_t *seq_data, score_t *score_matrix[], score_t minScore, int minSeparation, 
                    int **main_ends,int **match_ends,int **scores, int potential_count, int target, int idx)
{
#ifdef _OPENMP
#pragma omp parallel for shared(score_matrix,minScore, seq_data, main_ends, match_ends, scores, potential_count, target, minSeparation, indexes)
#endif
  for(int jdx=0; jdx < index_size; jdx++)
  {
    if(score_matrix[idx%3][indexes[0][jdx]] > minScore)
    {
#ifdef _OPENMP
#pragma omp critical
#endif
      potential_count = add_match(main_ends, match_ends, scores, potential_count,
                                  score_matrix[idx%3][indexes[0][jdx]],
                                  indexes[0][jdx],indexes[1][jdx]);
    }
  }

  return potential_count;
}

int check_closeness(int *best_mains, int *best_matches, int best_end, int best_length, int main, int match, int minSeperation)
{
  for(int idx=0; idx < best_end; idx++)
  {
    if(abs(best_mains[best_end-idx]-main) < minSeperation || abs(best_matches[best_end-idx]-match) < minSeperation)
      return 0;
  }

  return 1;
}

int extract_best(int *main_best, int *match_best, score_t *scores_best, 
                  int *main, int *match, score_t *scores, 
                  int max_matches, int total_matches, int minSeperation)
{
  int *sorted_scores = malloc(sizeof(score_t)*total_matches);
  int *sorted_mains = malloc(sizeof(int)*total_matches);
  int *sorted_matches = malloc(sizeof(int)*total_matches);
  int *indexes = malloc(sizeof(int)*total_matches);
  memcpy(sorted_scores, scores, sizeof(score_t) * total_matches);
  index_sort(sorted_scores, indexes, total_matches);
  int search_index = total_matches-1;
  int best_index = 0;

  for(int idx=0; idx < total_matches; idx++)
  {
    sorted_mains[idx] = main[indexes[idx]];
    sorted_matches[idx] = match[indexes[idx]];
  }
  
  while(best_index < max_matches && search_index >= 0)
  {
    if(check_closeness(sorted_mains, sorted_matches, best_index, total_matches,
                       main[indexes[search_index]], match[indexes[search_index]], minSeperation))
    {
      main_best[best_index] = main[indexes[search_index]];
      match_best[best_index] = match[indexes[search_index]];
      scores_best[best_index] = scores[indexes[search_index]];
      best_index++;
    }
    search_index--;
  }

  free(sorted_scores);
  free(sorted_mains);
  free(sorted_matches);
  free(indexes);

  return best_index;
}

good_match_t *pairwise_align(seq_data_t *seq_data, sim_matrix_t *sim_matrix, int minScore, int maxReports, int minSeparation)
{
  good_match_t *answer = malloc(sizeof(good_match_t));
  score_t *score_matrix[3];
  target = maxReports * 3;
  int *main_ends = malloc(sizeof(int)*target);
  int *match_ends = malloc(sizeof(int)*target);
  int *scores = malloc(sizeof(int)*target);
  int potential_count=0;

  score_matrix[0] = malloc(sizeof(score_t)*seq_data->mainLen);
  score_matrix[1] = malloc(sizeof(score_t)*seq_data->mainLen);
  score_matrix[2] = malloc(sizeof(score_t)*seq_data->mainLen);

  E[0] = malloc(sizeof(score_t)*seq_data->mainLen);
  E[1] = malloc(sizeof(score_t)*seq_data->mainLen);
  E[2] = malloc(sizeof(score_t)*seq_data->mainLen);

  F[0] = malloc(sizeof(score_t)*seq_data->mainLen);
  F[1] = malloc(sizeof(score_t)*seq_data->mainLen);
  F[2] = malloc(sizeof(score_t)*seq_data->mainLen);

  memset(score_matrix[0], '\0', sizeof(score_t)*seq_data->mainLen);
  memset(score_matrix[1], '\0', sizeof(score_t)*seq_data->mainLen);
  memset(score_matrix[2], '\0', sizeof(score_t)*seq_data->mainLen);

  memset(E[0], '\0', sizeof(score_t)*seq_data->mainLen);
  memset(E[1], '\0', sizeof(score_t)*seq_data->mainLen);
  memset(E[2], '\0', sizeof(score_t)*seq_data->mainLen);

  memset(F[0], '\0', sizeof(score_t)*seq_data->mainLen);
  memset(F[1], '\0', sizeof(score_t)*seq_data->mainLen);
  memset(F[2], '\0', sizeof(score_t)*seq_data->mainLen);

  indexes[0] = malloc(sizeof(int)*seq_data->mainLen);
  indexes[1] = malloc(sizeof(int)*seq_data->mainLen);

  answer->simMatrix = sim_matrix;
  answer->seqData = seq_data;

  for(int idx=2; idx < (seq_data->mainLen*2)-2; idx++)
  {
    score_diagnal(seq_data,sim_matrix,score_matrix,idx);
    potential_count = harvest_diagnal(seq_data, score_matrix, minScore, minSeparation, &main_ends, &match_ends, &scores, 
                                      potential_count, target, idx);
  }

  free(indexes[0]);
  free(indexes[1]);

  answer->goodScores=malloc(sizeof(score_t)*maxReports);
  answer->goodEnds[0]=malloc(sizeof(int)*maxReports);
  answer->goodEnds[1]=malloc(sizeof(int)*maxReports);

  answer->numReports = extract_best(answer->goodEnds[0],answer->goodEnds[1],answer->goodScores,
                                    main_ends,match_ends,scores,maxReports,potential_count,minSeparation);

  free(main_ends);
  free(match_ends);
  free(scores);

  free(score_matrix[0]);
  free(score_matrix[1]);
  free(score_matrix[2]);
  return answer;
}

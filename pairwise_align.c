/*
   This file is part of SSCA1.

   Copyright (C) 2008-2015, UT-Battelle, LLC.

   This product includes software produced by UT-Battelle, LLC under Contract No.
   DE-AC05-00OR22725 with the Department of Energy.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the New BSD 3-clause software license (LICENSE).

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   LICENSE for more details.

   For more information please contact the SSCA1 developers at:
   bakermb@ornl.gov
*/

#define _BSD_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <sort.h>
#include <pairwise_align.h>
#include <string.h>
#include <util.h>
#include <assert.h>
#include <unistd.h>

#include <sys/time.h>


#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

typedef struct
{
  index_t *goodEnds[2];
  score_t *goodScores;
  int report;
  int size;
  score_t min_score;
} current_ends_t;

static index_t unsigned_abs_diff(index_t A, index_t B) {
  return (A > B) ? (A - B) : (B - A);
}

void print_back_20(good_match_t *A, seq_t *in_seq, index_t end_index){
  index_t start = end_index - 20;
  char main_acid_chain[21];
  char main_codon_chain[61];
  seq_t this_seq;
  this_seq.sequence = in_seq->sequence+start;
  this_seq.length = 20;

  memset(main_acid_chain, '\0', 21);
  memset(main_codon_chain, '\0', 21);

  assemble_acid_chain(A, main_acid_chain, &this_seq, 20);

  assemble_codon_chain(A, main_codon_chain, &this_seq, 20);

  printf("%7ld  %s  %s  %7ld\n",
         start, main_acid_chain, main_codon_chain, end_index);
}

static void considerAdding(score_t score, int minSeparation, index_t main_index, index_t match_index,
                    int maxReports, current_ends_t *score_ends) {

  //first scan the list to see if there is a match already that is closer
  for(int idx=0; idx < score_ends->report; idx++){
    if(unsigned_abs_diff(score_ends->goodEnds[0][idx], main_index) < minSeparation || unsigned_abs_diff(score_ends->goodEnds[1][idx], match_index) < minSeparation) {
      if(score_ends->goodScores[idx] < score) {
        score_ends->goodEnds[0][idx] = main_index;
        score_ends->goodEnds[1][idx] = match_index;
        score_ends->goodScores[idx] = score;
        return;
      } else {
        return;
      }
    }
  }
  //enlarge if needed
  if(score_ends->report == score_ends->size) {
    index_t worst_keeper, new_best_index=0;
    index_t *index_array = NULL;
    score_t *sorted_array = NULL;
    index_t *best_index = NULL;

    if(index_array == NULL) index_array = (index_t *)malloc(score_ends->size*sizeof(index_t));
    if(sorted_array == NULL) sorted_array = (score_t *)malloc(score_ends->size*sizeof(score_t));
    if(best_index == NULL) best_index = (index_t *)malloc(score_ends->size*sizeof(index_t));

    memcpy(sorted_array, score_ends->goodScores, sizeof(score_t)*score_ends->report);
    index_sort(sorted_array, index_array, score_ends->report);

    worst_keeper = score_ends->size - maxReports;
    score_ends->min_score = score_ends->goodScores[best_index[worst_keeper]];

    for(int index_for_index=worst_keeper; index_for_index < score_ends->size; index_for_index++) {
      best_index[new_best_index] = index_array[index_for_index];
      new_best_index++;
    }

    sort(best_index, new_best_index);
    for(int idx=0; idx < new_best_index; idx++) {
      score_ends->goodScores[idx] =score_ends->goodScores[best_index[idx]];
      score_ends->goodEnds[0][idx]=score_ends->goodEnds[0][best_index[idx]];
      score_ends->goodEnds[1][idx]=score_ends->goodEnds[1][best_index[idx]];
    }
    score_ends->report = maxReports;
    free(index_array);
    free(sorted_array);
    free(best_index);
  }
  score_ends->goodEnds[0][score_ends->report] = main_index;
  score_ends->goodEnds[1][score_ends->report] = match_index;
  score_ends->goodScores[score_ends->report] = score;
  score_ends->report++;
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
    free_local_seq(doomed->bestSeqs[idx].main);
    free_local_seq(doomed->bestSeqs[idx].match);
  }
  free(doomed->bestSeqs);
  free(doomed);
}

//typedef score_t score_matrix_t;
typedef struct {
  score_t *scores;
  index_t length;
  index_t local_length;
} score_matrix_t;

typedef score_matrix_t gap_matrix_t;

static score_matrix_t *alloc_score_matrix(index_t matrix_length){
  score_matrix_t *new_alloc = (score_matrix_t *)malloc(sizeof(score_matrix_t));
  assert(new_alloc != NULL);
  new_alloc->length = matrix_length;
  new_alloc->local_length = (matrix_length / num_nodes);

  malloc_all(sizeof(score_t)*3*new_alloc->local_length, (void **)&new_alloc->scores);
  assert(new_alloc->scores != NULL);
  touch_memory(new_alloc->scores, sizeof(score_t)*3*new_alloc->local_length);
  return new_alloc;
}

static gap_matrix_t *alloc_gap_matrix(index_t matrix_length){
  gap_matrix_t *new_alloc = (gap_matrix_t *)malloc(sizeof(score_matrix_t));
  assert(new_alloc != NULL);
  new_alloc->length = matrix_length;
  new_alloc->local_length = (matrix_length / num_nodes);

  malloc_all(sizeof(score_t)*2*new_alloc->local_length, (void **)&new_alloc->scores);
  assert(new_alloc->scores != NULL);
  touch_memory(new_alloc->scores, sizeof(score_t)*2*new_alloc->local_length);
  return new_alloc;
}

static void free_score_matrix(score_matrix_t *doomed){
  FREE_ALL(doomed->scores);
  free(doomed);
}

static void free_gap_matrix(gap_matrix_t *doomed){
  FREE_ALL(doomed->scores);
  free(doomed);
}

#define index2d(x,y,stride) ((y) + ((x) * (stride)))

static void fetch_score(score_matrix_t *A, index_t m, index_t n, score_t *in){
#ifndef USE_NONE
  int target_ep = n / A->local_length;
  int local_index = n % A->local_length;
  SHORT_GET(in, &(A->scores[index2d(m%3,local_index,A->local_length)]), 1, target_ep);
#else
  *in = A->scores[index2d(m%3,n,A->local_length)];
#endif
}

static void fetch_gap(gap_matrix_t *A, index_t m, index_t n, score_t *in){
#ifndef USE_NONE
  int target_ep = n / A->local_length;
  int local_index = n %A->local_length;
  SHORT_GET(in, &(A->scores[index2d(m%2,local_index,A->local_length)]), 1, target_ep);
#else
  *in = A->scores[index2d(m%2,n,A->local_length)];
#endif
}

/* maybe useful for going to further extremes. Commented out because GCC complains of unused functions */
#if 0
static void fetch_score_nb(score_matrix_t *A, index_t m, index_t n, score_t *in){
  int target_ep = n / A->local_length;
  int local_index = n % A->local_length;

  SHORT_GET_NB((short*)in, &(A->scores[index2d(m%3,local_index,A->local_length)]), 1, target_ep);
}

static void fetch_gap_nb(gap_matrix_t *A, index_t m, index_t n, score_t *in){
  int target_ep = n / A->local_length;
  int local_index = n %A->local_length;

  SHORT_GET_NB((short*)in, &(A->scores[index2d(m%2,local_index,A->local_length)]), 1, target_ep);
}
#endif

static void assign_score(score_matrix_t *A, index_t m, index_t n, score_t new_value){
#ifndef USE_NONE
  int target_ep = n / A->local_length;
  int local_index = n %A->local_length;

  SHORT_PUT(&(A->scores[index2d(m%3,local_index,A->local_length)]), &new_value, 1, target_ep);
#else
  A->scores[index2d(m%3,n,A->local_length)] = new_value;
#endif    
}

static void assign_gap(score_matrix_t *A, index_t m, index_t n, score_t new_value){
#ifndef USE_NONE
  int target_ep = n / A->local_length;
  int local_index = n % A->local_length;

  SHORT_PUT(&(A->scores[index2d(m%2,local_index,A->local_length)]), &new_value, 1, target_ep);
#else
  A->scores[index2d(m%2,n,A->local_length)] = new_value;
#endif
}

#ifdef USE_SHMEM
long collect_pSync[_SHMEM_REDUCE_SYNC_SIZE];
int collect_pWrk[_SHMEM_REDUCE_MIN_WRKDATA_SIZE];
#endif

static void collect_thread_results(current_ends_t **good_ends, int max_reports)
{
    int max_threads = omp_get_max_threads();
    index_t copied = 0;
    sort_ends_t *sorted_list = malloc(sizeof(sort_ends_t)*max_threads * max_reports * 10);

    for(int idx=0; idx < max_threads; idx++)
    {
        for(int jdx=0; jdx < good_ends[idx]->report; jdx++)
        {
            sorted_list[copied].score = good_ends[idx]->goodScores[jdx];
            sorted_list[copied].main_end = good_ends[idx]->goodEnds[0][jdx];
            sorted_list[copied].match_end = good_ends[idx]->goodEnds[1][jdx];
            copied++;
        }
    }

    ends_sort(sorted_list, copied);
    good_ends[0]->report = copied > max_reports ? max_reports : copied;

    for(int idx=0; idx < good_ends[0]->report; idx++)
    {
        good_ends[0]->goodScores[idx] = sorted_list[idx].score;
        good_ends[0]->goodEnds[0][idx] = sorted_list[idx].main_end;
        good_ends[0]->goodEnds[1][idx] = sorted_list[idx].match_end;
    }

    free(sorted_list);
}

static void collect_best_results(current_ends_t **good_ends, int max_reports, int max_threads, good_match_t *answer){
  index_t copied=0;
  int max_values=0;
  current_ends_t collected_ends, current_end;

  collect_thread_results(good_ends, max_reports);
  BARRIER_ALL();
  if(rank != 0) return;

  memset(answer->goodEnds[0], 0, sizeof(index_t)*max_reports);
  memset(answer->goodEnds[1], 0, sizeof(index_t)*max_reports);
  memset(answer->goodScores, 0, sizeof(score_t)*max_reports);

  collected_ends.goodScores = malloc(sizeof(score_t)*good_ends[0]->size*num_nodes);
  collected_ends.goodEnds[0] = malloc(sizeof(index_t)*good_ends[0]->size*num_nodes);
  collected_ends.goodEnds[1] = malloc(sizeof(index_t)*good_ends[0]->size*num_nodes);

#ifndef USE_NONE
  for(int idx=0; idx < num_nodes; idx++){
    GETMEM(&current_end, good_ends[0], sizeof(current_ends_t), idx);
    SHORT_GET(&(collected_ends.goodScores[copied]), good_ends[0]->goodScores, current_end.report, idx);
    LONG_GET(&(collected_ends.goodEnds[0][copied]), good_ends[0]->goodEnds[0], current_end.report, idx);
    LONG_GET(&(collected_ends.goodEnds[1][copied]), good_ends[0]->goodEnds[1], current_end.report, idx);
    copied += current_end.report;
  }
#else
  memcpy(collected_ends.goodScores, good_ends[0]->goodScores, good_ends[0]->report*sizeof(score_t));
  memcpy(collected_ends.goodEnds[0], good_ends[0]->goodEnds[0],good_ends[0]->report*sizeof(index_t));
  memcpy(collected_ends.goodEnds[1], good_ends[0]->goodEnds[1],good_ends[0]->report*sizeof(index_t));
  copied = good_ends[0]->report;
#endif
  
  sort_ends_t *sorted_list = malloc(sizeof(sort_ends_t)*copied);

  for(int idx=0; idx < copied; idx++){
    sorted_list[idx].score = collected_ends.goodScores[idx];
    sorted_list[idx].main_end = collected_ends.goodEnds[0][idx];
    sorted_list[idx].match_end = collected_ends.goodEnds[1][idx];
  }

  ends_sort(sorted_list, copied);

  if(copied > max_reports){
    max_values = max_reports;
  } else {
    max_values = copied;
  }

  for(int idx=0; idx < max_values; idx++){
    answer->goodScores[idx] = sorted_list[idx].score;
    answer->goodEnds[0][idx] = sorted_list[idx].main_end;
    answer->goodEnds[1][idx] = sorted_list[idx].match_end;
  }

  free(sorted_list);

  answer->numReports = max_values;
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

good_match_t *pairwise_align(seq_data_t *seq_data, sim_matrix_t *sim_matrix, const int minScore, const int maxReports, const int minSeparation) {
    const int sortReports = maxReports * 10;
    const seq_t *main_seq = seq_data->main;
    const seq_t *match_seq = seq_data->match;
    const score_t gapExtend = sim_matrix->gapExtend;
    const score_t gapFirst = sim_matrix->gapStart + gapExtend;
    const index_t main_len = seq_data->main->length;
    codon_t current_main, current_match;
    codon_t next_main, next_match;
    const int max_threads = omp_get_max_threads();
/*
    for(int idx=0; idx < _SHMEM_REDUCE_SYNC_SIZE; idx++) {
        collect_pSync[idx] = _SHMEM_SYNC_VALUE;
    }
*/
/*
    int gogogo=0;
    while(gogogo==0);
    shmem_barrier_all();
*/
    current_ends_t **good_ends = (current_ends_t **)malloc(sizeof(current_ends_t *)*max_threads);
//#pragma omp parallel for
    for(int jdx=0; jdx < max_threads; jdx++) {
        int idx = idx;
        malloc_all(sizeof(current_ends_t), (void **)&good_ends[idx]);
        good_ends[idx]->size = sortReports;
        good_ends[idx]->report = 0;
        malloc_all(sizeof(score_t)*sortReports, (void **)&good_ends[idx]->goodScores);
        malloc_all(sizeof(index_t)*sortReports, (void **)&good_ends[idx]->goodEnds[0]);
        malloc_all(sizeof(index_t)*sortReports, (void **)&good_ends[idx]->goodEnds[1]);
        good_ends[idx]->min_score = minScore;
    }

    score_matrix_t *restrict score_matrix = alloc_score_matrix(seq_data->match->length);
    gap_matrix_t *restrict main_gap_matrix = alloc_gap_matrix(seq_data->match->length);
    gap_matrix_t *restrict match_gap_matrix = alloc_gap_matrix(seq_data->match->length);

    //shmem_barrier_all();
    //shmem_barrier(0,0,2,&collect_pSync[0]);

    index_t score_start, score_end;
    score_t G, W, E, F, cmp_a, cmp_b, cmp_c, new_score, next_G, next_F, next_E;
    index_t m, n;

    good_match_t *answer;
    codon_t main_codon;
    codon_t match_codon;

    index_t local_main_start = seq_data->main->local_size * rank;
    index_t local_main_end = seq_data->main->local_size * (rank+1) - 1;

    /*
      if(rank == 0){
      printf("Ready to debug on PID=%i\n", getpid());
      int gogogo=0;
      while(gogogo==0){}
      }
    */
    //First iteration, done by hand. Basically idx=0 in the big loop
    if(rank == 0){
        fetch_from_seq(main_seq,0,&main_codon);
        fetch_from_seq(match_seq,0,&match_codon);

        W = sim_matrix->similarity[main_codon][match_codon];
        assign_score(score_matrix,0,0,0 > W ? 0 : W);
        assign_gap(main_gap_matrix,0,0,-gapFirst + W);
        assign_gap(match_gap_matrix,0,0,-gapFirst + W);

        //idx=1 m=0,1 n =1,0
        fetch_from_seq(main_seq,0,&main_codon);
        fetch_from_seq(match_seq,1, &match_codon);

        W = sim_matrix->similarity[main_codon][match_codon];
        G = W;
        fetch_gap(main_gap_matrix,0,0,&E);
        cmp_a = 0 > E ? 0 : E;
        cmp_a = cmp_a > G ? cmp_a : G;
        assign_score(score_matrix,1,1,cmp_a);
        cmp_a = E - gapExtend;
        cmp_b = G - gapFirst;

        assign_gap(main_gap_matrix,1,0,cmp_a > cmp_b ? cmp_a : cmp_b);
        assign_gap(match_gap_matrix,1,0,-gapFirst > cmp_b ? -gapFirst : cmp_b);

        fetch_from_seq(main_seq,1,&main_codon);
        fetch_from_seq(match_seq,0,&match_codon);

        W = sim_matrix->similarity[main_codon][match_codon];
        G = W;
        fetch_gap(match_gap_matrix,0,0, &F);
        cmp_a = 0 > F ? 0 : F;
        cmp_a = cmp_a > G ? cmp_a : G;
        assign_score(score_matrix,1,0,cmp_a);
        cmp_a = F - gapExtend;
        cmp_b = G - gapFirst;
        assign_gap(main_gap_matrix,1,1,-gapFirst > cmp_b ? -gapFirst : cmp_b);
        assign_gap(match_gap_matrix,1,1,cmp_a > cmp_b ? cmp_a : cmp_b);
    }
    for(index_t idx=2; idx < seq_data->match->length * 2 - 1; idx++) {
        BARRIER_ALL();
        score_start = idx > (seq_data->match->length - 1) ? (idx-(seq_data->match->length-1)) : 0;
        score_end = idx < (seq_data->match->length-1) ? (idx) : (seq_data->match->length-1);

        if(idx < seq_data->match->length) {
            if(rank == 0){
                m = 0;
                n = idx;
                fetch_from_seq(main_seq,m,&main_codon);
                fetch_from_seq(match_seq,n,&match_codon);
                W = sim_matrix->similarity[main_codon][match_codon];
                G = W;
                fetch_gap(match_gap_matrix,idx-1,n-1,&F);
                cmp_a = F > 0 ? F : 0;
                cmp_a = cmp_a > G ? cmp_a : G;
                assign_score(score_matrix,idx,m,cmp_a);
                new_score= cmp_a;

                if((new_score > good_ends[omp_get_thread_num()]->min_score && W > 0 && new_score == G)){
                    if (m+1 == seq_data->main->length || n == 0) {
                        considerAdding(new_score, minSeparation, m, n, maxReports, good_ends[omp_get_thread_num()]);
                    }
                    else {
                        fetch_from_seq(main_seq, m+1, &next_main);
                        fetch_from_seq(match_seq, n-1, &next_match);

                        if((m == main_len - 1) || (n == 0) || sim_matrix->similarity[next_main][next_match] <= 0){
                            considerAdding(new_score, minSeparation, m, n, maxReports, good_ends[omp_get_thread_num()]);
                        }
                    }
                }

                cmp_a = F - gapExtend;
                cmp_b = G - gapFirst;
                assign_gap(match_gap_matrix,idx,n,cmp_a > cmp_b ? cmp_a : cmp_b);

                m = idx;
                n = 0;
                fetch_from_seq(main_seq,m,&main_codon);
                fetch_from_seq(match_seq,n,&match_codon);
                W = sim_matrix->similarity[main_codon][match_codon];
                G = W;
                fetch_gap(main_gap_matrix, idx-1, m-1, &E);
                cmp_a = E > 0 ? E : 0;
                cmp_a = cmp_a > G ? cmp_a : G;
                assign_score(score_matrix,idx,m,cmp_a);
                new_score = cmp_a;

                if((new_score > good_ends[omp_get_thread_num()]->min_score && W > 0 && new_score == G)){
                    if (m+1 == seq_data->main->length || n == 0) {
                        considerAdding(new_score, minSeparation, m, n, maxReports, good_ends[omp_get_thread_num()]);
                    } else {

                        fetch_from_seq(main_seq, m+1, &next_main);
                        fetch_from_seq(match_seq, n-1, &next_match);

                        if((m == main_len - 1) || (n == 0) || sim_matrix->similarity[next_main][next_match] <= 0){
                            considerAdding(new_score, minSeparation, m, n, maxReports, good_ends[omp_get_thread_num()]);
                        }
                    }
                }

                cmp_a = E - gapExtend;
                cmp_b = G - gapFirst;
                assign_gap(main_gap_matrix, idx, m, cmp_a > cmp_b ? cmp_a : cmp_b);
            }
            score_start++;
            score_end = score_end - 1;
        }

        index_t local_start, local_end;

        if(score_end < local_main_start && score_start > local_main_end){

            local_start = 1;
            local_end = 0;

        } else {

            if(local_main_start > score_start){
                local_start = local_main_start;
            } else {
                local_start = score_start;
            }

            if(score_end > local_main_end){
                local_end = local_main_end;
            } else {
                local_end = score_end;
            }

            fetch_from_seq(main_seq,local_start,&next_main);
            fetch_from_seq(match_seq,idx - score_start,&next_match);
            fetch_score(score_matrix, (idx-2)%3, local_start-1, &next_G);
            fetch_gap(match_gap_matrix, idx-1, idx - (score_start+1), &next_F);
            fetch_gap(main_gap_matrix, idx-1, local_start-1, &next_E);
        }

        //As a note, this loop is the program execution time. If you're looking to optimize this benchmark, this is all that counts.
#pragma omp parallel for                                                \
    default(none) \
    private(m,n,current_main, current_match, F, E, G, W, cmp_a, cmp_b, cmp_c, new_score, next_main, next_match, next_G, next_E, next_F) \
    shared(good_ends,idx, local_start, local_end, main_seq,match_seq,main_gap_matrix,match_gap_matrix,seq_data,score_matrix,sim_matrix) \
    schedule(static) ordered
        for(index_t antidiagonal = local_start; antidiagonal <= local_end; antidiagonal++) {
            m = antidiagonal;
            n = idx - m;
#ifdef USE_PREFETCH
            current_main = next_main;
            current_match = next_match;
            G = next_G;
            F = next_F;
            E = next_E;

            if (m < (seq_data->main->length-1))
                fetch_from_seq_nb(main_seq, m+1, &next_main);
            if (n > 0)
                fetch_from_seq_nb(match_seq, n-1, &next_match);
            if (n > 1)
                fetch_gap(match_gap_matrix, idx-1, n-2, &next_F);

            fetch_gap(main_gap_matrix, idx-1, m, &next_E);
            fetch_score(score_matrix, (idx-2)%3, m, &next_G);
#else
            fetch_from_seq(main_seq, m, &current_main);
            fetch_from_seq(match_seq, n, &current_match);

            fetch_gap(match_gap_matrix, idx-1, n-1, &F);
            fetch_gap(main_gap_matrix, idx-1, m-1, &E);
            fetch_score(score_matrix, (idx-2)%3, m-1, &G);
#endif
            cmp_a = 0;
            cmp_a = cmp_a > E ? cmp_a : E;
            cmp_a = cmp_a > F ? cmp_a : F;
            W = sim_matrix->similarity[current_main][current_match];
            G += W;
            new_score = cmp_a > G ? cmp_a : G;

            if((new_score > good_ends[omp_get_thread_num()]->min_score && W > 0 && new_score == G)){

#ifdef USE_PREFETCH
                if (m+1 == seq_data->main->length || n == 0) {
                    considerAdding(new_score, minSeparation, m, n, maxReports, good_ends[omp_get_thread_num()]);
                } else {
                    WAIT_NB();

                    if((m == main_len - 1) || (n == 0) || sim_matrix->similarity[next_main][next_match] <= 0){
                        considerAdding(new_score, minSeparation, m, n, maxReports, good_ends[omp_get_thread_num()]);
                    }
                }
#else
                if (m+1 == seq_data->main->length || n == 0) {
                    considerAdding(new_score, minSeparation, m, n, maxReports, good_ends[omp_get_thread_num()]);
                } else {
                    fetch_from_seq(main_seq, m+1, &next_main);
                    fetch_from_seq(match_seq, n-1, &next_match);
                    if((m == main_len - 1) || (n == 0) || sim_matrix->similarity[next_main][next_match] <= 0){
                        considerAdding(new_score, minSeparation, m, n, maxReports, good_ends[omp_get_thread_num()]);
                    }
                }
#endif
            }
            cmp_a = E - gapExtend;
            cmp_b = G - gapFirst;
            cmp_c = F - gapExtend;

#ifdef USE_PREFETCH
            WAIT_NB();
#endif
            assign_score(score_matrix,idx,m,new_score);
            assign_gap(main_gap_matrix, idx, m, cmp_a > cmp_b ? cmp_a : cmp_b);
            assign_gap(match_gap_matrix, idx, n, cmp_c > cmp_b ? cmp_c : cmp_b);
    
        }
    }

    answer = (good_match_t*)malloc(sizeof(good_match_t));
    answer->simMatrix = sim_matrix;
    answer->seqData = seq_data;
    answer->goodEnds[0] = (index_t*)malloc(sizeof(index_t)*maxReports);
    answer->goodEnds[1] = (index_t*)malloc(sizeof(index_t)*maxReports);
    answer->goodScores = (score_t*)malloc(sizeof(score_t)*maxReports);

    answer->bestEnds[0] = NULL;
    answer->bestStarts[0] = NULL;
    answer->bestEnds[1] = NULL;
    answer->bestStarts[1] = NULL;
    answer->bestSeqs = NULL;
    answer->bestScores = NULL;

    collect_best_results(good_ends, maxReports, max_threads, answer);

    free_score_matrix(score_matrix);
    free_gap_matrix(main_gap_matrix);
    free_gap_matrix(match_gap_matrix);

    for(int idx=0; idx < max_threads; idx++) {
        FREE_ALL(good_ends[idx]->goodScores);
        FREE_ALL(good_ends[idx]->goodEnds[0]);
        FREE_ALL(good_ends[idx]->goodEnds[1]);
        FREE_ALL(good_ends[idx]);
    }

    free(good_ends);

    return answer;
}

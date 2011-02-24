#ifndef _PAIRWISE_ALIGN_H
#define _PAIRWISE_ALIGN_H

#include <gen_sim_matrix.h>
#include <gen_scal_data.h>

#define index2d(x,y,stride) ((y) + ((x) * (stride)))

typedef struct _seq_t
{
  int *main;
  int *match;
  int length;
  int backing_memory; //NOTE: right now this is only used in multipleAlign, before that backing memory is the same as length
} seq_t;

// all pointers of of length numReports
typedef struct _good_match_t
{
  sim_matrix_t *simMatrix; // simMatrixed used to generate the matches
  seq_data_t *seqData; // sequences used to generate the matches
  int *goodEnds[2]; // end point for good sequences
  int *goodScores; // scores for the good sequences
  int numReports; // number of reports given back.
  int *bestStarts[2]; // location of the best starting points
  int *bestEnds[2]; // location of the best end points
  int *bestScores; // location of the best scores
  seq_t *bestSeqs; // list of the best sequences
  int bestLength;
} good_match_t;

good_match_t *pairwise_align(seq_data_t *seq_data, sim_matrix_t *sim_matrix, int K1_MIN_SCORE, int K1_MAX_REPORTS, int K1_MIN_SEPARATION, int threads);
void release_good_match(good_match_t *);

#endif

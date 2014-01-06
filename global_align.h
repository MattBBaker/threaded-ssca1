#ifndef __GLOBAL_ALIGN_H
#define __GLOBAL_ALIGN_H

#include <pairwise_align.h>
#include <types.h>

typedef struct _GA_struct
{
  index_t length;
  score_t **globalScores;
} ga_t;

ga_t *globalAlign(good_match_t *A, index_t size_s, good_match_t *S[], int misPenalty, int gapPenalty);
ga_t *test_ga();
void release_ga(ga_t *doomed, index_t size);
int verifyGlobal(ga_t G[], index_t size_g, int misPenalty, int gapPenalty, int maxDisplay);

#endif

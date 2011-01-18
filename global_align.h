#ifndef __GLOBAL_ALIGN_H
#define __GLOBAL_ALIGN_H

#include <pairwise_align.h>

typedef struct _GA_struct
{
  int length;
  short **globalScores;
} ga_t;

ga_t *globalAlign(good_match_t *A, int size_s, good_match_t *S[], int misPenalty, int gapPenalty);
ga_t *test_ga();
void release_ga(ga_t *doomed, int size);
int verifyGlobal(ga_t G[], int size_g, int misPenalty, int gapPenalty, int maxDisplay);

#endif

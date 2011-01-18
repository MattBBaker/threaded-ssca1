#ifndef __MULTIPLE_ALIGN_H
#define __MULTIPLE_ALIGN_H

#include <pairwise_align.h>
#include <global_align.h>

typedef struct _ma_t
{
  int *scores;
  seq_t *alignment;
  int length;
} ma_t;

ma_t *multipleAlign(good_match_t *A, int size_s, good_match_t *S[], ga_t *GA, int misPenalty, int gapPenalty);
void release_ma(ma_t *doomed, int length);
int verifyMultiple(ma_t MA[], int size_ma, int maxDisplay);

#endif

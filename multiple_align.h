#ifndef __MULTIPLE_ALIGN_H
#define __MULTIPLE_ALIGN_H

#include <pairwise_align.h>
#include <global_align.h>
#include <types.h>

typedef struct _ma_t
{
  score_t *scores;
  seq_data_t *alignment;
  index_t length;
} ma_t;

ma_t *multipleAlign(good_match_t *A, index_t size_s, good_match_t *S[], ga_t *GA, int misPenalty, int gapPenalty);
void release_ma(ma_t *doomed, index_t length);
int verifyMultiple(ma_t MA[], index_t size_ma, int maxDisplay);

#endif

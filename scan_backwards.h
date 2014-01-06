#ifndef _SCAN_BACKWARDS_H
#define _SCAN_BACKWARDS_H

#include <pairwise_align.h>
#include <types.h>

void scanBackward(good_match_t *A, int maxReports, int minSeparation);
int verify_alignment(good_match_t *A, int maxDisplay);

#endif

#ifndef _SORT_H
#define _SORT_H

#include <types.h>

void index_sort(score_t numbers[], index_t indexes[], index_t array_size);
void sort(index_t numbers[], index_t array_size);

typedef struct {
  score_t score;
  index_t main_end;
  index_t match_end;
} sort_ends_t;

void ends_sort(sort_ends_t *ends, index_t array_size);

#endif

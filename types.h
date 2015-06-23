#ifndef _TYPES_H
#define _TYPES_H

#include <stdint.h>

typedef uint64_t index_t;

extern int num_nodes;
extern int rank;

#if defined(USE_SHMEM) || defined(USE_MPI3)
typedef uint16_t codon_t;
typedef int16_t score_t;
#else
typedef int_fast8_t codon_t;
typedef int_fast16_t score_t;
#endif

typedef struct _sequence_t {
  codon_t *sequence;
  index_t length;
  index_t backing_memory;
  index_t local_size;
} seq_t;

typedef struct _seq_data_t {
  seq_t *main;
  seq_t *match;
  index_t max_validation;
} seq_data_t;

#endif

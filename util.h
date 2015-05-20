#ifndef __UTIL_H
#define __UTIL_H

#include <pairwise_align.h>
#include <types.h>

#ifdef SGI_SHMEM
#include <mpp/shmem.h>
#else
#include <shmem.h>
#endif

static inline int global_index_to_rank(seq_t *in, index_t codon_index){
  return codon_index / in->local_size;
}

static inline int global_index_to_local_index(seq_t *in, index_t codon_index){
  return codon_index % in->local_size;
}

codon_t _fetch_temp;

static inline codon_t fetch_from_seq(const seq_t *in, index_t codon_index){
  int target_pe = global_index_to_rank(in,codon_index);
  int local_index = global_index_to_local_index(in,codon_index);
  shmem_short_get(&_fetch_temp, &(in->sequence[local_index]), 1, target_pe);
  return _fetch_temp;
}

static inline void write_to_seq(const seq_t *in, index_t codon_index, codon_t data){
  int target_pe = global_index_to_rank(in,codon_index);
  int local_index = global_index_to_local_index(in,codon_index);
  shmem_short_put(&(in->sequence[local_index]), &data, 1, target_pe);
}

void distribute_rng_seed(unsigned int new_seed);
void seed_rng(int adjustment);
void touch_memory(void *mem, index_t size);
index_t scrub_hyphens(good_match_t *A, seq_t *dest, seq_t *source, index_t length);
void assemble_acid_chain(good_match_t *A, char *result, seq_t *chain, index_t length);
void assemble_codon_chain(good_match_t *A, char *result, seq_t *chain, index_t length);
score_t simple_score(good_match_t *A, seq_t *main, seq_t *match);
seq_t *alloc_seq(index_t seq_size);
void extend_seq(seq_t *extended, index_t extend_size);
void free_seq(seq_t *doomed);

#endif

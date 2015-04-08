#ifndef __UTIL_H
#define __UTIL_H

#include <pairwise_align.h>
#include <types.h>

void touch_memory(void *mem, index_t size);
index_t scrub_hyphens(good_match_t *A, seq_t *dest, seq_t *source, index_t length);
void assemble_acid_chain(good_match_t *A, char *result, seq_t *chain, index_t length);
void assemble_codon_chain(good_match_t *A, char *result, seq_t *chain, index_t length);
score_t simple_score(good_match_t *A, seq_t *main, seq_t *match);
seq_t *alloc_seq(index_t seq_size);
void extend_seq(seq_t *extended, index_t extend_size);
void free_seq(seq_t *doomed);

#endif

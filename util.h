#ifndef __UTIL_H
#define __UTIL_H

#include <pairwise_align.h>
#include <types.h>

void touch_memory(void *mem, index_t size);
index_t scrub_hyphens(good_match_t *A, codon_t *dest, codon_t *source, index_t length);
void assemble_acid_chain(good_match_t *A, char *result, codon_t *chain, index_t length);
void assemble_codon_chain(good_match_t *A, char *result, codon_t *chain, index_t length);
score_t simple_score(good_match_t *A, codon_t main[], codon_t match[], index_t length);

#endif

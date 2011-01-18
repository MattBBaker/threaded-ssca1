#ifndef __UTIL_H
#define __UTIL_H

#include <pairwise_align.h>

int scrub_hyphens(good_match_t *A, int *dest, int *source, int length);
void assemble_acid_chain(good_match_t *A, char *result, int *chain, int length);
void assemble_codon_chain(good_match_t *A, char *result, int *chain, int length);
int simple_score(good_match_t *A, int main[], int match[], int length);

#endif

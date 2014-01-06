#ifndef _SIM_MATRIX_H
#define _SIM_MATRIX_H

#include <types.h>

#define STAR ((int)48)
#define HYPHEN ((int)64)

typedef struct _sim_matrix_t
{
  score_t similarity[64][64]; /* - [2D array int] 1-based codon/codon similarity table */
  char aminoAcid[65];     /* - [1D char vector] 1-based codon to aminoAcid table */
  char bases[5];          /* - [1D char vector] 1-based encoding to base letter table */
  char codon[65][4];      /* - [64 x 3 char array] 1-based codon to base letters table */
  int encode[128];        /* - [int vector] aminoAcid character to last codon number */
  int hyphen;             /* - [int] encoding representing a hyphen (gap or space) */
  int star;               /* - [int] encoding representing a star */
  int exact;              /* - [integer] value for exactly matching codons */
  int similar;            /* - [integer] value for similar codons (same amino acid) */
  int dissimilar;         /* - [integer] value for all other codons */
  int gapStart;           /* - [integer] penalty to start gap (>=0) */
  int gapExtend;          /* - [integer] penalty to for each codon in the gap (>0) */
  int matchLimit;         /* - [integer] longest match including hyphens */
} sim_matrix_t;

sim_matrix_t *gen_sim_matrix(int exact, int similar, int dissimilar, int gapStart, int gapExtend, int matchLimit);
void release_sim_matrix(sim_matrix_t *doomed_matrix);

#endif

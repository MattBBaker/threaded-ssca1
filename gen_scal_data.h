#include <gen_sim_matrix.h>

#ifndef _GEN_SCAL_DATA
#define _GEN_SCAL_DATA

/*
 * seqData         - [structure] holds the generated sequences.
 *  main          - [1D uint8 array] generated main codon sequence (1-64).
 *  match         - [1D uint8 array] generated match codon sequence (1-64).
 *  maxValidation - [Integer] longest matching validation string.
 */

typedef struct _seq_data
{
  int *main;
  int *match;
  int maxValidation;
  int mainLen;
  int matchLen;
} seq_data_t;

seq_data_t *gen_scal_data( sim_matrix_t *simMatrix, int mainLen, int matchLen);
void release_scal_data(seq_data_t *doomed_scal_data);
void verifyData(sim_matrix_t *simMatrix, seq_data_t *seqData);

#endif

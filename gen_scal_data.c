#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sort.h>
#include <gen_scal_data.h>
#include <types.h>
#include <util.h>

#ifdef _OPENMP
#include <omp.h>
#endif

char validations[2][3][32] =
  {
    {"ACDEFG*IDENTICAL*HIKLMN", "ACDEFG*MISQRMATCHES*HIKLMN", "ACDEFG*STARTGAPMIDSTEND*HIKLMN"},
    {"MNLKIH*IDENTICAL*GFEDCA", "MNLKIH*MISRQMATCHES*GFEDCA", "MNLKIH*STARTMIDSTGAPEND*GFEDCA"}
  };

/* print the output */

void verifyData(sim_matrix_t *simMatrix, seq_data_t *seqData)
{
  printf("\n");
  printf("   Length of main sequence in codons: %lu\n", seqData->main->length);
  printf("  Length of match sequence in codons: %lu\n", seqData->match->length);
  printf("  Weight for exactly matching codons: %i\n", simMatrix->exact);
  printf("           Weight for similar codons: %i\n", simMatrix->similar);
  printf("        Weight for dissimilar codons: %i\n", simMatrix->dissimilar);
  printf("              Penalty to start a gap: %i\n", simMatrix->gapStart);
  printf("     Penalty for each codon in a gap: %i\n", simMatrix->gapExtend);
}

index_t *gen_indexes(int num_indexes, index_t rand_max, index_t min_seperation) {
  int not_done;
  index_t *indexes = (index_t*)malloc(sizeof(index_t) * num_indexes);
  do{
    not_done = 0;
    for(int idx=0; idx < num_indexes; idx++) {
      indexes[idx] = rand()%rand_max;
    }
    for(int idx=0; idx < num_indexes; idx++) {
      for(int jdx=idx+1; jdx < num_indexes; jdx++) {
        if(indexes[idx] - indexes[jdx] < min_seperation) not_done = 1;
      }
    }
  } while(not_done);
  return indexes;
}

void create_sequence(codon_t *sequence, char validations[][32], index_t sequence_length, int num_validations, sim_matrix_t *simMatrix){
  //void create_sequence(codon_t *sequence, index_t sequence_length){
  index_t total_length = sequence_length;
  index_t end;
  index_t *indexes = gen_indexes(num_validations, total_length-32, 32);

  for(index_t idx=0; idx < sequence_length; idx++) {
    sequence[idx] = rand()%64;
  }

  for(int idx=0; idx < num_validations; idx++) {
    end = strlen(validations[idx]);
    printf("Inserting sequence %s in location %lu\n", validations[idx], indexes[idx]);
    for(int jdx=0; jdx < end; jdx++){
      sequence[indexes[idx]+jdx] = simMatrix->encode[(int)validations[idx][jdx]];
    }
  }
  free(indexes);

}

//void insert_validation(codon_t 

seq_data_t *gen_scal_data( sim_matrix_t *simMatrix, index_t mainLen, index_t matchLen, int constant_rng) {
  seq_data_t *new_scal_data = (seq_data_t *)malloc(sizeof(seq_data_t));
  int validation_length, validation_size=0;

  memset(new_scal_data, '\0', sizeof(seq_data_t));

  for(int jdx=0; jdx < 3; jdx++){
    validation_length = strlen(validations[0][jdx]);
    validation_size += validation_length;
    if(validation_length > new_scal_data->max_validation)
     new_scal_data->max_validation = validation_length;
  }

  new_scal_data->max_validation -= 12;
  index_t main_size_with_validation = mainLen + validation_size;
  index_t match_size_with_validation = mainLen + validation_size;
  //new_scal_data->main->length = mainLen + validation_size;
  //new_scal_data->match->length = matchLen + validation_size;
  //Some fix ups to allow this to work distributed. Keep the amount of data the same across all nodes.
  if(main_size_with_validation % num_nodes != 0) {
    main_size_with_validation += (num_nodes - (main_size_with_validation % num_nodes));
  }

  if(match_size_with_validation % num_nodes!= 0) {
    match_size_with_validation += (num_nodes - (match_size_with_validation % num_nodes));
  }

  //new_scal_data->main->local_size = new_scal_data->main->length/num_nodes;
  //new_scal_data->match->local_size = new_scal_data->match->length/num_nodes;

  //new_scal_data->main->sequence = (codon_t*)ssca1_distributed_malloc(new_scal_data->main->local_size*sizeof(codon_t));
  //new_scal_data->match->sequence = (codon_t*)ssca1_distributed_malloc(new_scal_data->match->local_size*sizeof(codon_t));

  //new_scal_data->main = (codon_t*)ssca1_distributed_malloc(sizeof(codon_t)*new_scal_data->mainLen);
  //new_scal_data->match= (codon_t*)ssca1_distributed_malloc(sizeof(codon_t)*new_scal_data->matchLen);

  new_scal_data->main = alloc_seq(main_size_with_validation);
  new_scal_data->match = alloc_seq(match_size_with_validation);

  codon_t *gen_sequences[2] = {new_scal_data->main->sequence, new_scal_data->match->sequence};
  index_t seq_lengths[2] = {new_scal_data->main->local_size, new_scal_data->match->local_size};

#ifdef _OPENMP
  int thread_number = 2;
  if(constant_rng==1){
    thread_number = 1;
  }

#pragma omp parallel for num_threads(thread_number)
#endif
  for(int idx=0; idx < 2; idx++){
    touch_memory(gen_sequences[idx], sizeof(codon_t)*seq_lengths[idx]);
    //create_sequence(gen_sequences[idx], seq_lengths[idx]);
    create_sequence(gen_sequences[idx], validations[idx], seq_lengths[idx], 3, simMatrix);
  }

  return new_scal_data;
}

void release_scal_data(seq_data_t *doomed_seq) {
  //free((void *)(doomed_seq->main));
  //free((void *)(doomed_seq->match));
  free_seq(doomed_seq->main);
  free_seq(doomed_seq->match);
  free((void *)doomed_seq);
}

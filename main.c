#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <parameters.h>
#include <gen_sim_matrix.h>
#include <sys/time.h>
#include <time.h>
#include <gen_scal_data.h>
#include <pairwise_align.h>
#include <scan_backwards.h>
#include <locate_similar.h>
#include <global_align.h>
#include <multiple_align.h>
#include <limits.h>

unsigned int random_seed;

/* unsigned int get_dev_rand() - 
     routine to get an unsigned int from /dev/random.
   Intput-
     None
   Output-
     unsigned int- value from /dev/random
 */

unsigned int get_dev_rand()
{
  int rand_file = open("/dev/urandom", O_RDONLY);
  unsigned int value;
  ssize_t read_value = read(rand_file, &value, sizeof(unsigned int));
  close(rand_file);
  if(read_value != sizeof(unsigned int))
  {
    fprintf(stderr, "Error reading random values, aborting\n");
    abort();
  }
  return value;
}

/* void display_elapsed(struct timeval *start_time)
       Given a start time, display how much time has elapsed since thing
     Input- 
       struct timeval *start_time- pointer to the struct with the starting time
    Output-
       None
 */

void display_elapsed(struct timeval *start_time)
{
  struct timeval now;
  gettimeofday(&now, NULL);
  long hours_elapsed = 0;
  long minutes_elapsed = 0;
  time_t seconds_elapsed = now.tv_sec - start_time->tv_sec;
  long miliseconds_elapsed = 0;
  long u_elapsed = now.tv_usec - start_time->tv_usec;
  if(u_elapsed < 0)
  {
    seconds_elapsed--;
    u_elapsed = 1000000 + u_elapsed;
  }

  if(seconds_elapsed > 60)
  {
    minutes_elapsed = seconds_elapsed / 60;
    seconds_elapsed = seconds_elapsed % 60;
  }

  if(minutes_elapsed > 60)
  {
    hours_elapsed = minutes_elapsed / 60;
    minutes_elapsed = minutes_elapsed % 60;
  }

  if(u_elapsed > 1000)
  {
    miliseconds_elapsed = u_elapsed / 1000;
    u_elapsed = u_elapsed % 1000;
  }

  printf("\n\tElapsed time: %li hour(s), %li minute(s), %li second(s), %li milliseconds,  %li micro second(s).\n", hours_elapsed, minutes_elapsed, seconds_elapsed, miliseconds_elapsed, u_elapsed);
}

/* int main(int argc, char **argv)
     Entry routine.  Calls each kernel once, displaying elapsed time
 */

int main(int argc, char **argv)
{
  parameters_t global_parameters;
  sim_matrix_t *sim_matrix;
  seq_data_t *seq_data;
  struct timeval start_time;
  good_match_t *A;
  ga_t *GA;
  void *MA;

#ifdef _OPENMP
  printf("Running with OpenMP\n");
#endif

  init_parameters(&global_parameters);

  if(argc > 1 && !strcmp(argv[1],"--threads"))
  {
    global_parameters.threads = atoi(argv[2]);
  }
  else
  {
    global_parameters.threads = 1;
  }

  good_match_t *S[global_parameters.K2_MAX_REPORTS];
  memset(S, 0, sizeof(good_match_t *)*global_parameters.K2_MAX_REPORTS);

  printf("HPCS SSCA #1 Bioinformatics Sequence Alignment Executable Specification:\nRunning...\n");

  if(global_parameters.ENABLE_VERIF || global_parameters.CONSTANT_RNG)
  {
    printf("\n\tVerification run, using constant seed for RNG\n");
    // interesting values that have uncovered bugs in the past,
    // 2613174141 -- segfault caused by insert_validation producting two identical values
    // -550696422 -- segfault caused by the RNG producing 0.
    random_seed = (unsigned int)2613174141;
  }
  else
  {
    random_seed = (unsigned int)time(NULL); /* casting from time_t to unsigned int we can lose precision... no big deal here */
    random_seed += get_dev_rand();
  }

  printf("Using seed %u\n", random_seed);
  srand(random_seed);

  gettimeofday(&start_time, NULL);

  printf("\nScalable Data Generator - genScalData() beginning execution...\n");
  sim_matrix = gen_sim_matrix(global_parameters.SIM_EXACT, global_parameters.SIM_SIMILAR, global_parameters.SIM_DISSIMILAR, global_parameters.GAP_START, global_parameters.GAP_EXTEND, global_parameters.MATCH_LIMIT);

  seq_data = gen_scal_data(sim_matrix, global_parameters.MAIN_SEQ_LENGTH, global_parameters.MATCH_SEQ_LENGTH, global_parameters.CONSTANT_RNG);

  display_elapsed(&start_time);

  if(global_parameters.ENABLE_VERIF)
  {
    verifyData(sim_matrix, seq_data);
  }

  /* Kernel 1 run */

  printf("\nBegining Kernel 1 execution.\n");

  gettimeofday(&start_time, NULL);

  A=pairwise_align(seq_data, sim_matrix, global_parameters.K1_MIN_SCORE, global_parameters.K1_MAX_REPORTS, global_parameters.K1_MIN_SEPARATION);

  display_elapsed(&start_time);

  /* Kernel 2 run */

  printf("\nBegining Kernel 2 execution.\n");

  gettimeofday(&start_time, NULL);
 
  scanBackward(A, global_parameters.K2_MAX_REPORTS, global_parameters.K2_MIN_SEPARATION);

  display_elapsed(&start_time);

  if(global_parameters.ENABLE_VERIF)
  {
    verify_alignment(A, global_parameters.K2_DISPLAY);
  }

  /* Kernel 3 run */

  printf("\nBegining Kernel 3 execution.\n");

  gettimeofday(&start_time, NULL);

  locateSimilar(A,global_parameters.K3_MIN_SCORE, global_parameters.K3_MAX_REPORTS, global_parameters.K3_MIN_SEPARATION, global_parameters.K3_MAX_MATCH, S);

  display_elapsed(&start_time);

  if(global_parameters.ENABLE_VERIF)
    verify_similar(A, S, A->bestLength, global_parameters.K3_DISPLAY);

  /* Kernel 4 run */

  printf("\nBegining Kernel 4 execution.\n");

  gettimeofday(&start_time, NULL);

  GA = globalAlign(A, A->bestLength, S, global_parameters.MISMATCH_PENALTY, global_parameters.SPACE_PENALTY);

  display_elapsed(&start_time);

  if(global_parameters.ENABLE_VERIF)
    verifyGlobal(GA, A->bestLength, global_parameters.MISMATCH_PENALTY, global_parameters.SPACE_PENALTY, global_parameters.K4_DISPLAY);

  /* Kernel 5 run */

  printf("\nBegining Kernel 5 execution.\n");

  gettimeofday(&start_time, NULL);

  MA = multipleAlign(A, A->bestLength, S, GA, global_parameters.MISMATCH_PENALTY, global_parameters.SPACE_PENALTY);

  display_elapsed(&start_time);

  if(global_parameters.ENABLE_VERIF)
    verifyMultiple(MA, A->bestLength, global_parameters.K5_DISPLAY);

  for(int idx=0; idx < global_parameters.K3_MAX_REPORTS; idx++) release_good_match(S[idx]);
  release_ma(MA, A->bestLength);
  release_ga(GA, A->bestLength);
  release_good_match(A);
  release_sim_matrix(sim_matrix);
  release_scal_data(seq_data);

  return 0;
}

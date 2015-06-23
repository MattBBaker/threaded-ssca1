#ifndef __UTIL_H
#define __UTIL_H

#include <pairwise_align.h>
#include <types.h>
#include <stdio.h>

#ifdef USE_MPI3
#include <mpi.h>
extern MPI_Comm world;
extern MPI_Win window;
extern void *window_base;
extern size_t window_size;
extern void *next_window_address;
#else
#ifdef SGI_SHMEM
#include <mpp/shmem.h>
#else
#include <shmem.h>
#endif
#endif

#ifdef USE_MPI3
#define SHORT_GET(target, source, num_elems, rank)	MPI_Get(target, num_elems, MPI_SHORT, rank, (void *)source - window_base, num_elems, MPI_SHORT, window)
#else
#define SHORT_GET(target, source, num_elems, pe)	shmem_short_get(target, source, num_elems, pe)
#endif

#ifdef USE_MPI3
#define SHORT_GET_NB(target, source, num_elems, rank)	MPI_Get(target, num_elems, MPI_SHORT, rank, (void *)source - window_base, num_elems, MPI_SHORT, window) /* no non-blocking get available yet? */
#else
#define SHORT_GET_NB(target, source, num_elems, pe)	shmem_short_get_nb(target, source, num_elems, pe, NULL)
#endif

#ifdef USE_MPI3
#define LONG_GET(target, source, num_elems, rank)	MPI_Get(target, num_elems, MPI_LONG, rank, (void *)source - window_base, num_elems, MPI_LONG, window)
#else
#define LONG_GET(target, source, num_elems, pe)		shmem_long_get(target, source, num_elems, pe)
#endif

#ifdef USE_MPI3
#define GETMEM(target, source, length, rank)		MPI_Get(target, length, MPI_BYTE, rank, (void *)source - window_base, length, MPI_BYTE, window)
#else
#define GETMEM(target, source, length, pe)		shmem_getmem(target, source, length, pe)
#endif

#ifdef USE_MPI3
#define SHORT_PUT(target, source, num_elems, rank)	MPI_Put(source, num_elems, MPI_SHORT, rank, (void *)target - window_base, num_elems, MPI_SHORT, window)
#else
#define SHORT_PUT(target, source, num_elems, pe)	shmem_short_put(target, source, num_elems, pe)
#endif

#ifdef USE_MPI3
#define QUIET()		MPI_Win_flush_all(window)
#else
#define QUIET()		shmem_quiet()
#endif

#ifdef USE_MPI3
#define BARRIER_ALL()	QUIET(); MPI_Barrier(MPI_COMM_WORLD)
#else
#define BARRIER_ALL()	shmem_barrier_all()
#endif

#ifdef USE_MPI3
static int malloc_all(size_t size, void **address) {
  *address = next_window_address;
  next_window_address += size;
  MPI_Barrier(MPI_COMM_WORLD);
  if (next_window_address - window_base > window_size) {
    printf("ran out of memory!\n");
    return -1;
  } else
    return 0;
}
#else
static int malloc_all(size_t size, void **address) {
  *address = shmalloc(size);
  if (*address == NULL)
    return -1;
  else
    return 0;
}
#endif

#ifdef USE_MPI3
#define FREE_ALL(address) /* unable to free memory like this */
#else
#define FREE_ALL(address) shfree(address)
#endif

static inline int global_index_to_rank(const seq_t *in, const index_t codon_index){
  return codon_index / in->local_size;
}

static inline int global_index_to_local_index(const seq_t *in, const index_t codon_index){
  return codon_index % in->local_size;
}

static inline void fetch_from_seq(const seq_t *in, index_t const codon_index, codon_t *out){
  int target_ep = global_index_to_rank(in,codon_index);
  int local_index = global_index_to_local_index(in,codon_index);
  short *typed_seq = (short *)in->sequence;
  SHORT_GET((short *)out, &(typed_seq[local_index]), 1, target_ep);
}

static inline void fetch_from_seq_nb(const seq_t *in, index_t const codon_index, codon_t *out){
  int target_ep = global_index_to_rank(in,codon_index);
  int local_index = global_index_to_local_index(in,codon_index);
  short *typed_seq = (short *)in->sequence;
  SHORT_GET_NB((short *)out, &(typed_seq[local_index]), 1, target_ep);
}

static inline void write_to_seq(const seq_t *in, const index_t codon_index, codon_t data){
  int target_ep = global_index_to_rank(in,codon_index);
  int local_index = global_index_to_local_index(in,codon_index);
  short *typed_seq = (short *)in->sequence;
  short typed_data = (short)data;
  SHORT_PUT(&(typed_seq[local_index]), &typed_data, 1, target_ep);
}

#ifdef USE_MPI3
#define WAIT_NB() /* no non-blocking available yet? */
#else
#define WAIT_NB() QUIET()
#endif

void distribute_rng_seed(unsigned int new_seed);
void seed_rng(int adjustment);
void touch_memory(void *mem, index_t size);
index_t scrub_hyphens(good_match_t *A, seq_t *dest, seq_t *source, index_t length);
void assemble_acid_chain(good_match_t *A, char *result, seq_t *chain, index_t length);
void assemble_codon_chain(good_match_t *A, char *result, seq_t *chain, index_t length);
score_t simple_score(good_match_t *A, seq_t *main, seq_t *match);
seq_t *alloc_global_seq(index_t seq_size);
seq_t *alloc_local_seq(index_t seq_size);
void free_global_seq(seq_t *doomed);
void free_local_seq(seq_t *doomed);

#endif

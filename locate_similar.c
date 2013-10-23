#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gen_scal_data.h>
#include <gen_sim_matrix.h>
#include <pairwise_align.h>
#include <scan_backwards.h>
#include <sort.h>
#include <util.h>

int verify_similar(good_match_t *A, good_match_t *S[], int length_s, int maxDisplay)
{
  int display_count = maxDisplay < length_s ? maxDisplay : length_s;
  int score, retval=0;

  printf("\nFound %i acceptable sets of alignments.\n", length_s);
  if(display_count > 0)
  {
    if(display_count < length_s)
      printf("Displaying the first %i of them.\n", display_count);
    printf("\nStarting   Amino     Codon           Ending");
    printf("\nposition   acids     bases           position\n");
  }
  for(int b=0; b < length_s; b++)
  {
    for(int m=0; m < S[b]->bestLength; m++)
    {
      score = simple_score(A, S[b]->bestSeqs[m].main, S[b]->bestSeqs[m].match, S[b]->bestSeqs[m].length);
      if(score != S[b]->bestScores[m])
      {
        retval = 1;
        printf("\nverifySimilar %i/%i failed; reported %i vs actual %i:   ---------------------------\n",
               b, m, S[b]->bestScores[m], score);
      }
      else if(display_count > b)
      {
        printf("\nverifySimilar %i/%i, succeeded; score %i:\n",
               b, m, S[b]->bestScores[m]);
      }
      if(display_count > b)
      {
        char main_acid_chain[S[b]->bestSeqs[m].length+1];
        char match_acid_chain[S[b]->bestSeqs[m].length+1];
        char main_codon_chain[S[b]->bestSeqs[m].length*3+1];
        char match_codon_chain[S[b]->bestSeqs[m].length*3+1];

        memset(main_acid_chain, '\0', S[b]->bestSeqs[m].length+1);
        memset(match_acid_chain, '\0', S[b]->bestSeqs[m].length+1);
        memset(main_codon_chain, '\0', S[b]->bestSeqs[m].length*3+1);
        memset(match_codon_chain, '\0', S[b]->bestSeqs[m].length*3+1);

        assemble_acid_chain(A, main_acid_chain, S[b]->bestSeqs[m].main, S[b]->bestSeqs[m].length);
        assemble_acid_chain(A, match_acid_chain, S[b]->bestSeqs[m].match, S[b]->bestSeqs[m].length);

        //assemble_codon_chain(S[b]->bestSeqs[m].length, S[b]->bestSeqs[m].main, main_codon_chain, A);
        //assemble_codon_chain(S[b]->bestSeqs[m].length, S[b]->bestSeqs[m].match, match_codon_chain, A);
        assemble_codon_chain(A, main_codon_chain, S[b]->bestSeqs[m].main, S[b]->bestSeqs[m].length);
        assemble_codon_chain(A, match_codon_chain, S[b]->bestSeqs[m].match, S[b]->bestSeqs[m].length);

        printf("%7d  %s  %s  %7d\n%7d  %s  %s  %7d\n",
               S[b]->bestStarts[0][m], main_acid_chain, main_codon_chain, S[b]->bestEnds[0][m],
               S[b]->bestStarts[1][m], match_acid_chain, match_codon_chain, S[b]->bestEnds[1][m]);

      }
    }
  }
  return retval;
}

/*
 * Generate a typical matching pair.
 *
 * With a small change this recursive approach can be used to generate all
 * possible matches.  We only need to identify a single match.  C implementation
 * is limited in recursion only by memory.
 * slightly different from the one found in scan_backwards.c.  Maybe consider trying
 * to merge them?
 */

void tracepath_s(good_match_t *A, int maxResult, int sequence_length, short T[maxResult][sequence_length], int *mainSeq, int *matchSeq, int ci, int cj, int dir, int depth, seq_t *sequence, int rs[2])
{
  int Ci, Cj, iT;
  if(ci== -1 || cj == -1 || depth >= maxResult)  // if done
  {
    rs[0] = ci+1;
    rs[1] = cj+1;
    sequence->main = malloc(depth * sizeof(int));
    sequence->match = malloc(depth * sizeof(int));
    sequence->length = depth;
    return;
  }

  iT = ci % maxResult; // the traceback table wraps
  // bitand(dir, T(iT,cj)
  if(dir==0 || dir & T[iT][cj])
  {
    dir = T[iT][cj] >> 2;
  }

  if(dir==0)
  {
    rs[0] = -1;
    rs[1] = -1;
    return;
  }

  Ci = mainSeq[ci];
  Cj = matchSeq[cj];

  if(dir & 4) // match
  {
    tracepath_s(A, maxResult, sequence_length, T, mainSeq, matchSeq, ci-1, cj-1, 0, depth+1, sequence, rs);
    if(rs[0] != -1 && rs[1] != -1)
    {
      sequence->main[sequence->length-depth-1] = Ci;
      sequence->match[sequence->length-depth-1] = Cj;
      return;
    }
  }

  if(dir & 2) // down
  {
    tracepath_s(A, maxResult, sequence_length, T, mainSeq, matchSeq, ci-1, cj, 2, depth+1, sequence, rs);
    if(rs[0] != -1 && rs[1] != -1)
    {
      sequence->main[sequence->length-depth-1] = Ci;
      sequence->match[sequence->length-depth-1] = A->simMatrix->hyphen;
      return;
    }
  }

  if(dir & 1) // right
  {
    tracepath_s(A, maxResult, sequence_length, T, mainSeq, matchSeq, ci, cj-1, 1, depth+1, sequence, rs);
    if(rs[0] != -1 && rs[1] != -1)
    {
      sequence->main[sequence->length-depth-1] = A->simMatrix->hyphen;
      sequence->match[sequence->length-depth-1] = Cj;
      return;
    }
  }
  printf("Warning, reache the end of tracepath_s\n");
}

/*
 * Consider adding this end-point pair, and possibly deleting nearby end points.
 * This one is slightly different from the version found in pairwise align.  Maybe
 * try and merge them?  Looks non-trivial to do.
 */

void considerAdding_s(good_match_t *A, int *report, int sortReports, int maxReports, int *bestStarts[2], int *bestEnds[2], seq_t *bestSeqs, int *bestScores, int *minScore, int minSeparation, int *V, int i, int j, int maxResult, int sequence_length, short T[maxResult][sequence_length], int *mainSeq, int *matchSeq)
{
  int elements_to_copy;
  int rs[2];
  seq_t consider_seq;

  for(int r=(*report)-1; r >=0; r--)
  {
    if(i-bestEnds[0][r] >= minSeparation) break;
    if(bestScores[r] > V[j]) return;

    // discard point r
    elements_to_copy = (*report);

    free(bestSeqs[r].main);
    free(bestSeqs[r].match);

    for(int idx = r; idx < elements_to_copy; idx++)
    {
      bestScores[idx]=bestScores[idx+1];
      bestEnds[0][idx]=bestEnds[0][idx+1];
      bestEnds[1][idx]=bestEnds[1][idx+1];
      bestStarts[0][idx]=bestStarts[0][idx+1];
      bestStarts[0][idx]=bestStarts[1][idx+1];
      bestSeqs[idx].main = bestSeqs[idx+1].main;
      bestSeqs[idx].match = bestSeqs[idx+1].match;
      bestSeqs[idx].length = bestSeqs[idx+1].length;
    }

    bestSeqs[*report].main = NULL;
    bestSeqs[*report].match = NULL;
    (*report)--;
  }

  tracepath_s(A, maxResult, sequence_length, T, mainSeq, matchSeq, i, j, 0, 0, &consider_seq, rs);
  if(rs[1] != 0)
  {
    free(consider_seq.main);
    free(consider_seq.match);
    return; // drop the sequence if our result was truncated
  }

  for(int r=(*report)-1; r >=0; r--)
  {
    if(rs[0]-bestStarts[0][r] >= minSeparation) break;
    if(bestScores[r] > V[j])
    {
      free(consider_seq.main);
      free(consider_seq.match);
      return;
    }

    // discard point r
    elements_to_copy = (*report);

    free(bestSeqs[r].main);
    free(bestSeqs[r].match);

    for(int idx = r; idx < elements_to_copy; idx++)
    {
      bestScores[idx]=bestScores[idx+1];
      bestEnds[0][idx]=bestEnds[0][idx+1];
      bestEnds[1][idx]=bestEnds[1][idx+1];
      bestStarts[0][idx]=bestStarts[0][idx+1];
      bestStarts[0][idx]=bestStarts[1][idx+1];
      bestSeqs[idx].main = bestSeqs[idx+1].main;
      bestSeqs[idx].match = bestSeqs[idx+1].match;
      bestSeqs[idx].length = bestSeqs[idx+1].length;
    }
    bestSeqs[*report].main = NULL;
    bestSeqs[*report].match = NULL;
    (*report)--;
  }
  bestScores[*report] = V[j];
  bestStarts[0][*report] = rs[0];
  bestStarts[1][*report] = rs[1];
  bestEnds[0][*report] = i;
  bestEnds[1][*report] = j;
  bestSeqs[*report].main = consider_seq.main;
  bestSeqs[*report].match = consider_seq.match;
  bestSeqs[*report].length = consider_seq.length;
  (*report)++;

  if(*report == sortReports)
  {
    int index_array[sortReports];
    int sorted_array[sortReports];
    int best_index[sortReports];
    int worst_keeper;
    int new_best_index=0;

    sorted_array[0]=0;

    memcpy(sorted_array, bestScores, sizeof(int)*sortReports);
    index_sort(sorted_array, index_array, *report);

    worst_keeper = sortReports - maxReports;
    *minScore = sorted_array[worst_keeper] + 1;
    for(int index_for_index=worst_keeper; index_for_index < sortReports; index_for_index++)
    {
      best_index[new_best_index] = index_array[index_for_index];
      new_best_index++;
    }

    sort(best_index, new_best_index);
    for(int idx=0; idx < new_best_index; idx++)
    {
      bestScores[idx] =bestScores[best_index[idx]];
      bestEnds[0][idx]=bestEnds[0][best_index[idx]];
      bestEnds[1][idx]=bestEnds[1][best_index[idx]];
      bestStarts[0][idx]=bestStarts[0][best_index[idx]];
      bestStarts[1][idx]=bestStarts[1][best_index[idx]];
      bestSeqs[idx].main = bestSeqs[best_index[idx]].main;
      bestSeqs[idx].match = bestSeqs[best_index[idx]].match;
      bestSeqs[idx].length = bestSeqs[best_index[idx]].length;
    }
    for(int idx=new_best_index; idx < sortReports; idx++)
    {
      free(bestSeqs[best_index[idx]].main);
      free(bestSeqs[best_index[idx]].match);
    }
    *report = maxReports;
  }
}

/*
 * Scan the main sequence, locating promising alignments.
 * Report the score, start-pair, end-pair, and alignment for the best of these.
 */

void locateSeq(good_match_t *A, int report_number, int minScore, int maxReports, int minSeparation, int maxMatch, good_match_t *S)
{
  int report = 0;
  int sortReports = maxReports * 3;
  int gapFirst = A->simMatrix->gapExtend + A->simMatrix->gapStart;
  // allocating from the heap, because large values of sortReports causes the stack to explode
  int *working_scores = malloc(sortReports*sizeof(int));
  int *working_starts[2];
  int *working_ends[2];
  seq_t *working_seqs = malloc(sortReports*sizeof(seq_t));
  working_starts[0] = malloc(sortReports*sizeof(int));
  working_starts[1] = malloc(sortReports*sizeof(int));
  working_ends[0] = malloc(sortReports*sizeof(int));
  working_ends[1] = malloc(sortReports*sizeof(int));
  memset(working_scores, 0, sortReports*sizeof(int));
  memset(working_starts[0], 0, sortReports*sizeof(int));
  memset(working_ends[0], 0, sortReports*sizeof(int));
  memset(working_starts[1], 0, sortReports*sizeof(int));
  memset(working_ends[1], 0, sortReports*sizeof(int));
  memset(working_seqs, 0, sortReports*sizeof(seq_t));
  int matchSeq[A->bestSeqs[report_number].length];
  int *mainSeq = A->seqData->main;

  // sequence_length is the same as m
  int sequence_length = A->bestSeqs[report_number].length;
  int maxResult = sequence_length * maxMatch;
  short T[maxResult][sequence_length];
  memset(T, 0, maxResult*sequence_length*sizeof(short));

  int *V = malloc(sizeof(int)*sequence_length);
  int *F = malloc(sizeof(int)*sequence_length);
  int *G = malloc(sizeof(int)*sequence_length);
  int *E = malloc(sizeof(int)*sequence_length);
  int *Vg = malloc(sizeof(int)*sequence_length);
  int adder, iT;
  int worst;
  int *sorted_array;// = malloc(sizeof(int)*maxReports);
  int *index_array;// = malloc(sizeof(int)*maxReports);

  // previous kernels will insert hyphens into sequences to mark where a codon was inserted/deleted.  Clean those hyphens out
  int tidy_length = scrub_hyphens(A, matchSeq, A->bestSeqs[report_number].match, A->bestSeqs[report_number].length);

  for(int idx=0; idx < tidy_length; idx++)
  {
    V[idx] = -gapFirst - A->simMatrix->gapExtend*idx;
    F[idx] = V[idx];
  }

  /*
   *  Loop over each codon in the mainSeq sequence, matching it with each codon in
   *  the matchSeq subsequence, using a hybrid local/global affine-gap version of
   *  Smith-Waterman.
   */

  for(int idx=0; idx < A->seqData->mainLen; idx++)
  {
    for(int jdx=0; jdx < tidy_length; jdx++)
    {
      adder = (jdx == 0) ? 0 : V[jdx-1];      
      G[jdx] = A->simMatrix->similarity[mainSeq[idx]][matchSeq[jdx]] + adder;
    }

    for(int jdx=0; jdx < tidy_length; jdx++) 
    {
      V[jdx] = (F[jdx] > G[jdx]) ? F[jdx] : G[jdx];
    }

    E[0] = (-gapFirst > V[0] - gapFirst) ? -gapFirst : V[0] - gapFirst;

    for(int jdx=1; jdx < tidy_length; jdx++)
    {
      E[jdx] = (E[jdx-1] - A->simMatrix->gapExtend > V[jdx-1] - gapFirst) ? E[jdx-1] - A->simMatrix->gapExtend : V[jdx-1] - gapFirst;
      V[jdx] = (E[jdx] > V[jdx] ? E[jdx] : V[jdx]);
    }

    iT = idx % maxResult;

    for(int jdx=0; jdx < tidy_length; jdx++)
    {
      T[iT][jdx] = 4*(V[jdx] == E[jdx]) + 8*(V[jdx] == F[jdx]) + 16*(V[jdx] == G[jdx]);
      Vg[jdx] = V[jdx] - gapFirst;
      F[jdx] = (F[jdx] - A->simMatrix->gapExtend > Vg[jdx]) ? F[jdx]-A->simMatrix->gapExtend : Vg[jdx];
      adder = ((jdx+1 != tidy_length) ? E[jdx+1] : 0) == Vg[jdx];
      T[iT][jdx] = T[iT][jdx] + adder + 2 * (F[jdx] == Vg[jdx]);
    }

    // debug matrix dump
    /*
    if(idx==0)
    {
      printf("\nmaxResult=%i\n----------\n", maxResult);
      for(int adx=0; adx < maxResult; adx++)
      {
        for(int bdx=0; bdx < sequence_length; bdx++)
        printf("%hi ", T[adx][bdx]);
        printf("\n");
      }
    }
    */

    // if the match is good enough- consider adding it
    //printf("[%i]Thinking of adding, report=%i score=%i\n", report_number, report, V[tidy_length-1]);
    if(V[tidy_length-1] >= minScore)
    {
      //printf("[%i]Adding, report=%i, score=%i\n", report_number, report, V[tidy_length-1]);
      considerAdding_s(A, &report, sortReports, maxReports, working_starts, working_ends, working_seqs, working_scores, 
                       &minScore, minSeparation, V, idx, tidy_length-1, maxResult, sequence_length, T, mainSeq, matchSeq);
      //printf("[%i]Finished adding, report is now=%i\n", report_number, report);
    }
    //printf("[%i]Done thinking\n", report_number);
  }

  // once everything is done, add the sequences to the S array.

  index_array = malloc(sizeof(int)*report);
  sorted_array = malloc(sizeof(int)*report);

  if(report!=0)
  {
    memcpy(sorted_array, working_scores, report * sizeof(int));
    index_sort(sorted_array, index_array, report);
    worst = (report - maxReports) > 0 ? report - maxReports : 0;
    S->bestScores = malloc(report * sizeof(int));
    S->bestStarts[0] = malloc(report * sizeof(int));
    S->bestStarts[1] = malloc(report * sizeof(int));
    S->bestEnds[0] = malloc(report * sizeof(int));
    S->bestEnds[1] = malloc(report * sizeof(int));
    S->bestSeqs = malloc(report * sizeof(seq_t));

    for(int jdx=report-1; jdx >= 0; jdx--)
    {
      S->bestScores[report-jdx-1]=working_scores[index_array[jdx]];
      S->bestStarts[0][report-jdx-1]=working_starts[0][index_array[jdx]];
      S->bestStarts[1][report-jdx-1]=working_starts[1][index_array[jdx]];
      S->bestEnds[0][report-jdx-1]=working_ends[0][index_array[jdx]];
      S->bestEnds[1][report-jdx-1]=working_ends[1][index_array[jdx]];
      S->bestSeqs[report-jdx-1].main = working_seqs[index_array[jdx]].main;
      S->bestSeqs[report-jdx-1].match = working_seqs[index_array[jdx]].match;
      S->bestSeqs[report-jdx-1].length = working_seqs[index_array[jdx]].length;
    }
    S->bestLength=report;
  }

  free(working_seqs);
  free(working_ends[0]);
  free(working_starts[0]);
  free(working_ends[1]);
  free(working_starts[1]);
  free(working_scores);
  free(V);
  free(F);
  free(G);
  free(E);
  free(Vg);
  free(sorted_array);
  free(index_array);
}

/*
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Function locateSimilar() - Kernel 3 - Locate sets of similar sequences.
 * This function requires Matlab 7 due to its use of subfunctions.
 *
 * This function uses a variant of the Smith-Waterman dynamic programming
 * algorithm to match all subsequences of seqData.main against each of the
 * exact subsequences from seqData.match reported by Kernel 2.  This variant is
 * a sort of hybrid of local and global alignment -- it is local in main, but
 * global in match.  The match must involve the entire match subsequence.
 *
 * This kernel combines the functions of Kernel 1 and Kernel 2 in one pass by
 * building the trace table T during the forward pass.  Once a promising
 * end-point pair has been found, T is used to find the start-points and
 * generate an alignment.  For each of the A.bestSeqs input sequences, the best
 * of these are reported in the S output structure.
 *
 * This is possible since it is reasonable to limit the length of an alignment
 * in main via the maxMatch parameter, and the length of the match sequence is
 * known, so we can afford to generate the trace table T as we scan the main
 * sequence.  T is organized as a rectangular array operating as a ring buffer
 * of rows.  As previously, only the best match from a cluster of closely
 * related matches is reported.
 *
 * For a detailed description of the SSCA #1 Optimal Pattern Matching problem,
 * please see the SSCA #1 Written Specification.
 *
 * INPUT
 * good_match_t *A   - [structure pointer] structure containing the results of kernel 1 and kernel 2
 * int minScore      - [integer] minimum endpoint score
 * int maxReports    - [integer] maximum number of endpoints reported
 * int minSeparation - [integer] minimum endpoint separation in codons
 * int maxMatch      - [integer] factor which limits the match length
 *
 * OUTPUT
 * S                - [structure pointer array] similar sequence sets for each bestSeq
 *   bestLength     - [integer] The length for the array length (M below)
 *   bestScores     - [integer[M]] the scores for the bestSeqs
 *   bestStarts     - [integer[2][M]] main/match startpoints
 *   bestEnds       - [integer[2][M]] main/match endpoints
 *   bestSeqs       - [seq_t[M]] main/match sequences
 *      length      - [integer] length for the sequences (N below)
 *      main        - [integer[N]] array to the main sequence
 *      match       - [integer[N]] array to the match sequence
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */

void locateSimilar(good_match_t *A, int minScore, int maxReports, int minSeparation, int maxMatch, good_match_t *S[])
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int idx=0; idx < A->bestLength; idx++)
  {
    S[idx] = malloc(sizeof(good_match_t));
    memset(S[idx], 0, sizeof(good_match_t));
    locateSeq(A, idx, minScore, maxReports, minSeparation, maxMatch, S[idx]);
  }
}

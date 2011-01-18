#ifndef __LOCATE_SIMILAR_H
#define __LOCATE_SIMILAR_H

void locateSimilar(good_match_t *A, int minScore, int maxReports, int minSeparation, int maxMatch, good_match_t *S[]);
// throwing the debug function prototype too for good measure
good_match_t **test_similar_results();
int verify_similar(good_match_t *A, good_match_t *S[], int length_s, int maxDisplay);

#endif

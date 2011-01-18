#include <sort.h>
#include <stdlib.h>

/* function to be passed to to glibc's qsort function */

int sort_cmp(const void *first, const void *second)
{
  int *number_a = (int *)first;
  int *number_b = (int *)second;
  return *number_a - *number_b;
}

/* A wrapper around glibc's qsort.
     Input-
         int *numbers   - values to be sorted
         int array_size - size of the numbers array
     Output-
         int *numbers   - sorted values
*/

void sort(int numbers[], int array_size)
{
  qsort(numbers, array_size, sizeof(int), sort_cmp);
}

/* Sort the values and return a list of indexes of the sorted values
     Input-
         int *numbers   - values to be sorted
         int *indexes   - an int array size of array_size
         int array_size - size of the numbers array
     Output-
         int *numbers   - sorted values
         int *indexes   - the index where the value used to be before sorting
*/

void index_sort(int numbers[], int indexes[], int array_size)
{
  int big_index[array_size][2];
  for(int idx=0; idx < array_size; idx++)
  {
    big_index[idx][0] = numbers[idx];
    big_index[idx][1] = idx;
  }

  qsort(big_index, array_size, sizeof(int)*2, sort_cmp);

  for(int idx=0; idx < array_size; idx++)
  {
    numbers[idx] = big_index[idx][0];
    indexes[idx] = big_index[idx][1];
  }
}

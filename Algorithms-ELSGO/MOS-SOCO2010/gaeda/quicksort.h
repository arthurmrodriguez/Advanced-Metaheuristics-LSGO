#ifndef _QUICKSORT_H_
#define _QUICKSORT_H_

#ifdef __cplusplus
extern "C" {
#endif

void swap (int* keys, long double* values, int p1, int p2);

//This function orders the vector from the 0 to the IND_SIZE-1 elements
//Way to call it: quicksort(keys, values, 0, IND_SIZE-1);
void quicksort (int* keys, long double* values, int l, int r);

#ifdef __cplusplus
}
#endif

#endif

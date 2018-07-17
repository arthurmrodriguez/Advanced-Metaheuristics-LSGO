#include "quicksort.h"

void swap (int* keys, long double* values, int p1, int p2) {

   long double dDummy = values [p1];
   int    nDummy = keys   [p1];

   values [p1] = values [p2];
   values [p2] = dDummy;

   keys [p1] = keys [p2];
   keys [p2] = nDummy;

}

//This function orders the vector from the 0 to the IND_SIZE-1 elements
//Way to call it: quicksort(keys, values, 0, IND_SIZE-1);
void quicksort (int* keys, long double* values, int l, int r) {

   int i, j;
   long double v;

   if (r > l) {

      v = values [r];
      i = l - 1;
      j = r;

      for (;;) {

         while (values [++i] < v);
         while (values [--j] > v);

         if (i >= j)
            break;

         swap (keys, values, i, j);

      } // for

      swap (keys, values, i, r);

      quicksort (keys, values, l, i - 1);
      quicksort (keys, values, i + 1, r);

   } // if

}

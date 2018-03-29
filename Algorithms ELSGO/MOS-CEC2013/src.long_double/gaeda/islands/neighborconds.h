#ifndef NEIGHBORCONDS_H_
#define NEIGHBORCONDS_H_

#include "../genomes/GAGenome.h"
#include <vector>

using namespace std;

typedef long double (*neighbor_cond_func)(unsigned int, unsigned int, unsigned int, long double* dist_matrix, int num_rows); 

long double gabrielNeighborCond  (unsigned int pos_a, unsigned int pos_b, unsigned int pos_x, long double* dist_matrix, int num_rows);

long double relativeNeighborCond (unsigned int pos_a, unsigned int pos_b, unsigned int pos_x, long double* dist_matrix, int num_rows);

long double squareNeighborCond   (unsigned int pos_a, unsigned int pos_b, unsigned int pos_x, long double* dist_matrix, int num_rows);


#endif /*NEIGHBORCONDS_H_*/

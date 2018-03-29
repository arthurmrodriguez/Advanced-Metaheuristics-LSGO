#include "neighborconds.h"
#include <iostream>
#include <math.h>

using namespace std;

long double gabrielNeighborCond(unsigned int pos_a, unsigned int pos_b, unsigned int pos_x, long double* dist_matrix, int num_rows){
  long double dist_ab = dist_matrix[pos_a*num_rows+pos_b];
  long double dist_bx = dist_matrix[pos_b*num_rows+pos_x];
  long double dist_ax = dist_matrix[pos_a*num_rows+pos_x];
  return pow(dist_ab,2.0) < ( pow(dist_ax,2.0) + pow(dist_bx,2.0) );
}

long double relativeNeighborCond(unsigned int pos_a, unsigned int pos_b, unsigned int pos_x, long double* dist_matrix, int num_rows){
  long double dist_ab = dist_matrix[pos_a*num_rows+pos_b];
  long double dist_bx = dist_matrix[pos_b*num_rows+pos_x];
  long double dist_ax = dist_matrix[pos_a*num_rows+pos_x];
  long double max_dist_x = (dist_bx > dist_ax) ? dist_bx : dist_ax;
  return dist_ab < max_dist_x;
}

long double squareNeighborCond(unsigned int pos_a, unsigned int pos_b, unsigned int pos_x, long double* dist_matrix, int num_rows){
  //to be  completed
  return 0.0;
}

#include "neighborconds.h"
#include <iostream>
#include <math.h>

using namespace std;

double gabrielNeighborCond(unsigned int pos_a, unsigned int pos_b, unsigned int pos_x, double* dist_matrix, int num_rows){
  double dist_ab = dist_matrix[pos_a*num_rows+pos_b];
  double dist_bx = dist_matrix[pos_b*num_rows+pos_x];
  double dist_ax = dist_matrix[pos_a*num_rows+pos_x];
  return pow(dist_ab,2.0) < ( pow(dist_ax,2.0) + pow(dist_bx,2.0) );
}

double relativeNeighborCond(unsigned int pos_a, unsigned int pos_b, unsigned int pos_x, double* dist_matrix, int num_rows){
  double dist_ab = dist_matrix[pos_a*num_rows+pos_b];
  double dist_bx = dist_matrix[pos_b*num_rows+pos_x];
  double dist_ax = dist_matrix[pos_a*num_rows+pos_x];
  double max_dist_x = (dist_bx > dist_ax) ? dist_bx : dist_ax;
  return dist_ab < max_dist_x;
}

double squareNeighborCond(unsigned int pos_a, unsigned int pos_b, unsigned int pos_x, double* dist_matrix, int num_rows){
  //to be  completed
  return 0.0;
}

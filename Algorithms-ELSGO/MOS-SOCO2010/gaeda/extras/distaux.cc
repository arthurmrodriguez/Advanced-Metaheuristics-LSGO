#include "distaux.h"
#include <limits>
#include <math.h>

vector<long double> computeDistanceMatrix(vector<GAGenome*>& inds) {
  int num_inds = inds.size();

  vector<long double> dist_matrix( num_inds * num_inds, 0.0 );

  for (int i=0; i<num_inds; i++){
    dist_matrix[i*num_inds+i] = 0;

    for (int j=i+1; j<num_inds; j++){
      GAGenome& ind_i = *(inds[i]);
      GAGenome& ind_j = *(inds[j]);
      dist_matrix[i*num_inds+j] = dist_matrix[j*num_inds+i] = ind_i.compare(ind_j);
    }
  }
  return dist_matrix;
}

struct MinDist{
  int    other_pos;
  long double value;
};

MinDist minDistOfPos(int pos, vector<long double>& dist_matrix, vector<bool>* deleted_inds=NULL){
  if (deleted_inds != NULL) assert(deleted_inds->size() * deleted_inds->size() == dist_matrix.size() );

  int num_inds = (int) sqrt(dist_matrix.size());

  MinDist mindist;
  long double  tmp_dist;

  mindist.value = std::numeric_limits<long double>::max();

  for (int i=0; i<num_inds; i++) {
    if (pos == i || ( deleted_inds != NULL && (*deleted_inds)[i] ) ) continue;

    tmp_dist = dist_matrix[pos*num_inds+i];

    if (tmp_dist < mindist.value) {
      mindist.value     = tmp_dist;
      mindist.other_pos = i;
    }
  }

  return mindist;
}

vector<MinDist> computeMinDists(vector<GAGenome*>& inds, vector<long double>& dist_matrix) {
  int             num_inds = inds.size();
  vector<MinDist> mindists;

  for (int i=0; i<num_inds; i++) {
    mindists.push_back(minDistOfPos(i,dist_matrix));
  }
  return mindists;
}

int minDistPos(vector<MinDist>& mindists, vector<bool>& deleted_inds){
  assert(mindists.size() == deleted_inds.size());

  int pos = -1;
  MinDist mindist; mindist.value = std::numeric_limits<int>::max();

  for (unsigned i=0; i<mindists.size(); i++) {
    if (deleted_inds[i]) continue;

    if (mindists[i].value < mindist.value) {
      pos     = i;
      mindist = mindists[i];
    }
  }

  return pos;
}

void removeInd(int pos, vector<MinDist>& mindists, vector<long double>& dist_matrix, vector<bool>& deleted_inds) {
  assert(mindists.size() == deleted_inds.size());

  deleted_inds[pos] = true;

  int num_inds = mindists.size();
  for (int i=0; i<num_inds; i++) {
    if (i==pos || deleted_inds[i]) continue;

    if (mindists[i].other_pos == pos) {
      mindists[i] = minDistOfPos(i,dist_matrix,&deleted_inds);
    }
  }
}

/*
 * Computes the D2 Method from Glover, F., C.C., Kuo and K.S. Dhir, 1998 Heuristic
 * algorithms for the maximum diversity problem. Journal of Information and Optimization
 * Sciences. This is a different version from the one of the paper. Here, the function removes at each iteration
 * the point which has the minimum distance to the rest of the points. In the original version the the value used was the sum of the distances from the point
 * to the other points.
 */
void maxminD2Method (vector<GAGenome*>& all_inds, vector<GAGenome*>& slctd_inds) {
  assert(slctd_inds.size() <= all_inds.size());

  int num_inds         = all_inds.size();
  int num_inds_2delete = num_inds - slctd_inds.size();

  vector<bool>    deleted_inds(num_inds, false);

  if (slctd_inds.size() != all_inds.size() ) {

    vector<long double>  dist_matrix = computeDistanceMatrix(all_inds);
    vector<MinDist> mindists    = computeMinDists(all_inds, dist_matrix);

    for (int i=0; i<num_inds_2delete; i++) {                // Repeat num_inds_2delete times

     int mindist_pos = minDistPos(mindists,deleted_inds);

     removeInd(mindist_pos,mindists,dist_matrix,deleted_inds);
    }
  }

  unsigned slctd_pos=0;
  for (unsigned i=0; i<all_inds.size(); i++) {
    if (deleted_inds[i] == true) continue;

    assert(slctd_pos < slctd_inds.size());

    slctd_inds[slctd_pos]->copy(* all_inds[i]);
    slctd_pos++;
  }
}

long double computeDiversityAvg(vector<GAGenome*>& inds) {
  unsigned size    = inds.size();
  long double   div_sum = 0.0;

  for (unsigned i=0; i<size-1; i++) {
    for (unsigned j=i+1; j<size; j++) {
      div_sum += inds[i]->compare(* inds[j]);
    }
  }

  return div_sum/( size*(size-1)/2 );
}


// OLD: Used to represent the maxmin original problem
//
//vector<long double> computeSumDists(vector<GAGenome*>& inds, vector<long double>& dist_matrix) {
//  int num_inds = inds.size();
//
//  vector<long double> dist_sums (num_inds,0.0);
//  for (int i=0; i<num_inds; i++) {
//    long double sum_i = 0.0;
//    for (int j=0; j<num_inds; j++) {
//      if (i==j) continue;
//      sum_i += dist_matrix[i*num_inds+j];
//    }
//    dist_sums[i] = sum_i;
//  }
//
//  return dist_sums;
//}

//int minIndSumDistPos(vector<GAGenome*>& inds,
//                  vector<long double>&    dist_sums,
//                  vector<bool>&      deleted_inds) {
//
//  int    min_setdist_pos = -1;
//  long double min_setdist     = std::numeric_limits<int>::max();
//  long double tmp_setdist;
//
//  for (unsigned j=0; j<inds.size(); j++) {
//    if (deleted_inds[j]) continue;
//
//    tmp_setdist = dist_sums[j];
//    if (tmp_setdist < min_setdist) {
//      min_setdist     = tmp_setdist;
//      min_setdist_pos = j;
//    }
//  }
//  assert(min_setdist_pos>=0);
//  assert(deleted_inds[min_setdist_pos] == false);
//
//  return(min_setdist_pos);
//}

//void deleteInd(int             ind_pos,
//               vector<long double>& dist_matrix,
//               vector<long double>& dist_sums,
//               vector<bool>&   deleted_inds) {
//
//  assert(dist_sums.size() == deleted_inds.size());
//
//  int num_inds = dist_sums.size();
//
//  for (int i=0; i<num_inds; i++){
//    if (i==ind_pos) continue;
//
//    dist_sums[i] -= dist_matrix[ind_pos*num_inds+i];
//    deleted_inds[ind_pos] = true;
//  }
//}
//
//
//
// First version used to compute the minimum distance
//int minDistPos2(vector<GAGenome*>& inds,
//               vector<long double>&    dist_matrix,
//               vector<bool>&      deleted_inds) {
//  assert(inds.size() == deleted_inds.size() && deleted_inds.size() == (int) sqrt(dist_matrix.size()));
//  int    mindist_pos = -1;
//  long double mindist     = std::numeric_limits<int>::max();
//  int    num_inds     = inds.size();
//  long double tmp_dist;
//
//  for (int i=0; i<num_inds-1; i++) {
//    if (deleted_inds[i]) continue;
//    for (int j=i+1; j<num_inds; j++) {
//      if (deleted_inds[j]) continue;
//
//      tmp_dist = dist_matrix[i*num_inds+j];
//
//      if (tmp_dist < mindist) {
//        mindist     = tmp_dist;
//        mindist_pos = i;
//      }
//    }
//  }
//  assert(mindist_pos>=0);
//  assert(deleted_inds[mindist_pos] == false);
//
//  return(mindist_pos);
//}


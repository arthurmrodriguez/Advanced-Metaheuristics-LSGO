#include "RandomHighSeparatedCentGenerator.h"
#include "../logger/GALogger.h"
#include <math.h>
#include <assert.h>
#include <algorithm>

//Auxiliary functions (dont use class or instance variables)
void deleteCentroidsDistribution(vector<GAGenome*>* centroids);


//Methods
RandomHighSeparatedCentGenerator::RandomHighSeparatedCentGenerator(unsigned int cent_number, GAGenome& sample_gen, GAGenome::Initializer pinit_func, unsigned int pmax_tries) 
                                                  : RandomCentGenerator(cent_number,sample_gen, pinit_func) , max_tries(pmax_tries)  { }


vector<GAGenome*>* RandomHighSeparatedCentGenerator::generateCentroids() const {

  vector<GAGenome*>* centroids     = createARandomCentroidsDistribution();
  double             sum_cents     = sumMinDistancesBetweenCentroids(*centroids);

  for (unsigned int i=0; i<max_tries; i++){
    vector<GAGenome*>* centroids_temp = createARandomCentroidsDistribution();
    double             sum_cents_temp = sumMinDistancesBetweenCentroids(*centroids_temp);

    if (sum_cents_temp > sum_cents){                //We do the swap so centroids_temp at the end point to centroids
      vector<GAGenome*>* centroids_aux = centroids;
      centroids      = centroids_temp;
      sum_cents      = sum_cents_temp;              // We update the sum value
      centroids_temp = centroids_aux;
    }

    deleteCentroidsDistribution(centroids_temp);
  }
 return centroids; 
}

double RandomHighSeparatedCentGenerator::sumMinDistancesBetweenCentroids(vector<GAGenome*>& centroids) const{
  assert(centroids.size() == cent_number);
  vector<double> min_distances(cent_number);
  for (unsigned int i=0; i<cent_number; i++){
    min_distances[i] = findMinDistanceToOtherCentroids(centroids,i);
  }
 
  std::sort( min_distances.begin(), min_distances.end() );
//stringstream message; message << "y la suma de distancias ordenada es 0,1,2 " << min_distances[0] << ", " << min_distances[1] << ", " << min_distances[2] << endl;GAFileTracer::instance()->appendLogMessage("RandomHighSeparatedCentGenerator",message.str());

  double total_sum =0.0;
  unsigned int number_of_distances_to_sum = static_cast<unsigned int> ( round( cent_number * 0.2 ) );
  for (unsigned int i=0; i<number_of_distances_to_sum; i++){
    total_sum += min_distances[i];
  }
  return total_sum;
}
//vector<GAGenome*>* RandomHighSeparatedCentGenerator::generateCentroids() {
//  vector<GAGenome*>* centroids     = createARandomCentroidsDistribution();
//  double             sum_cents     = sumDistancesBetweenCentroids(*centroids);
//
//  for (unsigned int i=0; i<max_tries; i++){
//    vector<GAGenome*>* centroids_temp = createARandomCentroidsDistribution();
//    double             sum_cents_temp = sumDistancesBetweenCentroids(*centroids_temp);
//
//    if (sum_cents_temp > sum_cents){                //We do the swap so centroids_temp at the end point to centroids
//      vector<GAGenome*>* centroids_aux = centroids;
//      centroids      = centroids_temp;
//      sum_cents      = sum_cents_temp;              // We update the sum value
//      centroids_temp = centroids_aux;
//    }
//
//    deleteCentroidsDistribution(centroids_temp);
//  }
// return centroids; 
//}
RandomHighSeparatedCentGenerator::~RandomHighSeparatedCentGenerator(){}

void deleteCentroidsDistribution(vector<GAGenome*>* centroids){
  for (vector<GAGenome*>::iterator it=centroids->begin(); it!=centroids->end(); it++){
    delete *it;
  }
  delete centroids;
}

//double RandomHighSeparatedCentGenerator::sumDistancesBetweenCentroids(vector<GAGenome*>& centroids){
//  double total_sum = 0.0;
//
//  for(unsigned int i=0; i<cent_number; i++){
//    GAGenome& centroid_i = *centroids[i];
//    for (unsigned int j=i+1; j<cent_number; j++){
//      GAGenome& centroid_j = *centroids[j];
//      total_sum            += ( centroid_i.comparator() ) (centroid_i, centroid_j);
//    }
//  }
//  return total_sum;
//}

//double RandomHighSeparatedCentGenerator::sumDistancesBetweenCentroids(vector<GAGenome*>& centroids){
//  double total_sum = 0.0;
//
//  for(unsigned int i=0; i<cent_number; i++){
//    GAGenome& centroid_i = *centroids[i];
//    for (unsigned int j=i+1; j<cent_number; j++){
//      GAGenome& centroid_j = *centroids[j];
//      double distance       = ( centroid_i.comparator() ) (centroid_i, centroid_j);
//      if ( distance == 0) distance = 0.00001;
//      total_sum += (1/ (distance * distance * distance));
//    }
//  }
//  return total_sum;
//}

double RandomHighSeparatedCentGenerator::findMinDistanceToOtherCentroids(vector<GAGenome*>& centroids, unsigned int pos) const{
  GAGenome&    centroid      = *centroids[pos];
  unsigned int first_pos     = (pos == 0) ? 1 : 0;
  double       min_distance  = centroid.compare(*centroids[first_pos]);
  double       temp_distance = 0.0;
  
  for (unsigned int i=first_pos+1; i<centroids.size(); i++){
    if ( i==pos ) continue;
    temp_distance = centroid.compare(*centroids[i]);
    if (temp_distance < min_distance) min_distance = temp_distance;
  }
  return min_distance;
}

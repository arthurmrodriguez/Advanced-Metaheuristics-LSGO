#include "VoronoiIndInit.h"
#include "islandsutils.h"
#include "../genomes/GA1DArrayGenome.h"
#include "../logger/GALogger.h"
#include <iostream>
#include <stdexcept>
#include <list>
#include <assert.h>
#include <sstream>


using namespace std;

// File Static variables definitions
vector<GAGenome*>*    VoronoiIndInit::centroids__   = NULL;
GAGenome*             VoronoiIndInit::my_centroid__ = NULL;
GAGenome::Initializer VoronoiIndInit::init_func__   = NULL;

// Auxiliary functions
static vector<GAGenome*>* generateAndSendCentroidsToIslands(const CentroidsGenerator& cent_gen, const CommManager& comm_manager) {
  cerr << "Generating centroids..." << endl; 
 
  vector<GAGenome*>* centroids = cent_gen.generateCentroids();

  comm_manager.sendIndsToSlaveIslands(*centroids);

  return centroids;
}

static vector<GAGenome*>* getCentroids(const CentroidsGenerator& cent_gen, CommManager& comm_manager){
  vector<GAGenome*>* centroids = NULL;                // In order to isolate the main program, the centroids are not going to be deleted

  if (comm_manager.isIslandMaster() ){
    centroids = generateAndSendCentroidsToIslands(cent_gen,comm_manager);
  }
  else {
    GAPopulation* pop_received = comm_manager.receivePopFromMaster();
    centroids = convertPopToVector(pop_received );
    delete pop_received;
  }
    
  assert ((int)centroids->size() == comm_manager.getNumIslands());
  GALogger::instance()->appendPopulation( "Main", "Centroids generated",  *centroids  );

  return centroids;
}

// Class methods
void VoronoiIndInit::setVoronoiIndInitData(CommManager&              comm_manager,
                                           const CentroidsGenerator& cent_gen,
                                           GAGenome::Initializer     pinit_func){
  centroids__   = getCentroids(cent_gen,comm_manager);
  my_centroid__ = (*centroids__)[comm_manager.getMyRank()];
  init_func__   = pinit_func;
}

void VoronoiIndInit::initialize(GAGenome& gen){  
  assert(my_centroid__!= NULL);
  assert(centroids__->size() > 1);

  bool more_closer_to_centroid = false;
  gen.copy(*my_centroid__);                          // we copy the operators from the centroid
  do {
    init_func__(gen);                                // We initialize the values from the function defined by the user
    more_closer_to_centroid = isCloserToCentroid(gen); // Then we check that the new gene is closer to the centroid than to the rest of the centroids
  } while ( !more_closer_to_centroid );    
}


bool VoronoiIndInit::isCloserToCentroid(GAGenome& gen){
  bool   result             = true;
  double dist_to_mycentroid = gen.compare( *my_centroid__ );
  double dist_to_centroid_i = 0.0;
 
  for (unsigned int i=0; i<centroids__->size(); i++){
    if ( (*centroids__)[i] == my_centroid__) continue;

    GAGenome* centroid_i = (*centroids__)[i];
    dist_to_centroid_i   = gen.compare( *centroid_i );

    if (dist_to_centroid_i < dist_to_mycentroid) {
      result = false;
      break;
    }
  }
  return result;
} 


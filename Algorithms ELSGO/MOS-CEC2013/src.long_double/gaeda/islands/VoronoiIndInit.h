#ifndef VORONOIINDINIT_H_
#define VORONOIINDINIT_H_

#include "EAIslandsTopology.h"
#include "CommManager.h"
#include "CentroidsGenerator.h"
#include "../genomes/GAGenome.h"

#include <vector>

class VoronoiIndInit {
  static vector<GAGenome*>*    centroids__;
  static GAGenome*             my_centroid__;
  static GAGenome::Initializer init_func__;

  static bool isCloserToCentroid(GAGenome& gen);
public:
  static void setVoronoiIndInitData(CommManager&              comm_manager,
                                    const CentroidsGenerator& cent_gen, 
                                    GAGenome::Initializer     pinit_func);
  static void initialize(GAGenome& gen);
};

#endif 

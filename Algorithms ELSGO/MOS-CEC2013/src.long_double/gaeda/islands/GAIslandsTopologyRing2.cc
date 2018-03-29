#include <list>
#include "GAIslandsTopologyRing2.h"
#include <string>

using namespace std;

/**
 * Class constructor
 *
 * @param count Number of islands in the archipelago
 */
GAIslandsTopologyRing2::GAIslandsTopologyRing2( int count, int my_island_pos ) : EAIslandsTopology(count, my_island_pos){
  int j,k;
  for (int i=0; i<count; i++) {
    j = ( i+1         ) % count;
    k = ( i-1 + count ) % count;
    adjmatrix_[i][j] = 1;
    adjmatrix_[i][k] = 1;
  }  
}

string GAIslandsTopologyRing2::getTopologyName (){
  return "Ring Topology";
}



GAIslandsTopologyRing2::~GAIslandsTopologyRing2() {}

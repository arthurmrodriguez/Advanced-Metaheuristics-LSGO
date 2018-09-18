#include "GAIslandsTopologyRing.h"
#include <list>
#include <string>

using namespace std;

GAIslandsTopologyRing::GAIslandsTopologyRing( int count, int my_island_pos ) : 
                                           EAIslandsTopology(count, my_island_pos){
	int j;
	for (int i=0; i < count; i++) {
		j = (i+1) % count;
		adjmatrix_[i][j] = 1;
	}
}

GAIslandsTopologyRing::~GAIslandsTopologyRing() {}

string GAIslandsTopologyRing::getTopologyName (){
  return "Ring2 Topology";
}


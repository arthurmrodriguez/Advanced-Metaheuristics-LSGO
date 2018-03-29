#include <list>
#include "GAIslandsTopologyA2A.h"
#include <string>

using namespace std;

GAIslandsTopologyA2A::GAIslandsTopologyA2A( int count, int my_island_pos ) : EAIslandsTopology(count, my_island_pos){
  for (int i=0; i<count; i++) {
    for (int j=0; j<count; j++){
      if (i!=j) adjmatrix_[i][j] = 1;
    }
  }  
}

GAIslandsTopologyA2A::~GAIslandsTopologyA2A() {}

string GAIslandsTopologyA2A::getTopologyName (){
  return "All 2 all";
}

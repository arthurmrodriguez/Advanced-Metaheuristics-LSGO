#include "EAIslandsTopology.h"
#include <assert.h>

EAIslandsTopology::EAIslandsTopology( int count, int my_island_pos ): island_count_(count), my_island_pos_(my_island_pos) {
	int i,j;

	/* Create the adjacency matrix */
	adjmatrix_ = new short int*[count];
	for (i=0; i < count; i++) {
		adjmatrix_[i] = new short int[count];
	}

	/* Initialize the adjacency matrix */
	for (i=0; i < count; i++) {
		for (j=0; j < count; j++) {
			adjmatrix_[i][j] =0;
		}
	}
}

EAIslandsTopology::~EAIslandsTopology( ) {
	for (int i=0; i < island_count_; i++) delete []adjmatrix_[i];
	
	delete []adjmatrix_;
}

void EAIslandsTopology::addNeighbor(unsigned int pos){
  adjmatrix_[my_island_pos_][pos] = 1;
}

void EAIslandsTopology::removeNeighbor(unsigned int pos){
  adjmatrix_[my_island_pos_][pos] = 0;
}

void EAIslandsTopology::resetAdjMatrix(){
  int nislands = getIslandsCount();
  for (int i=0; i<nislands; i++)
    for (int j=0; j<nislands; j++)
      removeNeighbor(j);
}

int EAIslandsTopology::getIslandsCount( void ) const {
	return island_count_;
}

unsigned int EAIslandsTopology::getNumNeighbors() const {
  unsigned int ndest = 0;
  for (int i=0; i < island_count_; i++) {
    if (adjmatrix_[my_island_pos_][i])  ndest++;
  }
  return ndest;
}

bool EAIslandsTopology::areNeighbors(int i, int j) const{
  return adjmatrix_[i][j] == 1;
}

list<int>* EAIslandsTopology::getDestinations(){
  assert ( (int) my_island_pos_ < island_count_ );
  list<int> *res = new list<int>;

  for (int i=0; i < island_count_; i++) {
    if (adjmatrix_[my_island_pos_][i])
      res->push_back( i );
  }
  return res;
}

list<int>* EAIslandsTopology::getOrigins(){
  assert ( (int) my_island_pos_ < island_count_ );
  list<int> *res = new list<int>;

  for (int i=0; i < island_count_; i++) {
    if (adjmatrix_[i][my_island_pos_])
      res->push_back( i );
  }
  return res;
}

void EAIslandsTopology::generateNewNeighbors (const GAPopulation& pop){}

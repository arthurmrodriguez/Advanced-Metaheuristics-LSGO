#include "GAIslandsTopologyHyperCube.h"
#include "neighborconds.h"
#include "../logger/GALogger.h"
#include <list>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//Auxiliary functions
int                  powint(int base, int exp);
vector<unsigned int> intToBin(unsigned int number, unsigned int order);
unsigned int         binToInt(vector<unsigned int> bin);

GAIslandsTopologyHyperCube::GAIslandsTopologyHyperCube( int          count, 
                                                        int          my_island_pos, 
                                                        unsigned int conn_degree) :
                                                                               EAIslandsTopology(count,my_island_pos ) , 
                                                                               conn_degree_( conn_degree ){
  calculateNeighbors();
  assert ( powint(2,conn_degree_) <= count );
  if ( powint(2,conn_degree_) != count ) cerr << "GAIslandsTopologyHyperCube Warning: The number of islands " << count << " is not a power of 2 of the connection degree" << endl;
}

GAIslandsTopologyHyperCube::~GAIslandsTopologyHyperCube() { }

void GAIslandsTopologyHyperCube::calculateNeighbors(){
  if ( (int) my_island_pos_ < powint(2,conn_degree_) ) {
    vector<unsigned int> binary = intToBin(my_island_pos_, conn_degree_);

    for (unsigned int i=0; i<conn_degree_; i++){
      vector<unsigned int> tmp(binary);
      tmp[i] = ( tmp[i] + 1 ) % 2;
      addNeighbor( binToInt(tmp) );
    }
  }
  else{
    cout << "GAIslandsTopologyHyperCube Warning: The island " << my_island_pos_ << " will not participate in the hypercube of connectivity " << conn_degree_ << " topology " << endl;
  }
}

string GAIslandsTopologyHyperCube::getTopologyName ( ){return "Hypercube";}

vector<unsigned int> intToBin(unsigned int number, unsigned int order){
  assert( (int) number < powint(2,order) );
  vector<unsigned int> result(order);
  for (unsigned int i=0; i<order; i++){
    result[i] = number % 2;
    number = number >> 1;
  }
  return result;
}

unsigned int binToInt(vector<unsigned int> bin){
  unsigned int result = 0;
  for (unsigned int i=0; i<bin.size(); i++)  result += bin[i] * powint(2,i);
  return result;
}

int powint(int base, int exp){
  if      ( exp == 0 ) return 1;
  else if (exp == 1)   return base;
  else                 return base * powint(base,exp-1);
}


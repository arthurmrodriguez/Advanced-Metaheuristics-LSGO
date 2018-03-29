#ifndef GAISLANDSTOPOLOGY_H
#define GAISLANDSTOPOLOGY_H

#include "../GAPopulation.h"
#include <list>

using namespace std;

/**
 * @brief Island topology base class
 *
 * This is an abstract class that defines the interface that all 
 * topology subclasses must follow. Internally all topologies are
 * represented by an adajacency matrix (graph)
 */
class EAIslandsTopology  {
protected:
  short int**  adjmatrix_;   
  int          island_count_; 
  unsigned int my_island_pos_; /* The island pos is stored in this instance variable so this classes (and the derived ones)
                                * dont have to know the MPIIndividualTransmitter class or its equivalent (less dependencies
                                * and less calls to the function getMyRank */

  void addNeighbor    (unsigned int pos );
  void removeNeighbor (unsigned int pos );  
  void resetAdjMatrix (                 );  
  
public:

  EAIslandsTopology( int count, int my_island_pos );
		
  virtual ~EAIslandsTopology();
  
  virtual string getTopologyName()=0;  

	int                getIslandsCount      (                       ) const;
  unsigned int       getNumNeighbors      (                       ) const;
  bool               areNeighbors         (int i, int j           ) const;
	virtual list<int>* getDestinations      (                       );      
	virtual list<int>* getOrigins           (                       );
	virtual void       generateNewNeighbors (const GAPopulation& pop); // Only to be implemented by dynamic topologies
};


#endif

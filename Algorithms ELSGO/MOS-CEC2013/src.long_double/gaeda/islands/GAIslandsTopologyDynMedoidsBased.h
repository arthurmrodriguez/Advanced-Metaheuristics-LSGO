#include "GAIslandsTopologyDisc.h"
#include "neighborconds.h"
#include "CommManager.h"
#include <list>
#include <math.h>
#include <vector>

using namespace std; 

#ifndef GAISLANDSDYNAMICTOPOLOGY__H_
#define GAISLANDSDYNAMICTOPOLOGY__H_

/**
 * @brief Disconnected islands topology
 *
 * No edges on this topology
 */
class GAIslandsTopologyDynMedoidsBased : public EAIslandsTopology {
  CommManager&       comm_manager_;
protected:
  vector<GAGenome*>* medoids_;
  unsigned           num_islands_;
  long double*            dist_matrix_;

  void calculateDistanceMatrix   (vector<GAGenome*>& medoids) ; // This method is called from setNewMedoids and is the one that the child
                                                                                           // class must override

  virtual void regenerateNeihborhoods (vector<GAGenome*>& medoids );
public:
  GAIslandsTopologyDynMedoidsBased( int count, int my_island_pos, CommManager& comm_manager);

	virtual ~GAIslandsTopologyDynMedoidsBased ();	

  void generateNewNeighbors (const GAPopulation& pop);
  void setMedoids (vector<GAGenome*>* medoids );
};

#endif

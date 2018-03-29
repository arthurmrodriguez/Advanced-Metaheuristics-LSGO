#ifndef TSCORDEAU_H
#define TSCORDEAU_H

/*
 * Implementation of the Cordeau's paper
 * Cordeau, J.-F., & Laporte, G. (2003). A tabu search heuristic for the static multi-vehicle dial-a-ride problem.
 * Transportation Research Part B: Methodological, 37(6), 579-594. doi:10.1016/S0191-2615(02)00045-0
 */

#include "RoutingAlg.h"
#include "genomes/RoutingGenome.h"
#include "VNSOp.h"
#include <stdexcept>
#include <sstream>
#include <vector>
#include <list>
#include <algorithm>

using namespace std;

class TSCordeau : public RoutingAlg {
  int lsPeriod_;            // k in Cordeau's paper
  int paramsUpdatePeriod_;

public:

  TSCordeau(const  GAGenome& g, VNSOp *ls ); // takes ownership of passed values
  TSCordeau(const TSCordeau&);

  virtual ~TSCordeau() {}

  void        copy(const GAGeneticAlgorithm& other);
  RoutingAlg* clone();

  TSCordeau& operator=(const TSCordeau&);

  virtual void step();
};

#ifdef GALIB_USE_STREAMS
inline STD_OSTREAM& operator<< (STD_OSTREAM& os, TSCordeau & arg)
{ arg.write(os); return(os); }
inline STD_ISTREAM& operator>> (STD_ISTREAM& is, TSCordeau & arg)
{ arg.read(is); return(is); }
#endif

#endif

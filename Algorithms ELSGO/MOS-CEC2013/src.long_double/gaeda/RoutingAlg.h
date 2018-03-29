#ifndef ROUTINGALG_H
#define ROUTINGALG_H

#include "GAGeneticAlgorithm.h"
#include "genomes/RoutingGenome.h"
#include "VNSOp.h"
#include <stdexcept>
#include <sstream>
#include <vector>
#include <list>
#include <algorithm>

using namespace std;

class RoutingAlg : public GAGeneticAlgorithm {

protected:
  RoutingGenome* bestsol_;
  RoutingGenome* actualsol_;

  int maxRouteCalls_;

  VNSOp *localsearch_;

  // To allow a change for a different method in the future
  long double cost(const RoutingGenome& g) { return g.nonPenalizedScore(); }

public:

  static const int SBESTPOS = 0;
  static const int SPOS     = SBESTPOS + 1;

  RoutingAlg(const GAGenome& g, VNSOp *ls ); // takes ownership of passed values
  RoutingAlg(const RoutingAlg&);

  virtual ~RoutingAlg();

  RoutingAlg& operator=(const RoutingAlg&);

  virtual void        initialize();
  virtual void        copy(const GAGeneticAlgorithm&);
  virtual RoutingAlg* clone() = 0;

  virtual RoutingGenome& bestSol()   const { return *bestsol_; }
  virtual RoutingGenome& actualSol() const { return *actualsol_; }

  virtual void objectiveFunction(GAGenome::Evaluator f){ GAGeneticAlgorithm::objectiveFunction(f); }
  virtual void objectiveData(const GAEvalData& v)      { GAGeneticAlgorithm::objectiveData(v);     }

  virtual GAScalingScheme&   scaling (const GAScalingScheme&   s) { return GAGeneticAlgorithm::scaling(s);  }
  virtual GASelectionScheme& selector(const GASelectionScheme& s) { return GAGeneticAlgorithm::selector(s); }

  void offspring (GAPopulation*) {}

  long double getBestSolScore()   const { return (& (bestSol())   ) ? bestSol().score()   : 0.0;}
  long double getActualSolScore() const { return (& (actualSol()) ) ? actualSol().score() : 0.0;}

  virtual int evalRouteCalls() const { return actualSol().evalRouteCalls(); }

  virtual void maxEvalRouteCalls(int value) { maxRouteCalls_ = value; }

  static GABoolean TerminateUponRouteCalls (GAGeneticAlgorithm&);

  virtual void localSearch(RoutingGenome& g);
  virtual void updatePenalizations();

  virtual void setNewRoutes(GAGenome& inf);
  virtual void setNewRoutes(vector<GAGenome*>& inds);

};

#ifdef GALIB_USE_STREAMS
inline STD_OSTREAM& operator<< (STD_OSTREAM& os, RoutingAlg & arg)
{ arg.write(os); return(os); }
inline STD_ISTREAM& operator>> (STD_ISTREAM& is, RoutingAlg & arg)
{ arg.read(is); return(is); }
#endif

#endif

#ifndef VNS_H
#define VNS_H

#include "RoutingAlg.h"
#include "VNSOpAction.h"
#include "VNSOp.h"
#include "VNSShakerResults.h"
#include "genomes/RoutingGenome.h"
#include "extras/combinations.h"
#include "logger/GALogger.h"
#include <stdexcept>
#include <sstream>
#include <vector>
#include <list>
#include <algorithm>

using namespace std;

class VNS : public RoutingAlg {

  double   t_;
  double   tinit_;
  double   ratioPerEval_;
  unsigned shakerpos_;

  int initialEvalRouteCalls_; // Used for the update of the t_ parameter of the Simulated Anhealing behavior.

  // Set by user
  double tinitratio_;
  double tinitprob_;
  double ls_ratio_4improv1_;
  double ls_prob_;
  double ls_ratio_4improv2_;
  bool   update_penalizations_;
  int    maxevals_shaker_;
  bool   use_allshakers_in_eachstep_;
  int    max_shakers_res_;

  vector< VNSShaker* > shakers_;

  vector<const VNSShaker*> successfulshaker_; // Could be a combination of shakers

  bool stepCostCondition1(RoutingGenome& s1) { return cost(s1) <  ( ls_ratio_4improv1_ * cost(*actualsol_) );}
  bool stepCostCondition2(RoutingGenome& s1) { return cost(s1) >= ( ls_ratio_4improv2_ * cost(*actualsol_) );}

  RoutingGenome*    neighbor                   (const RoutingGenome& g);
  VNSShakerResults* bestGenFromAllShakers      (const RoutingGenome& gen);
  VNSShakerResults* bestGenFromShakerAttempts  (const RoutingGenome& gen, int shakerpos);
  RoutingGenome*    tryAllCombsAndCloneBestOpt (const RoutingGenome& gen, const VNSShakerResults& results);

  void resetShakerPos() { shakerpos_ = 0; }

  bool feasible(const RoutingGenome& g);

public:

  VNS(const  GAGenome& g,
      double tinitratio,           double tinitprob,
      double ls_ratio_4improv1,    double ls_prob, double ls_ratio4improv2,
      bool   update_penalizations,
      int    maxevalsshaker,
      bool   use_allshakers_in_eachstep_,
      int    num_shakres_stored_,
      vector< VNSShaker* > shakers,    // takes ownership of passed values
      VNSOp *ls );                     // takes ownership of passed values

  VNS(const VNS&);

  virtual ~VNS();

  VNS&                operator=(const VNS&);
  virtual void        copy(const VNS&);
  virtual RoutingAlg* clone();

  virtual void  initialize();
  virtual void  step();

  // For the logstats
  string successfulShaker() const;
  int    shakerPos()        const { return shakerpos_; }

  void setNewRoutes(GAGenome& inf);           // Overriden to reset the shaker pos
  void setNewRoutes(vector<GAGenome*>& inds); // Overriden to reset the shaker pos

};

#ifdef GALIB_USE_STREAMS
inline STD_OSTREAM& operator<< (STD_OSTREAM& os, VNS & arg)
{ arg.write(os); return(os); }
inline STD_ISTREAM& operator>> (STD_ISTREAM& is, VNS & arg)
{ arg.read(is); return(is); }
#endif

#endif

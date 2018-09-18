#ifndef VNSSHAKERRESULT
#define VNSSHAKERRESULT

#include "VNSOpAction.h"
#include "genomes/GAGenome.h"
#include "VNSOp.h"
#include <vector>
#include <list>

using namespace std;

// Only for managing the results of the shakers
struct VNSShakerResult {

  const GAGenome*    genome;
  const VNSOpAction* action;
  const VNSShaker*   shaker;

  VNSShakerResult(const GAGenome* genom, const VNSOpAction* act, const VNSShaker* shak ) : genome(genom->clone()),
                                                                                        action(act->clone()),
                                                                                        shaker(shak) {}
  ~VNSShakerResult() { delete genome; delete action; }

  VNSShakerResult* clone() { return new VNSShakerResult(genome,action,shaker); }
};

class VNSShakerResults {

  int                   max_size_;
  vector<VNSShakerResult*> results_;

  // Used for sorting the results. The best results should be placed first
  static bool shakComp(const VNSShakerResult* a, const VNSShakerResult* b);

  VNSShakerResult* last() { return results_[results_.size()-1]; }

  void addResult(VNSShakerResult* newres); // Takes ownership of the ShakerResult* passed

  bool isIncluded(VNSShakerResult& newres);

  const VNSOpAction& actionAt(int pos) const;

public:

  VNSShakerResults(int size) : max_size_(size) {}

  ~VNSShakerResults() {
    for (vector<VNSShakerResult*>::iterator it=results_.begin(); !(it==results_.end()); it++) delete *it;
  }

  void addResults(const VNSShakerResults& other);

  void addResult(const GAGenome* gen, const VNSOpAction* act, const VNSShaker* shak);

  int size() const { return results_.size(); }

  VNSShakerResult& best() const { return * results_.front(); }

  VNSShakerResult& resultOf(int pos) const { return * results_[pos]; }

  bool areSelectedActionsCompatible(const list<int>& sel_actions) const;

  GAGenome* newGenFromSelActionsExec(const GAGenome& gen, const list<int>& sel_comb) const;

};

#endif

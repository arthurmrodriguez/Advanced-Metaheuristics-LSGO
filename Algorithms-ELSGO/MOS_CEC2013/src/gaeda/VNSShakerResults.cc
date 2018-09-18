#include "VNSShakerResults.h"
#include "logger/GALogger.h"
#include <algorithm>

bool VNSShakerResults::shakComp(const VNSShakerResult* a, const VNSShakerResult* b) {
  return GAGenome::compareScores(a->genome->score(),b->genome->score()) == GAGenome::BETTER;
}

void VNSShakerResults::addResult(VNSShakerResult* newres) {
  if ( !isIncluded(*newres) and
      ( results_.size() < max_size_ or
        GAGenome::compareScores( newres->genome->score(), last()->genome->score() ) == GAGenome::BETTER) ) {

    results_.push_back(newres);
    sort(results_.begin(),results_.end(),shakComp);

    if (results_.size() > max_size_) {
      VNSShakerResult* tmp = results_[results_.size()-1];
      delete tmp;
      results_.resize(max_size_);
    }
  }
  else delete newres;
}

bool VNSShakerResults::isIncluded(VNSShakerResult& res) {
  bool included = false;
  // We start by the last score so whenever we find a better score we stop the search since the results should be
  // sorted by score
  for (int i=size()-1; i>=0; i--) {
    if (GAGenome::compareScores( res.genome->score(), resultOf(i).genome->score() ) == GAGenome::WORSE) break;
    else if (GAGenome::compareScores( res.genome->score(), resultOf(i).genome->score() ) == GAGenome::EQUAL) {
      if (res.action->equal(* resultOf(i).action) ) {
        included = true;
        break;
      }
    }
  }

  return included;
}

const VNSOpAction& VNSShakerResults::actionAt(int pos) const {
  assert(pos >=0 && pos < results_.size() );
  return * results_[pos]->action;
}

void VNSShakerResults::addResults(const VNSShakerResults& other) {
  for (int i=0; i<other.results_.size(); i++) addResult( other.results_[i]->clone() );
}

void VNSShakerResults::addResult(const GAGenome* gen, const VNSOpAction* act, const VNSShaker* shaker) {
  VNSShakerResult* newres = new VNSShakerResult(gen,act,shaker);

  addResult(newres);
}

bool VNSShakerResults::areSelectedActionsCompatible(const list<int>& sel_actions) const {
  for (list<int>::const_iterator it=sel_actions.begin(); it!= sel_actions.end(); it++) {
    for (list<int>::const_iterator itj=it;itj!=sel_actions.end();itj++) {
      if (it == itj) continue;
      if (! actionAt(*it).isCompatibleWith( actionAt(*itj) ) or
          ! actionAt(*itj).isCompatibleWith( actionAt(*it) ) ) {
        return false;
      }
    }
  }
  return true;
}

GAGenome* VNSShakerResults::newGenFromSelActionsExec(const GAGenome& gen, const list<int>& sel_comb) const {
  GAGenome* new_gen = gen.clone();

  for (list<int>::const_iterator it=sel_comb.begin(); it!=sel_comb.end(); it++) {
    actionAt(*it).executeOver(*new_gen);
  }

  return new_gen;
}




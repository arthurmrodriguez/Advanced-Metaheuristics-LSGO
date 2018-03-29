#include "math.h"
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "VNS.h"

#include "GAEDAConfig.h"
#include "garandom.h"
#include "genomes/GAGenome.h"
#include "logger/GALogger.h"
#include "islands/CommManager.h"

/* Original DARP VNS values
 *
 * tinitratio        = 0.005
 * tinitprob         = 0.2
 * ls_ratio_4improv1 = 1.02
 * ls_prob           = 0.01
 * ls_ratio_4improv2 = 1.05
 */

VNS::VNS(const  GAGenome& g,
         long double tinitratio,           long double tinitprob,
         long double ls_ratio_4improv1,    long double ls_prob, long double ls_ratio_4improv2,
         bool   update_penalizations,
         int    maxevalsshaker,
         bool   use_allshakers_in_eachstep,
         int    max_shakers_res,
         vector< VNSShaker* > shakers,
         VNSOp *ls ) :
                      RoutingAlg(g,ls),
                      t_(0),
                      tinit_(0),
                      ratioPerEval_(0),
                      tinitratio_(tinitratio),
                      tinitprob_(tinitprob),
                      ls_ratio_4improv1_(ls_ratio_4improv1),
                      ls_prob_(ls_prob),
                      ls_ratio_4improv2_(ls_ratio_4improv2),
                      update_penalizations_(update_penalizations),
                      maxevals_shaker_(maxevalsshaker),
                      use_allshakers_in_eachstep_(use_allshakers_in_eachstep),
                      max_shakers_res_(max_shakers_res),
                      shakers_(shakers) {
  resetShakerPos();
  // TODO: Llevar esto a alguna clase generica para todos los algoritmos (Ej clase Algorithm)
  if (CommManager::instance()->isIslandMaster()) {
    cout << endl;
    cout << "VNS algorithm parameter values " << endl;
    cout << "tinitratio                       = " << tinitratio_ << endl;
    cout << "tinitprob                        = " << tinitprob_ << endl;
    cout << "ls_ratio_4improv1                = " << ls_ratio_4improv1_ << endl;
    cout << "ls prob                          = " << ls_prob_ << endl;
    cout << "ls_ratio_4improv2                = " << ls_ratio_4improv2_ << endl;
    cout << "update penalizations             = "; if (update_penalizations_) cout << "true"; else cout << "false"; cout << endl;
    cout << "maxevals_shaker                  = " << maxevals_shaker_ << endl;
    cout << "use all shakers at each step     = "; if (use_allshakers_in_eachstep) cout << "true"; else cout << "false"; cout << endl;
    cout << "maximum number of stored results = " << max_shakers_res_ << endl;
    cout << "Shakers = "; cout << "[ "; for (int i=0; i<shakers_.size(); i++) cout << shakers_[i]->name() << ", "; cout << " ]" << endl;
    cout << "LocalSearch = " << localsearch_->name() << endl << endl;
  }

  assert(max_shakers_res_ > 0);
}

VNS::VNS(const VNS& alg) : RoutingAlg(alg){
  copy(alg);
}

VNS::~VNS(){
  for (int i=0; i<shakers_.size(); i++) delete shakers_[i];
}

VNS& VNS::operator=(const VNS& alg){
  if(&alg != this) copy(alg);
  return *this;
}

void VNS::copy(const VNS& alg){
  RoutingAlg::copy(alg);

  t_            = alg.t_;
  tinit_        = alg.tinit_;
  ratioPerEval_ = alg.ratioPerEval_;
  shakerpos_    = alg.shakerpos_;

  initialEvalRouteCalls_ = alg.initialEvalRouteCalls_;

  tinitratio_                 = alg.tinitratio_;
  tinitprob_                  = alg.tinitprob_;
  ls_ratio_4improv1_          = alg.ls_ratio_4improv1_;
  ls_prob_                    = alg.ls_prob_;
  ls_ratio_4improv2_          = alg.ls_ratio_4improv2_;
  update_penalizations_       = alg.update_penalizations_;
  maxevals_shaker_            = alg.maxevals_shaker_;
  use_allshakers_in_eachstep_ = alg.use_allshakers_in_eachstep_;
  max_shakers_res_            = alg.max_shakers_res_;

  for (int i=0; i<shakers_.size(); i++) delete shakers_[i];
  shakers_.resize( alg.shakers_.size() );
  for (int i=0; i<shakers_.size(); i++) shakers_[i] = dynamic_cast<VNSShaker*>(alg.shakers_[i]->clone());
}

RoutingAlg* VNS::clone() {
  RoutingAlg* alg = new VNS(*this);
  return alg;
}

/**
 * Initialize parameters (commentary needs to be completed).
 *
 * @param seed Random seed
 */
void VNS::initialize() {
  RoutingAlg::initialize();

  initialEvalRouteCalls_ = evalRouteCalls();

  // To initialize temperature if the initial solution is already
  // feasible (by chance or by the application of a heuristic)
  feasible (*actualsol_);

  successfulshaker_.clear();
}


bool VNS::feasible(const RoutingGenome& g) {
  bool f = g.feasible();

  // First time a feasible solution arrives, we should compute
  // the initial temperature
  if (f and tinit_ == 0) {
    t_ = tinit_ = ( tinitratio_ * bestsol_->score()) / log(1.0/tinitprob_ );
    ratioPerEval_ = t_ / (long double) GAEDAConfig::handle()->getEvals();
  }

  return f;
}

RoutingGenome* VNS::neighbor(const RoutingGenome& g) {

  if (g.numRoutes() == 1) return dynamic_cast<RoutingGenome*>(g.clone());
  // The shakers are intra-routes mutator so if only one route is found we
  // return the same genome just for the case where the local search is applied

  VNSShakerResults* shak_results = (use_allshakers_in_eachstep_) ? bestGenFromAllShakers(g) : bestGenFromShakerAttempts(g,shakerpos_);

  RoutingGenome* best_gen = tryAllCombsAndCloneBestOpt(g,*shak_results);

  best_gen->evaluate();

  delete shak_results;

  return best_gen;
}

VNSShakerResults* VNS::bestGenFromAllShakers(const RoutingGenome& gen) {
  VNSShakerResults* results = new VNSShakerResults(max_shakers_res_);

  for (int shakerpos=0; shakerpos<shakers_.size(); shakerpos++) {
    VNSShakerResults * tmp_res = bestGenFromShakerAttempts(gen,shakerpos);

    results->addResults(*tmp_res);

    delete tmp_res;
  }

  return results;
}

VNSShakerResults* VNS::bestGenFromShakerAttempts(const RoutingGenome& gen, int shakerpos) {
  VNSShakerResults* shaker_results = new VNSShakerResults(max_shakers_res_);

  int used_evals = 0;
  do {
    RoutingGenome* tmp_gen = dynamic_cast<RoutingGenome*>(gen.clone());
    assert(tmp_gen); assert(tmp_gen->nevals() == 0);

    VNSOpAction* act = (*shakers_[shakerpos])(*tmp_gen);
    tmp_gen->evaluate(); // So that the nevals are computed right afterwards (e.g. creating the ShakerResult would force
                         // the computation of extra nevals

    shaker_results->addResult(tmp_gen,act,shakers_[shakerpos]);

    bool score_improved = GAGenome::compareScores(tmp_gen->score(), gen.score()) == GAGenome::BETTER;
    bool cost_improved  = cost(*tmp_gen) < cost(gen);

    shakers_[shakerpos]->addResult(score_improved,cost_improved,tmp_gen->nevals());

    // Need to be placed after the call to tmp_gen->score()
    used_evals += tmp_gen->nevals();
    stats.nummut++;

    delete tmp_gen;
    delete act;

  } while (used_evals > 0 && used_evals < maxevals_shaker_ && (*shakers_[shakerpos]).isStochastic() );
  // if for any reason a shaker cannot be applied we stop trying to iterate over the same shaker

  stats.numeval += used_evals;

  return shaker_results;
}

RoutingGenome* VNS::tryAllCombsAndCloneBestOpt (const RoutingGenome& orig_gen, const VNSShakerResults& results) {
  vector<int> positions;
  for (int i=0; i<results.size(); i++) positions.push_back(i);

  list<int>            tmp_sol;
  vector< list<int>* > allcombs;
  //TODO: This should be analyzed and if it is inefficient, it should be replaced with an optimized method
  for (int combsize=2; combsize<=results.size(); combsize++) {
    tmp_sol.clear();
    getAllCombinations(positions,combsize,tmp_sol,0, allcombs);
  }

  RoutingGenome* best_gen = dynamic_cast<RoutingGenome*>(results.best().genome->clone());
  assert(best_gen);

  successfulshaker_.clear();
  if ( GAGenome::compareScores(best_gen->score(),orig_gen.score()) == GAGenome::BETTER ) {
    successfulshaker_.push_back(results.best().shaker);
  }

  RoutingGenome* tmp_gen  = 0;

  for (int i=0; i<allcombs.size(); i++) {
    assert( allcombs[i]->size() > 1 );
    if (results.areSelectedActionsCompatible( * allcombs[i] ) ) {

      tmp_gen =  dynamic_cast<RoutingGenome*>( results.newGenFromSelActionsExec(orig_gen,*allcombs[i]) );
      assert(tmp_gen);
      tmp_gen->evaluate(); // So that the nevals are computed right afterwatds
      stats.numeval += tmp_gen->nevals();

      if ( GAGenome::compareScores(tmp_gen->score(),best_gen->score()) == GAGenome::BETTER ) {
        delete best_gen;
        best_gen = tmp_gen;

        successfulshaker_.clear();
        for (list<int>::iterator it=allcombs[i]->begin(); it!=allcombs[i]->end(); it++) {
          successfulshaker_.push_back( results.resultOf(*it).shaker );
        }
      }
      else {
        delete tmp_gen;
      }
    }

  }

  for(int i=0; i<allcombs.size(); i++) delete allcombs[i];

  return best_gen;
}

void VNS::step(){
  /* LOG */ GALogger::instance()->appendPopulation("VNSStep::step", "Before doing the step", *pop);

  // Line 8: Randomly compute s1 in Nk(s) (SHAKING)
  RoutingGenome* s1 = neighbor(*actualsol_);

  /* LOG */ GALogger::instance()->appendOperatorResult("VNSStep::step: " + shakers_[shakerpos_]->name(), *actualsol_, *s1);

  // Line 10:: Apply Local Search to s1 yielding s2 (LOCAL SEARCH)
  long double prand = GARandomDouble(0, 1);
  if ( s1->score() != actualsol_->score() and ( stepCostCondition1(*s1) or prand < ls_prob_ ) ) {
    localSearch(*s1);

    /* LOG */ GALogger::instance()->appendOperatorResult("VNSStep::step: Local Search (1): ", *s1, *s1);
  }
  
  // Line 14: Move or not
  long double pSA = (t_ > 0) ? (1.0 / (exp((s1->score() - bestsol_->score()) / t_))) : -1;
//   cout << "t= " << t_ << " s2 - sbest = " << s2->score() - sbest_->score() << " dentro de la exp " << (s2->score() - sbest_->score()) / t_ << " exp = " << exp((s2->score() - sbest_->score()) / t_) << endl;
  prand = GARandomDouble(0,1);

  if ( s1->score() != actualsol_->score() and ( GAGenome::compareScores(s1->score(), actualsol_->score()) == GAGenome::BETTER or prand < pSA) ) {
    if ( stepCostCondition2(*s1 )) {
      localSearch(*s1);
      /* LOG */ GALogger::instance()->appendOperatorResult("VNSStep::step: Local Search (2): ", *s1, *s1);
      s1->evaluate();
    }

    actualsol_->copy(*s1);

    if (update_penalizations_) updatePenalizations();

    shakerpos_ = -1;  // -1 since later is increased by one unit
  }

  // Line 19: Best solution update
  if (feasible(*s1) and GAGenome::compareScores(s1->score(), bestsol_->score()) == GAGenome::BETTER ) {
    bestsol_->copy(*s1);
  }

  // Line 21: Update k
  shakerpos_ = (shakerpos_ + 1) % shakers_.size();

  // Linearly decrease temperature. We do not take into account the evaluations consumed in the initialization
  t_ = (t_ == 0) ? 0 : tinit_ - (ratioPerEval_ * (evalRouteCalls() - initialEvalRouteCalls_) );
//   cout << "t_ = " << t_ << "tinit_ " << tinit_ << " pSA antiguo " << pSA << " ratio per eval " << ratioPerEval_ << " evalRouteCalls= " << evalRouteCalls() << endl;

  delete s1;

//   s_->debugGenome();

  stats.update(*pop);     // update the statistics by one generation

  printStats("End of Step Stats");

  //GNAPA para imprimir las stats de los shakers  y de la ls al final
  /* LOG */ if (done()) {
  /* LOG */   stringstream msg; msg << endl << "Shakers and LS results" << endl;
  /* LOG */   for (int i=0; i<shakers_.size(); i++) msg << shakers_[i]->getStatsResults() << endl;
  /* LOG */   msg << localsearch_->getStatsResults() << endl;
  /* LOG */   GALogger::instance()->appendLogMessage("VNS:end step", msg.str(), GALogger::only_stats);
            }
}

string VNS::successfulShaker() const {
  switch (successfulshaker_.size()) {
  case 0:
    return "None";
  case 1:
    return successfulshaker_.front()->name();
  default:
    stringstream msg; msg << "Combination: ";
    for (int i=0; i<successfulshaker_.size(); i++) msg << successfulshaker_[i]->name() << "-";
    return msg.str();
  }
}

/*
 * Overriden to reset the shaker pos
 */
void VNS::setNewRoutes(GAGenome& inf) {
  RoutingAlg::setNewRoutes(inf);
  resetShakerPos();
}

/*
 * Overriden to reset the shaker pos
 */
void VNS::setNewRoutes(vector<GAGenome*>& inds) {
  RoutingAlg::setNewRoutes(inds);
  resetShakerPos();
}

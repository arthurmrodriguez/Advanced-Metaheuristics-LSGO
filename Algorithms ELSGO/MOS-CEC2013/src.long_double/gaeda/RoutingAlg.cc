#include "RoutingAlg.h"

#include "GAEDAConfig.h"
#include "garandom.h"
#include "genomes/GAGenome.h"
#include "logger/GALogger.h"
#include "islands/CommManager.h"

#include <math.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdlib.h>

RoutingAlg::RoutingAlg(const  GAGenome& g, VNSOp *ls ) : GAGeneticAlgorithm(g),
                                                         bestsol_(0),
                                                         actualsol_(0),
                                                         localsearch_(ls) {
  if (CommManager::instance()->isIslandMaster()) {
    cout << endl;
    cout << "RoutingAlg algorithm parameter values " << endl;
    if (ls) cout << "LocalSearch        = " << localsearch_->name() << endl << endl;
  }
}

RoutingAlg::RoutingAlg(const RoutingAlg& alg) : GAGeneticAlgorithm(alg){
  copy(alg);
}

RoutingAlg::~RoutingAlg(){
  delete localsearch_;
}

RoutingAlg& RoutingAlg::operator=(const RoutingAlg& alg){
  if(&alg != this) copy(alg);
  return *this;
}

void RoutingAlg::initialize() {
  // Needed here to override default population size configured in the main program
  // TODO: devolver a 1
  pop->size(2);

  bestsol_   = dynamic_cast<RoutingGenome*>( & pop->individual(SBESTPOS) );
  actualsol_ = dynamic_cast<RoutingGenome*>( & pop->individual(SPOS)     );
  assert(bestsol_ and actualsol_);

  actualsol_->initialize();
  actualsol_->evaluate(gaTrue);

  bestsol_->copy(*actualsol_);

  stats.reset(*pop);

  printStats("Initial Stats");
}

void RoutingAlg::copy(const GAGeneticAlgorithm& other){
  GAGeneticAlgorithm::copy(other);

  const RoutingAlg& alg = dynamic_cast<const RoutingAlg&>(other);
  assert(&alg);

  bestsol_   = dynamic_cast<RoutingGenome*>( & pop->individual(SBESTPOS) );
  actualsol_ = dynamic_cast<RoutingGenome*>( & pop->individual(SPOS)     );

  maxRouteCalls_ = alg.maxRouteCalls_;

  if (localsearch_ != 0) delete localsearch_;
  localsearch_ = alg.localsearch_->clone();

  assert(bestsol_ and actualsol_);
}

void RoutingAlg::localSearch(RoutingGenome& g) {
  int start_nevals = g.nevals();

  VNSOpAction* act = (*localsearch_)(g);
  delete act;

  int used_evals = g.nevals() - start_nevals;
  stats.numeval += used_evals;
  stats.nummut++;

  bool score_improved = GAGenome::compareScores(g.score(), g.score()) == GAGenome::BETTER;
  bool cost_improved  = cost(g) < cost(g);

  localsearch_->addResult(score_improved,cost_improved,used_evals);
}

GABoolean RoutingAlg::TerminateUponRouteCalls (GAGeneticAlgorithm& a) {
  RoutingAlg* alg = dynamic_cast<RoutingAlg*>(&a); assert(alg);

  GABoolean val = ( &(alg->actualSol()) != 0 and alg->evalRouteCalls() >= alg->maxRouteCalls_) ?
                  gaTrue : gaFalse;

  return val;
}

void RoutingAlg::updatePenalizations() {
  actualSol().updatePenalizations();

  // Update the values of the current solutions to the new penalizations.
  // Note that the objective function is going to be called but the computation is
  // going to be considerably less since only a sum is going to be carried out. However, the nevals
  // attribute is going to be incremented by the evaluate call.
  actualSol().evaluate(gaTrue);
  bestSol().evaluate(gaTrue);
}

/*
 * Sets the same routes as the genome passed. Updates both the bestSol and the actualSol
 */
void RoutingAlg::setNewRoutes(GAGenome& inf) {
  RoutingGenome& gen = dynamic_cast<RoutingGenome&>(inf); assert(&gen);
  bestSol().copy(gen); // OJO con que tenga el estado puesto a ya calculado
  actualSol().copy(bestSol());
}

/*
 * Similar method than the previous one. It updates both the bestSol and the actualSol.
 * In this case, both solutions get their flags set so that they do not use cached values.
 */
void RoutingAlg::setNewRoutes(vector<GAGenome*>& inds) {
  assert(inds.size() > 0);

  bestSol().emptyRoutes();

  for (int i=0; i<inds.size(); i++) bestSol().addRoutes(* inds[i] );

  actualSol().copy(bestSol());

  assert(&(bestSol()));
  assert(&(actualSol()));
}



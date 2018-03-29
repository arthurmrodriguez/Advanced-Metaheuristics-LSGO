#include "TSCordeau.h"

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

/* Original TSCordeau values
 *
 * localsearchfreq = 10
 *
 */

TSCordeau::TSCordeau(const  GAGenome& g, VNSOp *ls ) : RoutingAlg(g,ls),
                                                       lsPeriod_(10),
                                                       paramsUpdatePeriod_(10) {
  if (CommManager::instance()->isIslandMaster()) {
    cout << endl;
    cout << "TSCordeau algorithm parameter values " << endl;
    cout << "lsPeriod           = " << lsPeriod_ << endl;
    cout << "paramsupdatePeriod = " << paramsUpdatePeriod_ << endl;
    cout << "LocalSearch        = " << localsearch_->name() << endl << endl;
  }
}

TSCordeau::TSCordeau(const TSCordeau& alg) : RoutingAlg(alg){
  copy(alg);
}

TSCordeau& TSCordeau::operator=(const TSCordeau& alg){
  if(&alg != this) copy(alg);
  return *this;
}

void TSCordeau::copy(const GAGeneticAlgorithm& other) {
  RoutingAlg::copy(other);
  const TSCordeau& alg = dynamic_cast<const TSCordeau&>(other);

  lsPeriod_           = alg.lsPeriod_;
  paramsUpdatePeriod_ = alg.paramsUpdatePeriod_;
}

RoutingAlg* TSCordeau::clone() {
  TSCordeau* alg = new TSCordeau(*this);
  return alg;
}

void TSCordeau::step(){
  /* LOG */ GALogger::instance()->appendPopulation("TSCordeauStep::step", "Before doing the step", *pop);

  // bestFromNeighborhood returns a solution that is better from the actual or worst but compatible with the TABU
  // constraints. Therefore, if a 0 is returned, that would imply that no better solution or compatible with the
  // TABU memory was found.

  pair<RoutingGenome*,int> neighborsol = actualsol_->bestFromNeighborhood();
  RoutingGenome* bestneighbor = neighborsol.first;
  stats.numeval += neighborsol.second; assert(neighborsol.second>0);

  if (bestneighbor) {
    actualsol_->copy(*bestneighbor);
    assert(bestneighbor->score() == actualsol_->score());

    /* LOG */ GALogger::instance()->appendOperatorResult("TSCordeauStep::step: actual sol + neighbor " , *actualsol_, *bestneighbor);

    if (bestneighbor->feasible() and GAGenome::compareScores(bestneighbor->score(),bestsol_->score()) == GAGenome::BETTER ) {
      localSearch(*actualsol_);
      if (actualsol_->feasible()) bestsol_->copy(*actualsol_); // could be unfeasible after the local search
    }

    delete bestneighbor;
  }

  if ( stats.generation() > 1 and stats.generation() % lsPeriod_ == 0) localSearch(*actualsol_);

  /* LOG */ GALogger::instance()->appendOperatorResult("TSCordeauStep::step: After Possible Local Search: ", *actualsol_, *actualsol_);

  if ( stats.generation() > 1 and stats.generation() % paramsUpdatePeriod_ == 0) actualsol_->updateDynParameters();

  updatePenalizations();

  stats.update(*pop);     // update the statistics by one generation

  printStats("End of Step Stats");

  //GNAPA para imprimir las stats de los shakers  y de la ls al final
  /* LOG */ if (done()) {
  /* LOG */   stringstream msg; msg << endl << "LS results" << endl;
  /* LOG */   msg << localsearch_->getStatsResults() << endl;
  /* LOG */   msg << "Frequency Memory (only values > 1 are displayed): " << endl << actualsol_->memoryDataInformation() << endl;
  /* LOG */   GALogger::instance()->appendLogMessage("TSCordeau:end step", msg.str(), GALogger::only_stats);
  /* LOG */ }

}

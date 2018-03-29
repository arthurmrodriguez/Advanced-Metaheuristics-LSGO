#include "MOSTechniqueCMAES.h"

#include "GAEDAConfig.h"
#include "MOSTechniqueSet.h"
#include "PopElitism.h"
#include "logger/GALogger.h"
#include "genomes/MOSGenome.h"
#include "genomes/GA1DArrayGenome.h"

#include "cmaes_interface.h"

#include <vector>

/**
 * MOSTechniqueCMAES Constructor
 */
MOSTechniqueCMAES::MOSTechniqueCMAES(techIdType id, std::string description, GAGenome::Initializer init,
                                     GAGenome::Evaluator evaluator, encodingType encoding, GAGenome* genomeBase,
                                     GASelectionScheme* selector) {

  _id          = id;
  _description = description;

  _encoding    = encoding;
  _genomeBase  = genomeBase;

  _quality     = 0.0;
  _partRatio   = 0.0;

  _initializer = init;
  _evaluator   = evaluator;
  _selector    = selector;

  _recombinator = new PopElitism();

  unsigned psize = GAEDAConfig::handle()->getProblemSize();
  long double stddevs [psize];
  long double initialX[psize];

  std::fill_n(stddevs,  psize, 0.3);
  std::fill_n(initialX, psize, 0.5);

  GA1DArrayAlleleGenome<long double>& tmp_gen = dynamic_cast<GA1DArrayAlleleGenome<long double>&>(*genomeBase);
  for (unsigned i = 0; i < psize; i++)
    initialX[i] = tmp_gen.gene(i);

  cmaes_output_dir (&_cmaes_evo, (GAEDAConfig::handle()->getBBOBOutputPath() + "/" + GAEDAConfig::handle()->getBBOBConfigName() + "/config/").c_str());
  _fitvals = cmaes_init (&_cmaes_evo, psize, initialX, stddevs, 0, GAEDAConfig::handle()->getPopSize(), GAEDAConfig::handle()->getCMAESInitials().c_str());

  _old_best = 0.0;
  _old_worst = 0.0;

  return;
}


bool MOSTechniqueCMAES::restartRequired() {
  return (((cmaes_TestForTermination(&_cmaes_evo))) ? true : false);
}


bool MOSTechniqueCMAES::restartInnerData(GAPopulation* pop) {

  unsigned psize = GAEDAConfig::handle()->getProblemSize();
  long double stddevs [psize];
  long double initialX[psize];
  std::fill_n(stddevs,  psize, 0.3);

  GA1DArrayAlleleGenome<long double>& tmp_gen = dynamic_cast<GA1DArrayAlleleGenome<long double>&>(*dynamic_cast<MOSGenome&>(pop->individual(GARandomInt(0, pop->size()-1))).getGenome(GAID::RealEncoding));
  for (unsigned i = 0; i < psize; i++)
    initialX[i] = tmp_gen.gene(i);

  cmaes_exit(&_cmaes_evo);
  _fitvals = cmaes_init(&_cmaes_evo, psize, initialX, stddevs, 0, pop->size(), GAEDAConfig::handle()->getCMAESInitials().c_str());
  _old_best = 0.0;
  _old_worst = 0.0;

  return true;

}


unsigned MOSTechniqueCMAES::offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged) {

  assert(pop->size() == auxPop->size());

  _times_improved = 0;

  long double*const* cmaes_pop;

  unsigned old_lambda = _cmaes_evo.sp.lambda;

  if (cmaes_Get(&_cmaes_evo, "lambda") > maxEvals)
    _cmaes_evo.sp.lambda = maxEvals;

  long double fits[(int)cmaes_Get(&_cmaes_evo, "lambda")];

  cmaes_pop = cmaes_SamplePopulation(&_cmaes_evo);

  for (unsigned i = 0; i < cmaes_Get(&_cmaes_evo, "lambda"); i++) {
    MOSGenome& gen = dynamic_cast<MOSGenome&>(auxPop->individual(i));
    GA1DArrayAlleleGenome<long double>& g = dynamic_cast<GA1DArrayAlleleGenome<long double>&> (*gen.getGenome(GAID::RealEncoding));

    // Clean old encodings (except that used by this technique)
    gen.purgeGenome(this);

    for (unsigned j = 0; j < cmaes_Get(&_cmaes_evo, "dim"); j++)
      g.gene(j, cmaes_pop[i][j]);

    g.evaluate();
    fits[i] = gen.score(g.score());

    if (gen.precissionReached()) converged = true;

    gen.setFitnessIncrement(gen.computeFitnessIncrement(pop->ave()));

    gen.mustComputeQuality(true);

    if (GAGenome::compareScores(gen.score(), pop->ave()) == GAGenome::BETTER)
      _times_improved++;
  }

  cmaes_UpdateDistribution(&_cmaes_evo, fits);

  unsigned newInds = (maxEvals < cmaes_Get(&_cmaes_evo, "lambda")) ? maxEvals : cmaes_Get(&_cmaes_evo, "lambda");

  _cmaes_evo.sp.lambda = old_lambda;

  usedEvals      += newInds;
  _stats.numrep  += newInds;
  _stats.numeval += newInds;

  _old_best = auxPop->best().score();
  _old_worst = auxPop->worst().score();

  return newInds;

}


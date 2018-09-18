#include "MOSTechniqueAdapDE.h"

// #include "DEElitism.h"
// #include "MOSConversion.h"
#include "MOSTechniqueSet.h"
#include "GAEDAConfig.h"
#include "GAGenealogy.h"
#include "GAGenealogyMemory.h"
#include "logger/GALogger.h"
#include "genomes/MOSGenome.h"
#include "genomes/GA1DArrayGenome.h"

// #include <vector>

MOSTechniqueAdapDE::MOSTechniqueAdapDE(techIdType id, std::string description, GAGenome::Initializer init,
                                       GAGenome::Evaluator evaluator, encodingType encoding, GAGenome* genomeBase,
                                       GAGenome::DECrossover crossover, double F, double CR, GASelectionScheme* selector )
              : MOSTechniqueDE(id, description, init, evaluator, encoding, genomeBase, crossover, F, CR, selector ) {}

void MOSTechniqueAdapDE::offspring(MOSGenome& dad, MOSGenome& mom, GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset) {
  selfAdjustDEParams(*origPop);
  MOSTechniqueAdapDE::offspring(dad, mom, origPop, destPop, size, offset);
}

void MOSTechniqueAdapDE::offspring(GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset) {
  selfAdjustDEParams(*destPop);
  MOSTechniqueDE::offspring(origPop, destPop, size, offset);
}

unsigned MOSTechniqueAdapDE::offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged) {
  selfAdjustDEParams(*pop);
  return MOSTechniqueDE::offspring (maxEvals, usedEvals, pop, auxPop, converged);
}

void MOSTechniqueAdapDE::selfAdjustDEParams(const GAPopulation& pop) {
  static bool first_call = true;

  GAEDAConfig* cfg = GAEDAConfig::handle();

  double tau1 = cfg->getAdapDETauF ();
  double tau2 = cfg->getAdapDETauCR();
  double Fl   = cfg->getAdapDEFl   ();
  double Fu   = cfg->getAdapDEFu   ();

  unsigned sz = pop.size();
  if ( _Fv.size() != sz)  _Fv.resize(sz);
  if (_CRv.size() != sz) _CRv.resize(sz);

  // If it is the first call to the method, we initialize
  // both vectors to initial values defined in the technique.
  if (first_call) {
    for (unsigned i = 0; i < sz; i++) {
      _Fv[i]=_F;
      _CRv[i]=_CR;
    }
    first_call=false;
  }

  for (unsigned i = 0; i < sz; i++) {
    double r1 = GARandomDouble(0, 1);
    double r2 = GARandomDouble(0, 1);
    double r3 = GARandomDouble(0, 1);
    double r4 = GARandomDouble(0, 1);

     _Fv[i] = (r2 < tau1) ? (Fl + (r1 * Fu)) :  _Fv[i];
    _CRv[i] = (r4 < tau2) ? (r3)             : _CRv[i];
  }
}

void MOSTechniqueAdapDE::offspring_internal(MOSGenome& genomeM, MOSGenome& newGenomeM, MOSGenome& genome1M, MOSGenome& genome2M, MOSGenome& genome3M) {
  GA1DArrayAlleleGenome<double>* genome    = dynamic_cast <GA1DArrayAlleleGenome<double>*> (   genomeM.getGenome(_encoding));
  GA1DArrayAlleleGenome<double>* newGenome = dynamic_cast <GA1DArrayAlleleGenome<double>*> (newGenomeM.getGenome(_encoding));
  GA1DArrayAlleleGenome<double>* genome1   = dynamic_cast <GA1DArrayAlleleGenome<double>*> (  genome1M.getGenome(_encoding));
  GA1DArrayAlleleGenome<double>* genome2   = dynamic_cast <GA1DArrayAlleleGenome<double>*> (  genome2M.getGenome(_encoding));
  GA1DArrayAlleleGenome<double>* genome3   = dynamic_cast <GA1DArrayAlleleGenome<double>*> (  genome3M.getGenome(_encoding));

  static std::vector<GAGenome*> parents (4, (GAGenome*)0);
  parents[0]= genomeM.getGenome(_encoding);
  parents[1]=genome1M.getGenome(_encoding);
  parents[2]=genome2M.getGenome(_encoding);
  parents[3]=genome3M.getGenome(_encoding);

  double bestParentFit = genomeM.score();
  if (GAGenome::compareScores(genome1M.score(), bestParentFit) == GAGenome::BETTER) bestParentFit = genome1M.score();
  if (GAGenome::compareScores(genome2M.score(), bestParentFit) == GAGenome::BETTER) bestParentFit = genome2M.score();
  if (GAGenome::compareScores(genome3M.score(), bestParentFit) == GAGenome::BETTER) bestParentFit = genome3M.score();

  // BEGIN: Update strategies
  MOSTechniqueSet::handle()->updateStrategies(parents, newGenomeM.getGenome(_encoding));
  // END: Update strategies

  unsigned dim = genome->length();

  // Clean old encodings (except that used by this technique)
  newGenomeM.purgeGenome(this);

  // Mutation
  for (unsigned i = 0; i < dim; i++) {
    double new_val = genome1->gene(i) + _Fv[_curTarget] * (genome2->gene(i) - genome3->gene(i));

    if (new_val > newGenome->alleleset(i).upper())
      new_val = newGenome->alleleset(i).upper();
    else if (new_val < newGenome->alleleset(i).lower())
      new_val = newGenome->alleleset(i).lower();

    newGenome->gene (i, new_val);

  }

  GALogger::instance()->appendOperatorResult("Mutation result: ", *genome, *newGenome);

  // Crossover
  _stats.numcro += (*_crossover) (*genome, *newGenome, newGenome, _CRv[_curTarget]);

  GALogger::instance()->appendOperatorResult("Crossover result: ", *genome, *newGenome);

  // BEGIN: Genealogy
  // Put crossover children in the genealogy
  if(GAGenealogy::handle() != NULL)
    if (GAGenealogy::handle()->isGenealogyMemory())
      GAGenealogy::handle()->familyRelationship(genomeM, genomeM, newGenomeM, newGenomeM);
    else {

      std::vector<GAGenome*> children, parents;

      children.push_back(&newGenomeM);
      parents.push_back(&genomeM);
      parents.push_back(&genome1M);
      parents.push_back(&genome2M);
      parents.push_back(&genome3M);

      GAGenealogy::handle()->familyRelationship(parents, children);

    }
  // END: Genealogy

  _stats.numsel  += 4;
  _stats.nummut  += 1;
  _stats.numrep  += 1;
  _stats.numeval += 1;

  double newScore = newGenomeM.score();

  newGenomeM.setFitnessIncrement(newGenomeM.computeFitnessIncrement(bestParentFit));
  newGenomeM.mustComputeQuality(true);
  if (GAGenome::compareScores(newScore, bestParentFit) == GAGenome::BETTER)
    _times_improved++;

  // No need to mark children individuals to be evaluated. purgeGenome did that for us...
  return;
}

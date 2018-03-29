/**
 * @file
 * @brief EvolutionStrategy class impl.
 *
 * Implementation of the EvolutionStrategy class
 */

#include "EvolutionStrategy.h"

#include "normal.h"
#include "garandom.h"
#include "GAEDAConfig.h"
#include "logger/GALogger.h"
#include "genomes/GA1DArrayGenome.h"

GAParameterList& EvolutionStrategy::registerDefaultParameters(GAParameterList& p) {
  GAGeneticAlgorithm::registerDefaultParameters(p);

  p.add(gaNelitism, gaSNelitism, GAParameter::BOOLEAN, &DefaultParameters::elitism);
  return p;
}


EvolutionStrategy::EvolutionStrategy(const GAGenome& c) : GAGeneticAlgorithm(c),
							  _mu(GAEDAConfig::handle()->getMu()),
							  _ro(GAEDAConfig::handle()->getRo()),
							  _lambda(GAEDAConfig::handle()->getLambda())
{
  oldPop = pop->clone();
}


EvolutionStrategy::EvolutionStrategy(const GAPopulation& p) : GAGeneticAlgorithm(p),
							      _mu(GAEDAConfig::handle()->getMu()),
							      _ro(GAEDAConfig::handle()->getRo()),
							      _lambda(GAEDAConfig::handle()->getLambda())
{
  oldPop = pop->clone();
}


EvolutionStrategy::EvolutionStrategy(const EvolutionStrategy& ga) : GAGeneticAlgorithm(ga) {
  oldPop = (GAPopulation*) 0;
  copy(ga);
}


EvolutionStrategy::~EvolutionStrategy() {
  delete oldPop;
}


EvolutionStrategy& EvolutionStrategy::operator=(const EvolutionStrategy& ga) {
  if(&ga != this)
    copy(ga);
  return *this;
}

void EvolutionStrategy::copy(const GAGeneticAlgorithm& g) {
  GAGeneticAlgorithm::copy(g);

  const EvolutionStrategy& ga = DYN_CAST(const EvolutionStrategy&, g);

  if(oldPop)
    oldPop->copy(*(ga.oldPop));
  else
    oldPop = ga.oldPop->clone();

  oldPop->geneticAlgorithm(*this);

  _mu = ga._mu;
  _ro = ga._ro;
  _lambda = ga._lambda;
  _tau0 = ga._tau0;
  _tau = ga._tau;

}


int EvolutionStrategy::setptr(const char* name, const void* value) {
  int status = GAGeneticAlgorithm::setptr(name, value);

  if(strcmp(name, gaNelitism) == 0 || strcmp(name, gaSNelitism) == 0) {
    status = 0;
  }
  return status;
}


int EvolutionStrategy::get(const char* name, void* value) const {
  int status = GAGeneticAlgorithm::get(name, value);

  if(strcmp(name, gaNelitism) == 0 || strcmp(name, gaSNelitism) == 0) {
    status = 0;
  }
  return status;
}


void EvolutionStrategy::objectiveFunction(GAGenome::Evaluator f) {
  GAGeneticAlgorithm::objectiveFunction(f);
  for(unsigned i=0; i<pop->size(); i++)
    oldPop->individual(i).evaluator(f);
}


void EvolutionStrategy::objectiveData(const GAEvalData& v) {
  GAGeneticAlgorithm::objectiveData(v);
  for(unsigned i=0; i<pop->size(); i++)
    oldPop->individual(i).evalData(v);
}


const GAPopulation& EvolutionStrategy::population(const GAPopulation& p) {
  GAGeneticAlgorithm::population(p);

  GAPopulation* tmpPop = pop->clone();
  oldPop->copy(*tmpPop);
  delete tmpPop;

  oldPop->geneticAlgorithm(*this);

  return *pop;
}


unsigned EvolutionStrategy::populationSize(unsigned value) {
  GAGeneticAlgorithm::populationSize(value);
  oldPop->size(value);
  return value;
}


/**
 * Initialize the population, do a few stupidity
 * checks, reset the stats.  We must initialize the old pop because there is no
 * guarantee that each individual will get initialized during the course of our
 * operator++ operations.  We do not evaluate the old pop because that will
 * happen as-needed later on.
 *
 */
void EvolutionStrategy::initialize () {
  unsigned N = GAEDAConfig::handle()->getProblemSize();

  // Initialization of Tau vars
  _tau0 = 1 / (sqrt(2*N));
  _tau  = 1 / (sqrt(2*sqrt(N)));

  // Hard-coded initialization
  for (unsigned i = 0; i < pop->size(); i++) {

    GA1DArrayAlleleGenome<long double>& g = DYN_CAST(GA1DArrayAlleleGenome<long double>&, pop->individual(i));

    for (int j = 0; j < g.size(); j++) {
      g.stdDev(j, 0.3);
    }

    //g.initialize();

    for (int j = 0; j < g.size(); j++)
      g.gene(j, 0.5*normal_random(0,0.3*0.3));


  }

  //  pop->initialize();
  pop->evaluate(gaTrue); // the old pop will get it when the pops switch
  //  oldPop->initialize();

  stats.reset(*pop);

  if(!scross)
    GAErr(GA_LOC, className(), "initialize", gaErrNoSexualMating);

   printStats("Initial Stats");

}


void EvolutionStrategy::offspring (GAPopulation* offpop){

  offpop->size(_lambda);

  for (unsigned i = 0; i < _lambda; i++) {

    GA1DArrayAlleleGenome<long double>& g = DYN_CAST(GA1DArrayAlleleGenome<long double>&, offpop->individual(i));

    // Selection of parents for recombination
    static std::vector<GA1DArrayAlleleGenome<long double>*> parents (_ro, (GA1DArrayAlleleGenome<long double>*)0);

    for (unsigned j = 0; j < _ro; j++) {
      int indiv = GARandomInt(0, pop->size()-1);
      parents[j] = DYN_CAST(GA1DArrayAlleleGenome<long double>*, &(pop->individual(indiv)));
    }

    // Recombine both endogenous parameters and variables
    //IntermediateRecombination(parents, g);

    // Recombine both endogenous parameters and variables
    for (int j = 0; j < g.size(); j++) {
      long double sum = 0.0;

      // First, genes
      for (unsigned k = 0; k < _ro; k++)
	sum += (1.0/_ro) * parents[k]->gene(j);

      g.gene(j, sum);

      sum = 0.0;

      // Then, Std Devs
      for (unsigned k = 0; k < _ro; k++)
	sum += (1.0/_ro) * parents[k]->stdDev(j);

      g.stdDev(j, sum);

    }

    // Update mutation rates
    long double F = exp(_tau0*normal_random(0,1));

    for (int j = 0; j < g.size(); j++)
      g.stdDev(j, F*g.stdDev(j)*exp(_tau*normal_random(0,1)));

    // Mutate variables
    for (int j = 0; j < g.size(); j++)
      g.gene(j, g.gene(j)+(g.stdDev(j)*normal_random(0,1)));

    stats.numsel += _ro;
    stats.numcro++;
    stats.nummut++;
    stats.numeval++;
    stats.numrep += _lambda;

  }

}


/**
 * Evolve a new generation of genomes.  When we start this routine, pop
 * contains the current generation.  When we finish, pop contains the new
 * generation and oldPop contains the (no longer) current generation.
 * The previous old generation is lost. We don't deallocate any memory, we
 * just reset the contents of the genomes. The selection routine must return
 * a pointer to a genome from the old population. Needs to be refactored
 * into smaller functions!!!!
 */
void EvolutionStrategy::step() {

  /*LOG*/ GALogger::instance()->appendPopulation("EvolutionStrategyStep::step", "Before doing the step", *pop);

  offspring(oldPop);

  oldPop->evaluate(gaTrue);
  oldPop->sort();

  GALogger::instance()->appendPopulation("The generated children is: ", "", *oldPop);

  /*
  std::cout << "Best genome from oldPop: " << std::endl;
  std::cout << oldPop->best() << " score: " << oldPop->best().score() << std::endl;

  std::cout << "Meto de oldPop en pop: " << std::endl;
  */
  for (unsigned i = 0; i < _mu; i++) {
    pop->add(oldPop->individual(i));
    //std::cout << oldPop->individual(i) << " score: " << oldPop->individual(i).score() << std::endl;
  }

  pop->scale();
  pop->sort();

  //  std::cout << "Quito de pop: " << std::endl;
  for (unsigned i = 0; i < _mu; i++) {
    // Memory leak here...
    GAGenome* gen = pop->remove(GAPopulation::WORST, GAPopulation::SCALED);
    //std::cout << *gen << " score: " << gen->score() << std::endl;
    delete gen;
  }

  GALogger::instance()->appendPopulation("After recombination the population is: ", "", *pop);

  stats.update(*pop);   // update the statistics by one generation
  printStats("End of Step Stats");
}


/**
 * Some syntactic sugar for the step() method
 */
EvolutionStrategy & EvolutionStrategy::operator++() {
  step();
  return *this;
}


/**
 * Sets the scaling scheme
 *
 * @param s New scaling scheme
 */
GAScalingScheme& EvolutionStrategy::scaling(const GAScalingScheme& s) {
  oldPop->scaling(s);
  return GAGeneticAlgorithm::scaling(s);
}


/**
 * Sets the selection scheme
 *
 * @param s New selection scheme
 */
GASelectionScheme& EvolutionStrategy::selector(const GASelectionScheme& s) {
  oldPop->selector(s);
  return GAGeneticAlgorithm::selector(s);
}


int EvolutionStrategy::IntermediateRecombination(const std::vector<GA1DArrayAlleleGenome<long double>*>& parents,
						 GA1DArrayAlleleGenome<long double>& child) const {

  // Recombine both endogenous parameters and variables
  for (int j = 0; j < child.size(); j++) {
    long double sum = 0.0;

    // First, genes
    for (unsigned k = 0; k < _ro; k++)
      sum += (1.0/_ro) * parents[k]->gene(j);

    child.gene(j, sum);

    sum = 0.0;

    // Then, Std Devs
    for (unsigned k = 0; k < _ro; k++)
      sum += (1.0/_ro) * parents[k]->stdDev(j);

    child.stdDev(j, sum);

  }

  return 1;

}

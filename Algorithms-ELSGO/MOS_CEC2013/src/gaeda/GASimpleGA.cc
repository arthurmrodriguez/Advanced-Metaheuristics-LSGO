/**
 * @file
 * @brief GASimpleGA class impl.
 *
 * Implementation of the GASimpleGA class
 */

#include "GASimpleGA.h"

#include "garandom.h"
#include "Recombinator.h"
#include "logger/GALogger.h"

GAParameterList& GASimpleGA::registerDefaultParameters(GAParameterList& p) {
	GAGeneticAlgorithm::registerDefaultParameters(p);

	p.add(gaNelitism, gaSNelitism,
			GAParameter::BOOLEAN, &DefaultParameters::elitism);

	return p;
}



GASimpleGA::GASimpleGA(const GAGenome& c) : GAGeneticAlgorithm(c){
	oldPop = pop->clone();
}

GASimpleGA::GASimpleGA(const GAPopulation& p) : GAGeneticAlgorithm(p){
	oldPop = pop->clone();
}

GASimpleGA::GASimpleGA(const GASimpleGA& ga) : GAGeneticAlgorithm(ga){
	oldPop = (GAPopulation *)0;
	copy(ga);
}
GASimpleGA::~GASimpleGA(){
	delete oldPop;
}
GASimpleGA&
GASimpleGA::operator=(const GASimpleGA& ga){
	if(&ga != this) copy(ga);
	return *this;
}
void
GASimpleGA::copy(const GAGeneticAlgorithm & g){
	GAGeneticAlgorithm::copy(g);
	const GASimpleGA& ga = DYN_CAST(const GASimpleGA&,g);

	if(oldPop) oldPop->copy(*(ga.oldPop));
	else oldPop = ga.oldPop->clone();
	oldPop->geneticAlgorithm(*this);
}


int
GASimpleGA::setptr(const char* name, const void* value){
	int status = GAGeneticAlgorithm::setptr(name, value);

	if(strcmp(name, gaNelitism) == 0 || strcmp(name, gaSNelitism) == 0){
		status = 0;
	}
	return status;
}

int
GASimpleGA::get(const char* name, void* value) const {
	int status = GAGeneticAlgorithm::get(name, value);

	if(strcmp(name, gaNelitism) == 0 || strcmp(name, gaSNelitism) == 0){
		status = 0;
	}
	return status;
}

void
GASimpleGA::objectiveFunction(GAGenome::Evaluator f){
	GAGeneticAlgorithm::objectiveFunction(f);
	for(unsigned i=0; i<pop->size(); i++)
		oldPop->individual(i).evaluator(f);
}

void
GASimpleGA::objectiveData(const GAEvalData& v){
	GAGeneticAlgorithm::objectiveData(v);
	for(unsigned i=0; i<pop->size(); i++)
		oldPop->individual(i).evalData(v);
}

const GAPopulation&
GASimpleGA::population(const GAPopulation& p) {
	GAGeneticAlgorithm::population(p);

	GAPopulation* tmpPop = pop->clone();
	oldPop->copy(*tmpPop);
	delete tmpPop;

	oldPop->geneticAlgorithm(*this);

	return *pop;
}

unsigned int
GASimpleGA::populationSize(unsigned int value) {
	GAGeneticAlgorithm::populationSize(value);
	oldPop->size(value);
	return value;
}

/**
 * Initialize the population, set the random seed as needed, do a few stupidity
 * checks, reset the stats.  We must initialize the old pop because there is no
 * guarantee that each individual will get initialized during the course of our
 * operator++ operations.  We do not evaluate the old pop because that will
 * happen as-needed later on.
 *
 * @param seed Random seed
 */
void GASimpleGA::initialize() {
	pop->initialize();
	pop->evaluate(gaTrue);	// the old pop will get it when the pops switch
	//  oldPop->initialize();

	stats.reset(*pop);

	if(!scross) GAErr(GA_LOC, className(), "initialize", gaErrNoSexualMating);

   printStats("Initial Stats");

}



void GASimpleGA::offspring (GAPopulation* offpop){

  unsigned i, mut, c1, c2;
  GAGenome *mom    = NULL;
  GAGenome *dad    = NULL;
  GAGenome *child1 = NULL;
  GAGenome *child2 = NULL;  // tmp holders for selected genomes

  // For updating the stats, we need to save the relationship between which parents generated which children
  vector<GAGenome*>  selected_parents  ( (offpop->size()%2 != 0) ? offpop->size()+1 : offpop->size() );

  // Generate the individuals in the temporary population from individuals in the main population.
  for(i=0; i<offpop->size()-1; i+=2){         // takes care of odd population

    mom    = &( pop->select()     );
    dad    = &( pop->select()     );

    stats.numsel += 2;    // keep track of number of selections

    child1 = &( offpop->individual( i ) );
    child2 = &( offpop->individual(i+1) );

    selected_parents[i]     = mom;
    selected_parents[i+1]   = dad;

    c1 = c2 = 0;

    if(GAFlipCoin(pCrossover())){
      stats.numcro += (*scross)( *mom, *dad, child1, child2 );
      c1 = c2 = 1;
    }
    else{
      child1->copy(*mom);
      child2->copy(*dad);
    }

    /*LOG*/ GALogger::instance()->appendOperatorResult("GASimpleGA: Crossovers (stats as age is not update yet)",*dad, *mom, *child1,*child2);

    /*LOG*/ GALogger::instance()->appendOperatorResult("GASimpleGA: Before applying mutation to these inds (stats as age is not update yet)", *child1, *child2);

    stats.nummut += ( mut = child1->mutate( pMutation()) );

    if(mut > 0) c1 = 1;

    stats.nummut += ( mut = child2->mutate( pMutation()) );

    if(mut > 0) c2 = 1;

    /*LOG*/ GALogger::instance()->appendOperatorResult("GASimpleGA: After applying mutation to these inds (stats as age is not update yet)",*child1,*child2);

    stats.numeval += c1 + c2;
  }

  if(offpop->size() % 2 != 0){ // do the remaining population member
    mom          =  &( pop->select()   );
    dad          =  &( pop->select()   );

    child1       =  &( offpop->individual(i) );
    stats.numsel += 2;                   // keep track of number of selections

    selected_parents[offpop->size()-1]   = mom; // In this case, the selected_parents population is bigger than the children one
    selected_parents[offpop->size()]     = dad;

    c1 = 0;
    if(GAFlipCoin(pCrossover())){
      stats.numcro += (*scross)(*mom, *dad, child1, (GAGenome*)0);
      c1 = 1;
      /*Cross LOG*/ GALogger::instance()->appendOperatorResult("GASimpleGA: Crossovers (stats as age is not update yet)",*mom,*dad,*child1);
    }
    else{
      if(GARandomBit()) {
        child1->copy(*mom);
      }
      else {
        child1->copy(*dad);
      }
    }

    GALogger::instance()->appendOperatorResult("GASimpleGA: Crossovers (stats as age is not update yet) ",
                                               *dad,*mom,
                                               *child1);

    /*Mutator LOG*/ GALogger::instance()->appendInd("GASimpleGA: Before applying mutation to the ind(stats as age is not update yet) ",*child1);

    stats.nummut += (mut = child1->mutate(pMutation()));
    if(mut > 0) c1 = 1;

    /*Mutator LOG*/ GALogger::instance()->appendInd("GASimpleGA: After applying mutation to the ind(stats as age is not update yet) ",*child1);

    stats.numeval += c1;
  }

  stats.numrep += pop->size();

  setBirthData(*offpop);

  stats.updateFitnessIncrement(*pop,selected_parents, *offpop);
}

/**
 * Evolve a new generation of genomes.  When we start this routine, pop contains the current generation.  When we finish,
 * pop contains the new generation and oldPop contains the (no longer) current generation.  The previous old generation is lost.
 * We don't deallocate any memory, we just reset the contents of the genomes. The selection routine must return a pointer to a
 * genome from the old population. Needs to be refactored into smaller functions!!!!
 */
void GASimpleGA::step(){

  /*LOG*/ GALogger::instance()->appendPopulation("GASimpleGAStep::step", "Before doing the step", *pop);

  GAPopulation *tmppop;		// Swap the old population with the new pop.
  // Generate the individuals in the temporary population from individuals in
  // the main population.

  offspring(oldPop);
  tmppop = oldPop;		// When we finish the ++ we want the newly
  oldPop = pop;			// generated population to be current (for
  pop = tmppop;			// references to it from member functions).


  GALogger::instance()->appendPopulation("The generated children is: ", "", *pop);

  recombinator_->recombine(*oldPop,*pop);

  pop->evaluate(gaTrue);  // get info about current pop for next time

  GALogger::instance()->appendPopulation("After recombination the population is: ", "", *pop);

  stats.update(*pop);   // update the statistics by one generation
  printStats("End of Step Stats");
}

/**
 * Some syntactic sugar for the step() method
 */
GASimpleGA & GASimpleGA::operator++() {
	step();
	return *this;
}

/**
 * Sets the scaling scheme
 *
 * @param s New scaling scheme
 */
GAScalingScheme& GASimpleGA::scaling(const GAScalingScheme & s) {
	oldPop->scaling(s);
	return GAGeneticAlgorithm::scaling(s);
}

/**
 * Sets the selection scheme
 *
 * @param s New selection scheme
 */
GASelectionScheme& GASimpleGA::selector(const GASelectionScheme& s) {
	oldPop->selector(s);
	return GAGeneticAlgorithm::selector(s);
}

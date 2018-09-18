#include "PopElitism.h"

#include "GAPopulation.h"
#include "GAGenealogy.h"

PopElitism::PopElitism() : Recombinator() {}

PopElitism::~PopElitism(){}

/*
 * It adds the first n individuals from the old population to the new population and then it resizes
 * to the population size. It assumes that the old population is sorted in descendant order being the
 * first individual the best one
 */
void PopElitism::recombine (const GAPopulation& old_pop, GAPopulation& new_pop){

  for (unsigned i = 0; i < old_pop.size(); i++)
    new_pop.add(old_pop.individual(i));

  new_pop.sort();

  for( unsigned i = 0; i < old_pop.size(); i++) {
    GAGenome* gen = new_pop.remove(GAPopulation::WORST, GAPopulation::RAW);
    // BEGIN: Genealogy
    if (GAGenealogy::handle())
      GAGenealogy::handle()->deceased(*gen);
    // END: Genealogy
    delete gen;
  }

}

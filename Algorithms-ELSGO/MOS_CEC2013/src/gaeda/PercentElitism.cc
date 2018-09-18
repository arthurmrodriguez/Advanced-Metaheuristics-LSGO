#include "PercentElitism.h"

#include "GAPopulation.h"
#include "GAGenealogy.h"

PercentElitism::PercentElitism(double percent) : Recombinator(), _percentage(percent) {}

PercentElitism::~PercentElitism(){}


void PercentElitism::recombine (const GAPopulation& old_pop, /*inout*/ GAPopulation& new_pop) {

  unsigned numInds = (unsigned) (old_pop.size() * _percentage);

  new_pop.sort();

  for( unsigned i = 0; i < numInds; i++)
    new_pop.add(old_pop.individual(i));

  new_pop.sort();

  for( unsigned i = 0; i < numInds; i++) {
    GAGenome* gen = new_pop.remove(GAPopulation::WORST, GAPopulation::RAW);
    // BEGIN: Genealogy
    if (GAGenealogy::handle())
      GAGenealogy::handle()->deceased(*gen);
    // END: Genealogy
    delete gen;
  }

}

#include "IncElitism.h"

#include "GAPopulation.h"

IncElitism::IncElitism() : Recombinator() {}

IncElitism::~IncElitism(){}

/*
 * It adds the first n individuals from the old population to the new population and then it resizes
 * to the population size. It assumes that the old population is sorted in descendant order being the
 * first individual the best one
 */
void IncElitism::recombine (const GAPopulation& old_pop, GAPopulation& new_pop){
  new_pop.copy(old_pop);
}

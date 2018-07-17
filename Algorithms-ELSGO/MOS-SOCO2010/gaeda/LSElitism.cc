#include "LSElitism.h"

#include "GAPopulation.h"

LSElitism::LSElitism() : Recombinator() {}

LSElitism::~LSElitism(){}

void LSElitism::recombine (const GAPopulation& old_pop, GAPopulation& new_pop){

  // We copy all the individuals from the original ppulation, except the first one,
  // which is the one that has been modified by the LS
  for (unsigned i = 1; i < new_pop.size(); i++)
    new_pop.individual(i).copy(old_pop.individual(i));

  new_pop.sort(gaTrue);

}

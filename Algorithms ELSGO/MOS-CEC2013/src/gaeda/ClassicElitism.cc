#include "ClassicElitism.h"

#include "GAPopulation.h"
#include "GAGenealogy.h"

ClassicElitism::ClassicElitism() : Recombinator() {}

ClassicElitism::~ClassicElitism(){}

void ClassicElitism::recombine (const GAPopulation& old_pop, /*inout*/ GAPopulation& new_pop){

  new_pop.sort();

  // BEGIN: Genealogy
  if (GAGenealogy::handle())
    GAGenealogy::handle()->deceased(new_pop.worst());
  // END: Genealogy

  new_pop.worst().copy(old_pop.best());
  new_pop.sort(gaTrue);

}

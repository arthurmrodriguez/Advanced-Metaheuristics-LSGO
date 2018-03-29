#include "GenerationalElitism.h"

#include "GAPopulation.h"
#include "GAGenealogy.h"

GenerationalElitism::GenerationalElitism() : Recombinator() {}

GenerationalElitism::~GenerationalElitism(){}

/*
 * It exchanges old and new populations.
 */
void GenerationalElitism::recombine (const GAPopulation& old_pop, GAPopulation& new_pop){

  for (unsigned i = 0; i < old_pop.size(); i++) {
    // BEGIN: Genealogy
    if (GAGenealogy::handle())
      GAGenealogy::handle()->deceased(old_pop.individual(i));
    // END: Genealogy
  }

}

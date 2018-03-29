#include "DEElitism.h"

#include "GAPopulation.h"

DEElitism::DEElitism() : Recombinator() {}

DEElitism::~DEElitism(){}

/*
 * 1 on 1 elitism. This kind of elitism is used with the differential evolution algorithm. An individual from the
 * old population replaces the one from the new population if is in the same position and has greater fitness.
 */
void DEElitism::recombine (const GAPopulation& old_pop, GAPopulation& new_pop){

  GAGenome* x_i   = NULL;
  GAGenome* x_new = NULL;
  for (unsigned i=0; i<new_pop.size(); i++){
    x_i   = (& old_pop.individual( i ) );
    x_new = (& new_pop.individual( i ) );
    if (x_i->hasBetterScoreThan(*x_new)) x_new->copy(*x_i);
  }
  new_pop.touch();
}

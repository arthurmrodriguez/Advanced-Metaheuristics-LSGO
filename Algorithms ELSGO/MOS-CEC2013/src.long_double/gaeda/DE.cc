/**
 * @file
 * @brief GASimpleGA class impl.
 *
 * Implementation of the DE class
 */

#include "math.h"
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "DE.h"

#include "garandom.h"
#include "Recombinator.h"
#include "logger/GALogger.h"
#include "genomes/GA1DArrayGenome.h"

//GAParameterList& DE::registerDefaultParameters(GAParameterList& p) {
//  GAGeneticAlgorithm::registerDefaultParameters(p);
//
//  return p;
//}

DE::DE(const GAGenome& c, long double F, long double CR) : GAGeneticAlgorithm(c){
  F_ = F;
  CR_ = CR;
  oldPop = pop->clone();
  deCross_ = c.de ();
}

DE::DE(const GAPopulation& p, long double F, long double CR) : GAGeneticAlgorithm(p){
  F_ = F;
  CR_ = CR;
  oldPop = pop->clone();
  deCross_ = p.individual (0).de ();
}

DE::DE(const DE& ga) : GAGeneticAlgorithm(ga){
  copy(ga);
}

DE::~DE(){
  delete oldPop;
}

DE& DE::operator=(const DE& alg){
  if(&alg != this) copy(alg);
  return *this;
}

void DE::copy(const DE& de){
  GAGeneticAlgorithm::copy(de);

  F_ = de.F_;
  CR_ = de.CR_;
  deCross_ = de.deCross_;

  if(oldPop) oldPop->copy(*(de.oldPop));
  else       oldPop = de.oldPop->clone();

  oldPop->geneticAlgorithm(*this);
}

void DE::objectiveFunction(GAGenome::Evaluator f){
  GAGeneticAlgorithm::objectiveFunction(f);
  for(unsigned i=0; i<pop->size(); i++) oldPop->individual(i).evaluator(f);
}

void DE::objectiveData(const GAEvalData& v){
  GAGeneticAlgorithm::objectiveData(v);
  for(unsigned i=0; i<pop->size(); i++) oldPop->individual(i).evalData(v);
}

const GAPopulation& DE::population(const GAPopulation& p) {
  GAGeneticAlgorithm::population(p);

  GAPopulation* tmpPop = pop->clone();
  oldPop->copy(*tmpPop);
  delete tmpPop;

  oldPop->geneticAlgorithm(*this);

  return *pop;
}

unsigned int DE::populationSize(unsigned int value) {
  GAGeneticAlgorithm::populationSize(value);
  oldPop->size(value);
  return value;
}

/**
 * Initialize parameters (commentary needs to be completed).
 *
 * @param seed Random seed
 */
void DE::initialize() {
  pop->initialize();
  pop->evaluate(gaTrue);

  stats.reset(*pop);

  printStats("Initial Stats");
}

void DE::step(){
  /*LOG*/ GALogger::instance()->appendPopulation("DEStep::step", "Before doing the step", *pop);

  GAPopulation *tmppop;   // Swap the old population with the new pop.


  tmppop = oldPop;        // When we finish the ++ we want the newly
  oldPop = pop;           // generated population to be current (for
  pop    = tmppop;        // references to it from member functions).
  offspring(pop);


  recombinator_->recombine(*oldPop,*pop);

  GALogger::instance()->appendPopulation("After recombination the population is: ", "", *pop);

  pop->evaluate(gaTrue);  // get info about current pop for next time
  stats.update(*pop);     // update the statistics by one generation

  printStats("End of Step Stats");
}

void DE::offspring(GAPopulation* pop){
  int           dim    = dynamic_cast<GA1DArrayAlleleGenome<long double>&>(pop->individual(0)).length();
  GA1DArrayAlleleGenome<long double>* t1     = NULL;
  GA1DArrayAlleleGenome<long double>* t2     = NULL;
  GA1DArrayAlleleGenome<long double>* t3     = NULL;
  GA1DArrayAlleleGenome<long double>* x_i    = NULL;
  GA1DArrayAlleleGenome<long double>* x_new  = NULL;

  //GARealGenome Vdonor    = new GARealGenome(dynamic_cast <GARealGenome&> (pop->individual(i)));
  //GARealGenome U         = new GARealGenome(dynamic_cast <GARealGenome&> (pop->individual(i)));

  for(unsigned ind_pos=0; ind_pos<pop->size(); ind_pos++){
    x_i   = dynamic_cast <GA1DArrayAlleleGenome<long double>*> (& oldPop->individual( ind_pos ) );
    x_new = dynamic_cast <GA1DArrayAlleleGenome<long double>*> (&    pop->individual( ind_pos ) );

    /*
     * TODO esto es una chapuza para arreglar error que hay en GASelector.
     * se debe corregir!!
     */
    unsigned attempts=0;
    do {
      t1 = dynamic_cast <GA1DArrayAlleleGenome<long double>*> (& oldPop->select() );
      attempts += 1;
      if (attempts >= pop->size() ){
        attempts=0;
        break;
      }
    }while (t1 == x_i);

    do {
       t1 = dynamic_cast <GA1DArrayAlleleGenome<long double>*> (& oldPop->individual(GARandomInt(0,pop->size()-1)) );
    }while (t1 == x_i);

    do{
      t2 = dynamic_cast <GA1DArrayAlleleGenome<long double>*> (& oldPop->select() );
      attempts+=1;
      if (attempts >= pop->size() ){
           attempts=0;
           break;
      }
    }while (t2 == x_i or t2 == t1);

    do{
       t2 = dynamic_cast <GA1DArrayAlleleGenome<long double>*> (& oldPop->individual(GARandomInt(0,pop->size()-1)) );
    }while (t2 == x_i or t2 == t1);

    do{
      t3 = dynamic_cast <GA1DArrayAlleleGenome<long double>*> (& oldPop->select() );
      attempts+=1;
      if (attempts >= pop->size() ){
           attempts=0;
           break;
         }
    }while (t3 == x_i or t3 == t2 or t3 == t1);

    do{
        t3 = dynamic_cast <GA1DArrayAlleleGenome<long double>*> (& oldPop->individual(GARandomInt(0,pop->size()-1)) );
    }while (t3 == x_i or t3 == t2 or t3 == t1);

    for(int gene_pos=0; gene_pos<dim; gene_pos++){
      long double new_val = t1->gene(gene_pos) + F_ * (t2->gene(gene_pos) - t3->gene(gene_pos));

      if (new_val > x_new->alleleset(gene_pos).upper())
        new_val = x_new->alleleset(gene_pos).upper();
      else if (new_val < x_new->alleleset(gene_pos).lower())
        new_val = x_new->alleleset(gene_pos).lower();

      x_new->gene (gene_pos, new_val);
    }

    stats.numsel += 4;
    stats.nummut += 1;
    stats.numcro += (deCross_)( *x_i, *x_new, x_new, CR_ );

  }

  stats.numrep  += pop->size();
  stats.numeval += pop->size();
}


/**
 * Sets the scaling scheme
 *
 * @param s New scaling scheme
 */
GAScalingScheme& DE::scaling(const GAScalingScheme & s) {
  oldPop->scaling(s);
  return GAGeneticAlgorithm::scaling(s);
}

/**
 * Sets the selection scheme
 *
 * @param s New selection scheme
 */
GASelectionScheme& DE::selector(const GASelectionScheme& s) {
  oldPop->selector(s);
  return GAGeneticAlgorithm::selector(s);
}




// $Header: /home/cvs/galib/ga/GAPopulation.C,v 1.5 2004/12/29 16:25:42 mwall Exp $
/* ----------------------------------------------------------------------------
   population.C
   mbwall 11aug94
   Copyright (c) 1995 Massachusetts Institute of Technology
   all rights reserved
   ---------------------------------------------------------------------------- */
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "GAPopulation.h"
#include "GASelector.h"
#include "garandom.h"
#include "GAGeneticAlgorithm.h"		// for the sake of flaky g++ compiler

#include <sstream>
#include "logger/GALogger.h"
#include <stdexcept>
#include "genomes/GA1DBinStrGenome.h"
#include "genomes/GARealGenome.h"

#include "GAGenealogy.h"
#include "GAGenealogyMemory.h"

// windows is promiscuous in its use of min/max, and that causes us grief.  so
// turn of the use of min/max macros in this file.   thanks nick wienholt
#if !defined(NOMINMAX)
#define NOMINMAX
#endif

// This is the default population initializer.  It simply calls the initializer
// for each member of the population.  Then we touch the population to tell it
// that it needs to update stats and/or sort (but we don't actually force
// either one to occur.
//   The population object takes care of setting/unsetting the status flags.
void
GAPopulation::DefaultInitializer(GAPopulation & p, double perc){
	for(unsigned int i=0; i<p.size(); i++) {
		p.individual(i).initialize();
	}
}

// Percentage Initializer: it only initializes a percentage of the population
void GAPopulation::PercInitializer(GAPopulation & p, double perc){
  // We must evaluate each genome one by one, as they have the score
  // value cached, their flag evaluated is gaTrue and p.evaluate(gaTrue);
  // does not enforce the evaluation of each genome under those circumstances
  for (unsigned i = 0; i < p.size(); i++)
    p.individual(i).evaluate(gaTrue);

  p.sort();

  unsigned initSize = p.size() * perc;
  int diff = p.size() - initSize;

  std::cout << "initSize: " << initSize << std::endl;

	for(int i=p.size() - 1; i>=diff; i--) {
		p.individual(i).initialize();
	}
}

//  The default evaluator simply calls the evaluate member of each genome in
// the population.  The population object takes care of setting/unsetting the
// status flags for indicating when the population needs to be updated again.
void
GAPopulation::DefaultEvaluator(GAPopulation & p){
	for(unsigned int i=0; i<p.size(); i++)
		p.individual(i).evaluate();
}



#define GA_POP_CHUNKSIZE 10	// allocate chrom ptrs in chunks of this many


/* ----------------------------------------------------------------------------
   Population

   The population class is basically just a holder for the genomes.  We also
   keep track of statistics about the fitness of our genomes.  We don't care
   what kind of genomes we get.  To create the population we call the clone
   method of the genome we're given.
   By default we do not calculate the population's diversity, so we set the
   div matrix to NULL.
   ---------------------------------------------------------------------------- */
GAPopulation::GAPopulation() {
	csz = N = GA_POP_CHUNKSIZE;
	n = 0;
	while(N < n) N += csz;

	rind = new GAGenome * [N];
	sind = new GAGenome * [N];
	memset(rind, 0, N * sizeof(GAGenome*));
	memset(sind, 0, N * sizeof(GAGenome*));
	//  indDiv = new double[N*N];
	indDiv = 0;

	neval = 0;
	rawSum = rawAve = rawDev = rawVar = rawMax = rawMin = 0.0;
	fitSum = fitAve = fitDev = fitVar = fitMax = fitMin = 0.0;
	popDiv = -1.0;
	rsorted = ssorted = evaluated = gaFalse;
	scaled = statted = divved = selectready = gaFalse;
	init = DefaultInitializer;
	eval = DefaultEvaluator;
	slct = new DEFAULT_SELECTOR;
	slct->assign(*this);
	sclscm = new DEFAULT_SCALING;
	evaldata = (GAEvalData*)0;
	ga = (GAGeneticAlgorithm*)0;

  initPerc = 1.0;
}

GAPopulation::GAPopulation(const GAGenome & c, unsigned int popsize) {

	csz = N = GA_POP_CHUNKSIZE;
	n = (popsize < 1 ? 1 : popsize);
	while(N < n) N += csz;

	rind = new GAGenome * [N];
	sind = new GAGenome * [N];

	for(unsigned int i=0; i<n; i++)
		rind[i] = c.clone(GAGenome::ATTRIBUTES);

	memcpy(sind, rind, N * sizeof(GAGenome*));
	//  indDiv = new double[N*N];
	indDiv = 0;

	neval = 0;
	rawSum = rawAve = rawDev = rawVar = rawMax = rawMin = 0.0;
	fitSum = fitAve = fitDev = fitVar = fitMax = fitMin = 0.0;
	popDiv = -1.0;
	rsorted = ssorted = evaluated = gaFalse;
	scaled = statted = divved = selectready = gaFalse;
	init = DefaultInitializer;
	eval = DefaultEvaluator;
	slct = new DEFAULT_SELECTOR;
	slct->assign(*this);
	sclscm = new DEFAULT_SCALING;
	evaldata = (GAEvalData*)0;
	ga = (GAGeneticAlgorithm*)0;

  initPerc = 1.0;
}

GAPopulation::GAPopulation(const GAPopulation & orig){
	n = N = 0;
	rind = sind = (GAGenome**)0;
	indDiv = (double*)0;
	sclscm = (GAScalingScheme*)0;
	slct = (GASelectionScheme*)0;
	evaldata = (GAEvalData*)0;
	copy(orig);
}

GAPopulation::GAPopulation(const vector<GAGenome*>& genomes){
  csz = N = GA_POP_CHUNKSIZE;
  n = ( genomes.size() < 1 ? 1 : genomes.size() );
  while(N < n) N += csz;

  rind = new GAGenome * [N];
  sind = new GAGenome * [N];

  for(unsigned int i=0; i<n; i++) rind[i] = genomes[i]->clone();

  memcpy(sind, rind, N * sizeof(GAGenome*));

  indDiv = 0;

  neval   = 0;
  rawSum = rawAve = rawDev = rawVar = rawMax = rawMin = 0.0;
  fitSum  = fitAve = fitDev = fitVar = fitMax = fitMin = 0.0;
  popDiv  = -1.0;
  rsorted = ssorted = evaluated = gaFalse;
  scaled  = statted = divved = selectready = gaFalse;
  init = DefaultInitializer;
  eval = DefaultEvaluator;
  slct = new DEFAULT_SELECTOR;
  slct->assign(*this);
  sclscm = new DEFAULT_SCALING;
  evaldata = (GAEvalData*)0;
  ga = (GAGeneticAlgorithm*)0;

  initPerc = 1.0;
}


GAPopulation::~GAPopulation(){
	for(unsigned int i=0; i<n; i++)
		delete rind[i];
	delete [] rind;
	delete [] sind;
	delete [] indDiv;
	delete sclscm;
	delete slct;
	delete evaldata;
}


// Make a complete copy of the original population.  This is a deep copy of
// the population object - we clone everything in the genomes and copy all of
// the population's information.
	void
GAPopulation::copy(const GAPopulation & arg)
{
	unsigned int i;
	for(i=0; i<n; i++)
		delete rind[i];
	delete [] rind;
	delete [] sind;
	delete [] indDiv;
	delete sclscm;
	delete slct;
	delete evaldata;

	csz = arg.csz; N = arg.N; n = arg.n;
	rind = new GAGenome * [N];
	for(i=0; i<n; i++)
		rind[i] = arg.rind[i]->clone();
	sind = new GAGenome * [N];
	memcpy(sind, rind, N * sizeof(GAGenome*));

	if(arg.indDiv) {
		indDiv = new double[N*N];
		memcpy(indDiv, arg.indDiv, (N*N*sizeof(double)));
	}
	else {
		indDiv = 0;
	}

	sclscm = arg.sclscm->clone();
	scaled = gaFalse;
	if(arg.scaled == gaTrue) scale();

	slct = arg.slct->clone();
	slct->assign(*this);
	selectready = gaFalse;
	if(arg.selectready == gaTrue) prepselect();

	if(arg.evaldata) evaldata = arg.evaldata->clone();
	else evaldata = (GAEvalData*)0;

	neval = 0;			// don't copy the evaluation count!
	rawSum = arg.rawSum; rawAve = arg.rawAve;
	rawMax = arg.rawMax; rawMin = arg.rawMin;
	rawVar = arg.rawVar; rawDev = arg.rawDev;
	popDiv = arg.popDiv;

	fitSum = arg.fitSum; fitAve = arg.fitAve;
	fitMax = arg.fitMax; fitMin = arg.fitMin;
	fitVar = arg.fitVar; fitDev = arg.fitDev;

	rsorted = arg.rsorted;
	ssorted = gaFalse;		// we must sort at some later point
	statted = arg.statted;
	evaluated = arg.evaluated;
	divved = arg.divved;

	init = arg.init;
	eval = arg.eval;
	ud = arg.ud;
	ga = arg.ga;

  initPerc = arg.initPerc;
}


// Resize the population.  If we shrink, we delete the extra genomes.  If
// we grow, we clone new ones (and we DO NOT initialize them!!!).  When we
// trash the genomes, we delete the worst of the population!  We do not
// free up the space used by the array of pointers, but we do free up the
// space used by the genomes.
//   We do a clone of the genome contents so that we don't have to initialize
// the new ones (what if the population has a custom initilizer?).  We randomly
// pick which ones to clone from the existing individuals.  If the population
// contains no genomes, then we post an error message (since there are no
// individuals from which to clone the new ones).
//   If the population was evaluated, then we evaluate the new genomes.  We
// do not sort nor restat the population, and we tag the statted and sorted
// flags to reflect the fact that they are no longer valid.
//   Resizing to a bigger size is the same as a batch 'add'
unsigned int
GAPopulation::size(unsigned int popsize){
	if(popsize == n) return n;
	if(n == 0 && popsize > 0) {
		GAErr(GA_LOC, "GAPopuluation", "size", gaErrNoIndividuals);
		return n;
	}

	if(popsize > n){
		grow(popsize);
		for(unsigned int i=n; i<popsize; i++)
			rind[i] = rind[GARandomInt(0,n-1)]->clone(GAGenome::CONTENTS);
		rsorted = gaFalse;
	}
	else{
		for(unsigned int i=popsize; i<n; i++) {// trash the worst ones (if sorted)
		  // BEGIN: Genealogy
		  // IMPORTANT: Decease the genomes that will be trashed
		  if (GAGenealogy::handle() != NULL)
		    GAGenealogy::handle()->deceased(*rind[i]);
		  // END: Genealogy
		  delete rind[i];			  // may not be sorted!!!!

		}
	}

	memcpy(sind, rind, N * sizeof(GAGenome*));
	ssorted = scaled = statted = divved = selectready = gaFalse;
	n = popsize;

	if(evaluated == gaTrue) evaluate(gaTrue);

	return n;
}


// This is a private method for adjusting the size of the arrays used by the
// population object.  Unlike the size method, this method does not allocate
// more genomes (but it will delete genomes if the specified size is smaller
// than the current size).
//   This maintains the integrity of the diversity scores (but the new ones
// will not have been set yet).
//   We return the total amount allocated (not the amount used).
int
GAPopulation::grow(unsigned int s) {
	if(s <= N) return N;

	int oldsize = N;
	while(N < s) N += csz;

	GAGenome ** tmp;

	tmp = rind;
	rind = new GAGenome * [N];
	memcpy(rind, tmp, oldsize*sizeof(GAGenome *));
	delete [] tmp;
	tmp = sind;
	sind = new GAGenome * [N];
	memcpy(sind, tmp, oldsize*sizeof(GAGenome *));
	delete [] tmp;

	if(indDiv) {
		double *tmpd = indDiv;
		indDiv = new double[N*N];
		for(int i=0; i<oldsize; i++)
			memcpy(&(indDiv[i*N]), &(tmpd[i*oldsize]), oldsize*sizeof(double));
		delete [] tmpd;
	}

	return N;
}


// Get rid of 'extra' memory that we have allocated.  We just trash the
// diversity matrix and flag it as being invalid.  Return the amount
// allocated (which is also the amount used).
	int
GAPopulation::compact()
{
	if(n == N) return N;

	GAGenome ** tmp;

	tmp = rind;
	rind = new GAGenome * [n];
	memcpy(rind, tmp, n*sizeof(GAGenome *));
	delete [] tmp;
	tmp = sind;
	sind = new GAGenome * [n];
	memcpy(sind, tmp, n*sizeof(GAGenome *));
	delete [] tmp;

	//  if(indDiv) {
	//    double *tmpd = indDiv;
	//    indDiv = new double [n*n];
	//    for(unsigned int i=0; i<n; i++)
	//      memcpy(&(indDiv[i*n]), &(tmpd[i*N]), n*sizeof(double));
	//    delete [] tmpd;
	//  }

	if(indDiv) {
		delete [] indDiv;
		indDiv = 0;
	}

	return N = n;
}


// Sort using the quicksort method.  The sort order depends on whether a high
// number means 'best' or a low number means 'best'.  Individual 0 is always
// the 'best' individual, Individual n-1 is always the 'worst'.
//   We may sort either array of individuals - the array sorted by raw scores
// or the array sorted by scaled scores.
void
GAPopulation::sort(GABoolean flag, SortBasis basis) const {
	GAPopulation * This = (GAPopulation *)this;
	if(basis == RAW){
		if(rsorted == gaFalse || flag == gaTrue){
         if(GAGenome::optCriterion() == GAGenome::MINIMIZATION)
				GAPopulation::QuickSortAscendingRaw(This->rind, 0, n-1);
			else
				GAPopulation::QuickSortDescendingRaw(This->rind, 0, n-1);
			This->selectready = gaFalse;
		}
		This->rsorted = gaTrue;
	}
	else if(basis == SCALED){
		if(ssorted == gaFalse || flag == gaTrue){
         if(GAGenome::optCriterion() == GAGenome::MINIMIZATION)
				GAPopulation::QuickSortAscendingScaled(This->sind, 0, n-1);
			else
				GAPopulation::QuickSortDescendingScaled(This->sind, 0, n-1);
			This->selectready = gaFalse;
		}
		This->ssorted = gaTrue;
	}
}


// Evaluate each member of the population and store basic population statistics
// in the member variables.  It is OK to run this on a const object - it
// changes to physical state of the population, but not the logical state.
//   The partial sums are normalized to the range [0,1] so that they can be
// used whether the population is sorted as low-is-best or high-is-best.
// Individual 0 is always the best individual, and the partial sums are
// calculated so that the worst individual has the smallest partial sum.  All
// of the partial sums add to 1.0.
void
GAPopulation::statistics(GABoolean flag) const {
	if(statted == gaTrue && flag != gaTrue) return;
	GAPopulation * This = (GAPopulation *)this;

	if(n > 0) {
		double tmpsum;
		This->rawMin = This->rawMax = tmpsum = rind[0]->score();

		unsigned int i;
		for(i=1; i<n; i++){
			double scr = rind[i]->score();
            tmpsum += scr;
            This->rawMax = GAMax(rawMax, scr);
            This->rawMin = GAMin(rawMin, scr);
		}
		double tmpave = tmpsum / n;
		This->rawAve = tmpave;
		This->rawSum = tmpsum;	// if scores are huge we'll lose data here

		double tmpvar = 0.0;
		if(n > 1){
			for(i=0; i<n; i++){
				double s = rind[i]->score() - This->rawAve;
				s *= s;
				tmpvar += s;
			}
			tmpvar /= (n-1);
		}
		This->rawDev = (double)sqrt(tmpvar);
		This->rawVar = tmpvar;	// could lose data if huge variance
	}
	else {
		This->rawMin = This->rawMax = This->rawSum = 0.0;
		This->rawDev = This->rawVar = 0.0;
	}

	This->statted = gaTrue;
}


// Do the scaling on the population.  Like the statistics and diversity, this
// method does not change the contents of the population, but it does change
// the values of the status members of the object.  So we allow it to work on
// a const population.
void
GAPopulation::scale(GABoolean flag) const {
	if(scaled == gaTrue && flag != gaTrue) return;
	GAPopulation* This = (GAPopulation*)this;

	if(n > 0) {
		sclscm->evaluate(*This);

		double tmpsum;
		This->fitMin = This->fitMax = tmpsum = sind[0]->fitness();

		unsigned int i;
		for(i=1; i<n; i++){
			tmpsum += sind[i]->fitness();
			This->fitMax = GAMax(fitMax, sind[i]->fitness());
			This->fitMin = GAMin(fitMin, sind[i]->fitness());
		}
		double tmpave = tmpsum / n;
		This->fitAve = tmpave;
		This->fitSum = tmpsum;	// if scores are huge we'll lose data here

		double tmpvar = 0.0;
		if(n > 1){
			for(i=0; i<n; i++){
				double s = sind[i]->fitness() - This->fitAve;
				s *= s;
				tmpvar += s;
			}
			tmpvar /= (n-1);
		}
		This->fitDev = (double)sqrt(tmpvar);
		This->fitVar = tmpvar;	// could lose data if huge variance
	}
	else {
		This->fitMin = This->fitMax = This->fitSum = 0.0;
		This->fitVar = This->fitDev = 0.0;
	}

	This->scaled = gaTrue;
	This->ssorted = gaFalse;
}


// Calculate the population's diversity score.  The matrix is triangular and
// we don't have to calculate the diagonals.  This assumes that div(i,j) is
// the same as div(j,i) (for our purposes this will always be true, but it is
// possible for someone to override some of the individuals in the population
// and not others).
//   For now we keep twice as many diversity numbers as we need.  We need only
// n*(n-1)/2, but I can't seem to figure out an efficient way to map i,j to the
// reduced n*(n-1)/2 set (remember that the diagonals are always 0.0).
//   The diversity of the entire population is just the average of all the
// individual diversities.  So if every individual is completely different from
// all of the others, the population diversity is > 0.  If they are all the
// same, the diversity is 0.0.  We don't count the diagonals for the population
// diversity measure.  0 means minimal diversity means all the same.
void
GAPopulation::diversity(GABoolean flag) const {
	if(divved == gaTrue && flag != gaTrue) return;
	GAPopulation* This = (GAPopulation*)this;

	if(n > 1) {
		if(This->indDiv == 0) This->indDiv = new double[N*N];

		This->popDiv = 0.0;
		for(unsigned int i=0; i<n; i++){
			This->indDiv[i*n+i] = 0.0;
			for(unsigned int j=i+1; j<n; j++){
				This->indDiv[j*n+i] = This->indDiv[i*n+j] =
					individual(i).compare(individual(j));
				This->popDiv += indDiv[i*n+j];
			}
		}
		This->popDiv /= n*(n-1)/2;
	}
	else {
		This->popDiv = 0.0;
	}

	This->divved = gaTrue;
}


void
GAPopulation::prepselect(GABoolean flag) const {
	if(selectready == gaTrue && flag != gaTrue) return;
	GAPopulation* This = (GAPopulation*)this;
	This->slct->update();
	This->selectready = gaTrue;
}


// Return a reference to the scaling object.
GAScalingScheme &
GAPopulation::scaling(const GAScalingScheme& s){
	delete sclscm;
	sclscm = s.clone();
	scaled = gaFalse;
	return *sclscm;
}


// Return a reference to the selection object.
GASelectionScheme&
GAPopulation::selector(const GASelectionScheme& s) {
	delete slct;
	slct = s.clone();
	slct->assign(*this);
	selectready = gaFalse;
	return *slct;
}



// Replace the specified genome with the one that is passed to us then
// return the one that got replaced.  Use the replacement flags to determine
// which genome will be replaced.  If we get a genome as the second
// argument, then replace that one.  If we get a NULL genome, then we
// return a NULL and don't do anything.
//   If the population is sorted, then we maintain the sort by doing a smart
// replacement.
//   If the population is not sorted, then we just do the replacement without
// worrying about the sort.  Replace best and worst both require that we know
// which chromsomes are which, so we do a sort before we do the replacement,
// then we do a smart replacement.
//   In both cases we flag the stats as out-of-date, but we do not update the
// stats.  Let that happen when it needs to happen.
//   If which is < 0 then it is a flag that tells us to do a certain kind of
// replacement.  Anything non-negative is assumed to be an index to a
// genome in the population.
//   This does not affect the state of the evaluated member - it assumes that
// the individual genome has a valid number for its score.
	GAGenome *
GAPopulation::replace(GAGenome * repl, int which, SortBasis basis)
{
	int i=-1;
	GAGenome * orig=(GAGenome *)0;
	if(repl == (GAGenome *)0) return orig;

	switch(which){
		case BEST:
			sort(gaFalse, basis);
			i = 0;
			break;

		case WORST:
			sort(gaFalse, basis);
			i = n-1;
			break;

		case RANDOM:
			i = GARandomInt(0, n-1);
			break;

		default:
			if(0 <= which && which < (int)n)
				i = which;
			break;
	}

	if(i >= 0){
		// We could insert this properly if the population is sorted, but that would
		// require us to evaluate the genome, and we don't want to do that 'cause that
		// will screw up any parallel implementations.  So we just stick it in the
		// population and let the sort take care of it at a later time as needed.
		if(basis == RAW){
			orig = rind[i];		// keep the original to return at the end
			rind[i] = repl;
			memcpy(sind, rind, N * sizeof(GAGenome*));
		}
		else{
			orig = sind[i];		// keep the original to return at the end
			sind[i] = repl;
			memcpy(rind, sind, N * sizeof(GAGenome*));
		}
		rsorted = ssorted = gaFalse;	// must sort again
		// flag for recalculate stats
		statted = gaFalse;
		// Must flag for a new evaluation.
		evaluated = gaFalse;
		// No way to do incremental update of scaling info since we don't know what the
		// scaling object will do.
		scaled = gaFalse;
		// *** should do an incremental update of the diversity here so we don't
		// recalculate all of the diversities when only one is updated
		divved = gaFalse;
		// selector needs update
		selectready = gaFalse;

		// make sure the genome has the correct genetic algorithm pointer
		if(ga) repl->geneticAlgorithm(*ga);
	}

	return orig;
}


// Replace the genome o in the population with the genome r.  Return a
// pointer to the original genome, o.  This assumes that o exists in the
// population.   If it does not, we return a NULL.  If the genomes are the
// same, do nothing and return a pointer to the genome.
	GAGenome *
GAPopulation::replace(GAGenome * r, GAGenome * o)
{
	GAGenome * orig=(GAGenome *)0;
	if(r == (GAGenome *)0 || o == (GAGenome *)0) return orig;
	if(r == o) return r;
	unsigned int i;
	for(i=0; i<n && rind[i] != o; i++);
	if(i < n) orig = replace(r, i, RAW);
	return orig;
}


//   Remove the xth genome from the population.  If index is out of bounds, we
// return NULL.  Otherwise we return a pointer to the genome that was
// removed.  The population is now no longer responsible for freeing the
// memory used by that genome.
//   We don't touch the sorted flag for the array we modify - a remove will not
// affect the sort order.
	GAGenome *
GAPopulation::remove(int i, SortBasis basis)
{
	GAGenome * removed=(GAGenome *)0;
	if(i == BEST) { sort(gaFalse, basis); i = 0; }
	else if(i == WORST) { sort(gaFalse, basis); i = n-1; }
	else if(i == RANDOM) i = GARandomInt(0,n-1);
	else if(i < 0 || i >= (int)n) return removed;

	if(basis == RAW){
		removed = rind[i];
		memmove(&(rind[i]), &(rind[i+1]), (n-i-1)*sizeof(GAGenome *));
		memcpy(sind, rind, N * sizeof(GAGenome*));
		ssorted = gaFalse;
	}
	else if(basis == SCALED){
		removed = sind[i];
		memmove(&(sind[i]), &(sind[i+1]), (n-i-1)*sizeof(GAGenome *));
		memcpy(rind, sind, N * sizeof(GAGenome*));
		rsorted = gaFalse;
	}
	else return removed;

	n--;
	evaluated = gaFalse;

	// *** should be smart about these and do incremental update?
	scaled = statted = divved = selectready = gaFalse;

	return removed;
}


//   Remove the specified genome from the population.  If the genome is
// not in the population, we return NULL.  We do a linear search here (yuk for
// large pops, but little else we can do).  The memory used by the genome is
// now the responsibility of the caller.
	GAGenome *
GAPopulation::remove(GAGenome * r)
{
	GAGenome * removed=(GAGenome *)0;
	if(r == (GAGenome *)0) return removed;
	unsigned int i;
	for(i=0; i<n && rind[i] != r; i++);
	if(i < n) removed = remove(i, RAW);
	return removed;
}


// Add the specified individual to the population.  We don't update the stats
// or sort - let those get updated next time they are needed.
//   Notice that it is possible to add individuals to the population that are
// not the same type as the other genomes in the population.  Eventually we
// probably won't allow this (or at least we'll have to fix things so that the
// behaviour is completely defined).
//   If you invoke the add with a genome reference, the population will make
// a clone of the genome then it owns it from then on.  If you invoke add with
// a genome pointer, then the population does not allocate any memory - it uses
// the memory pointed to by the argument.  So don't trash the genome without
// first letting the population know about the change.
GAGenome* GAPopulation::add(const GAGenome& g)
{
  return GAPopulation::add(g.clone());
}


double GAPopulation::SumOfDistanceToPop(unsigned ind_pos) const {
  double result = 0.0;
  for (unsigned i=0; i< (unsigned) size(); i++){
    if ( i==ind_pos )   continue;
    result += indDiv[ind_pos*n+i];
  }
  return result;
}


// This one does *not* allocate space for the genome - it uses the one that
// was passed to us.  So the caller should not free it up or leave it dangling!
// We own it from now on (unless remove is called on it), and the population
// will destroy it when the population destructor is invoked.
	GAGenome*
GAPopulation::add(GAGenome* c)
{
	if(c == (GAGenome *)0) return c;
	grow(n+1);
	rind[n] = sind[n] = c;
	if(ga) rind[n]->geneticAlgorithm(*ga);
	n++;

	rsorted = ssorted = gaFalse;	// may or may not be true, but must be sure
	evaluated = scaled = statted = divved = selectready = gaFalse;

	return c;
}

void GAPopulation::add(const GAGenome& g, unsigned size)
{
  for (unsigned i=0; i<size; i++)
    GAPopulation::add(g.clone());
}



GAGeneticAlgorithm *
GAPopulation::geneticAlgorithm(GAGeneticAlgorithm& g){
	for(unsigned int i=0; i<n; i++)
		rind[i]->geneticAlgorithm(g);
	return(ga = &g);
}


#ifdef GALIB_USE_STREAMS
void
GAPopulation::write(STD_OSTREAM & os, SortBasis basis) const {
	for(unsigned int i=0; i<n; i++){
		if(basis == RAW)
			os << *rind[i] << "\n";
		else
			os << *sind[i] << "\n";
	}
	os << "\n";
}
#endif




/*
bool goesFirstBeforeSecond(const GAPopulation::SortOrder sort_order,double v1, double v2){
  if (sort_order == GAPopulation::LOW_IS_BEST) return v1 < v2;
  else                                         return v1 > v2;
}
*/


void
GAPopulation::QuickSortAscendingRaw(GAGenome **c, int l, int r) {
	int i,j; double v; GAGenome *t;
	if(r > l){
		v = c[r]->score(); i = l-1; j = r;
		for(;;){
			while(c[++i]->score() < v && i <= r);
			while(c[--j]->score() > v && j > 0);
			if(i >= j) break;
			t = c[i]; c[i] = c[j]; c[j] = t;
		}
		t = c[i]; c[i] = c[r]; c[r] = t;
		GAPopulation::QuickSortAscendingRaw(c,l,i-1);
		GAPopulation::QuickSortAscendingRaw(c,i+1,r);
	}
}
void
GAPopulation::QuickSortDescendingRaw(GAGenome **c, int l, int r) {
	int i,j; double v; GAGenome *t;
	if(r > l){
		v = c[r]->score(); i = l-1; j = r;
		for(;;){
			while(c[++i]->score() > v && i <= r);
			while(c[--j]->score() < v && j > 0);
			if(i >= j) break;
			t = c[i]; c[i] = c[j]; c[j] = t;
		}
		t = c[i]; c[i] = c[r]; c[r] = t;
		GAPopulation::QuickSortDescendingRaw(c,l,i-1);
		GAPopulation::QuickSortDescendingRaw(c,i+1,r);
	}
}

void
GAPopulation::QuickSortAscendingScaled(GAGenome **c, int l, int r) {
	int i,j; double v; GAGenome *t;
	if(r > l){
		v = c[r]->fitness(); i = l-1; j = r;
		for(;;){
			while(c[++i]->fitness() < v && i <= r);
			while(c[--j]->fitness() > v && j > 0);
			if(i >= j) break;
			t = c[i]; c[i] = c[j]; c[j] = t;
		}
		t = c[i]; c[i] = c[r]; c[r] = t;
		GAPopulation::QuickSortAscendingScaled(c,l,i-1);
		GAPopulation::QuickSortAscendingScaled(c,i+1,r);
	}
}
void
GAPopulation::QuickSortDescendingScaled(GAGenome **c, int l, int r) {
	int i,j; double v; GAGenome *t;
	if(r > l){
		v = c[r]->fitness(); i = l-1; j = r;
		for(;;){
			while(c[++i]->fitness() > v && i <= r);
			while(c[--j]->fitness() < v && j > 0);
			if(i >= j) break;
			t = c[i]; c[i] = c[j]; c[j] = t;
		}
		t = c[i]; c[i] = c[r]; c[r] = t;
		GAPopulation::QuickSortDescendingScaled(c,l,i-1);
		GAPopulation::QuickSortDescendingScaled(c,i+1,r);
	}
}

/*
void GAPopulation::quickSortPopRec(const SortOrder sort_order, GAGenome** all_inds, int left_pos, int right_pos){
  int        i,j;
  GAGenome*  tmp_ind;
  double     pivot_value;
  double     pivot_age;

  if(right_pos > left_pos){
    pivot_value = all_inds[right_pos]->score();
    pivot_age   = all_inds[right_pos]->age();
    i           = left_pos;
    j           = right_pos-1;

     while(true){
       while ( i <= right_pos && ( goesFirstBeforeSecond(sort_order,all_inds[i]->score(),pivot_value) ||
                                   ( all_inds[i]->score() == pivot_value &&
                                     goesFirstBeforeSecond(sort_order,all_inds[i]->age(),pivot_age) ) ) ) {
         i++;
       }
       while ( j>0 && ( goesFirstBeforeSecond(sort_order,pivot_value,all_inds[j]->score()) ||
                        (all_inds[j]->score() == pivot_value &&
                         goesFirstBeforeSecond(sort_order,pivot_age,all_inds[j]->age()) ) ) ) {
         j--;
       }
       if ( i>=j ) break;

       tmp_ind = all_inds[i]; all_inds[i] = all_inds[j]; all_inds[j] = tmp_ind; // Swap between i and j
       i++; j--;
     }

     tmp_ind = all_inds[i]; all_inds[i] = all_inds[right_pos]; all_inds[right_pos] = tmp_ind;

     GAPopulation::quickSortPopRec(sort_order,all_inds,left_pos,i-1);
     GAPopulation::quickSortPopRec(sort_order,all_inds,i+1,right_pos);
   }
}
*/


GAGenome& GAPopulation::medoid() const {
  assert( size() > 0 );

  diversity();

  // We calculate the individual which has the lesser sum of distances to the rest of the population
  GAGenome* medoid          = &( individual(0) );                                   // We initialize our initial value with the first individual
  double    medoid_distance = SumOfDistanceToPop(0);
  double    tmp_distance    = 0.0;

  for (unsigned i=1; i< (unsigned) size(); i++){
    tmp_distance = SumOfDistanceToPop(i);
    if (tmp_distance < medoid_distance){
      medoid_distance = tmp_distance;
      medoid          = &( individual(i) );
    }
  }

  { //LOG
    stringstream log_messg;
    log_messg << "The new medoid obtained from the population is " << *medoid << endl;
    GALogger::instance()->appendLogMessage("Sync Dynamic Neighbor Model", log_messg.str());
  }
  return *medoid;
}


double GAPopulation::entropy() const{
  GA1DBinaryStringGenome* gen_tmp = dynamic_cast<GA1DBinaryStringGenome*>( &( best() ) );

  if (!gen_tmp) throw runtime_error("In EntropyLogStat, genomes are not binary");

  unsigned pop_size = size();

  double tmp_num_ones, tmp_num_zeros, entropy = 0.0;

  int num_genes = gen_tmp->length();
  for (int gen_pos=0; gen_pos<num_genes; gen_pos++){
    tmp_num_ones = 0;
    for (unsigned ind_pos=0; ind_pos<pop_size; ind_pos++){
      gen_tmp = dynamic_cast<GA1DBinaryStringGenome*>( &( individual(ind_pos) ) );
      if (gen_tmp->gene(gen_pos) == 1 ) tmp_num_ones++;
    }
    tmp_num_zeros = pop_size - tmp_num_ones;

    // with no ones or zeros the logarithm returns nan so we dont consider this case
    if (tmp_num_zeros != 0 && tmp_num_ones != 0) {
      entropy += -tmp_num_ones/pop_size * log2(tmp_num_ones/pop_size) - tmp_num_zeros/pop_size * log2(tmp_num_zeros/pop_size);
    }
  }
  entropy /= num_genes;

  assert(!isnan(entropy));

  return entropy;
}

double GAPopulation::distAvgPoint() const{

  // The dynamic casting part was initially done right with the dynamic_cast operator. I dont know why but with the
  // dlopen libraries the dynamic_cast operator does not identigy the GARealGenome instances. Aa hack has been used
  // to fix this problem: we only check
  // if the genome is from a binary class, otherwise we assume that it derives the GA1DArrayGenome so therefore we can
  // do a casting that derives from this class (eg. GARealGenome) and call the apropiate method
  GA1DArrayGenome<double>* gen = dynamic_cast<GA1DArrayGenome<double>* >(& best());
  unsigned num_genes = gen->length();

  unsigned pop_size  = size();

  double sum, sqr_sum, avg, var, devs_sum, ind_gene;

  devs_sum = 0.0;

  GA1DBinaryStringGenome* gen_binary_tmp;

  for (unsigned gen_pos=0; gen_pos<num_genes; gen_pos++){
    sum      = 0.0;
    sqr_sum  = 0.0;
    avg      = 0.0;

    for (unsigned ind_pos=0; ind_pos<pop_size; ind_pos++){
      gen_binary_tmp = dynamic_cast<GA1DBinaryStringGenome*>( &( individual(ind_pos) ) );

      if      ( gen_binary_tmp ) ind_gene = gen_binary_tmp->bit(gen_pos);
      else {
        GA1DArrayGenome<double>* tmp_ind = dynamic_cast<GA1DArrayGenome<double>* >(& individual(ind_pos) );
        ind_gene = tmp_ind->gene(gen_pos);
      }

      sum     += ind_gene;
      sqr_sum += pow(ind_gene,2.0);
    }

    avg = sum / pop_size;
    var = (sqr_sum /( pop_size - 1.0) ) - ( ( (double)pop_size/(pop_size - 1.0) ) * pow(avg,2.0) );

    if (var < 0.0 && var > -0.0000001 ) var = 0.0;

    assert(!isnan(var)); assert(var >= 0 );

    devs_sum += sqrt(var);
  }
  return devs_sum/num_genes;
}

double GAPopulation::grefenstetteBias() const {

  GA1DBinaryStringGenome* gen_tmp = dynamic_cast<GA1DBinaryStringGenome*>( &( best() ) );

  if (!gen_tmp) throw runtime_error("In GrefenstetteLogStat, genomes are not binary");

  unsigned pop_size = size();

  double sum_genes, sum_compl_genes, gref_bias = 0.0;

  int num_genes = gen_tmp->length();
  for (int gen_pos=0; gen_pos<num_genes; gen_pos++){
    sum_genes = sum_compl_genes = 0;
    for (unsigned ind_pos=0; ind_pos<pop_size; ind_pos++){
      gen_tmp = dynamic_cast<GA1DBinaryStringGenome*>( &( individual(ind_pos) ) );

      sum_genes       += gen_tmp->gene(gen_pos);
      sum_compl_genes += 1-gen_tmp->gene(gen_pos);
    }
    gref_bias += (sum_genes > sum_compl_genes) ? sum_genes : sum_compl_genes;
  }
  gref_bias /= (num_genes * pop_size);

  gref_bias = 2*(1-gref_bias); // Corrected version so its range is in (0,1)

  assert(!isnan(gref_bias));
  assert(gref_bias <= 1.0 && gref_bias >= 0.0);

  return gref_bias;
}

double GAPopulation::nativePrcnt(unsigned rank) const {
  double   value = 0.0;
  unsigned pop_size = size();

  for (unsigned i=0; i<pop_size; i++) value  += ( individual(i).origin() == (int) rank ) ? 1.0 : 0.0;

  value  /= pop_size;

  return value;
}


double GAPopulation::avgGenDiv() const {
    diversity();                            // Computes (if not already) the distance between the individuals
      double* gen_dists = (const_cast<GAPopulation* > (this))->genDistsMatrix();

        double sum = 0.0;

          double   tmp;
            unsigned pop_size = size();

              for (unsigned i=0; i<pop_size; i++){
                    for (unsigned j=i+1; j<pop_size; j++){
                            tmp     =  gen_dists[i*pop_size+j];
                                  sum     += tmp;
                                      }
                      }

                int ncomparations = pop_size*(pop_size-1)/2;

                  return sum / ncomparations;
}


// BEGIN: Genealogy
void GAPopulation::computeAge() {

  int actgen, i, npop = (int) n;

  if(GAGenealogy::handle() != NULL)
    if(GAGenealogy::handle()->isGenealogyMemory()) {

      GAGenealogyMemory *genealmem = DYN_CAST(GAGenealogyMemory*, GAGenealogy::handle());

      actgen = genealmem->getGeneration();

      for(i=0; i<npop; i++) {
	GAGenomeNode *genNode = genealmem->getGeneNode(individual(i));
	individual(i).setAge(actgen - genNode->getFirstG());
      }

    }

}
// END: Genealogy

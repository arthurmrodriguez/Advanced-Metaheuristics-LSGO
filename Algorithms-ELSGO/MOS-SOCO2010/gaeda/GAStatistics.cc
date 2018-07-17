// $Header: /home/cvs/galib/ga/GAStatistics.C,v 1.2 2004/12/28 00:12:12 mwall Exp $
/* ----------------------------------------------------------------------------
  statistics.C
  mbwall 28jul94
  Copyright (c) 1995 Massachusetts Institute of Technology
                     all rights reserved

 DESCRIPTION:
  Definition of the statistics object.
---------------------------------------------------------------------------- */
#include <string.h>

#include "GAStatistics.h"

#include "gaerror.h"
#include "logger/GALogger.h"

const unsigned HFAME_INITIAL_SIZE = 500;   // Chunk size of the hall of fame population. We set a high number

// Default settings and their names.
int  gaDefNumBestGenomes   = 1;
int  gaDefScoreFrequency1  = 1;
int  gaDefScoreFrequency2  = 1; //TODO: Why this was 100???
int  gaDefFlushFrequency   = 0;
char gaDefScoreFilename[]  = "generations.dat";


GAStatistics::GAStatistics() {
  curgen = 0;
  numsel = numcro = nummut = numrep = numeval = numpeval = 0;
  maxever = minever = 0.0;
  on = offmax = offmin = 0.0;
  aveInit = maxInit = minInit = devInit = 0.0;
  divInit = -1.0;

  //score_avg_cur_ = score_max_cur_ = score_min_cur_ = score_dev_cur_ = 0.0;

  aveCur = maxCur = minCur = devCur = 0.0;
  divCur = -1.0;

  scoreFreq = gaDefScoreFrequency1;
  dodiv = gaFalse;		// default is do not calculate diversity

  nconv=0;
  Nconv = 10;
  cscore = new long double[Nconv]; memset(cscore, 0, Nconv*sizeof(long double));

  nscrs=0;
  Nscrs = gaDefFlushFrequency;
  gen = new int[Nscrs]; memset(gen, 0, Nscrs*sizeof(int));
  aveScore = new long double[Nscrs]; memset(aveScore, 0, Nscrs*sizeof(long double));
  maxScore = new long double[Nscrs]; memset(maxScore, 0, Nscrs*sizeof(long double));
  minScore = new long double[Nscrs]; memset(minScore, 0, Nscrs*sizeof(long double));
  devScore = new long double[Nscrs]; memset(devScore, 0, Nscrs*sizeof(long double));
  divScore = new long double[Nscrs]; memset(divScore, 0, Nscrs*sizeof(long double));
  scorefile = new char[strlen(gaDefScoreFilename)+1];
  strcpy(scorefile, gaDefScoreFilename);
  which = Maximum;

  fit_inc_       = 0.0;
  num_new_child_ = 0;

  boa = (GAPopulation *)0;

  hall_of_fame_ = new GAPopulation();
  hall_of_fame_->chunksize(HFAME_INITIAL_SIZE);
  use_hall_of_fame_ = false;
}

GAStatistics::GAStatistics(const GAStatistics & orig){
  cscore=(long double *)0;
  gen=(int *)0;
  aveScore=(long double *)0;
  maxScore=(long double *)0;
  minScore=(long double *)0;
  devScore=(long double *)0;
  divScore=(long double *)0;
  scorefile=(char *)0;
  boa=(GAPopulation *)0;

  hall_of_fame_ = new GAPopulation();
  hall_of_fame_->chunksize(HFAME_INITIAL_SIZE);
  use_hall_of_fame_ = false;

  copy(orig);
}
GAStatistics::~GAStatistics(){
  delete [] cscore;
  delete [] gen;
  delete [] aveScore;
  delete [] maxScore;
  delete [] minScore;
  delete [] devScore;
  delete [] divScore;
  delete [] scorefile;
  delete hall_of_fame_;
  delete boa;
}
void
GAStatistics::copy(const GAStatistics & orig){
  if(&orig == this) return;

  curgen = orig.curgen;
  numsel = orig.numsel;
  numcro = orig.numcro;
  nummut = orig.nummut;
  numrep = orig.numrep;
  numeval = orig.numeval;
  numpeval = orig.numpeval;
  maxever = orig.maxever;
  minever = orig.minever;
  on = orig.on;
  offmax = orig.offmax;
  offmin = orig.offmin;
  aveInit = orig.aveInit;
  maxInit = orig.maxInit;
  minInit = orig.minInit;
  devInit = orig.devInit;
  divInit = orig.divInit;
  aveCur = orig.aveCur;
  maxCur = orig.maxCur;
  minCur = orig.minCur;
  devCur = orig.devCur;
  divCur = orig.divCur;

  scoreFreq = orig.scoreFreq;
  dodiv = orig.dodiv;

  nconv=orig.nconv; Nconv=orig.Nconv;
  delete [] cscore;
  cscore = new long double [Nconv]; memcpy(cscore, orig.cscore, Nconv*sizeof(long double));

  nscrs=orig.nscrs; Nscrs=orig.Nscrs;
  delete [] gen;
  gen = new int [Nscrs];
  memcpy(gen, orig.gen, Nscrs*sizeof(int));
  delete [] aveScore;
  aveScore = new long double [Nscrs];
  memcpy(aveScore, orig.aveScore, Nscrs*sizeof(long double));
  delete [] maxScore;
  maxScore = new long double [Nscrs];
  memcpy(maxScore, orig.maxScore, Nscrs*sizeof(long double));
  delete [] minScore;
  minScore = new long double [Nscrs];
  memcpy(minScore, orig.minScore, Nscrs*sizeof(long double));
  delete [] devScore;
  devScore = new long double [Nscrs];
  memcpy(devScore, orig.devScore, Nscrs*sizeof(long double));
  delete [] divScore;
  divScore = new long double [Nscrs];
  memcpy(divScore, orig.divScore, Nscrs*sizeof(long double));

  delete [] scorefile;
  if(orig.scorefile){
    scorefile = new char [strlen(orig.scorefile)+1];
    strcpy(scorefile, orig.scorefile);
  }
  else scorefile = (char*)0;

  which = orig.which;

  delete boa;
  if(orig.boa) boa = orig.boa->clone();

  delete hall_of_fame_;
  if (orig.hall_of_fame_) hall_of_fame_ = orig.hall_of_fame_->clone();
  use_hall_of_fame_ = orig.use_hall_of_fame_;

  fit_inc_       = orig.fit_inc_;
  num_new_child_ = orig.num_new_child_;
}


// Update the genomes in the 'best of all' population to reflect any
// changes made to the current population.  We just grab the genomes with
// the highest scores from the current population, and if they are higher than
// those of the genomes in the boa population, they get copied.  Note that
// the bigger the boa array, the bigger your running performance hit because
// we have to look through all of the boa to figure out which are better than
// those in the population.  The fastest way to use the boa is to keep only
// one genome in the boa population.  A flag of 'True' will reset the boa
// population so that it is filled with the best of the current population.
//   Unfortunately it could take a long time to update the boa array using the
// copy method.  We'd like to simply keep pointers to the best genomes, but
// the genomes change from generation to generation, so we can't depend on
// that.
//   Notice that keeping boa is useful even for overlapping populations.  The
// boa keeps individuals that are different from each other - the overlapping
// population may not.  However, keeping boa is most useful for populations
// with little overlap.
//   When we check to see if a potentially better member is already in our
// best-of-all population, we use the operator== comparator not the genome
// comparator to do the comparison.
void
GAStatistics::updateBestIndividual(const GAPopulation & pop, GABoolean flag){
  if(boa == (GAPopulation *)0 || boa->size() == 0) return; // do nothing

  if(flag == gaTrue){		// reset the BOA array
    unsigned int j=0;
    for(unsigned int i=0; i<boa->size(); i++){
      boa->best(i).copy(pop.best(j));
      if(j < pop.size()-1) j++;
    }
    return;
  }

  if(boa->size() == 1){		// there's only one boa so replace it with bop
   if(GAGenome::optCriterion() == GAGenome::MAXIMIZATION &&
       pop.best().score() > boa->best().score())
      boa->best().copy(pop.best());
   if(GAGenome::optCriterion() == GAGenome::MINIMIZATION &&
       pop.best().score() < boa->best().score())
      boa->best().copy(pop.best());
  }
  else{
    unsigned int i=0, j, k;
   if(GAGenome::optCriterion() == GAGenome::MAXIMIZATION) {
      while(i < pop.size() && pop.best(i).score() > boa->worst().score()){
	for(k=0;
	    pop.best(i).score() < boa->best(k).score() && k < boa->size();
	    k++);
	for(j=k; j<boa->size(); j++){
	  if(pop.best(i) == boa->best(j)) break;
	  if(pop.best(i).score() > boa->best(j).score()){
	    boa->worst().copy(pop.best(i));        // replace worst individual
	    boa->sort(gaTrue, GAPopulation::RAW);  // re-sort the population
	    break;
	  }
	}
	i++;
      }
    }
   if(GAGenome::optCriterion() == GAGenome::MINIMIZATION) {
      while(i < pop.size() && pop.best(i).score() < boa->worst().score()){
	for(k=0;
	    pop.best(i).score() > boa->best(k).score() && k < boa->size();
	    k++);
	for(j=k; j<boa->size(); j++){
	  if(pop.best(i) == boa->best(j)) break;
	  if(pop.best(i).score() < boa->best(j).score()){
	    boa->worst().copy(pop.best(i));        // replace worst individual
	    boa->sort(gaTrue, GAPopulation::RAW);  // re-sort the population
	    break;
	  }
	}
	i++;
      }
    }
  }
  return;
}


unsigned getNumChildrenInserted(const GAPopulation& pop) {
  unsigned num_child_inserted = 0;
  for(unsigned int i=0; i<pop.size(); i++){
    if (pop.individual(i).age() == 0) num_child_inserted++;
  }
  return num_child_inserted;
}


// Use this method to update the statistics to account for the current
// population.  This routine increments the generation counter and assumes that
// the population that gets passed is the current population.
//   If we are supposed to flush the scores, then we dump them to the specified
// file.  If no flushing frequency has been specified then we don't record.

void GAStatistics::update(const GAPopulation & pop){
  ++curgen;			// must do this first so no
                                // divide-by-zero
  for (unsigned int i=0; i<pop.size(); i++) pop.individual(i).incAge(); // We increment each genome age

  if(scoreFreq > 0 && (curgen % scoreFreq == 0)) setScore(pop);
  if(Nscrs > 0 && nscrs >= Nscrs) flushScores();
  maxever = (pop.max() > maxever) ? pop.max() : maxever;
  minever = (pop.min() < minever) ? pop.min() : minever;
  long double tmpval;
  tmpval = (on*(curgen-1) + pop.ave()) / curgen;
  on = tmpval;
  tmpval = (offmax*(curgen-1) + pop.max()) / curgen;
  offmax = tmpval;
  tmpval = (offmin*(curgen-1) + pop.min()) / curgen;
  offmin = tmpval;
  setConvergence((GAGenome::optCriterion() == GAGenome::MAXIMIZATION) ? pop.max() : pop.min());

  updateHallOfFame(pop.best());

  updateBestIndividual(pop);
  numpeval = pop.nevals();
}

GAPopulation* GAStatistics::hallOfFame() const{
  assert(use_hall_of_fame_);
  return hall_of_fame_;
}

void GAStatistics::updateHallOfFame(const GAGenome& best_ind){
  if (use_hall_of_fame_) {
    // An elitist method is assumed so the score from each generation is always better or the same than the
    // previous generation. This assumption can improve consideraby the efficiency because we
    // only need to traverse the array until an individual with less fitness is found
    int      i              = hall_of_fame_->size()-1;
    long double   best_ind_score = best_ind.score();
    bool     already_in_hof = false;

    while ( i>=0 && best_ind_score >= hall_of_fame_->individual(i).score() ) {
      if ( hall_of_fame_->individual(i).score() == best_ind_score ) {
        if ( best_ind.compare( hall_of_fame_->individual(i) ) == 0.0 ) {
          already_in_hof = true;
          break;
        }
      }
      i--;
    }
    if (!already_in_hof) hall_of_fame_->add(best_ind);     // We update the hall of fame

    GALogger::instance()->appendPopulation("Hall of Fame","",*hall_of_fame_);
  }
}

void GAStatistics::updateNoStepsStats(const GAPopulation& pop) {
  update(pop);
  fit_inc_ = 0;
}


// Reset the GA's statistics based on the population.  To do this right you
// should initialize the population before you pass it to this routine.  If you
// don't, the stats will be based on a non-initialized population.
void
GAStatistics::reset(const GAPopulation & pop){
  curgen = 0;
  numsel = numcro = nummut = numrep = numeval = numpeval = 0;

  memset(gen, 0, Nscrs*sizeof(int));
  memset(aveScore, 0, Nscrs*sizeof(long double));
  memset(maxScore, 0, Nscrs*sizeof(long double));
  memset(minScore, 0, Nscrs*sizeof(long double));
  memset(devScore, 0, Nscrs*sizeof(long double));
  memset(divScore, 0, Nscrs*sizeof(long double));
  nscrs = 0;
  setScore(pop);
  if(Nscrs > 0) flushScores();

  memset(cscore, 0, Nconv*sizeof(long double));
  nconv = 0;			// should set to -1 then call setConv
  cscore[0] = ((GAGenome::optCriterion() == GAGenome::MAXIMIZATION) ? pop.max() : pop.min());

//  cscore[0] = pop.max();
//  setConvergence(maxScore[0]);

  updateBestIndividual(pop, gaTrue);
  aveCur = aveInit = pop.ave();
  maxCur = maxInit = maxever = pop.max();
  minCur = minInit = minever = pop.min();
  devCur = devInit = pop.dev();
  divCur = divInit = ((dodiv == gaTrue) ? pop.div() : (long double)-1.0);

  on = pop.ave();
  offmax = pop.max();
  offmin = pop.min();
  numpeval = pop.nevals();
  for(unsigned int i=0; i<pop.size(); i++)
    numeval += pop.individual(i).nevals();
}

void
GAStatistics::flushScores(){
  if(nscrs == 0) return;
  writeScores();
  memset(gen, 0, Nscrs*sizeof(int));
  memset(aveScore, 0, Nscrs*sizeof(long double));
  memset(maxScore, 0, Nscrs*sizeof(long double));
  memset(minScore, 0, Nscrs*sizeof(long double));
  memset(devScore, 0, Nscrs*sizeof(long double));
  memset(divScore, 0, Nscrs*sizeof(long double));
  nscrs = 0;
}


// Set the score info to the appropriate values.  Update the score count.
void
GAStatistics::setScore(const GAPopulation& pop){
  aveCur = pop.ave();
  maxCur = pop.max();
  minCur = pop.min();
  devCur = pop.dev();
  divCur = ((dodiv == gaTrue) ? pop.div() : (long double)-1.0);

  if(Nscrs == 0) return;
  gen[nscrs] = curgen;
  aveScore[nscrs] = aveCur;
  maxScore[nscrs] = maxCur;
  minScore[nscrs] = minCur;
  devScore[nscrs] = devCur;
  divScore[nscrs] = divCur;
  nscrs++;
}


// For recording the convergence we have to keep a running list of the past N
// best scores.  We just keep looping around and around the array of past
// scores.  nconv keeps track of which one is the current one.  The current
// item is thus nconv%Nconv.  The oldest is nconv%Nconv+1 or 0.
void
GAStatistics::setConvergence(long double s){
  nconv++;
  cscore[nconv%Nconv] = s;
}


// When a new number of gens to conv is specified, keep all the data that we
// can in the transfer.  Make sure we do it in the right order!  Then just
// continue on as before.
//   If someone passes us a zero then we set to 1.
int
GAStatistics::nConvergence(unsigned int n){
  if(n == 0) n = 1;
  long double *tmp = cscore;
  cscore = new long double[n];
  if(Nconv < n){
    if(nconv < Nconv){
      memcpy(cscore, tmp, (nconv+1) * sizeof(long double));
    }
    else{
      memcpy(&(cscore[Nconv-(nconv%Nconv)-1]), tmp,
	     ((nconv%Nconv)+1) * sizeof(long double));
      memcpy(cscore, &(tmp[(nconv%Nconv)+1]),
	     (Nconv-(nconv%Nconv)-1) * sizeof(long double));
    }
  }
  else{
    if(nconv < n){
      memcpy(cscore, tmp, (nconv+1) * sizeof(long double));
    }
    else{
      if((nconv % Nconv) + 1 < n){
	memcpy(&(cscore[n-(nconv%Nconv)-1]), tmp,
	       ((nconv%Nconv)+1) * sizeof(long double));
	memcpy(cscore, &(tmp[Nconv-(1+n-(nconv%Nconv))]), sizeof(long double));
      }
      else{
	memcpy(cscore, &(tmp[1+(nconv%Nconv)-n]), n * sizeof(long double));
      }
    }
  }

  Nconv = n;
  delete [] tmp;
  return Nconv;
}


int
GAStatistics::nBestGenomes(const GAGenome& genome, unsigned int n){
  if(n == 0){
    delete boa;
    boa = (GAPopulation*)0;
  }
  else if(boa == (GAPopulation*)0){
    boa = new GAPopulation(genome, n);
  }
  else {
    boa->size(n);
  }
  return n;
}

const GAGenome &
GAStatistics::bestIndividual(unsigned int n) const {
  if(boa == 0 || n >= boa->size()){
    GAErr(GA_LOC, "GAStatistics", "bestIndividual", gaErrBadPopIndex);
    n = 0;
  }
  return boa->best(n);		// this will crash if no boa
}

// Adjust the scores buffers to match the specified amount.  If someone
// specifies zero then we don't keep the scores, so set all to NULL.
int
GAStatistics::flushFrequency(unsigned int freq){
  if(freq == 0){
    if(nscrs > 0) flushScores();
    resizeScores(freq);
  }
  else if(freq > Nscrs){
    resizeScores(freq);
  }
  else if(freq < Nscrs){
    if(nscrs > freq) flushScores();
    resizeScores(freq);
  }
  Nscrs = freq;
  return freq;
}


// Resize the scores vectors to the specified amount.  Copy any scores that
// exist.
void
GAStatistics::resizeScores(unsigned int n){
  int *tmpi;
  long double *tmpf;

  if(n == 0){
    delete [] gen;
    gen = (int*)0;
    delete [] aveScore;
    aveScore = (long double*)0;
    delete [] maxScore;
    maxScore = (long double*)0;
    delete [] minScore;
    minScore = (long double*)0;
    delete [] devScore;
    devScore = (long double*)0;
    delete [] divScore;
    divScore = (long double*)0;

    nscrs = n;
  }
  else{
    tmpi = gen;
    gen = new int [n];
    memcpy(gen, tmpi, (n < Nscrs ? n : Nscrs)*sizeof(int));
    delete [] tmpi;

    tmpf = aveScore;
    aveScore = new long double [n];
    memcpy(aveScore, tmpf, (n < Nscrs ? n : Nscrs)*sizeof(long double));
    delete [] tmpf;

    tmpf = maxScore;
    maxScore = new long double [n];
    memcpy(maxScore, tmpf, (n < Nscrs ? n : Nscrs)*sizeof(long double));
    delete [] tmpf;

    tmpf = minScore;
    minScore = new long double [n];
    memcpy(minScore, tmpf, (n < Nscrs ? n : Nscrs)*sizeof(long double));
    delete [] tmpf;

    tmpf = devScore;
    devScore = new long double [n];
    memcpy(devScore, tmpf, (n < Nscrs ? n : Nscrs)*sizeof(long double));
    delete [] tmpf;

    tmpf = divScore;
    divScore = new long double [n];
    memcpy(divScore, tmpf, (n < Nscrs ? n : Nscrs)*sizeof(long double));
    delete [] tmpf;

    if(nscrs > n) nscrs = n;
  }
  Nscrs = n;
}


// Write the current scores to file.  If this is the first chunk (ie gen[0]
// is 0) then we create a new file.  Otherwise we append to an existing file.
// We give no notice that we're overwriting the existing file!!
void
GAStatistics::writeScores(){
  if(!scorefile) return;
#ifdef GALIB_USE_STREAMS
  STD_OFSTREAM outfile(scorefile, ((gen[0] == 0) ?
			       (STD_IOS_OUT | STD_IOS_TRUNC) :
			       (STD_IOS_OUT | STD_IOS_APP)));
// should be done this way, but SGI systems (and others?) don't do it right...
//  if(! outfile.is_open()){
  if(outfile.fail()){
    GAErr(GA_LOC, "GAStatistics", "writeScores", gaErrWriteError, scorefile);
    return;
  }
  scores(outfile, which);
  outfile.close();
#endif
}


#ifdef GALIB_USE_STREAMS
int
GAStatistics::write(const char* filename) const {
  STD_OFSTREAM outfile(filename, (STD_IOS_OUT | STD_IOS_TRUNC));
// should be done this way, but SGI systems (and others?) don't do it right...
//  if(! outfile.is_open()){
  if(outfile.fail()){
    GAErr(GA_LOC, "GAStatistics", "write", gaErrWriteError, filename);
    return 1;
  }
  write(outfile);
  outfile.close();
  return 0;
}

int
GAStatistics::write(STD_OSTREAM & os) const {
  os << curgen << "\t# current generation\n";
  os << convergence() << "\t# current convergence\n";
  os << numsel << "\t# number of selections since initialization\n";
  os << numcro << "\t# number of crossovers since initialization\n";
  os << nummut << "\t# number of mutations since initialization\n";
  os << numrep << "\t# number of replacements since initialization\n";
  os << numeval << "\t# number of genome evaluations since initialization\n";
  os << numpeval << "\t# number of population evaluations since initialization\n";
  os << maxever << "\t# maximum score since initialization\n";
  os << minever << "\t# minimum score since initialization\n";
  os << on << "\t# average of all scores ('on-line' performance)\n";
  os << offmax << "\t# average of maximum scores ('off-line' performance)\n";
  os << offmin << "\t# average of minimum scores ('off-line' performance)\n";
  os << "\n";
  os << aveInit << "\t# mean score in initial population\n";
  os << maxInit << "\t# maximum score in initial population\n";
  os << minInit << "\t# minimum score in initial population\n";
  os << devInit << "\t# standard deviation of initial population\n";
  os <<divInit<<"\t# diversity of initial population (0=identical,-1=unset)\n";
  os << "\n";
  os << aveCur << "\t# mean score in current population\n";
  os << maxCur << "\t# maximum score in current population\n";
  os << minCur << "\t# minimum score in current population\n";
  os << devCur << "\t# standard deviation of current population\n";
  os << divCur<<"\t# diversity of current population (0=identical,-1=unset)\n";
  os << "\n";
  os << Nconv << "\t# how far back to look for convergence\n";
  os << scoreFreq << "\t# how often to record scores\n";
  os << Nscrs << "\t# how often to write scores to file\n";
  os << scorefile << "\t# name of file to which scores are written\n";
  return 0;
}


// You can specify the data that you want to dump out when you call this
// routine, or you can just let it use the selection from the object.  If you
// specify a data set, that will be used rather than the 'which' in the object.
int
GAStatistics::scores(const char* filename, int w){
  STD_OFSTREAM outfile(filename, (STD_IOS_OUT | STD_IOS_TRUNC));
// should be done this way, but SGI systems (and others?) don't do it right...
//  if(! outfile.is_open()){
  if(outfile.fail()){
    GAErr(GA_LOC, "GAStatistics", "scores", gaErrWriteError, filename);
    return 1;
  }
  scores(outfile, w);
  outfile.close();
  return 0;
}


int
GAStatistics::scores(STD_OSTREAM & os, int w){
  if(w == NoScores) w = which;

  for(unsigned int i=0; i<nscrs; i++){
    os << gen[i];
    if(w & Mean) os << "\t" << aveScore[i];
    if(w & Maximum)  os << "\t" << maxScore[i];
    if(w & Minimum)  os << "\t" << minScore[i];
    if(w & Deviation)  os << "\t" << devScore[i];
    if(w & Diversity)  os << "\t" << divScore[i];
    os << "\n";
  }
  return 0;
}

// TODO: Moved here to gain compatibility with gcc versions prior to 4.X
// and linking as a shared library (this symbol was discarded by the
// linker when it was an inline method in the header file)
long double GAStatistics::current(int w) const {
  long double val = 0.0;
  switch(w){
  case Mean:      val = aveCur; break;
  case Maximum:   val = maxCur; break;
  case Minimum:   val = minCur; break;
  case Deviation: val = devCur; break;
  case Diversity: val = divCur; break;
  default: break;
  }
  return val;
}

void GAStatistics::updateFitnessIncrement(const GAPopulation&      old_pop,
                                          const vector<GAGenome*>& parents,
                                          const GAPopulation& children){
  assert(parents.size() == children.size() || parents.size() == (children.size() + 1) );

  old_pop.sort();
  long double    max_score          = old_pop.best().score();
  long double    min_score          = old_pop.worst().score();
  long double    best_parent_score, best_child_score;
  long double    fit_inc_sum_tmp    = 0.0;

  if ( max_score == min_score ) {
    fit_inc_ = 0.0;
  }
  else {
    GAGenome* parent1;
    GAGenome* parent2;
    GAGenome* child1;
    GAGenome* child2;

    for (unsigned int i=0; i<parents.size()-2; i++){
      parent1 = parents[i];
      parent2 = parents[i+1];
      child1  = & children.individual(i);
      // When the population is odd, the parents pop is greater than the children so in order to have a cleaner code,
      // child2 is assigned child1 so the rest of the code remains the same.
      if ( ( parents.size() > children.size() ) && ( i==children.size()-1 ) ) child2 = child1;
      else                                                                    child2 = & children.individual(i+1);

      best_parent_score =  ( parent1->score() > parent2->score() ) ? parent1->score() : parent2->score();
      best_child_score  =  ( child1->score()  > child2->score()  ) ? child1->score()  : child2->score();
      fit_inc_sum_tmp   += getScaledFitness(max_score,min_score, best_child_score) -
                           getScaledFitness(max_score,min_score, best_parent_score);
    }
    fit_inc_ = fit_inc_sum_tmp / children.size(); // children pop size not parent since, since the latter is bigger


    assert(!isnan(fit_inc_));

  }
}

long double GAStatistics::getScaledFitness(long double max_score, long double min_score, long double genome_score) {
  long double result;
  if      (genome_score == max_score)      result = 1.0;  // Needed for representation errors
  else if (genome_score == min_score)      result = 0.0;
  else result = ( ( 1/(max_score - min_score) ) * genome_score + 1 - ( max_score/(max_score - min_score) ) );
  return result;
}

long double GAStatistics::fitInc() const { return fit_inc_;}

unsigned GAStatistics::numNewChild() const { return num_new_child_; }

#endif

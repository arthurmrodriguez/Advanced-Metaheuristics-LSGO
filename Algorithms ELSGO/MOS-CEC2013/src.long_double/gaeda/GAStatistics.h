// $Header: /home/cvs/galib/ga/GAStatistics.h,v 1.3 2004/12/28 14:38:44 mwall Exp $
/* ----------------------------------------------------------------------------
  statistics.h
  mbwall 14jul95
  Copyright (c) 1995 Massachusetts Institute of Technology
                   - all rights reserved

 DESCRIPTION:
  Header for the statistics object used by the GA objects.
---------------------------------------------------------------------------- */
#ifndef GASTATISTICS_H
#define GASTATISTICS_H

/* INCLUDES */
#include <math.h>
#include <string.h>

#include "gatypes.h"
#include "gaconfig.h"
#include "genomes/GAGenome.h"
#include "GAPopulation.h"

// Default settings and their names.
extern int  gaDefNumBestGenomes;
extern int  gaDefScoreFrequency1;
extern int  gaDefScoreFrequency2;
extern int  gaDefFlushFrequency;
extern char gaDefScoreFilename[];

/* ----------------------------------------------------------------------------
Statistics class
  We define this class as a storage object for the current state of the GA.
Whereas the parameters object keeps track of the user-definable settings for
the GA, the statistics object keeps track of the data that the GA generates
along the way.
---------------------------------------------------------------------------- */
class GAStatistics {
public:
  enum {
    NoScores=0x00,
    Mean=0x01,
    Maximum=0x02,
    Minimum=0x04,
    Deviation=0x08,
    Diversity=0x10,
    AllScores=0xff
    };

  GAStatistics();
  GAStatistics(const GAStatistics&);
  GAStatistics& operator=(const GAStatistics& orig){copy(orig); return *this;}
  virtual ~GAStatistics();
  void copy(const GAStatistics &);

  long double online() const {return on;}
  long double offlineMax() const {return offmax;}
  long double offlineMin() const {return offmin;}
  long double initial(int w=Maximum) const;
  long double current(int w=Maximum) const;
  long double maxEver() const {return maxever;}
  long double minEver() const {return minever;}

  int generation() const {return curgen;}
  unsigned long int selections() const {return numsel;}
  unsigned long int crossovers() const {return numcro;}
  unsigned long int mutations() const {return nummut;}
  unsigned long int replacements() const {return numrep;}
  unsigned long int indEvals() const {return numeval;}
  unsigned long int popEvals() const {return numpeval;}
  long double convergence() const;

  int nConvergence() const { return Nconv; }
  int nConvergence(unsigned int);
  int nBestGenomes(const GAGenome&, unsigned int);
  int nBestGenomes() const { return(boa ? boa->size() : 0); }
  int scoreFrequency(unsigned int x) { return(scoreFreq = x); }
  int scoreFrequency() const { return scoreFreq; }
  int flushFrequency(unsigned int x);
  int flushFrequency() const { return Nscrs; }
  const char* scoreFilename(const char *filename);
  const char* scoreFilename() const { return scorefile; }
  int selectScores(int w){ return which = w; }
  int selectScores() const { return which; }
  GABoolean recordDiversity(GABoolean flag){ return dodiv=flag; }
  GABoolean recordDiversity() const { return dodiv; }
  void flushScores();

  // Estos métodos vienen del desarrollo de smuelas
  void updateFitnessIncrement  (const GAPopulation& oldpop,
				 const vector<GAGenome*>& parents,
				 const GAPopulation& children );

  void updateHallOfFame(const GAGenome& best_ind);//

  void          updateNoStepsStats      (const GAPopulation& pop );
  long double        fitInc                  (                        ) const;
  unsigned      numNewChild             (                        ) const;
  GAPopulation* hallOfFame              (                        ) const;
  void          storeHallOfFame           (                        ) {use_hall_of_fame_ = true;}


  void update(const GAPopulation & pop);
  void reset(const GAPopulation & pop);
  const GAPopulation & bestPopulation() const { return *boa; }
  const GAGenome & bestIndividual(unsigned int n=0) const;
  long double getScaledFitness(long double max_score, long double min_score, long double genome_score);

#ifdef GALIB_USE_STREAMS // NO_STREAMS
  int scores(const char* filename, int which=NoScores);
  int scores(STD_OSTREAM& os, int which=NoScores);
  int write(const char* filename) const;
  int write(STD_OSTREAM& os) const;
#endif

// These should be protected (accessible only to the GA class) but for now they
// are publicly accessible.  Do not try to set these unless you know what you
// are doing!!
  unsigned long int numsel;     // number of selections since reset
  unsigned long int numcro;	// number of crossovers since reset
  unsigned long int nummut;	// number of mutations since reset
  unsigned long int numrep;	// number of replacements since reset
  unsigned long int numeval;	// number of individual evaluations since reset
  unsigned long int numpeval;	// number of population evals since
                                // reset
  void setScore(const GAPopulation&);

protected:
  unsigned int curgen;		// current generation number
  unsigned int scoreFreq;	// how often (in generations) to record scores
  GABoolean dodiv;		// should we record diversity?

  long double maxever;		// maximum score since initialization
  long double minever;		// minimum score since initialization
  long double on;			// "on-line" performance (ave of all scores)
  long double offmax;			// "off-line" performance (ave of maximum)
  long double offmin;			// "off-line" performance (ave of minimum)

  long double aveInit;		// stats from the initial population
  long double maxInit;
  long double minInit;
  long double devInit;
  long double divInit;

  long double aveCur;			// stats from the current population
  long double maxCur;
  long double minCur;
  long double devCur;
  long double divCur;

  unsigned int nconv, Nconv;	// how many scores we're recording (flushFreq)
  long double * cscore;		// best score of last n generations

  unsigned int nscrs, Nscrs;	// how many scores do we have?
  int * gen;			// generation number corresponding to scores
  long double * aveScore;		// average scores of each generation
  long double * maxScore;		// best scores of each generation
  long double * minScore;		// worst scores of each generation
  long double * devScore;		// stddev of each generation
  long double * divScore;		// diversity of each generation
  char * scorefile;		// name of file to which scores get written
  int which;			// which data to write to file
  GAPopulation * boa;		// keep a copy of the best genomes

  unsigned      num_new_child_;// keeps the number of new child accepted into the population in the last generation
  GAPopulation* hall_of_fame_; // keeps the best genomes from the previous generations
  bool          use_hall_of_fame_;

  void setConvergence(long double);

  void updateBestIndividual(const GAPopulation&, GABoolean flag=gaFalse);
  void writeScores();
  void resizeScores(unsigned int);

  long double fit_inc_;

  friend class GA;
};



inline const char* GAStatistics::scoreFilename(const char* filename){
  delete [] scorefile;
  scorefile = 0;
  if(filename){
    scorefile = new char[strlen(filename)+1];
    strcpy(scorefile, filename);
  }
  return scorefile;
}

inline long double GAStatistics::convergence() const {
  long double cnv = 0.0;
  if(nconv >= Nconv-1 && cscore[nconv%Nconv] != 0)
    cnv = (long double)(cscore[(nconv+1)%Nconv]) / (long double)(cscore[nconv%Nconv]);
  return (long double)cnv;
}

inline long double GAStatistics::initial(int w) const {
  long double val = 0.0;
  switch(w){
  case Mean:      val = aveInit; break;
  case Maximum:   val = maxInit; break;
  case Minimum:   val = minInit; break;
  case Deviation: val = devInit; break;
  case Diversity: val = divInit; break;
  default: break;
  }
  return val;
}


#ifdef GALIB_USE_STREAMS
inline STD_OSTREAM& operator<< (STD_OSTREAM& os, const GAStatistics& s)
{ s.write(os); return os; }
#endif

#endif

#include "GARealOps.h"

#include "normal.h"
#include "GACommonOps.h"
#include "GAEDAConfig.h"
#include "GAPopulation.h"
#include "ImprovPerDimManager.h"
#include "solis.h"
#include "genomes/GAGenome.h"
#include "genomes/GA1DArrayGenome.h"
#include "genomes/MOSGenome.h"

#include <vector>

// Initializers

void RealUniformInitializer (GAGenome& g) {

   GA1DArrayAlleleGenome<double>& gen = DYN_CAST (GA1DArrayAlleleGenome<double>&, g);

   for (int i = 0; i < gen.length (); i++) {

      double upper = gen.alleleset (i).upper ();
      double lower = gen.alleleset (i).lower ();

      gen.gene (i, GARandomDouble (lower, upper));

   }

   return;

}

 void RealOrderedInitializer (GAGenome& g) {
   return OrderedInitializer<double> (g);
 }



// Comparators

double RealElementComparator (const GAGenome& g1, const GAGenome& g2) {
  return ElementComparator<double> (g1, g2);
}


double RealEuclideanComparator (const GAGenome& a, const GAGenome& b) {

   const GA1DArrayGenome<double>& sis = DYN_CAST (const GA1DArrayGenome<double>&, a);
   const GA1DArrayGenome<double>& bro = DYN_CAST (const GA1DArrayGenome<double>&, b);

   if (sis.length () != bro.length ()) return -1;
   if (sis.length () == 0            ) return  0;

   double xn, yn, tmp;
   double result = 0.0;

   for (int i=0; i < sis.length(); i++) {

      xn      = sis.gene (i);
      yn      = bro.gene (i);
      tmp     = xn  -  yn;
      result += tmp * tmp;

   }

   return sqrt (result);

}


double RealEuclideanLinearizedComparator (const GAGenome& g1, const GAGenome& g2) {

   const GA1DArrayAlleleGenome<double>& genome = dynamic_cast <const GA1DArrayAlleleGenome<double>&> (g1);

   double max_diff  = genome.alleleset(1).upper () - genome.alleleset(1).lower ();
   double max_value = sqrt (max_diff * max_diff) * genome.length ();

   assert(max_value != 0.0);

   double return_value = (RealEuclideanComparator (g1, g2) + max_value) / (2.0 * max_value);

   assert (return_value <= 1.0 && return_value >= 0.0);

   return return_value;

}


double RealChebyshevComparator (const GAGenome& g1, const GAGenome& g2) {

   const GA1DArrayGenome<double>& sis = DYN_CAST (const GA1DArrayGenome<double>&, g1);
   const GA1DArrayGenome<double>& bro = DYN_CAST (const GA1DArrayGenome<double>&, g2);

	double xn, yn, tmp;
   double max_distance = 0;

   for (int i = 0; i < sis.length (); i++) {

   	xn = sis.gene (i);
		yn = bro.gene (i);
		tmp = fabs (xn - yn);

		if (tmp > max_distance)
         max_distance = tmp;

	}

	return max_distance;

}


double RealManhattanComparator (const GAGenome& g1, const GAGenome& g2) {

   const GA1DArrayGenome<double>& sis = DYN_CAST (const GA1DArrayGenome<double>&, g1);
   const GA1DArrayGenome<double>& bro = DYN_CAST (const GA1DArrayGenome<double>&, g2);

	double xn, yn, tmp;
	double sum = 0;

   for (int i = 0; i < sis.length (); i++) {

		xn = sis.gene (i);
		yn = bro.gene (i);
		tmp = fabs (xn - yn);
		sum +=  tmp;

	}

	return sum;

}


// Crossovers

int RealOnePointCrossover     (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return OnePointCrossover<double>     (p1, p2, c1, c2);}
int RealTwoPointsCrossover    (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return TwoPointsCrossover<double>    (p1, p2, c1, c2);}
int RealEvenOddCrossover      (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return EvenOddCrossover<double>      (p1, p2, c1, c2);}
int RealPartialMatchCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return PartialMatchCrossover<double> (p1, p2, c1, c2);}
int RealUniformCrossover      (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return UniformCrossover<double>      (p1, p2, c1, c2);}
int RealOrderCrossover        (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return OrderCrossover<double>        (p1, p2, c1, c2);}
int RealCycleCrossover        (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return CycleCrossover<double>        (p1, p2, c1, c2);}

int RealAlternativeOrderCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return AlternativeOrderCrossover<double> (p1, p2, c1, c2);}
int RealAlternativeCycleCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return AlternativeCycleCrossover<double> (p1, p2, c1, c2);}


// Blend crossover generates a new value based on the interval between parents.
// We generate a uniform distribution based on the distance between parent
// values, then choose the child value based upon that distribution.
int RealBlendCrossover(const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {

  const GA1DArrayAlleleGenome<double> &mom = DYN_CAST(const GA1DArrayAlleleGenome<double> &, p1);
  const GA1DArrayAlleleGenome<double> &dad = DYN_CAST(const GA1DArrayAlleleGenome<double> &, p2);

  int n=0;

  if(c1 && c2){
    GA1DArrayGenome<double> &sis=DYN_CAST(GA1DArrayGenome<double> &, *c1);
    GA1DArrayGenome<double> &bro=DYN_CAST(GA1DArrayGenome<double> &, *c2);

    int len = GAMax(mom.length(), dad.length());
    for(int i=0; i<len; i++) {
      double dist = 0;
      if(mom.gene(i) > dad.gene(i))
        dist = mom.gene(i) - dad.gene(i);
      else
        dist = dad.gene(i) - mom.gene(i);

      double lo = (GAMin(mom.gene(i), dad.gene(i))) - 0.5*dist;
      if (lo < mom.alleleset(i).lower())
        lo = mom.alleleset(i).lower();    // these were added to assure that the limits are not passed

      double hi = (GAMax(mom.gene(i), dad.gene(i))) + 0.5*dist;
      if (hi > mom.alleleset(i).upper())
        hi = mom.alleleset(i).upper();// these were added to assure that the limits are not passed

      sis.gene(i, GARandomDouble(lo, hi));
      bro.gene(i, GARandomDouble(lo, hi));
    }
    n = 2;
  }
  else if(c1 || c2){
    GA1DArrayGenome<double> &sis = (c1 ?
				   DYN_CAST(GA1DArrayGenome<double> &, *c1) :
				   DYN_CAST(GA1DArrayGenome<double> &, *c2));

    int len = GAMax(mom.length(), dad.length());
    for(int i=0; i<len; i++) {
      double dist = 0;
      if(mom.gene(i) > dad.gene(i)) 	dist = mom.gene(i) - dad.gene(i);
      else                          	dist = dad.gene(i) - mom.gene(i);

      double lo = (GAMin(mom.gene(i), dad.gene(i))) - 0.5*dist;
      if (lo < mom.alleleset(i).lower())
        lo = mom.alleleset(i).lower();   // these were added to assure that the limits are not passed
      double hi = (GAMax(mom.gene(i), dad.gene(i))) + 0.5*dist;
      if (hi > mom.alleleset(i).upper())
        hi = mom.alleleset(i).upper();   // these were added to assure that the limits are not passed
      sis.gene(i, GARandomDouble(lo, hi));
    }
    n = 1;
  }
  return n;
}


// Arithmetic crossover generates a new value that is the average of the parent
// values.  Note that this means both children in a sexual crossover will be
// identical.  If parents are not the same length, the extra elements are not
// set!  You might want to add some noise to this so that both children are not
// the same...
int RealArithmeticCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {

  const GA1DArrayGenome<double> &mom=
    DYN_CAST(const GA1DArrayGenome<double> &, p1);
  const GA1DArrayGenome<double> &dad=
    DYN_CAST(const GA1DArrayGenome<double> &, p2);

  int n=0;

  if(c1 && c2){
    GA1DArrayGenome<double> &sis=DYN_CAST(GA1DArrayGenome<double> &, *c1);
    GA1DArrayGenome<double> &bro=DYN_CAST(GA1DArrayGenome<double> &, *c2);

    int len = GAMax(mom.length(), dad.length());
    for(int i=0; i<len; i++) {
      sis.gene(i, 0.5 * (mom.gene(i) + dad.gene(i)));
      bro.gene(i, 0.5 * (mom.gene(i) + dad.gene(i)));
    }
    n = 2;
  }
  else if(c1 || c2){
    GA1DArrayGenome<double> &sis = (c1 ?
				   DYN_CAST(GA1DArrayGenome<double> &, *c1) :
				   DYN_CAST(GA1DArrayGenome<double> &, *c2));

    int len = GAMax(mom.length(), dad.length());
    for(int i=0; i<len; i++) {
      sis.gene(i, 0.5 * (mom.gene(i) + dad.gene(i)));
    }
    n = 1;
  }

  return n;
}



// DE Crossover

int RealExponentialCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, const double prob) {return ExponentialCrossover<double> (p1, p2, c1, prob);}
int RealBinomialCrossover    (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, const double prob) {return BinomialCrossover<double>    (p1, p2, c1, prob);}



// Mutators

int RealSwapMutator             (GAGenome& g, double pmut) {return SwapMutator<double>             (g, pmut);}
int RealFlipMutator             (GAGenome& g, double pmut) {return FlipMutator<double>             (g, pmut);}
int RealRepeatedExchangeMutator (GAGenome& g, double pmut) {return RepeatedExchangeMutator<double> (g, pmut);}
int RealSimpleInversionMutator  (GAGenome& g, double pmut) {return SimpleInversionMutator<double>  (g, pmut);}

int RealUniformMutator (GAGenome& g, double pmut) {

   int nMut = 0;
   GA1DArrayAlleleGenome<double>& gen = DYN_CAST (GA1DArrayAlleleGenome<double>&, g);

   if (pmut <= 0.0)
      return 0;

   for (int i = 0; i < gen.length (); i++) {

      if (GAFlipCoin (pmut)) {

         double upper = gen.alleleset (i).upper ();
         double lower = gen.alleleset (i).lower ();

         gen.gene (i, GARandomDouble (lower, upper));

         nMut++;

      }

   }

   return nMut;

}

int RealGaussianMutator (GAGenome& g, double pmut) {

  GA1DArrayAlleleGenome<double>& child = DYN_CAST (GA1DArrayAlleleGenome<double>&, g);
  //register int i,n;
  if(pmut <= 0.0) return(0);

  double nMut       = pmut * (double)(child.length());
  int    length     = child.length()-1;
  double new_value  = 0.0;

//  if (nMut < 1.0){		// we have to do a flip test on each element
//    nMut = 0;
//    for(i=length; i>=0; i--){
//      double gene_value = child.gene(i);
//
//      if(GAFlipCoin(pmut)){
//	      if (child.alleleset(i).type() == GAAllele::ENUMERATED ||
//	          child.alleleset(i).type() == GAAllele::DISCRETIZED) {
//	        new_value = child.alleleset(i).allele();
//	      }
//	      else if(child.alleleset(i).type() == GAAllele::BOUNDED){
//	        double allele_region = child.alleleset(i).upper() - child.alleleset(i).lower();
//	        do {
//	          new_value = gene_value + GAUnitGaussian()*allele_region;
//	        } while (new_value > child.alleleset(i).upper() || new_value < child.alleleset(i).lower());
//	      }
//	      child.gene(i, new_value);
//	      nMut++;
//      }
//
//    }
//  }
//  else{				        // only mutate the ones we need to
//  }
  if (GAFlipCoin(nMut - (int)nMut) ) nMut++;

  int    idx;
  double gene_value, allele_region;

  for(int n=0; n<(int)nMut; n++){
    idx        = GARandomInt(0,length);
    gene_value = child.gene(idx);

    if (child.alleleset(idx).type() == GAAllele::ENUMERATED ||
	      child.alleleset(idx).type() == GAAllele::DISCRETIZED) {
	    new_value = child.alleleset(idx).allele();
    }
    else if(child.alleleset(idx).type() == GAAllele::BOUNDED){
      allele_region = child.alleleset(idx).upper() - child.alleleset(idx).lower();
      do {
        new_value = gene_value + GAUnitGaussian()*allele_region;
      } while (new_value > child.alleleset(idx).upper() || new_value < child.alleleset(idx).lower());
    }
    child.gene(idx, new_value);
  }
  return((int)nMut);
}


// ES Crossovers

int RealIntermediateCrossover (const std::vector<GAGenome*> parents, GAGenome* child) {

  GA1DArrayGenome<double>& g = DYN_CAST(GA1DArrayGenome<double>&, *child);

  unsigned ro = parents.size();

  std::vector<GA1DArrayAlleleGenome<double>*> pars (ro, (GA1DArrayAlleleGenome<double>*)0);

  for (unsigned i = 0; i < ro; i++)
    pars [i] = DYN_CAST(GA1DArrayAlleleGenome<double>*, parents[i]);

  // Recombine both endogenous parameters and variables
  for (int i = 0; i < g.size(); i++) {

    double sum = 0.0;

    // First, genes
    for (unsigned j = 0; j < ro; j++)
      sum += (1.0/ro) * pars[j]->gene(i);

    g.gene(i, sum);

    sum = 0.0;

    // Then, Std Devs
    for (unsigned j = 0; j < ro; j++)
      sum += (1.0/ro) * pars[j]->stdDev(i);

    g.stdDev(i, sum);

  }

  return 1;

}


int RealIntermediateCrossoverUpdateOnly (const std::vector<GAGenome*> parents, GAGenome* child) {

  GA1DArrayGenome<double>& g = DYN_CAST(GA1DArrayGenome<double>&, *child);

  unsigned ro = parents.size();

  std::vector<GA1DArrayAlleleGenome<double>*> pars (ro, (GA1DArrayAlleleGenome<double>*)0);

  for (unsigned i = 0; i < ro; i++)
    pars [i] = DYN_CAST(GA1DArrayAlleleGenome<double>*, parents[i]);

  // Recombine both endogenous parameters and variables
  for (int i = 0; i < g.size(); i++) {

    double sum = 0.0;

    // Then, Std Devs
    for (unsigned j = 0; j < ro; j++)
      sum += (1.0/ro) * pars[j]->stdDev(i);

    g.stdDev(i, sum);

  }

  return 1;

}


int RealDominantCrossover (const std::vector<GAGenome*> parents, GAGenome* child) {

  GA1DArrayGenome<double>& g = DYN_CAST(GA1DArrayGenome<double>&, *child);

  unsigned ro = parents.size();

  std::vector<GA1DArrayAlleleGenome<double>*> pars (ro, (GA1DArrayAlleleGenome<double>*)0);

  for (unsigned i = 0; i < ro; i++)
    pars [i] = DYN_CAST(GA1DArrayAlleleGenome<double>*, parents[i]);

  // Recombine both endogenous parameters and variables
  for (int i = 0; i < g.size(); i++) {

    unsigned p = GARandomInt (0, ro-1);

    g.gene  (i, pars[p]->gene  (i));
    g.stdDev(i, pars[p]->stdDev(i));

  }

  return 1;

}


int RealDominantCrossoverUpdateOnly (const std::vector<GAGenome*> parents, GAGenome* child) {

  GA1DArrayGenome<double>& g = DYN_CAST(GA1DArrayGenome<double>&, *child);

  unsigned ro = parents.size();

  std::vector<GA1DArrayAlleleGenome<double>*> pars (ro, (GA1DArrayAlleleGenome<double>*)0);

  for (unsigned i = 0; i < ro; i++)
    pars [i] = DYN_CAST(GA1DArrayAlleleGenome<double>*, parents[i]);

  // Recombine both endogenous parameters and variables
  for (int i = 0; i < g.size(); i++) {

    unsigned p = GARandomInt (0, ro-1);

    g.stdDev(i, pars[p]->stdDev(i));

  }

  return 1;

}


// ES Mutators

int RealIsotropicMutator (GAGenome& g) {

  GA1DArrayGenome<double>& gen = DYN_CAST(GA1DArrayGenome<double>&, g);

  // Initialization of Tau vars
  unsigned N    = GAEDAConfig::handle()->getProblemSize();
  double   tau  = 1 / (sqrt(2*sqrt(N)));

  // Update mutation rates
  for (int i = 0; i < gen.size(); i++)
    gen.stdDev(i, gen.stdDev(i) * exp(tau * normal_random(0,1)));

  // Mutate variables
  for (int i = 0; i < gen.size(); i++)
    gen.gene(i, gen.gene(i) + (gen.stdDev(0) * normal_random(0, 1)));

  return 1;

}


int RealIsotropicMutatorUpdateOnly (GAGenome& g) {

  GA1DArrayGenome<double>& gen = DYN_CAST(GA1DArrayGenome<double>&, g);

  // Initialization of Tau vars
  unsigned N    = GAEDAConfig::handle()->getProblemSize();
  double   tau  = 1 / (sqrt(2*sqrt(N)));

  // Update mutation rates
  for (int i = 0; i < gen.size(); i++)
    gen.stdDev(i, gen.stdDev(i) * exp(tau * normal_random(0,1)));

  return 1;

}


int RealNonIsotropicMutator (GAGenome& g) {

  GA1DArrayGenome<double>& gen = DYN_CAST(GA1DArrayGenome<double>&, g);

  // Initialization of Tau vars
  unsigned N    = GAEDAConfig::handle()->getProblemSize();
  double   tau0 = 1 / (sqrt(2*N));
  double   tau  = 1 / (sqrt(2*sqrt(N)));

  // Update mutation rates
  double F = exp(tau0*normal_random(0,1));

  for (int i = 0; i < gen.size(); i++)
    gen.stdDev(i, F * gen.stdDev(i) * exp(tau * normal_random(0,1)));

  // Mutate variables
  for (int i = 0; i < gen.size(); i++)
    gen.gene(i, gen.gene(i) + (gen.stdDev(i) * normal_random(0, 1)));

  return 1;

}


int RealNonIsotropicMutatorUpdateOnly (GAGenome& g) {

  GA1DArrayGenome<double>& gen = DYN_CAST(GA1DArrayGenome<double>&, g);

  // Initialization of Tau vars
  unsigned N    = GAEDAConfig::handle()->getProblemSize();
  double   tau0 = 1 / (sqrt(2*N));
  double   tau  = 1 / (sqrt(2*sqrt(N)));

  // Update mutation rates
  double F = exp(tau0*normal_random(0,1));

  for (int i = 0; i < gen.size(); i++)
    gen.stdDev(i, F * gen.stdDev(i) * exp(tau * normal_random(0,1)));

  return 1;

}


// ES Initializers

void RealESUniformInitializer (GAGenome& g) {

   GA1DArrayAlleleGenome<double>& gen = DYN_CAST (GA1DArrayAlleleGenome<double>&, g);

   for (int i = 0; i < gen.length (); i++) {

      double upper = gen.alleleset (i).upper ();
      double lower = gen.alleleset (i).lower ();

      gen.gene (i, GARandomDouble (lower, upper));

   }

   for (int i = 0; i < gen.length (); i++)
      gen.stdDev (i, 0.3);

   return;

}


// MTS Local Searches

//#define DEBUG_MTS
   #include <iomanip>
#ifdef DEBUG_MTS
   #include <iomanip>
#endif

void computeDiffInd(GA1DArrayAlleleGenome<double>& main_ind, GA1DArrayAlleleGenome<double>& other_ind, double step_size, /* inout */ unsigned& used_evals) {
  double new_value;
  for (unsigned dim_pos = 0; dim_pos < main_ind.length(); dim_pos++) {
    new_value = other_ind.gene(dim_pos) + (main_ind.gene(dim_pos) - other_ind.gene(dim_pos)) * step_size;
    main_ind.gene(dim_pos, new_value);
  }
  used_evals++;
}

bool DiffsLS (GAPopulation& pop, unsigned& used_evals, double& fit_inc_acum, unsigned& improvements, unsigned max_evals) {

  GA1DArrayAlleleGenome<double>& best_ind = dynamic_cast<GA1DArrayAlleleGenome<double>&>(*dynamic_cast<MOSGenome&>(pop.best()).getGenome(GAID::RealEncoding));
  GA1DArrayAlleleGenome<double>* tmp_ind  = dynamic_cast<GA1DArrayAlleleGenome<double>*>(best_ind.clone());

  double SRLSDiff = best_ind.SR();
  bool   improved = false;

  for (unsigned i = 0; i < pop.size() && used_evals < max_evals; i++) {

    GA1DArrayAlleleGenome<double>& pop_ind = dynamic_cast<GA1DArrayAlleleGenome<double>&>(*dynamic_cast<MOSGenome&>(pop.individual(i)).getGenome(GAID::RealEncoding));
    if (&pop_ind == &best_ind) continue;

    double step_size = SRLSDiff * GARandomDouble(0.1, 1.0);

    computeDiffInd(*tmp_ind, pop_ind, step_size, used_evals);

    if (GAGenome::compareScores(tmp_ind->score(), best_ind.score()) == GAGenome::WORSE) {

      tmp_ind->copy(best_ind); // restore

      computeDiffInd(*tmp_ind,pop_ind,-step_size,used_evals);

      if (GAGenome::compareScores(tmp_ind->score(), best_ind.score()) == GAGenome::WORSE)
        tmp_ind->copy(best_ind);

    }

    if (GAGenome::compareScores(tmp_ind->score(), best_ind.score()) == GAGenome::BETTER) {
      fit_inc_acum += tmp_ind->computeFitnessIncrement(best_ind.score());
      best_ind.copy(*tmp_ind);
      improvements++;
      improved = true;
    }

  }

  if (!improved) {
    SRLSDiff = SRLSDiff / 2.0;
    if (SRLSDiff < 1e-13) SRLSDiff = 1.0;
    best_ind.SR(SRLSDiff);
  }

  delete tmp_ind;

  return improved;

}


double MTS_LS1 (GAGenome& Xk, double& fitBest, unsigned& evals, double& fit_inc_acum, unsigned& improvements, double& SR, unsigned maxIters) {

  static const double adjustFailed = GAEDAConfig::handle()->getMTSAdjustFailed();
  static const double adjustMin    = GAEDAConfig::handle()->getMTSAdjustMin   ();
  static const double moveLeft     = GAEDAConfig::handle()->getMTSMoveLeft    ();
  static const double moveRight    = GAEDAConfig::handle()->getMTSMoveRight   ();

  GA1DArrayAlleleGenome<double>& g = dynamic_cast <GA1DArrayAlleleGenome<double>&> (Xk);

  unsigned sz  = g.size();
  double grade = 0;
  // double SR    = g.SR();
  fit_inc_acum = 0.0;

  double allele_min = g.alleleset(0).lower();
  double allele_max = g.alleleset(0).upper();

  // If this solution was not improved during the last COMPLETE iteration,
  // update the SR interval
  if (!g.improve() && g.lastIterPos() == 0 && g.step() == 0) {
    SR *= adjustFailed;

    if (SR < 1e-14)
      SR = (allele_max - allele_min) * adjustMin;

    g.SR(SR);
  }

  // Mark solution as not improved (only at the beginning of a complete iteration)
  if (g.lastIterPos() == 0 && g.step() == 0)
    g.improve(0);

  // Boolean value which says if the LS was interrupted in the middle of
  // an iteration in the second step of one dimension
  bool ignoreFirstStep = (g.step() == 1) ? true : false;

  // For each dimension, perform the LS
  for (unsigned i = g.lastIterPos(); i < sz; i++) {

    if (evals == maxIters) {
      g.lastIterPos(i);
      g.step(0);
      return grade;
    }

    double old_gene  = g.gene(i);
    double old_score = g.score();
    double new_score;

    bool isgene_allele_min = old_gene == allele_min;

    // Only enter this branch if the search was not interrupted during the las call
    if (!isgene_allele_min && !ignoreFirstStep) {
      double new_gene = old_gene - (moveLeft * SR);

      if (new_gene < allele_min)
        new_gene = allele_min;

      g.gene(i, new_gene);

      new_score = g.evaluate();
      evals++;
      fit_inc_acum += g.computeFitnessIncrement(old_score);
      if (GAGenome::compareScores(new_score, old_score) == GAGenome::BETTER)
        improvements++;

      if (GAGenome::compareScores(new_score, fitBest) == GAGenome::BETTER) {
        grade += BONUS1;
        fitBest = new_score;
      }

    }

    if (!isgene_allele_min && !ignoreFirstStep && GAGenome::compareScores(new_score, old_score) == GAGenome::EQUAL) {
      g.gene(i, old_gene);
      g.score(old_score);
    }
    else {
      if (isgene_allele_min || ignoreFirstStep || GAGenome::compareScores(new_score, old_score) == GAGenome::WORSE) {

        g.gene(i, old_gene);
        g.score(old_score);

        if (ignoreFirstStep) {
          g.step(0);
          ignoreFirstStep = false;
        }

        if (evals == maxIters) {
          g.step(1);
          g.lastIterPos(i);
          return grade;
        }

        if (old_gene != allele_max) {
          double new_gene = old_gene + (moveRight * SR);
          if (new_gene > allele_max)
            new_gene = allele_max;

          assert(!isnan(new_gene));

          g.gene(i, new_gene);
          new_score = g.evaluate();

          evals++;
          fit_inc_acum += g.computeFitnessIncrement(old_score);
          if (GAGenome::compareScores(new_score, old_score) == GAGenome::BETTER)
            improvements++;

          if (GAGenome::compareScores(new_score, fitBest  ) == GAGenome::BETTER) {
            grade += BONUS1;
            fitBest = new_score;
          }

          if (GAGenome::compareScores(new_score, old_score) <  GAGenome::BETTER) {
            g.gene(i, old_gene);
            g.score(old_score);
          }
          else {
            grade += BONUS2;
            g.improve(1);
          }
        }
      }
      else {
        grade += BONUS2;
        g.improve(1);
      }
    }

  }

  g.lastIterPos(0);
  g.step(0);

  assert(GAGenome::compareScores(g.score(), fitBest) == GAGenome::EQUAL  ||
         GAGenome::compareScores(g.score(), fitBest) == GAGenome::BETTER    );

  return grade;

}


double MTS_LS1_Reduced_Dims (GAGenome& Xk, double& fitBest, unsigned& evals, double& fit_inc_acum, unsigned& improvements, unsigned maxIters) {
  static const double adjustFailed = GAEDAConfig::handle()->getMTSAdjustFailed();
  static const double adjustMin    = GAEDAConfig::handle()->getMTSAdjustMin   ();
  static const double moveLeft     = GAEDAConfig::handle()->getMTSMoveLeft    ();
  static const double moveRight    = GAEDAConfig::handle()->getMTSMoveRight   ();

  GA1DArrayAlleleGenome<double>& g = dynamic_cast <GA1DArrayAlleleGenome<double>&> (Xk);

  unsigned sz  = g.size();
  double grade = 0;
  double SR    = g.SR();
  fit_inc_acum = 0.0;

  double allele_min = g.alleleset(0).lower();
  double allele_max = g.alleleset(0).upper();

  // If this solution was not improved during the last COMPLETE iteration,
  // update the SR interval
  if (!g.improve()) {
    SR *= adjustFailed;

    if (SR < 1e-14)
      SR = (allele_max - allele_min) * adjustMin;

    g.SR(SR);
  }

  // Mark solution as not improved (only at the beginning of a complete iteration)
  g.improve(0);

  // BEGIN: Cambios Santi
  std::vector<double> improvements_per_dim(sz);

  std::vector<int> dimensions;
  static const double searchProb=GAEDAConfig::handle()->getMTSRedSearchProb();
  static const double minProb   =GAEDAConfig::handle()->getMTSRedMinProb   ();
  ImprovPerDimManager::instance()->getDimsWhichProbSumAThreshold(dimensions, searchProb, sz*minProb);
  std::sort(dimensions.begin(), dimensions.end()); // For some problems, exploring the dims in the ascending order is beneficial
  // END: Cambios Santi

  // For each dimension, perform the LS
  for (unsigned i = 0, dim; i < dimensions.size(); i++) {
    dim = dimensions[i];

    if (evals == maxIters) {
      return grade;
    }

    double old_gene  = g.gene(dim);
    double old_score = g.score();
    double new_score;

    bool isgene_allele_min = old_gene == allele_min;

    // Only enter this branch if the search was not interrupted during the las call
    if (!isgene_allele_min) {
      double new_gene = old_gene - (moveLeft * SR);

      if (isnan(new_gene)) {
        std::cout << "SR=" << SR << std::endl;
        std::cout << "new gene" << new_gene << std::endl;
        exit(-1);
      }

      if (new_gene < allele_min)
        new_gene = allele_min;

      assert(!isnan(new_gene));
      g.gene(dim, new_gene);

      new_score = g.evaluate();
      evals++;
      fit_inc_acum += g.computeFitnessIncrement(old_score);
      if (GAGenome::compareScores(new_score, old_score) == GAGenome::BETTER)
        improvements++;

      if (GAGenome::compareScores(new_score, fitBest) == GAGenome::BETTER) {
        grade += BONUS1;
        fitBest = new_score;
      }
    }

    if (!isgene_allele_min && GAGenome::compareScores(new_score, old_score) == GAGenome::EQUAL) {
      g.gene(dim, old_gene);
      g.score(old_score);
    }
    else {
      if (isgene_allele_min || GAGenome::compareScores(new_score, old_score) == GAGenome::WORSE) {
        g.gene(dim, old_gene);
        g.score(old_score);

        if (evals == maxIters) {
          return grade;
        }

        if (old_gene != allele_max) {

          double new_gene = old_gene + (moveRight * SR);
          if (new_gene > allele_max)
            new_gene = allele_max;

          assert(!isnan(new_gene));

          g.gene(dim, new_gene);
          new_score = g.evaluate();

          evals++;
          fit_inc_acum += g.computeFitnessIncrement(old_score);
          if (GAGenome::compareScores(new_score, old_score) == GAGenome::BETTER)
            improvements++;

          if (GAGenome::compareScores(new_score, fitBest) == GAGenome::BETTER) {
            grade += BONUS1;
            fitBest = new_score;
          }

          if (GAGenome::compareScores(new_score, old_score) < GAGenome::BETTER) {
            g.gene(dim, old_gene);
            g.score(old_score);
          }
          else {
            grade += BONUS2;
            g.improve(1);
          }
        }
      }
      else {
        grade += BONUS2;
        g.improve(1);
      }
    }

    improvements_per_dim[dim] = (GAGenome::compareScores(new_score, old_score) == GAGenome::BETTER) ? fabs(new_score - old_score) : 0.0;
    ImprovPerDimManager::instance()->updateValue(dim, improvements_per_dim[dim]);

  }

  assert(GAGenome::compareScores(g.score(), fitBest) == GAGenome::EQUAL  ||
         GAGenome::compareScores(g.score(), fitBest) == GAGenome::BETTER    );

  return grade;
}


double MTS_LS2 (GAGenome& Xk, double& fitBest, unsigned& evals, double& fit_inc_acum, unsigned maxIters) {

  GA1DArrayAlleleGenome<double>& g = dynamic_cast <GA1DArrayAlleleGenome<double>&> (Xk);

  unsigned sz = g.size();
  double grade = 0;

  double SR = g.SR();

  fit_inc_acum = 0.0;

/*
  std::cout << "    => IMPROVE: ";
  if (g.improve())
     std::cout << "TRUE" << std::endl;
  else
     std::cout << "FALSE" << std::endl;
*/

  // If this solution was not improved the last time, change
  // the SR interval
  if (!g.improve()) {

    SR /= 2;

    //std::cout << "       => NEW SR: " << SR << std::endl;

    if (SR < 1e-14) {
      double upper = g.alleleset (0).upper ();
      double lower = g.alleleset (0).lower ();

      SR = (upper - lower) * 0.4;
      //std::cout << "       => NEW SR: " << SR << std::endl;
    }

    g.SR(SR);

  }

  // Mark solution as not improved
  g.improve(0);

  // For each dimension, perform the LS
  for (unsigned i = 0; i < sz; i++) {

    if (evals == maxIters)
      return grade;

    double old_score = g.score();

    std::vector<std::pair<unsigned,double> > old_genes;

    std::vector<bool> r (sz);
    std::vector<bool> D (sz);

    for (unsigned j = 0; j < sz; j++) {

      if (GARandomDouble(0, 1) < 0.25) {

        r[j] = true;

        if (GARandomInt(0,1))
          D[j] = 1;
        else
          D[j] = -1;

        old_genes.push_back(make_pair(j, g.gene(j)));
        g.gene(j, g.gene(j) - SR * D[j]);

      }
      else
        r[j] = false;

    }

    unsigned nChanges = old_genes.size();

    double new_score = g.evaluate();
    evals++;

    //std::cout << std::endl;
    //std::cout << "    => OLD SCORE: " << old_score << std::endl;
    //std::cout << "    => NEW SCORE: " << new_score << std::endl;

    if (GAGenome::compareScores(new_score, fitBest) == GAGenome::BETTER) {
      grade += BONUS1;
      fitBest = new_score;
      //std::cout << "       => CHANGE ACCEPTED." << std::endl;
      //std::cout << "       => NEW GRADE: " << grade << std::endl;
      //std::cout << "       => NEW BEST FITNESS: " << fitBest << std::endl;
    }

    if (GAGenome::compareScores(new_score, old_score) == GAGenome::EQUAL) {

      for (unsigned j = 0; j < nChanges; j++)
        g.gene(old_genes[j].first, old_genes[j].second);

      //std::cout << "       => CHANGE NOT ACCEPTED. REVERTING..." << std::endl;

    }
    else {

      if (GAGenome::compareScores(new_score, old_score) == GAGenome::WORSE) {

        for (unsigned j = 0; j < nChanges; j++)
          g.gene(old_genes[j].first, old_genes[j].second);

        if (evals == maxIters)
          return grade;

        for (unsigned j = 0; j < nChanges; j++)
          g.gene(old_genes[j].first, g.gene(old_genes[j].first) + 0.5 * SR * D[old_genes[j].first]);

        new_score = g.evaluate();
        evals++;

        //std::cout << "    => NEW SCORE: " << new_score << std::endl;

        if (GAGenome::compareScores(new_score, fitBest) == GAGenome::BETTER) {
          grade += BONUS1;
          fitBest = new_score;
          //std::cout << "       => CHANGE ACCEPTED." << std::endl;
          //std::cout << "       => NEW GRADE: " << grade << std::endl;
          //std::cout << "       => NEW BEST FITNESS: " << fitBest << std::endl;
        }

        if (GAGenome::compareScores(new_score, old_score) < GAGenome::BETTER) {
          for (unsigned j = 0; j < nChanges; j++)
            g.gene(old_genes[j].first, old_genes[j].second);
          g.score(old_score);
          //std::cout << "       => CHANGE NOT ACCEPTED. REVERTING..." << std::endl;
        }
        else {
          grade += BONUS2;
          g.improve(1);
          //std::cout << "       => CHANGE ACCEPTED." << std::endl;
          //std::cout << "       => NEW GRADE: " << grade << std::endl;
        }

      }
      else {
        grade += BONUS2;
        g.improve(1);
        //std::cout << "       => CHANGE ACCEPTED." << std::endl;
        //std::cout << "       => NEW GRADE: " << grade << std::endl;
      }

    }

  }

  return grade;

}

double MTS_LS3 (GAGenome& Xk, double& fitBest, unsigned& evals, double& fit_inc_acum, unsigned maxIters) {

  GA1DArrayAlleleGenome<double>& g = dynamic_cast <GA1DArrayAlleleGenome<double>&> (Xk);

  // Store a copy of the original individual to restore it if not improvement
  // has been obtaine at the end of the LS
  GA1DArrayAlleleGenome<double>* X_copy = dynamic_cast <GA1DArrayAlleleGenome<double>*> (g.clone());

  unsigned sz = g.size();
  double grade = 0.0;
  double old_score = g.score();
  double X_score = old_score;

  fit_inc_acum = 0.0;

  for (unsigned i = 0; i < sz; i++) {

    if (evals == maxIters)
      return grade;

    double old_gene = g.gene(i);

    // X1
    g.gene(i, old_gene + 0.1);
    double X1_score = X_score = g.evaluate();
    evals++;

    if (evals == maxIters)
      return grade;

    // Y1
    g.gene(i, old_gene - 0.1);
    double Y1_score = X_score = g.evaluate();
    evals++;

    if (evals == maxIters)
      return grade;

    // X2
    g.gene(i, old_gene + 0.2);
    double X2_score = X_score = g.evaluate();
    evals++;

    // Check BONUS1 for X1
    if (GAGenome::compareScores(X1_score, fitBest) == GAGenome::BETTER) {

      grade += BONUS1;
      fitBest = X1_score;

    }

    // Check BONUS1 for Y1
    if (GAGenome::compareScores(Y1_score, fitBest) == GAGenome::BETTER) {

      grade += BONUS1;
      fitBest = Y1_score;

    }

    // Check BONUS1 for X2
    if (GAGenome::compareScores(X2_score, fitBest) == GAGenome::BETTER) {

      grade += BONUS1;
      fitBest = X2_score;

    }

    // Check BONUS2 for X1
    double D1 = fabs(X_score - X1_score);
    if (GAGenome::compareScores(X1_score, X_score) == GAGenome::BETTER) // X1 is better than X
      grade += BONUS2;
    else
      D1 *= -1; // This way we respect the original behavior of the LS, regardless the problem is minimization or maximization

    // Check BONUS2 for Y1
    double D2 = fabs(X_score - Y1_score);
    if (GAGenome::compareScores(Y1_score, X_score) == GAGenome::BETTER) // Y1 is better than X
      grade += BONUS2;
    else
      D1 *= -1; // This way we respect the original behavior of the LS, regardless the problem is minimization or maximization

    // Check BONUS2 for X2
    double D3 = fabs(X_score - X2_score);
    if (GAGenome::compareScores(X2_score, X_score) == GAGenome::BETTER) // X2 is better than X
      grade += BONUS2;
    else
      D1 *= -1; // This way we respect the original behavior of the LS, regardless the problem is minimization or maximization

    double a = GARandomDouble(0.4, 0.5);
    double b = GARandomDouble(0.1, 0.3);
    double c = GARandomDouble(0, 1);

    if (evals == maxIters)
      return grade;

    g.gene(i, old_gene + (a * (D1 - D2)) + (b * (D3 - 2 * D1)) + c);

    X_score = g.evaluate();
    evals++;

  }

  if (GAGenome::compareScores(X_score, old_score) < GAGenome::BETTER)
    g.copy(*X_copy);
  else
    grade+=BONUS2;

  delete X_copy;

  return grade;

}


double RandomGreedyLS (GAGenome& Xk, double& fitBest, unsigned& evals, double& fit_inc_acum, unsigned& improvements, double& SR, unsigned maxIters) {

  GA1DArrayAlleleGenome<double>& g = dynamic_cast<GA1DArrayAlleleGenome<double>&>(Xk);

  int      numdims   = g.length();

  fit_inc_acum = 0.0;

  while (evals < maxIters) {

    double old_value, new_value, old_score;

    int pos = GARandomInt(0, numdims-1);

    double min_value = g.alleleset(pos).lower();
    double max_value = g.alleleset(pos).upper();

    old_value = g.gene(pos);
    old_score = g.score();

    new_value = GARandomDouble(min_value, max_value);
    g.gene(pos , new_value);

    double new_score = g.score();
    fit_inc_acum += g.computeFitnessIncrement(old_score);

    if (GAGenome::compareScores(new_score, old_score) == GAGenome::BETTER)
      improvements++;
    else {
      g.gene(pos, old_value);
      g.score(old_score);
    }

    evals++;

  }

  return fit_inc_acum;

}

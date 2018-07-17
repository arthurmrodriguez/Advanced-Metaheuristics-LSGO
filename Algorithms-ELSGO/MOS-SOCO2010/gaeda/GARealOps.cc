#include "GARealOps.h"

#include "GACommonOps.h"
#include "genomes/GAGenome.h"
#include "genomes/GA1DArrayGenome.h"

// Initializers
void RealUniformInitializer (GAGenome& g) {

   GA1DArrayAlleleGenome<long double>& gen = DYN_CAST (GA1DArrayAlleleGenome<long double>&, g);

   for (int i = 0; i < gen.length (); i++) {

      long double upper = gen.alleleset (i).upper ();
      long double lower = gen.alleleset (i).lower ();

      gen.gene (i, GARandomDouble (lower, upper));

   }

   return;

}

// Comparators
long double RealEuclideanComparator (const GAGenome& a, const GAGenome& b) {

   const GA1DArrayGenome<long double>& sis = DYN_CAST (const GA1DArrayGenome<long double>&, a);
   const GA1DArrayGenome<long double>& bro = DYN_CAST (const GA1DArrayGenome<long double>&, b);

   if (sis.length () != bro.length ()) return -1;
   if (sis.length () == 0            ) return  0;

   long double xn, yn, tmp;
   long double result = 0.0;

   for (int i=0; i < sis.length(); i++) {

      xn      = sis.gene (i);
      yn      = bro.gene (i);
      tmp     = xn  -  yn;
      result += tmp * tmp;

   }

   return sqrt (result);

}

// GA Crossovers
int RealBlendCrossover(const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {

  const GA1DArrayAlleleGenome<long double> &mom = DYN_CAST(const GA1DArrayAlleleGenome<long double> &, p1);
  const GA1DArrayAlleleGenome<long double> &dad = DYN_CAST(const GA1DArrayAlleleGenome<long double> &, p2);

  int n=0;

  if(c1 && c2){
    GA1DArrayGenome<long double> &sis=DYN_CAST(GA1DArrayGenome<long double> &, *c1);
    GA1DArrayGenome<long double> &bro=DYN_CAST(GA1DArrayGenome<long double> &, *c2);

    int len = GAMax(mom.length(), dad.length());
    for(int i=0; i<len; i++) {
      long double dist = 0;
      if(mom.gene(i) > dad.gene(i))
        dist = mom.gene(i) - dad.gene(i);
      else
        dist = dad.gene(i) - mom.gene(i);

      long double lo = (GAMin(mom.gene(i), dad.gene(i))) - 0.5*dist;
      if (lo < mom.alleleset(i).lower())
        lo = mom.alleleset(i).lower();    // these were added to assure that the limits are not passed

      long double hi = (GAMax(mom.gene(i), dad.gene(i))) + 0.5*dist;
      if (hi > mom.alleleset(i).upper())
        hi = mom.alleleset(i).upper();// these were added to assure that the limits are not passed

      sis.gene(i, GARandomDouble(lo, hi));
      bro.gene(i, GARandomDouble(lo, hi));
    }
    n = 2;
  }
  else if(c1 || c2){
    GA1DArrayGenome<long double> &sis = (c1 ?
				   DYN_CAST(GA1DArrayGenome<long double> &, *c1) :
				   DYN_CAST(GA1DArrayGenome<long double> &, *c2));

    int len = GAMax(mom.length(), dad.length());
    for(int i=0; i<len; i++) {
      long double dist = 0;
      if(mom.gene(i) > dad.gene(i)) 	dist = mom.gene(i) - dad.gene(i);
      else                          	dist = dad.gene(i) - mom.gene(i);

      long double lo = (GAMin(mom.gene(i), dad.gene(i))) - 0.5*dist;
      if (lo < mom.alleleset(i).lower())
        lo = mom.alleleset(i).lower();   // these were added to assure that the limits are not passed
      long double hi = (GAMax(mom.gene(i), dad.gene(i))) + 0.5*dist;
      if (hi > mom.alleleset(i).upper())
        hi = mom.alleleset(i).upper();   // these were added to assure that the limits are not passed
      sis.gene(i, GARandomDouble(lo, hi));
    }
    n = 1;
  }
  return n;
}

// GA Mutators
int RealGaussianMutator (GAGenome& g, long double pmut) {

  GA1DArrayAlleleGenome<long double>& child = DYN_CAST (GA1DArrayAlleleGenome<long double>&, g);
  //register int i,n;
  if(pmut <= 0.0) return(0);

  long double nMut       = pmut * (long double)(child.length());
  int    length     = child.length()-1;
  long double new_value  = 0.0;

  if (GAFlipCoin(nMut - (int)nMut) ) nMut++;

  int    idx;
  long double gene_value, allele_region;

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


// DE Crossover
int RealExponentialCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, const long double prob) {return ExponentialCrossover<long double> (p1, p2, c1, prob);}

// MTS Local Searches
long double MTS_LS1 (GAGenome& Xk, long double& fitBest, unsigned& evals, long double& fit_inc_acum, unsigned& improvements, long double& SR, unsigned maxIters) {

  GA1DArrayAlleleGenome<long double>& g = dynamic_cast <GA1DArrayAlleleGenome<long double>&> (Xk);

  unsigned sz  = g.size();
  long double grade = 0;
  fit_inc_acum = 0.0;

  // If this solution was not improved during the last COMPLETE iteration,
  // update the SR interval
  if (!g.improve() && g.lastIterPos() == 0 && g.step() == 0) {

    SR /= 2;

    if (SR < 1e-14) {
      long double upper = g.alleleset (0).upper ();
      long double lower = g.alleleset (0).lower ();

      SR = (upper - lower) * 0.4;
    }

    g.SR(SR);

  }

  // Mark solution as not improved (only at the beginning of a complete iteration)
  if (g.lastIterPos() == 0 && g.step() == 0)
    g.improve(0);

  long double allele_min = g.alleleset(0).lower();
  long double allele_max = g.alleleset(0).upper();

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

    long double old_gene  = g.gene(i);
    long double old_score = g.score();
    long double new_score;

    bool isgene_allele_min = old_gene == allele_min;

    // Only enter this branch if the search was not interrupted during the las call
    if (!isgene_allele_min && !ignoreFirstStep) {

      long double new_gene = old_gene - SR;

      if (new_gene < allele_min)
        new_gene = allele_min;

      g.gene(i, new_gene);

      new_score = g.evaluate();
      evals++;
      fit_inc_acum += g.computeFitnessIncrement(old_score);
      if (GAGenome::compareScores(new_score, old_score) == GAGenome::BETTER) {
        improvements++;
      }

      if (GAGenome::compareScores(new_score, fitBest) == GAGenome::BETTER) {
        grade += BONUS1;
        fitBest = new_score;
      }

    }

    if (!isgene_allele_min && GAGenome::compareScores(new_score, old_score) == GAGenome::EQUAL && !ignoreFirstStep) {

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

          long double new_gene = old_gene + (0.5 * SR);
          if (new_gene > allele_max)
            new_gene = allele_max;

          assert(!isnan(new_gene));

          g.gene(i, new_gene);
          new_score = g.evaluate();

          evals++;
          fit_inc_acum += g.computeFitnessIncrement(old_score);
          if (GAGenome::compareScores(new_score, old_score) == GAGenome::BETTER) {
            improvements++;
          }

          if (GAGenome::compareScores(new_score, fitBest) == GAGenome::BETTER) {
            grade += BONUS1;
            fitBest = new_score;
          }

          if (GAGenome::compareScores(new_score, old_score) < GAGenome::BETTER) {

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

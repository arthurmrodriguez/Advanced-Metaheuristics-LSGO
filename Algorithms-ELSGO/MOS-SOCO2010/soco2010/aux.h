#ifndef AUX_H_
#define AUX_H_

#include <GARealOps.h>
#include <MOSGenomeFactory.h>
#include <GAEDAConfig.h>
#include <GAGeneticAlgorithm.h>
#include <genomes/GA1DArrayGenome.h>

void   setInitValues(GAGenome& ind, long double param_m);
long double composeFuncs (GAGenome& ind, long double (*f1)(GAGenome&), long double bias1, long double (*f2)(GAGenome&), long double bias2, long double* trans);

#define CREATEHYBRIDFUNC( NAME, FUNC1, FUNC2, M, MIN_VALUE, MAX_VALUE, TRANS) \
\
   GA1DArrayAlleleGenome<long double>* optimum_genome = NULL;\
\
   long double NAME (GAGenome& g) {\
      return composeFuncs(g, FUNC1, FUNC1##_BIAS, FUNC2, FUNC2##_BIAS, TRANS);\
   }\
\
   extern "C" long double objective (GAGenome& g) {\
      return NAME(g);\
   }\
\
   extern "C" void individualInit(GAGenome & g) {\
      return RealUniformInitializer (g);\
   }\
\
   extern "C" GAGenome* defineProblem (int size, char *data) {\
      GAEDAConfig* cfg = GAEDAConfig::handle();\
\
      GAAlleleSet<long double> alleles (MIN_VALUE, MAX_VALUE);\
      GA1DArrayAlleleGenome<long double>* genome = new GA1DArrayAlleleGenome<long double> (cfg->getProblemSize (), alleles, objective);\
\
      genome->initializer (RealUniformInitializer);\
      genome->comparator  (RealEuclideanComparator);\
\
      genome->crossover   (RealBlendCrossover);\
      genome->mutator     (RealGaussianMutator);\
\
      genome->crossover   (RealExponentialCrossover);\
\
      MOSGenomeFactory::handle()->registerGenome (GAID::RealEncoding, genome);\
\
      optimum_genome = new GA1DArrayAlleleGenome<long double> (cfg->getProblemSize (), alleles, objective);\
      optimum_genome->comparator (RealEuclideanComparator);\
      optimum_genome->crossover  (RealBlendCrossover);\
\
      for (int i = 0; i < optimum_genome->length (); i++)\
         optimum_genome->gene (i, 0);\
\
      setInitValues(*genome, M);\
\
      return genome;\
   }\
\
   extern "C" bool postprocess (GAPopulation* pop, int rank) {\
      delete optimum_genome;\
      return true;\
   }\
\
   extern "C" const char* describeProblem (void) {\
      char* msg = new char[100];\
      sprintf (msg, "%s Hybrid from %s and %s", #NAME, #FUNC1, #FUNC2);\
      return msg;\
   }\
\
   extern "C" GAGenome* getOptimumGenome () {\
      return optimum_genome;\
   }\
\
   extern "C" GAGenome::OptCriterion optCriterion(){\
      return GAGenome::MINIMIZATION;\
   }

#endif /* AUX_H_ */

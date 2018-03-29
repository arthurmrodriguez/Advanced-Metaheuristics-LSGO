#include "QualityLogStat.h"

#include "../GAGeneticAlgorithm.h"
#include "../MOSTechniqueSet.h"

QualityLogStat::QualityLogStat(const Algorithm& alg) : SetLogStat("quality",alg) {
   _setSize = MOSTechniqueSet::handle()->nTechniques();
}

void QualityLogStat::computeValues(vector<long double>* v) {
   const GAGeneticAlgorithm& alg = DYN_CAST(const GAGeneticAlgorithm&, alg_);

   if (alg.generation () == 0)
      for (unsigned i=0; i<_setSize; i++)
         (*v)[i]=0;
   else
      MOSTechniqueSet::handle()->getQualityVector(v);
}

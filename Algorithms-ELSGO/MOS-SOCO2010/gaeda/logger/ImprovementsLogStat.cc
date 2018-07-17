#include "ImprovementsLogStat.h"

#include "../GAGeneticAlgorithm.h"
#include "../MOSTechniqueSet.h"
#include "../MOSTechnique.h"

ImprovementsLogStat::ImprovementsLogStat(const Algorithm& alg) : SetLogStat("improvements",alg) {
   _setSize = MOSTechniqueSet::handle()->nTechniques();
}

void ImprovementsLogStat::computeValues(vector<long double>* v) {
   const GAGeneticAlgorithm& alg = DYN_CAST(const GAGeneticAlgorithm&, alg_);

   if (alg.generation () == 0)
      for (unsigned i=0; i<_setSize; i++)
         (*v)[i]=0;
   else {
     MOSTechniqueSet::MOSTechniqueSetIterator it = MOSTechniqueSet::handle()->begin();
     unsigned i = 0;
     for (;it != MOSTechniqueSet::handle()->end(); it++, i++)
       (*v)[i] = it->second->getRatioImproved();
   }
}

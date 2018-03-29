#include "ParticipationLogStat.h"

#include "../GAGeneticAlgorithm.h"
#include "../MOSTechniqueSet.h"

ParticipationLogStat::ParticipationLogStat(const Algorithm& alg) : SetLogStat("participation",alg) {
   _setSize = MOSTechniqueSet::handle()->nTechniques();
}

void ParticipationLogStat::computeValues(vector<double>* v) {
   const GAGeneticAlgorithm& alg = DYN_CAST(const GAGeneticAlgorithm&, alg_);

   if (alg.generation () == 0)
      for (unsigned i=0; i<_setSize; i++)
         (*v)[i]=0;
   else
      MOSTechniqueSet::handle()->getPartRatiosVector(v);
}

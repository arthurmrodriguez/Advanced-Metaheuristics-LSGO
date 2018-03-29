#include "FSSLogStat.h"

#include "../GAPopulation.h"
#include "../genomes/MOSGenome.h"
#include "../genomes/FSSGenome.h"

FSSLogStat::FSSLogStat(const Algorithm& alg, int type) : SetLogStat(type == AVG ? "fss-avg" : (type == BEST ? "fss-best" : "fss-worst"), alg) {
   _setSize = 5;
   _type = type;
}

void FSSLogStat::computeValues(vector<long double>* v) {

   const GAPopulation& pop = alg_.population();

   if (_type == AVG) {

      unsigned size = pop.size();

      for (unsigned i = 0; i < size; i++) {

         const MOSGenome&      mosGen = DYN_CAST (const MOSGenome&, pop.individual(i));
         const FSSGenome<int>&    gen = DYN_CAST (const FSSGenome<int>&, *(mosGen.getDefaultGenome()));

         (*v)[0] += gen.getAUC();
         (*v)[1] += gen.getLL ();
         (*v)[2] += gen.getLLS();
         (*v)[3] += gen.getBS ();
         (*v)[4] += gen.getACC();

      }

      (*v)[0] /= size;
      (*v)[1] /= size;
      (*v)[2] /= size;
      (*v)[3] /= size;
      (*v)[4] /= size;

   }
   else if (_type == BEST) {

      const MOSGenome&      mosGen = DYN_CAST (const MOSGenome&, pop.best());
      const FSSGenome<int>&    gen = DYN_CAST (const FSSGenome<int>&, *(mosGen.getDefaultGenome()));

      (*v)[0] += gen.getAUC();
      (*v)[1] += gen.getLL ();
      (*v)[2] += gen.getLLS();
      (*v)[3] += gen.getBS ();
      (*v)[4] += gen.getACC();

   }
   else if (_type == WORST) {

      const MOSGenome&      mosGen = DYN_CAST (const MOSGenome&, pop.worst());
      const FSSGenome<int>&    gen = DYN_CAST (const FSSGenome<int>&, *(mosGen.getDefaultGenome()));

      (*v)[0] += gen.getAUC();
      (*v)[1] += gen.getLL ();
      (*v)[2] += gen.getLLS();
      (*v)[3] += gen.getBS ();
      (*v)[4] += gen.getACC();

   }

}

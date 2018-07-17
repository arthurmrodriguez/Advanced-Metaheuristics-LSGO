#ifndef MOSDEFINEPROBLEM_H
#define MOSDEFINEPROBLEM_H

#include <GAPopulation.h>
#include <genomes/MOSGenome.h>
#include <GAEDAConfig.h>
#include <MOSEA.h>

MOSEA* genericDefineMOSProblemCentral (unsigned popSize) {

   MOSTechniqueSet* techniqueSet = MOSTechniqueSet::handle ();
   unsigned nTechs = techniqueSet->nTechniques ();

   // Check for a good set of techniques to use
   if (nTechs < 1) {
      std::cerr << "Error: wrong selection of techniques. Aborting." << std::endl;
      exit (-1);
   }

   GAPopulation* pop = new GAPopulation ();

   unsigned sizePerTech = popSize / nTechs;
   unsigned extraIndivs = popSize - (sizePerTech * nTechs);

   // Inicializar individuos de la poblacion por tecnica
   MOSTechniqueSet::MOSTechniqueSetIterator it;
   for (it=techniqueSet->begin(); it!=techniqueSet->end(); it++) {
      MOSGenome* g = new MOSGenome(it->second);
      pop->add (*g, sizePerTech);
      delete g;
   }

   // Individuos extra por el redondeo
   unsigned i = 0;
   for (it=techniqueSet->begin(); it!=techniqueSet->end() && i<extraIndivs; it++,i++) {
      MOSGenome* g = new MOSGenome(it->second);
      pop->add(g);
   }

   GAEDAConfig* cfg = GAEDAConfig::handle();

   MOSParticipation* part = new DynamicParticipation();
   MOSQuality*       qual = cfg->getQualityFunction();

   MOSEA* ga = new MOSEA(*pop, part, qual);
   ga->populationSize(popSize);

   delete pop;

   return ga;

}

#endif

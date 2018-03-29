#ifndef MOSDEFINEPROBLEM_H
#define MOSDEFINEPROBLEM_H

#include <GAPopulation.h>
#include <MOSTechnique.h>
#include <genomes/MOSGenome.h>
#include <MOSParticipationFunc.h>
#include <MOSParticipationFunction.h>
#include <GAEDAConfig.h>
#include <MOSEAMultiDeme.h>
#include <MOSEA2.h>
#include <MOSEARL.h>

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

   if (cfg->getPopInitFunction())
      pop->initializer(cfg->getPopInitFunction());

   MOSEA* ga = new MOSEA(*pop, cfg->getElisitsmPercent(), cfg->getMinPart());
   ga->populationSize(popSize);
   ga->setParticipationFunction(cfg->getPartFunction());

   delete pop;

   return ga;

}

MOSEA2* genericDefineMOS2ProblemCentral (unsigned popSize) {

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

   if (cfg->getPopInitFunction())
      pop->initializer(cfg->getPopInitFunction());

   MOSParticipation* part = cfg->getParticipationFunction();
   MOSQuality*       qual = cfg->getQualityFunction();

   MOSEA2* ga = new MOSEA2(*pop, part, qual);
   ga->populationSize(popSize);

   delete pop;

   return ga;

}

MOSEA2* genericDefineMOS2ProblemCentralReusing (unsigned popSize, GAPopulation* pop) {

   MOSTechniqueSet* techniqueSet = MOSTechniqueSet::handle ();
   unsigned nTechs = techniqueSet->nTechniques ();

   // Check for a good set of techniques to use
   if (nTechs < 1) {
      std::cerr << "Error: wrong selection of techniques. Aborting." << std::endl;
      exit (-1);
   }

   GAEDAConfig* cfg = GAEDAConfig::handle();

   pop->initializer(GAPopulation::PercInitializer);
   pop->initPercentage(cfg->getRerunPercentage());

   if (cfg->getPopInitFunction())
      pop->initializer(cfg->getPopInitFunction());

   MOSParticipation* part = cfg->getParticipationFunction();
   MOSQuality*       qual = cfg->getQualityFunction();

   MOSEA2* ga = new MOSEA2(*pop, part, qual);
   ga->populationSize(popSize);

   if (pop->size() > popSize)
     std::cerr << "Warning: population shrinking after re-executing the algorithm..." << std::endl;
   else if (pop->size() < popSize)
     std::cerr << "Warning: population growing after re-executing the algorithm..." << std::endl;

   return ga;

}

MOSEA* genericDefineMOSProblemAutonomic (unsigned popSize) {

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
   for (it = techniqueSet->begin (); it != techniqueSet->end (); it++) {
      MOSGenome* g = new MOSGenome (it->second);
      pop->add (*g, sizePerTech);
      delete g;
   }

   // Individuos extra por el redondeo
   unsigned i = 0;
   for (it = techniqueSet->begin (); it != techniqueSet->end () && i < extraIndivs; it++, i++) {
      MOSGenome* g = new MOSGenome (it->second);
      pop->add (g);
   }

   GAEDAConfig* cfg = GAEDAConfig::handle();

   if (cfg->getPopInitFunction())
      pop->initializer(cfg->getPopInitFunction());

   MOSEA* ga = new MOSEA(*pop, cfg->getElisitsmPercent(), cfg->getMinPart(), AutonomicEvolution);

   ga->populationSize(popSize);

   delete pop;

   return ga;

}

MOSEAMultiDeme* genericDefineMOSMultiDemeProblemCentral (unsigned popSize) {

   MOSTechniqueSet* techniqueSet = MOSTechniqueSet::handle ();
   unsigned nTechs = techniqueSet->nTechniques ();

   std::vector<MOSGenome*> genomes (nTechs, (MOSGenome*)0);

   // Check for a good set of techniques to use
   if (nTechs < 1) {
      std::cerr << "Error: wrong selection of techniques. Aborting." << std::endl;
      exit (-1);
   }

   // Inicializar individuos por tecnica
   MOSTechniqueSet::MOSTechniqueSetIterator it;
   for (it=techniqueSet->begin(); it!=techniqueSet->end(); it++)
     genomes[it->first] = new MOSGenome(it->second);

   GAEDAConfig* cfg = GAEDAConfig::handle();

   MOSEAMultiDeme* ga = new MOSEAMultiDeme(genomes, cfg->getElisitsmPercent(), cfg->getMinPart());

   ga->populationSize(popSize);
   ga->setParticipationFunction(cfg->getPartFunction());

   for (unsigned i = 0; i < genomes.size(); i++)
     delete genomes[i];

   return ga;

}

MOSEARL* genericDefineMOSProblemRL (unsigned popSize) {

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

   if (cfg->getPopInitFunction())
      pop->initializer(cfg->getPopInitFunction());

   MOSEARL* ga = new MOSEARL(*pop, cfg->getElisitsmPercent());
   ga->populationSize(popSize);

   delete pop;

   return ga;

}

#endif

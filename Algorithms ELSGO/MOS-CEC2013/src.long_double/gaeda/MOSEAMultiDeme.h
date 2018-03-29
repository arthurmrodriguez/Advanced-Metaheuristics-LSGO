#ifndef MOSEAMultiDeme_H
#define MOSEAMultiDeme_H

#include <map>

#include "GAGeneticAlgorithm.h"
#include "MOSConfig.h"
#include "MOSEA.h"

/**
 * Class: MOSEAMultiDeme
 *
 * @brief Description of MOSEAMultiDeme class
 */

class NSC;
class MOSGenome;
class MOSConversion;
class MOSTechnique;
class MOSTechniqueSet;

class MOSEAMultiDeme : public MOSEA {

 public:

  GADefineIdentity("MOSEAMultiDeme", GAID::MOSMultiDeme);

  // Constructors and destructor
  MOSEAMultiDeme(const std::vector<MOSGenome*>& genome, const long double elitismPercent, const long double minPart, const evolutionType evolutiveApproach = CentralEvolution);

  ~MOSEAMultiDeme();

  // Initialization of the algorithm
  void initialize();

  void step();

 protected:

  // Attributes
  std::vector<GAPopulation*> _pops;
  std::vector<GAPopulation*> _auxPops;

};

#endif

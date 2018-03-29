#ifndef MOSEA2_H
#define MOSEA2_H

#include "GAGeneticAlgorithm.h"
#include "MOSConfig.h"

#include <vector>

/**
 * Class: MOSEA2
 *
 * @brief Description of MOSEA2 class
 */

class MOSParticipation;
class MOSTechnique;

#include "MOSParticipationFunction.h"

class MOSEA2 : public GAGeneticAlgorithm {

 public:

  GADefineIdentity("MOSEA2", GAID::MOS2);

  // Constructors and destructor
  MOSEA2(const GAGenome& genome,         MOSParticipation* part, MOSQuality* qual);
  MOSEA2(const GAPopulation& population, MOSParticipation* part, MOSQuality* qual);
  MOSEA2(const GAGeneticAlgorithm& alg,  MOSParticipation* part, MOSQuality* qual);

  ~MOSEA2();

  // Initialization of the algorithm
  void initialize();

  // Carry out a step of the algorithm
  void step();

  // Creates a new population (not used but must be defined as it is abstract in base class)
  void offspring(GAPopulation* offpop) {};

  // Population size
  virtual unsigned populationSize(unsigned size)       {return GAGeneticAlgorithm::populationSize(size);}
  virtual unsigned populationSize(             ) const {return pop->size();                             }

  // Set and get PF
  MOSParticipation* setParticipationFunction(MOSParticipation* partFunction)       {return _partFunction = partFunction;}
  MOSParticipation* getParticipationFunction(                              ) const {return _partFunction;               }

  // Set and get QF
  MOSQuality* setQualityFunction(MOSQuality* qualFunction)       {return _qualFunction = qualFunction;}
  MOSQuality* getQualityFunction(                        ) const {return _qualFunction;               }

  GAScalingScheme&   scaling  (const GAScalingScheme&   s) {return GAGeneticAlgorithm::scaling(s); }
  GASelectionScheme& selector (const GASelectionScheme& s) {return GAGeneticAlgorithm::selector(s);}

 protected:

  // Attributes
  MOSParticipation* _partFunction; // Participation function used
  MOSQuality*       _qualFunction; // Quality function used

  std::vector<MOSTechnique*> _additionalLS;

};

#endif

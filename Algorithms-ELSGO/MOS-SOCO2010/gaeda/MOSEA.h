#ifndef MOSEA_H
#define MOSEA_H

#include "GAGeneticAlgorithm.h"
#include "MOSConfig.h"

#include <vector>

/**
 * Class: MOSEA
 *
 * @brief Description of MOSEA class
 */

class MOSParticipation;
class MOSTechnique;

#include "MOSParticipationFunction.h"

class MOSEA : public GAGeneticAlgorithm {

 public:

  GADefineIdentity("MOSEA", GAID::MOS);

  // Constructors and destructor
  MOSEA(const GAGenome& genome,         MOSParticipation* part, MOSQuality* qual);
  MOSEA(const GAPopulation& population, MOSParticipation* part, MOSQuality* qual);
  MOSEA(const GAGeneticAlgorithm& alg,  MOSParticipation* part, MOSQuality* qual);

  ~MOSEA();

  // Initialization of the algorithm
  void initialize(unsigned seed = 0);

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

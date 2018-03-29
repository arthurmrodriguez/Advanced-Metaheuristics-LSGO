#ifndef MOSEA_H
#define MOSEA_H

#include <map>

#include "GAGeneticAlgorithm.h"
#include "MOSConfig.h"

/**
 * Class: MOSEA
 *
 * @brief Description of MOSEA class
 */

class NSC;
class MOSConversion;
class MOSTechnique;
class MOSTechniqueSet;

class MOSEA : public GAGeneticAlgorithm {

 public:

  GADefineIdentity("MOSEA", GAID::MOS);

  // Constructors and destructor
  MOSEA(const GAGenome& genome, const double elitismPercent, const double minPart, const evolutionType evolutiveApproach = CentralEvolution);
  MOSEA(const GAPopulation& population, const double elitismPercent, const double minPart, const evolutionType evolutiveApproach = CentralEvolution);
  MOSEA(const GAGeneticAlgorithm& alg, const double elitismPercent, const double minPart, const evolutionType evolutiveApproach = CentralEvolution);

  ~MOSEA();

  // Initialization of the algorithm
  void initialize();

  // Wrapper that actually calls stepCentral() or stepAutonomic()
  void step();

  // Creates a new population (not used but must be defined as it is abstract in base class)
  void offspring(GAPopulation* offpop) {};

  // Returns a reference to the set of techniques used
  MOSTechniqueSet& getTechSet() const {return (MOSTechniqueSet&)*_techniqueSet;}

  // Population size
  virtual unsigned populationSize(unsigned size);
  virtual unsigned populationSize () const;

  // Set and get PF data
  void  setParticipationFunctionData(void* data) {_PFData = data;}
  void* getParticipationFunctionData() {return _PFData;}

  // Set and get PF
  void setParticipationFunction(ParticipationFunction partFunction) {_partFunction = partFunction; return;}
  ParticipationFunction getParticipationFunction() {return _partFunction;}

  // Restore Generation
  unsigned getParticipationRestoreGen() {return _partRestoreGen;}
  void setParticipationRestoreGen(unsigned int generation) {_partRestoreGen = generation;}

  // Actual population size
  unsigned getActualPopSize() {return _actualPopSize;}
  void setActualPopSize(unsigned int popsize) {_actualPopSize = popsize;}

  // Minimum participation
  double getMinParticipation(void) {return _minPart;}

  const std::string& getQualityMeasure (void) const {return qualityMeasure;}

  unsigned offPopulationSize (unsigned value) {return pop->size (value);}

  GAScalingScheme&   scaling  (const GAScalingScheme&   s);
  GASelectionScheme& selector (const GASelectionScheme& s);

 protected:

  // Evolves a generation using the central approach
  GAPopulation* stepCentral(unsigned maxEvals);

  // Evolves a generation using the autonomic approach
  GAPopulation* stepAutonomic(unsigned maxEvals);

  // Update participation ratios
  void updatePartRatios();

  // Evaluate quality of techniques (for central approach only)
  void evalTechQuality(GAPopulation*);

  void fitnessAverageQuality         (GAPopulation* evalPop);
  void fitnessIncrementAverageQuality(GAPopulation* evalPop);
  void NSCQuality                    (GAPopulation* evalPop);
  void compassQuality                (GAPopulation* evalPop);
  void diversityAvgQuality           (GAPopulation* evalPop);

  std::map<techIdType,double> computeDiversitiesAvg(GAPopulation* evalPop);

  // Select a technique based on the probability vector passed (autonomic approach only)
  MOSTechnique* selectTechnique(MOSProbVector& probVector);

  // Attributes
  double _bonus;                       // Coefficient to adapt the participation of each technique (in central approach only)
  GAPopulation* _auxPop;               // Aux population to store new offspring
  MOSTechniqueSet* _techniqueSet;      // Set of techniques used in the process
  evolutionType _evolutiveApproach;    // Evolutive approach used (central or autonomic)
  MOSConversion* _conversion;          // Converts among different encodings
  unsigned _partRestoreGen;            // Participation restore generation
  ParticipationFunction _partFunction; // Participation function used
  void* _PFData;                       // Data used by the participation function
  double _minPart;                     // Minimal participation for a technique
  double _elitismPercent;              // The percent of the old Pop that are brought into the new Pop
  std::string qualityMeasure;
  NSC* _nsc;                           // Object to compute and store NSC information
  unsigned _actualPopSize;

};

#endif

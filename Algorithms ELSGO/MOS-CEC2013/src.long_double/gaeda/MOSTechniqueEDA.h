/**
 * @file
 * @brief MOSTechniqueEDA class hdr.
 *
 */

#ifndef MOSTechniqueEDA_H
#define MOSTechniqueEDA_H


#include "MOSTechnique.h"
#include "GAGraphModel.h"
#include "GABayesianNetwork.h"
#include "GAGaussianNetwork.h"

/**
 * @brief Clase que representa las técnicas de tipo EDA
 */
class MOSTechniqueEDA : public MOSTechnique {

 public:

  GADefineIdentity("MOSTechniqueEDA", GAID::TechniqueEDA);

  MOSTechniqueEDA(techIdType id, std::string description,
		    long double selPerc, GAGraphModel::ModelType networkType,
		    GABayesianNetwork::LearningMethod learnMethod,
		    GABayesianNetwork::EBNALocalScoring localScoring,
		    GABayesianNetwork::SimulationMethod simMethod,
		    encodingType encoding, GAGenome* genomeBase,
		    GAGenome::Initializer init, GAGenome::Evaluator evaluator);

  MOSTechniqueEDA(techIdType id, std::string description,
		    long double selPerc, GAGraphModel::ModelType networkType,
		    GAGaussianNetwork::LearningMethod learnMethod,
		    GAGaussianNetwork::ScoreMethod scoreMethod,
		    encodingType encoding, GAGenome* genomeBase,
		    GAGenome::Initializer init, GAGenome::Evaluator evaluator);

  ~MOSTechniqueEDA();

  // Produces a new offspring from a population
  void offspring (GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset);

  unsigned offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged) {return 0;}

  // Produces a new child from learnt distribution
  // De momento no se implementa en EDA's
  //bool offspring (const MOSGenome& dad, const MOSGenome& mom, MOSGenome* child) const;

  // For learning a distribution
  //virtual bool scan (const GAPopulation& pop) const;

 protected:

  // The Bayesian network estimated from the selected individuals.
  GAGraphModel* _network;
  long double _selPerc;
};

#endif

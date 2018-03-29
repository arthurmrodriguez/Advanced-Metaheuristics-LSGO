/**
 * @file
 * @brief MOSEARL class impl.
 *
 * Implementation of the MOSEARL class
 */

#include <stdlib.h>
#include <algorithm>
#include <iomanip>

#include "MOSEARL.h"

#include "garandom.h"
#include "GAEDAConfig.h"
#include "GAPopulation.h"
#include "Recombinator.h"
#include "MOSConversion.h"
#include "MOSTechnique.h"
#include "MOSTechniqueSet.h"
#include "genomes/MOSGenome.h"

/**
 * MOSEARL Constructors
 */
MOSEARL::MOSEARL(const GAGenome& genome, const double elitismPercent)
  : GAGeneticAlgorithm(genome),
    _auxPop(NULL),
    _techniqueSet(MOSTechniqueSet::handle()),
    _conversion(MOSConversion::handle()),
    _elitismPercent(elitismPercent),
    _parts(_techniqueSet->nTechniques(), 10/_techniqueSet->nTechniques()),
    _Inds(_techniqueSet->nTechniques(), 0),
    _partialInds(_techniqueSet->nTechniques(), 0),
    _avgFit(_techniqueSet->nTechniques(), 0.0),
    _totalInds(0),
    _state(-1),
    _stateGen(0),
    _policy(GAEDAConfig::handle()->getMOSRLPolicy()),
    ALFA (GAEDAConfig::handle()->getMOSRLAlfa()),
    BETA (GAEDAConfig::handle()->getMOSRLBeta()),
    GAMMA (GAEDAConfig::handle()->getMOSRLGamma()),
    PMAX (GAEDAConfig::handle()->getMOSQLPMax()),
    DELTA (GAEDAConfig::handle()->getMOSRLDelta()),
    DELTA_L (GAEDAConfig::handle()->getMOSRLDeltaL()),
    DELTA_W (GAEDAConfig::handle()->getMOSRLDeltaW())
{
  _auxPop = new GAPopulation(*pop);
}


MOSEARL::MOSEARL(const GAPopulation& population, const double elitismPercent)
  : GAGeneticAlgorithm(population),
    _auxPop(NULL),
    _techniqueSet(MOSTechniqueSet::handle()),
    _conversion(MOSConversion::handle()),
    _elitismPercent(elitismPercent),
    _parts(_techniqueSet->nTechniques(), 10/_techniqueSet->nTechniques()),
    _Inds(_techniqueSet->nTechniques(), 0),
    _partialInds(_techniqueSet->nTechniques(), 0),
    _avgFit(_techniqueSet->nTechniques(), 0.0),
    _totalInds(0),
    _state(-1),
    _stateGen(0),
    _policy(GAEDAConfig::handle()->getMOSRLPolicy()),
    ALFA (GAEDAConfig::handle()->getMOSRLAlfa()),
    BETA (GAEDAConfig::handle()->getMOSRLBeta()),
    GAMMA (GAEDAConfig::handle()->getMOSRLGamma()),
    PMAX (GAEDAConfig::handle()->getMOSQLPMax()),
    DELTA (GAEDAConfig::handle()->getMOSRLDelta()),
    DELTA_L (GAEDAConfig::handle()->getMOSRLDeltaL()),
    DELTA_W (GAEDAConfig::handle()->getMOSRLDeltaW())
{
  _auxPop = new GAPopulation(*pop);
}


MOSEARL::MOSEARL(const GAGeneticAlgorithm& alg, const double elitismPercent)
  : GAGeneticAlgorithm(alg),
    _auxPop(NULL),
    _techniqueSet(MOSTechniqueSet::handle()),
    _conversion(MOSConversion::handle()),
    _elitismPercent(elitismPercent),
    _parts(_techniqueSet->nTechniques(), 10/_techniqueSet->nTechniques()),
    _Inds(_techniqueSet->nTechniques(), 0),
    _partialInds(_techniqueSet->nTechniques(), 0),
    _avgFit(_techniqueSet->nTechniques(), 0.0),
    _totalInds(0),
    _state(-1),
    _stateGen(0),
    _policy(GAEDAConfig::handle()->getMOSRLPolicy()),
    ALFA (GAEDAConfig::handle()->getMOSRLAlfa()),
    BETA (GAEDAConfig::handle()->getMOSRLBeta()),
    GAMMA (GAEDAConfig::handle()->getMOSRLGamma()),
    PMAX (GAEDAConfig::handle()->getMOSQLPMax()),
    DELTA (GAEDAConfig::handle()->getMOSRLDelta()),
    DELTA_L (GAEDAConfig::handle()->getMOSRLDeltaL()),
    DELTA_W (GAEDAConfig::handle()->getMOSRLDeltaW())
{
  _auxPop = new GAPopulation(*pop);
}


/**
 * MOSEARL Destructor
 */
MOSEARL::~MOSEARL() {

  printQLMatrix();
  writeRLFile();

  if (_auxPop) {
    delete _auxPop;
  }

}


/**
 * Initializes the algorithm
 * @param seed Random seed
 */
void MOSEARL::initialize() {
  // Loads previous QL matrix and action probabilities matrix
  bool previousRL = readRLFile();

  // TODO: Not really necessary
  _techniqueSet->initPartRatios();

  // Initialize and evaluate initial population
  pop->initialize();
  pop->evaluate(gaTrue);
  pop->scale();

  // Reset statistics
  stats.scoreFrequency(gaDefScoreFrequency2);
  stats.reset(*pop);

  unsigned nTechs    = _techniqueSet->nTechniques();
  unsigned extraPart = 10 % nTechs;

  // Add extra participation to techniques in order (techniques are
  // initialized in the constructor with equal participation)
  for (unsigned i = 0; i < extraPart; i++)
    _parts[i]++;

  // Shake vector to assign extra participation to random techniques
  for (unsigned i = 0; i < nTechs; i++)
    std::swap(_parts[GARandomInt(0, nTechs-1)], _parts[GARandomInt(0, nTechs-1)]);

  // Register initial state
  _state = registerState (_parts);

  // Retrieve pointers to techniques for easy access
  unsigned tech = 0;
  for (MOSTechniqueSet::MOSTechniqueSetIterator it = _techniqueSet->begin();
       it != _techniqueSet->end(); it++, tech++)
    _techs [tech] = it->second;

  printStats("Initial Stats");

  switch (_policy) {
  case GAEDAConfig::PMAX:
    std::cout << "Using policy PMax." << std::endl;
    break;
  case GAEDAConfig::PHC:
    std::cout << "Using policy PHC." << std::endl;
    break;
  case GAEDAConfig::WOLF:
    std::cout << "Using policy WoLF." << std::endl;
    break;
  default:
    std::cout << "Error: Wrong policy. This should not happen..." << std::endl;
    exit(-1);
    break;
  }

  return;

}


/**
 * Evolve one generation
 */
void MOSEARL::step() {

  // Update generation state if necessary (every 100 generations)
  if (stats.generation() != 0 && stats.generation() % 100 == 0) {
    _stateGen++;
    _state = registerState(_parts);
  }

  // TODO: Fix techniques to allow the generation of a single offspring individual
  for (unsigned i = 0; i < pop->size(); i+=2) {

    // Select the technique used to produce new offspring. Retrieve both
    // pointer to the technique and action ID associated
    unsigned action;

    // Select a technique. Different methods are used
    // depending on the selected policy
    MOSTechnique* tech = NULL;

    switch (_policy) {
    case GAEDAConfig::PMAX:
      tech = selectTechniquePMax(action);
      break;
    case GAEDAConfig::PHC:
    case GAEDAConfig::WOLF:
      tech = selectTechniquePHCWoLF(action);
      break;
    default:
      std::cout << "Error: Wrong policy. This should not happen..." << std::endl;
      exit(-1);
      break;
    }

    // Select parents from current population
    GAGenome* newGenome;
    MOSGenome& mom = static_cast<MOSGenome&> (pop->select());
    MOSGenome& dad = static_cast<MOSGenome&> (pop->select());

    int encoding = tech->getEncoding();

    // Convert parents to the encoding of the selected technique (if necessary)
    if (!dad.existEncoding(encoding)) {
      newGenome = tech->getGenome();
      dad.addEncoding(encoding, newGenome);
      _conversion->convertGenome(dad.getDefaultEncoding(), encoding,
                                 dad.getDefaultGenome  (), newGenome);
    }

    if (!mom.existEncoding(encoding)) {
      newGenome = tech->getGenome();
      mom.addEncoding(encoding, newGenome);
      _conversion->convertGenome(mom.getDefaultEncoding(), encoding,
                                 mom.getDefaultGenome  (), newGenome);
    }

    // Offspring with selected technique
    tech->offspring(dad, mom, _auxPop, _auxPop->size(), i);

    // TODO: make techniques return the number of individuals generated to avoid
    // checking it here all the time
    unsigned inc = 0;

    if (i == pop->size() - 1)
      inc = 1;
    else
      inc = 2;

    _Inds[action]+=inc;
    _partialInds[action]+=inc;
    _totalInds+=inc;

    // Check for state transition
    double fracPart, intPart;
    fracPart = modf (((double)_Inds[action]/(double)_totalInds) * 10.0, &intPart);

    unsigned discretePart = (int) intPart;

    // The new state is, initially, the same as the current one
    int newState = _state;

    // If there's been a state transition, we should recompute the discrete
    // participations for all the techniques and change to the new state
    if (discretePart > _parts[action] && stats.generation() > 25) {

      std::vector< std::pair<unsigned, double> > decimals (_parts.size());
      unsigned totalDisc = 0;

      for (unsigned j = 0; j < _parts.size(); j++) {

         fracPart = modf (((double)_Inds[j]/(double)_totalInds) * 10.0, &intPart);
         _parts[j] = (int) intPart;

         decimals[j].first = j;
         decimals[j].second = fracPart;

         totalDisc += (int) intPart;

      }

      std::sort (decimals.begin(), decimals.end(), pairCmp);

      for (unsigned j = 0; j < 10 - totalDisc; j++)
         _parts[decimals[decimals.size()-j-1].first]++;

      newState = registerState (_parts);

    }

    // Define vars to update Q(s,a)
    double popAvgFit = stats.current(GAStatistics::Mean);
    double popMaxFit = stats.current(GAStatistics::Maximum);
    double fitThreshold = (popAvgFit+popMaxFit)/2;

    double fit, r, oldQ, maxQ;

    // Update the value of Q(s,a) twice, as we have generated two children
    for (unsigned j = 0; j < 2; j++) {

      // maxQ is the maximum value for any action in the current state (if there's not
      // been a transition) or of the new state (in case a transition happened)
      maxQ = *max_element (_QLearner[newState].begin(), _QLearner[newState].end());

      fit = _auxPop->individual(i+j).score();

      // Only update quality if the individual is good enough
      //if (fit > fitThreshold)
         r = (scaleFitness(fit) - scaleFitness(popAvgFit)) * (GAGenome::compareScores(popAvgFit, fit) < GAGenome::BETTER ? BETA : (1 - BETA));
      //else
         //r = 0;

      // Old Q(s,a) for current state and action
      oldQ = _QLearner[_state][action];

      // Update Q(s,a) for current state and action
      _QLearner[_state][action] = oldQ + ALFA * (r + GAMMA * maxQ - oldQ);

      // Update counter
      _QLearnerCounter[_state][action]++;

    }

    // Update action probabilities matrix. Different methods are use
    // depending on the selected policy.
    switch (_policy) {
    case GAEDAConfig::PMAX:
      // Nothing to do here...
      break;
    case GAEDAConfig::PHC:
      updateActionProbabilityPHC(action);
      break;
    case GAEDAConfig::WOLF:
      updateAveragePolicyEstimateWoLF();
      updateActionProbabilityWoLF(action);
      break;
    default:
      std::cout << "Error: Wrong policy. This should not happen..." << std::endl;
      exit(-1);
      break;
    }

    // Update current state (if a transition happened)
    _state = newState;

  }

  // Evaluate fitness of new population
  _auxPop->evaluate();

  // Scale fitness of new population
  _auxPop->scale();

  // Store average fitness of each technique (for quality statistics)
  for (unsigned i = 0; i < _auxPop->size(); i++) {
    MOSGenome& gen = dynamic_cast<MOSGenome&>(_auxPop->individual(i));
    unsigned id = gen.getTechniqueId();
    _avgFit[id] += (gen.score() / (double) _partialInds[id]);
  }

  // Recombine old and new populations
  recombinator_->recombine(*pop, *_auxPop);
  std::swap(pop, _auxPop);

  // TODO: Is this necessary?
  pop->evaluate();

  // Scale fitness values. TODO: Is this necessary
  pop->scale();

  updateQualityAndParticipation();

  // Combine statistics from different techniques
  stats.numsel  = _techniqueSet->getAllSelections  ();
  stats.numcro  = _techniqueSet->getAllCrossovers  ();
  stats.nummut  = _techniqueSet->getAllMutations   ();
  stats.numrep  = _techniqueSet->getAllReplacements();
  stats.numeval = _techniqueSet->getAllEvals       ();

  // Update statistics of the algorithm
  stats.update(*pop);

  printStats("End of Step Stats");

  return;

}


/**
 * Method to adjust the size of the populations
 */
unsigned MOSEARL::populationSize(unsigned size) {
  GAGeneticAlgorithm::populationSize(size);
  _auxPop->size(size);
  return size;
}


/**
 * Registers a new state
 */
unsigned MOSEARL::registerState(const std::vector<unsigned>& parts) {

  unsigned nTechs = _techniqueSet->nTechniques();

  // Convert the vector into a string to use as an identifier
  std::string state;

  for (unsigned i = 0; i < nTechs; i++) {
    if (parts[i] == 10)
      state+='A';
    else {
      char buff [10];
      sprintf (buff, "%d", parts[i]);
      state += buff;
    }
  }

  std::cout << "New state: " << state << std::endl;

  MapKey key = {state, _stateGen};

  if (_QLStates.find(key) != _QLStates.end())
    return _QLStates[key];
  else {

    std::vector<double> v (nTechs, 0.0);
    _QLearner.push_back (v);
    _averagePolicyEstimate.push_back (v);

    std::vector<unsigned> vc (nTechs, 0);
    _QLearnerCounter.push_back (vc);

    std::vector<double> va (nTechs, 1.0 / (double) nTechs);
    _actionProb.push_back (va);

    _counter.push_back (0);

    QLStatesT::const_iterator it;

    _QLStates.insert (pair<MapKey, unsigned> (key, _QLearner.size() - 1));
    return _QLearner.size()-1;

  }

}


/**
 * Selects the technique for the next offspring based on the QL matrix
 */
MOSTechnique* MOSEARL::selectTechniquePMax(unsigned& action) const {

  // Random number to see if we must use the technique with best Q(s,a)
  // or a random one among the rest
  double prob = GARandomDouble (0.0, 1.0);

  unsigned maxPos = bestTechnique(_state);

  // Select one of the other techniques proportionally to their Q(s,a)
  if (prob > PMAX) {

    std::vector<double> probs = _QLearner[_state];
    std::vector<unsigned> pos (probs.size(), 0);

    for (unsigned i = 0; i < pos.size(); i++)
      pos[i] = i;

    for (unsigned i = 0; i < pos.size(); i++)
      for (unsigned j = i; j < pos.size(); j++)
        if (probs[pos[i]] > probs[pos[j]])
          std::swap (pos[i], pos[j]);

    double suma = (double)(probs.size() * (probs.size() + 1))/2.0;

    for (unsigned i = 0; i< probs.size();) {
      unsigned j=i+1;
      while (j < probs.size() && probs[pos[j]]==probs[pos[i]]) {
        j++;
      }
        unsigned base=i+1;
        unsigned end=j+1;
        for(;i<j;i++)
          probs[pos[i]] = (((double)(end - 1 + base)/2.0))/suma;

    }

    double total = 1.0;

    // Erase best technique from this aux var (to avoid selecting it)
    //probs.erase(std::find(probs.begin(), probs.end(), maxVal));

    // Adjust probs vector with the proportional value for each technique
    for (unsigned i = 0; i < probs.size(); i++)
      probs[i] = (probs[i] / total) + ((i == 0) ? 0 : probs [i-1]);

    // Construct a vector of pairs (technique, Q(s,a)) to be able to
    // select a technique
    std::vector< std::pair<unsigned, double> > probsPairs;
    for (unsigned i = 0; i < probs.size(); i++)
      probsPairs.push_back(std::pair<unsigned, double> (i, probs[i]));

    // Sort the vector
    std::sort(probsPairs.begin(), probsPairs.end(), pairCmp);

    // New random number to select a technique
    double rnd = GARandomDouble(0.0, 1.0);

    for (unsigned i = 0; i < probsPairs.size(); i++)
      if (rnd <= probsPairs[i].second) {
        action = probsPairs[i].first;
        return _techs.find(action)->second;
      }

  }

  // Select best technique
  action = maxPos;
  return _techs.find(maxPos)->second;

}


/**
 * Selects the technique for the next offspring based on the Action
 * Probability matrix
 */
MOSTechnique* MOSEARL::selectTechniquePHCWoLF(unsigned& action) const {

  // Random number to see if we must use the technique with best Q(s,a)
  // or a random one among the rest
  double prob = GARandomDouble (0.0, 1.0);

  unsigned nTechs = _techs.size();
  double acum = 0.0;

  for (unsigned i = 0; i < nTechs; i++) {
    acum += _actionProb[_state][i];
    if (prob <= acum) {
      action = i;
      return _techs.find(action)->second;
    }
  }

  std::cerr << "Error: Exiting from selectTechniquePHC without selecting a technique." << std::endl;
  action = 0;
  return _techs.find(action)->second;

}


/**
 * Selects the best action (technique) in the given state
 */
unsigned MOSEARL::bestTechnique(unsigned state) const {

  double      max = -999999999.0;
  unsigned posMax = 0;

  for (unsigned i = 0; i < _QLearner[state].size(); i++) {

    if (_QLearner[state][i] > max) {
      max = _QLearner[state][i];
      posMax = i;
    }
    else if (_QLearner[state][i] == max && GAFlipCoin(0.5)) {
      max = _QLearner[state][i];
      posMax = i;
    }

  }

  return posMax;

}


/**
 * Selects the worst action (technique) in the given state
 */
unsigned MOSEARL::worstTechnique(unsigned state) const {

  double      min = 999999999.0;
  unsigned posMin = 0;

  for (unsigned i = 0; i < _QLearner[state].size(); i++) {

    if (_QLearner[state][i] < min) {
      min = _QLearner[state][i];
      posMin = i;
    }
    else if (_QLearner[state][i] == min && GAFlipCoin(0.5)) {
      min = _QLearner[state][i];
      posMin = i;
    }

  }

  return posMin;

}


/*
 * Updates quality and participation
 */
void MOSEARL::updateQualityAndParticipation() {

   std::map<unsigned, MOSTechnique*>::const_iterator it;

   for (it = _techs.begin(); it != _techs.end(); it++) {
      it->second->setQuality(_avgFit[it->first]);
      it->second->setPartRatio((double) _partialInds[it->first] / (double) pop->size());
      _partialInds[it->first] = 0;
      _avgFit[it->first] = 0;
   }

}


/**
 * Scales given fitness within maximum and minimum fitness
 */
double MOSEARL::scaleFitness(double fit) const {

   double popMaxFit = stats.current(GAStatistics::Maximum);
   double popMinFit = stats.current(GAStatistics::Minimum);

   if (popMaxFit == popMinFit) {
      if (popMaxFit == 0.0)
         return fit;
      else
         return fit / popMaxFit;
   }

   return (fit - popMinFit) / (popMaxFit - popMinFit);

}


/**
 * Updates Action Probability for the PHC policy
 */
void MOSEARL::updateActionProbabilityPHC(unsigned action) {

  // Index of the action with highest Q-value
  unsigned maxIndex = bestTechnique(_state);

  // Index of the action with lowest Q-value
  unsigned minIndex = worstTechnique(_state);

  unsigned nTechs = _techs.size();

  if (action == maxIndex)
    _actionProb[_state][action] += DELTA;
  else {
    _actionProb[_state][action] -= DELTA / (nTechs - 1);
    if (_actionProb[_state][action] < 0.0)
      _actionProb[_state][action] = 0.0;
  }

  double total = 0.0;
  for (unsigned i = 0; i < nTechs; i++)
    total += _actionProb[_state][i];

  for (unsigned i = 0; i < nTechs; i++)
    _actionProb[_state][i] /= total;

}


/**
 * Updates Action Probability for the WoLF policy
 */
void MOSEARL::updateActionProbabilityWoLF(unsigned action) {

  // Index of the action with highest Q-value
  unsigned maxIndex = bestTechnique(_state);

  // Index of the action with lowest Q-value
  unsigned minIndex = worstTechnique(_state);

  double currentSum = 0.0, averageSum = 0.0;

  for (unsigned i = 0; i < _QLearner[_state].size(); i++)
    currentSum += _actionProb[_state][i] * _QLearner[_state][i];

  for (unsigned i = 0; i < _QLearner[_state].size(); i++)
    averageSum += _averagePolicyEstimate[_state][i] * _QLearner[_state][i];

  double delta = (currentSum > averageSum) ? DELTA_W : DELTA_L;

  unsigned nTechs = _techs.size();

  if (action == maxIndex)
    _actionProb[_state][action] += delta;
  else {
    _actionProb[_state][action] -= delta / (nTechs - 1);
    if (_actionProb[_state][action] < 0.0)
      _actionProb[_state][action] = 0.0;
  }

  double total = 0.0;
  for (unsigned i = 0; i < nTechs; i++)
    total += _actionProb[_state][i];

  for (unsigned i = 0; i < nTechs; i++)
    _actionProb[_state][i] /= total;

}


/**
 * Updates Average Policy Estimate for the WoLF policy
 */
void MOSEARL::updateAveragePolicyEstimateWoLF() {

  _counter[_state]++;

  assert(_counter[_state] > 0);

  for (unsigned i = 0; i < _averagePolicyEstimate[_state].size(); i++) {
    _averagePolicyEstimate[_state][i] += (1.0 / (double) _counter[_state]) *
              (_actionProb[_state][i] - _averagePolicyEstimate[_state][i]);
  }

}


/**
 * Prints the QL matrix
 */
void MOSEARL::printQLMatrix() const {

  QLStatesT::const_iterator it;

  for (it = _QLStates.begin(); it != _QLStates.end(); it++) {

    std::cout << "State: ";

    for (unsigned i = 0; i < it->first.state.size(); i++)
      std::cout << it->first.state[i] << " ";

    std::cout << ", Gen: " << it->first.gen << " ==> ";

    for (unsigned i = 0; i < _QLearner[it->second].size(); i++)
      std::cout << _QLearner[it->second][i] << "\t";

    std::cout << std::endl;

  }

  std::cout << std::endl;

  for (it = _QLStates.begin(); it != _QLStates.end(); it++)
    for (unsigned i = 0; i < _QLearnerCounter[it->second].size(); i++) {

      std::cout << "State: ";

      for (unsigned j = 0; j < it->first.state.size(); j++)
        std::cout << it->first.state[j] << " ";

      std::cout << ", Gen: " << it->first.gen << " ";

      std::cout << ", Action: " << i << " => " << _QLearnerCounter[it->second][i] << std::endl;

    }

}


/**
 * Reads the QL matrix from a file
 */
bool MOSEARL::readRLFile() {

  std::ifstream f ("rl.dat", std::ios::binary);

  if (f.good()) {

    unsigned nStates, size;

    f.read ((char*) &nStates, sizeof(unsigned));

    for (unsigned i = 0; i < nStates; i++) {

      f.read ((char*) &size, sizeof(unsigned));
      char buff [size+1];

      f.read (buff, size*sizeof(char));
      buff[size]='\0';
      std::string state(buff);

      unsigned gen, id;
      f.read ((char*) &gen, sizeof(unsigned));
      f.read ((char*) &id, sizeof(unsigned));

      MapKey key = {state, gen};
      _QLStates.insert (pair<MapKey, unsigned> (key, id));

    }

    unsigned maxi, maxj;

    // Read the QL matrix
    f.read ((char*) &maxi, sizeof(unsigned));
    f.read ((char*) &maxj, sizeof(unsigned));

    std::vector<double> v (maxj, 0.0);
    _QLearner.insert (_QLearner.begin(), maxi, v);

    std::vector<unsigned> vc (maxj, 0);
    _QLearnerCounter.insert (_QLearnerCounter.begin(), maxi, vc);

    for (unsigned i = 0; i < maxi; i++)
      for (unsigned j = 0; j < maxj; j++) {
         double val;
         f.read ((char*) &val, sizeof(double));
         _QLearner[i][j] = val;
      }

    // Read the Action Probability matrix
    f.read ((char*) &maxi, sizeof(unsigned));
    f.read ((char*) &maxj, sizeof(unsigned));

    _actionProb.insert (_actionProb.begin(), maxi, v);

    for (unsigned i = 0; i < maxi; i++)
      for (unsigned j = 0; j < maxj; j++) {
         double val;
         f.read ((char*) &val, sizeof(double));
         _actionProb[i][j] = val;
      }

    // Read the Average Policy Estimate matrix
    f.read ((char*) &maxi, sizeof(unsigned));
    f.read ((char*) &maxj, sizeof(unsigned));

    _averagePolicyEstimate.insert (_averagePolicyEstimate.begin(), maxi, v);

    for (unsigned i = 0; i < maxi; i++)
      for (unsigned j = 0; j < maxj; j++) {
         double val;
         f.read ((char*) &val, sizeof(double));
         _averagePolicyEstimate[i][j] = val;
      }

    // Read the Counter vector
    f.read ((char*) &maxi, sizeof(unsigned));

    _counter.insert (_counter.begin(), maxi, 0);

    for (unsigned i = 0; i < maxi; i++) {
      unsigned val;
      f.read ((char*) &val, sizeof(unsigned));
      _counter[i] = val;
    }

    //    printQLMatrix();

    return true;

  }

  return false;

}


/**
 * Writes the QL matrix to a file
 */
bool MOSEARL::writeRLFile() const {

  std::ofstream f ("rl.dat", std::ios::binary);

  if (f.good()) {

    QLStatesT::const_iterator it;

    unsigned size = _QLStates.size();

    f.write ((char*) &size, sizeof(unsigned));

    for (it = _QLStates.begin(); it != _QLStates.end(); it++) {
      size = it->first.state.size();
      f.write ((char*) &size, sizeof(unsigned));
      f.write ((char*) (it->first.state.c_str()), size*sizeof(char));
      f.write ((char*) &(it->first.gen), sizeof(unsigned));
      f.write ((char*) &(it->second), sizeof(unsigned));
    }

    // Write the Q-Learner matrix
    size = _QLearner.size();
    f.write ((char*) &size, sizeof(unsigned));

    size = _QLearner[0].size();
    f.write ((char*) &size, sizeof(unsigned));

    for (unsigned i = 0; i < _QLearner.size(); i++)
      for (unsigned j = 0; j < _QLearner[i].size(); j++) {
        double val = _QLearner[i][j];
        f.write ((char*) &val, sizeof(double));
      }

    // Write the Action Probability matrix
    size = _actionProb.size();
    f.write ((char*) &size, sizeof(unsigned));

    size = _actionProb[0].size();
    f.write ((char*) &size, sizeof(unsigned));

    for (unsigned i = 0; i < _actionProb.size(); i++)
      for (unsigned j = 0; j < _actionProb[i].size(); j++) {
        double val = _actionProb[i][j];
        f.write ((char*) &val, sizeof(double));
      }

    // Write the Average Policy Estimate matrix
    size = _averagePolicyEstimate.size();
    f.write ((char*) &size, sizeof(unsigned));

    size = _averagePolicyEstimate[0].size();
    f.write ((char*) &size, sizeof(unsigned));

    for (unsigned i = 0; i < _averagePolicyEstimate.size(); i++)
      for (unsigned j = 0; j < _averagePolicyEstimate[i].size(); j++) {
        double val = _averagePolicyEstimate[i][j];
        f.write ((char*) &val, sizeof(double));
      }

    // Write the Counter vector
    size = _counter.size();
    f.write ((char*) &size, sizeof(unsigned));

    for (unsigned i = 0; i < _counter.size(); i++) {
        unsigned val = _counter[i];
        f.write ((char*) &val, sizeof(unsigned));
      }

    return true;

  }

  return false;

}

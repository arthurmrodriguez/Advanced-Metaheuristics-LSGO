#include "MOSParticipationFunc.h"

#include "MOSEA.h"
#include "MOSTechnique.h"
#include "MOSTechniqueSet.h"
#include "logger/GALogger.h"

#include <algorithm>
#include <iomanip>
#include <sstream>

/**
 * Constant participation function
 */
extern "C" void constantPF (MOSEA& algorithm) {
  return;
}


/**
 * Alternating participation function which distributes a given number of
 * iterations among the available techniques.
 */
extern "C" void alternatingQualityPF (MOSEA& alg) {

  // Base fitness for the Average Fitness QM
  static double base = 0.0;

  // Iteration of the current technique
  static unsigned gen = 0;

  // Current technique
  static unsigned currentTech = 0;

  // Reference to the Technique Set
  static MOSTechniqueSet& techSet = alg.getTechSet();

  // Number of techniques in the Techniques Set
  static unsigned nTechs = techSet.nTechniques();

  // Cache for storing the number of generations allowed for each technique
  static unsigned gensPerTech [10];

  // Caches for storing continuous values for the number of generations per technique
  static double gensPerTechCont [10];

  // Cache for storing the cumulated quality of each technique
  static double qualPerTech [10];

  // Total number of iterations to distribute among the techniques
  static unsigned sharedGens = 50;

  // Bool array that says if minimum ratio holds for a particular technique or not
  static bool canDie [10];

  // Current generation of the main algorithm
  int curGen = alg.generation();

  // Calculate the base value if we are in the first generation.
  if (curGen == 0) {
    if (alg.getQualityMeasure() == "fAvg")
       base = alg.population().worst().score();
  }

  // Init part ratios and distribute initial iterations
  if (curGen == 0) {
    unsigned total = 0;
    unsigned gens = sharedGens / nTechs;

    // Distribute the iterations uniformly
    for (unsigned i = 0; i < nTechs; i++) {
      gensPerTechCont [i] = (double) gens / (double) sharedGens;
      gensPerTech [i] = gens;
      total += gens;
      if (i == 0)
        techSet.setPartRatio(i, 1.0);
      else
        techSet.setPartRatio(i, 0.0);

      if (techSet.getTechnique(i)->classID() == GAID::TechniqueSTSDE)
        canDie[i] = true;
      else
        canDie[i] = false;
    }

    // Distribute the remaining iterations in the case that the division is not exact
    while (total < sharedGens) {
      (*(std::min_element(gensPerTech, gensPerTech + nTechs)))++;
      total++;
    }

    // Start with first technique
    currentTech = 0;

    stringstream message;
    message << curGen << " alt_part=[";
    for (unsigned i = 0; i < nTechs; i++)
      message << gensPerTech[i] << ((i != nTechs - 1) ? "," : "");
    message << "]";
    GALogger::instance()->appendLogMessage("", message.str(), GALogger::only_stats);

    return;
  }

  // Iterate through the techniques to retrieve quality values of the last generation
  for (unsigned i = 0; i < nTechs; i++)
    qualPerTech [i] += techSet.getTechnique(i)->getQuality();

  // Increment the number of iterations of the current technique
  gen++;

  // If the maximum numbers of iterations for that technique has been reached
  if (gen == gensPerTech [currentTech]) {

    while (++currentTech < nTechs && gensPerTech[currentTech] <= 0) ;

    // If this is not the last technique
    if (currentTech < nTechs) {
      techSet.setPartRatio(currentTech - 1, 0.0);
      techSet.setPartRatio(currentTech    , 1.0);
      gen = 0;
      return;
    }

  }
  else
    return;

  // If we have finished with the last technique, compute average quality during
  // the last set of iterations for each technique and readjust participation ratios

  double avgQual [nTechs];

  // Compute the average quality obtained by each technique in its assigned iterations
  for (unsigned i = 0; i < nTechs; i++)
    avgQual [i] = (gensPerTech[i] > 0) ? qualPerTech [i] / (double) gensPerTech [i] : 0.0;

  // Retrieve the index of the best technique
  unsigned bestTechnique = std::max_element(avgQual,avgQual+nTechs) - avgQual;

  techSet.setPartRatio(0, 1.0);

  for (unsigned i = 1; i < nTechs; i++)
    techSet.setPartRatio(i, 0.0);

  gen = 0;
  currentTech = 0;

  if (avgQual [bestTechnique] == 0) {
    stringstream message;
    message << curGen << " alt_part=[";
    for (unsigned i = 0; i < nTechs; i++)
      message << gensPerTech[i] << ((i != nTechs - 1) ? "," : "");
    message << "]";
    GALogger::instance()->appendLogMessage("", message.str(), GALogger::only_stats);
    return;
  }

  // Get the adjusting value
  void* tmp = alg.getParticipationFunctionData();
  double adjust = (tmp == NULL) ? 0.05 : *((double*) tmp);

  if (adjust <= 0 || adjust > 1)
    adjust = 0.05;

  // Get the min participation value
  double minGens = 0.10;

  double bestQuality = avgQual [bestTechnique];

  for (unsigned i = 0; i < nTechs; i++) {

    double diff = (bestQuality == 0.0) ? 1.0 : (bestQuality - avgQual [i]) / (bestQuality - base);
    double gensInc = gensPerTechCont [i] * adjust * diff;

    if (((gensPerTechCont [i] - gensInc) > minGens) || canDie[i]) {
       gensPerTechCont [bestTechnique] += gensInc;
       gensPerTechCont [i]             -= gensInc;
    }
    else if (gensPerTechCont [i] > minGens) {
       gensPerTechCont [bestTechnique] += gensPerTechCont [i] - minGens;
       gensPerTechCont [i]              = minGens;
    }

    qualPerTech [i] = 0;
    gensPerTech [i] = gensPerTechCont [i] * sharedGens;

  }

  gensPerTech [bestTechnique] = gensPerTechCont [bestTechnique] * sharedGens;

  stringstream message;
  message << curGen << " alt_part=[";
  for (unsigned i = 0; i < nTechs; i++)
    message << gensPerTech[i] << ((i != nTechs - 1) ? "," : "");
  message << "]";
  GALogger::instance()->appendLogMessage("", message.str(), GALogger::only_stats);

  return;

}

/**
 * Dynamic participation function. To calculate averages, it uses the population
 * quality.
 */
extern "C" void dynQualityMOSPF (MOSEA& alg) {

  static double base = 0.0;

  // Calculate stats so we can decide
  int curGen = alg.generation();
  unsigned int restoreGen = alg.getParticipationRestoreGen();
  MOSTechniqueSet& techSet = alg.getTechSet();

  // Calculate the base value if we are in the first generation.
  if (curGen == 0) {
    if (alg.getQualityMeasure() == "fAvg")
       base = alg.population().worst().score();
  }

  // Init part ratios
  if (curGen == 0 || (restoreGen != 0 && (curGen % restoreGen) == 0)) {
    techSet.initPartRatios();
    return;
  }

  techIdType bestTechId = techSet.getBestTechniqueId();

  if (bestTechId == -1)
    return;

  // Get the adjusting value
  void* tmp = alg.getParticipationFunctionData();
  double adjust;

  if (tmp == NULL)
    adjust = 0.05;
  else
    adjust = *((double*) tmp);

  if (adjust <= 0 or adjust > 1)
    adjust = 0.05;

  // Get the min participation value
  double minPart = alg.getMinParticipation();

  // Iterate through the algorithms and recalculate percents
  MOSTechniqueSet::MOSTechniqueSetIterator it;

  for (it=techSet.begin(); it!=techSet.end(); it++) {
    if (it->first != bestTechId) {
      double bestQuality = techSet.getTechnique(bestTechId)->getQuality();
      double diff;
      if (bestQuality == 0.0)
         diff = 1.0;
      else
         diff = (bestQuality - it->second->getQuality()) / (bestQuality - base);

      double currentTechPart = techSet.getPartRatio(it->first);
      double ratioInc = currentTechPart * adjust * diff;

      if ((currentTechPart - ratioInc) > minPart) {
         techSet.setPartRatio(bestTechId, techSet.getPartRatio(bestTechId) + ratioInc);
         techSet.setPartRatio(it->first, techSet.getPartRatio(it->first) - ratioInc);
      }
      else if (currentTechPart > minPart) {
         techSet.setPartRatio(bestTechId, techSet.getPartRatio(bestTechId) + (techSet.getPartRatio(it->first) - minPart));
         techSet.setPartRatio(it->first, minPart);
      }

      currentTechPart = techSet.getPartRatio(it->first);

      if (currentTechPart < 0.0)
         techSet.setPartRatio(it->first, 0.0);

      if (techSet.getPartRatio(bestTechId) > 1.0)
         techSet.setPartRatio(bestTechId, 1.0);
    }
  }
}

/**
 * Dynamic participation function with maximum participation and decreasing
 * population size. To calculate averages, it uses the population quality.
 */
extern "C" void dynQualityMaxPartPF (MOSEA& alg) {

  static unsigned origPopSize = alg.population().size();

  // Calculate stats so we can decide
  int curGen = alg.generation();
  MOSTechniqueSet& techSet = alg.getTechSet();

  static double maxPart = 1.0 / (double) techSet.nTechniques();

  // Init part ratios
  if (curGen == 0) {
    techSet.initPartRatios();
    return;
  }

  techIdType bestTechId    = techSet.getBestTechniqueId();
  double     bestQuality   = techSet.getTechnique(bestTechId)->getQuality();
  double     bestPartRatio = techSet.getPartRatio(bestTechId);

  // Get the min participation value
  double minPart = alg.getMinParticipation();

  // Iterate through the algorithms and recalculate percents
  MOSTechniqueSet::MOSTechniqueSetIterator it;

  double sumRatios = 0.0;
  double diff = 0.0;

  for (it=techSet.begin(); it!=techSet.end(); it++) {

    double alfa = it->second->getQuality() / bestQuality;

    double currentTechPart = techSet.getPartRatio(it->first);
    double newRatio = alfa * currentTechPart;

    if (newRatio < minPart)
      newRatio = minPart;

    sumRatios += newRatio;

    if (it->first != bestTechId)
      diff += (currentTechPart - newRatio);

    techSet.setPartRatio(it->first, newRatio);

  }

  if (bestPartRatio + diff <= maxPart) {
    techSet.setPartRatio(bestTechId, bestPartRatio + diff);
    sumRatios += diff;
  }
  else {
    techSet.setPartRatio(bestTechId, maxPart);
    sumRatios += (maxPart - bestQuality);
  }

  int newPopSize = 0;
  for (it=techSet.begin(); it!=techSet.end(); it++) {
    int size = (int) (techSet.getPartRatio(it->first) * origPopSize);
    newPopSize += size;
  }

  if (newPopSize > origPopSize)
    newPopSize = origPopSize;

  alg.offPopulationSize(newPopSize);

}


/**
 * Dynamic participation function with maximum participation and decreasing
 * population size. To calculate averages, it uses the population quality.
 */
extern "C" void dynQualityMOSPF2 (MOSEA& alg) {

  static unsigned origPopSize = alg.population().size();

  // Calculate stats so we can decide
  int curGen = alg.generation();
  MOSTechniqueSet& techSet = alg.getTechSet();

  // Init part ratios
  if (curGen == 0) {
    techSet.initPartRatios();
    return;
  }

  unsigned nTechs = techSet.nTechniques();

  std::vector<double> quals (nTechs, 0);
  std::vector<unsigned> pos (nTechs, 0);

  for (unsigned i = 0; i < nTechs; i++) {
    pos[i] = i;
    quals[i] = techSet.getTechnique(i)->getQuality();
  }

  for (unsigned i = 0; i < nTechs; i++)
    for (unsigned j = i; j < nTechs; j++)
      if (quals[pos[i]] > quals[pos[j]])
        std::swap (pos[i], pos[j]);

  // Get the min participation value
  double minPart = alg.getMinParticipation();

  unsigned nLose = nTechs / 2;

  double bestQuality = techSet.getTechnique(pos[nTechs-1])->getQuality();

  double diff = 0.0;

  for (unsigned i = 0; i < nLose; i++) {

    double alfa = techSet.getTechnique(pos[i])->getQuality() / bestQuality * 0.05;

    double currentTechPart = techSet.getPartRatio(pos[i]);
    double newRatio = currentTechPart - alfa;

    if (newRatio < minPart)
      newRatio = minPart;

    diff += (currentTechPart - newRatio);

    techSet.setPartRatio(pos[i], newRatio);

  }

  double sumQualities = 0.0;

  for (unsigned i = nLose; i < nTechs; i++)
    sumQualities += techSet.getTechnique(pos[i])->getQuality();

  for (unsigned i = nLose; i < nTechs; i++) {

    double currentTechPart = techSet.getPartRatio(pos[i]);
    double inc = techSet.getTechnique(pos[i])->getQuality() / sumQualities * diff;

    techSet.setPartRatio(pos[i], currentTechPart + inc);

  }

}

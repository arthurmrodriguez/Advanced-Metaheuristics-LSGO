#include "MOSParticipationFunction.h"

#include "MOSEA2.h"
#include "MOSTechnique.h"
#include "MOSTechniqueSet.h"
#include "MOSQualityFunction.h"
//#include "logger/GALogger.h"

#include <algorithm>
/*
#include <iomanip>
#include <sstream>
*/

/**
 * Constant participation function.
 */

ConstantParticipation::ConstantParticipation(long double minPart, unsigned restoreGen, void* pfData)
 : MOSParticipation(minPart, restoreGen, pfData) {

  MOSTechniqueSet::handle()->initPartRatios();

}

/**
 * Dynamic participation function.
 */

DynamicParticipation::DynamicParticipation(long double adjust, long double minPart, unsigned restoreGen, void* pfData)
 : MOSParticipation(minPart, restoreGen, pfData), _adjust(adjust) {

  MOSTechniqueSet::handle()->initPartRatios();

  if (_adjust < 0.0 || _adjust > 1.0)
    _adjust = 0.05;

}

DynamicParticipation::~DynamicParticipation() {
  if (_pfData)
    free(_pfData);
}

long double DynamicParticipation::setBaseQual (MOSEA2& alg) {
  _base = alg.getQualityFunction()->classID() == GAID::QualityFit ? alg.population().worst().score() : 0.0;
  return _base;
}

void DynamicParticipation::update (MOSEA2& alg) {

  int curGen = alg.generation();
  MOSTechniqueSet& techSet = *MOSTechniqueSet::handle();
  unsigned nTechs = techSet.nTechniques();

  // Reset participation ratios and return
  if (_restoreGen != 0 && (curGen % _restoreGen) == 0) {
    techSet.initPartRatios();
    return;
  }

  // Retrieve the set of the best techniques
  std::vector<techIdType> bestTechIds = techSet.getBestTechniqueIds(alg.generation());

  // If no best technique or all the techniques are best techniques,
  // quality is zero or equal for every technique => No adjustment needed
  if (bestTechIds.size() == 0 || bestTechIds.size() == nTechs)
    return;

  // Quality of the best techniques (it should be the same for all of them)
  long double bestQuality = techSet.getTechnique(bestTechIds[0])->getQuality();

  for (MOSTechniqueSet::MOSTechniqueSetIterator it =techSet.begin(); it != techSet.end(); it++) {
    if (std::find(bestTechIds.begin(), bestTechIds.end(), it->first) == bestTechIds.end()) {
      long double diff = (bestQuality - it->second->getQuality()) / (bestQuality - _base);

      long double currentTechPart = techSet.getPartRatio(it->first);
      long double ratioInc = currentTechPart * _adjust * diff;

      if ((currentTechPart - ratioInc) > _minPart) {
        long double sharedRatio = ratioInc / (long double) bestTechIds.size();
        for (unsigned i = 0; i < bestTechIds.size(); i++)
          techSet.setPartRatio(bestTechIds[i], techSet.getPartRatio(bestTechIds[i]) + sharedRatio);

        techSet.setPartRatio(it->first, techSet.getPartRatio(it->first) - ratioInc);
      }
      else if (currentTechPart > _minPart) {
        long double sharedRatio = (techSet.getPartRatio(it->first) - _minPart) / (long double) bestTechIds.size();
        for (unsigned i = 0; i < bestTechIds.size(); i++)
          techSet.setPartRatio(bestTechIds[i], techSet.getPartRatio(bestTechIds[i]) + sharedRatio);

        techSet.setPartRatio(it->first, _minPart);
      }

    }

  }

  // Correct possibly wrong values due to floating operations
  for (MOSTechniqueSet::MOSTechniqueSetIterator it =techSet.begin(); it != techSet.end(); it++) {
    long double currentTechPart = techSet.getPartRatio(it->first);
    if (currentTechPart < 0.0)
       techSet.setPartRatio(it->first, 0.0);
    else if (currentTechPart > 1.0)
       techSet.setPartRatio(it->first, 1.0);
  }

  MOSTechnique::improvement_override = false;

}

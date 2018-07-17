/**
 * @file
 * @brief MOSTechniqueSet class impl.
 *
 */

#include <algorithm>

#include "MOSTechniqueSet.h"

#include "GAEDAConfig.h"
#include "garandom.h"
#include "GAPopulation.h"
#include "MOSTechnique.h"
#include "MOSConfig.h"

MOSTechniqueSet* MOSTechniqueSet::pSelf = NULL;

/**
 * Manejador Singleton
 */
MOSTechniqueSet* MOSTechniqueSet::handle() {
  if( !pSelf ) {
    pSelf = new MOSTechniqueSet;
  }
  return pSelf;
}

void MOSTechniqueSet::destroy() {
  if(pSelf)
    delete pSelf;
  pSelf=NULL;
}

/**
 * Constructor
 */
MOSTechniqueSet::MOSTechniqueSet() {
  _bonus = GAEDAConfig::handle()->getBonus();
}

/**
 * Destructor
 */
MOSTechniqueSet::~MOSTechniqueSet () {
  TechniqueSet::iterator it;
  for (it = _techniqueSet.begin(); it != _techniqueSet.end(); it++)
    delete it->second;
  _techniqueSet.clear();
}

/**
 * Registrar una técnica
 *
 * @param technique Referencia a MOSTechnique
 * @return bool Cierto si técnica registrada correctamente
 */
bool MOSTechniqueSet::registerTechnique (MOSTechnique* technique) {
  techIdType id = technique->getId();
  TechniqueSet::iterator iterTech = _techniqueSet.find(id);
  if (iterTech == _techniqueSet.end()) {
    technique->setPartRatio(0.0);
    _techniqueSet[id] = technique;

    return true;
  }
  else {
    return false;
  }
}

/**
 * Dar de baja una técnica
 *
 * @param techniqueID id de la tecnica
 * @return bool Cierto si técnica dada de baja correctamente
 */
bool MOSTechniqueSet::unregisterTechnique (const techIdType techniqueId) {
  TechniqueSet::iterator iter = _techniqueSet.find(techniqueId);

   if (iter == _techniqueSet.end()) {
    return false;
   }
}

/**
 * Obtener referencia a técnica
 *
 */
MOSTechnique* MOSTechniqueSet::getTechnique(const techIdType techniqueId) {
  TechniqueSet::iterator iter = _techniqueSet.find(techniqueId);
  if (iter == _techniqueSet.end()) {
    return NULL;
  } else {
    return _techniqueSet[techniqueId];
  }
}

/**
 * Contar técnicas registradas
 *
 * @return unsigned Numero de técnicas registradas
 */
unsigned MOSTechniqueSet::nTechniques () const {
    return _techniqueSet.size();
}

/**
 * Fijar ratio de participación de una técnica
 *
 * @param techniqueID string que identifica la técnica
 * @param value ratio de participación
 * @return bool Cierto si ratio fijado correctamente
 */
bool MOSTechniqueSet::setPartRatio (const techIdType techniqueId, const long double partRatio) {
  TechniqueSet::iterator it = _techniqueSet.find(techniqueId);
  if (it == _techniqueSet.end()) {
    return false;
  }
  else {
    /*
    assert (partRatio <= 1.0);
    assert (partRatio >= 0.0);
    */
    // TODO: Corregir esto
    if (partRatio >= 1.0)
      _techniqueSet[techniqueId]->setPartRatio(1.0);
    else if (partRatio <= 0.0)
      _techniqueSet[techniqueId]->setPartRatio(0.0);

    _techniqueSet[techniqueId]->setPartRatio(partRatio);
    return true;
  }
}

/**
 * Consultar ratio de participación de una técnica
 *
 * @param techniqueID string que identifica la técnica
 * @return long double ratio de participación
 */
long double MOSTechniqueSet::getPartRatio (const techIdType techniqueId) const {
  TechniqueSet::const_iterator it = _techniqueSet.find(techniqueId);
  if (it == _techniqueSet.end()) {
    return 0.0;
  } else {
    return it->second->getPartRatio();
  }
}

/**
 * Recorrer las técnicas y acumular sus ratios de participación
 */
long double MOSTechniqueSet::sumPartRatios() {
  long double ratiosAcum = 0.0;
  TechniqueSet::iterator it;
  for (it=_techniqueSet.begin(); it!=_techniqueSet.end(); it++) {
    ratiosAcum += getPartRatio(it->first);
  }
  return ratiosAcum;
}

/**
 * Resetear las calidades acumuladas de las técnicas
 */
void MOSTechniqueSet::resetTechQuality() {
  TechniqueSet::iterator it;
  for (it=_techniqueSet.begin(); it!=_techniqueSet.end (); it++) {
    it->second->setQuality(0.0);
  }
  return;
}

/**
 * Obtienes la calidad de una tecnica
 * @param techniqueId Identificador de la técnica
 */
long double MOSTechniqueSet::getTechQuality(techIdType techniqueId) {
  MOSTechnique* technique = _techniqueSet[techniqueId];
  return technique->getQuality();
}

/**
 * Fijar la calidad de las técnicas
 * @param techniqueId Identificador de la técnica
 * @param quality Calidad de la técnica
 */
void MOSTechniqueSet::setTechQuality(techIdType techniqueId, long double quality) {

  MOSTechnique* technique = _techniqueSet[techniqueId];
  technique->setQuality(quality);
  return;
}

/**
 * Agregar estadísticas de seleccionados
 */
unsigned long int MOSTechniqueSet::getAllSelections() {
  unsigned long int numsel = 0;
  TechniqueSet::iterator it;
  for (it=_techniqueSet.begin(); it!=_techniqueSet.end(); it++) {
    numsel += it->second->getSelections();
  }
  return numsel;
}

/**
 * Agregar estadisticas de cruces
 */
unsigned long int MOSTechniqueSet::getAllCrossovers() {
  unsigned long int numcross = 0;
  TechniqueSet::iterator it;
  for (it=_techniqueSet.begin(); it!=_techniqueSet.end(); it++) {
    numcross += it->second->getCrossovers();
  }
  return numcross;
}

/**
 * Agregar estadisticas de mutaciones
 */
unsigned long int MOSTechniqueSet::getAllMutations() {
  unsigned long int nummut = 0;
  TechniqueSet::iterator it;
  for (it=_techniqueSet.begin(); it!=_techniqueSet.end(); it++) {
    nummut += it->second->getMutations();
  }
  return nummut;
}

/**
 * Agregar estadisticas de reemplazos
 */
unsigned long int MOSTechniqueSet::getAllReplacements() {
  unsigned long int numrep = 0;
  TechniqueSet::iterator it;
  for (it=_techniqueSet.begin(); it!=_techniqueSet.end(); it++) {
    numrep += it->second->getReplacements();
  }
  return numrep;
}

/**
 * Agregar estadisticas de evaluaciones
 */
unsigned long int MOSTechniqueSet::getAllEvals() {
  unsigned long int numeval = 0;
  TechniqueSet::iterator it;
  for (it=_techniqueSet.begin(); it!=_techniqueSet.end(); it++) {
    numeval += it->second->getEvals();
  }
  return numeval;
}

/**
 * Inicia los ratios de participacion
 */
void MOSTechniqueSet::initPartRatios() {
  TechniqueSet::iterator it;
  long double partRatio = 1.0 / _techniqueSet.size();

  for (it=_techniqueSet.begin(); it!=_techniqueSet.end(); it++) {
    it->second->setPartRatio(partRatio);
  }
  return;
}

/**
 * Recupera los ratios de participación en un vector
 * PRE: Presuponemos que el vector que nos pasan es del mismo tamano que el numero de tecnicas
 * @return std::vector<long double>* Vector con los ratios de participacion
 */
std::vector<long double>* MOSTechniqueSet::getPartRatiosVector(std::vector<long double>* ratiosVector) const {
  TechniqueSet::const_iterator it;
  unsigned i = 0;

  for (it=_techniqueSet.begin(); it!=_techniqueSet.end(); it++) {
    (*ratiosVector)[i]=(it->second->getPartRatio());
    i++;
  }
  return ratiosVector;
}

/**
 * Recupera la calidad en un vector
 * PRE: Presuponemos que el vector que nos pasan es del mismo tamano que el numero de tecnicas
 * @return std::vector<long double>* Vector con la calidad
 */
std::vector<long double>* MOSTechniqueSet::getQualityVector(std::vector<long double>* qualityVector) const {
  TechniqueSet::const_iterator it;
  unsigned i = 0;

  for (it=_techniqueSet.begin(); it!=_techniqueSet.end(); it++) {
    (*qualityVector)[i]=(it->second->getQuality(true));
    i++;
  }
  return qualityVector;
}

std::vector<techIdType> MOSTechniqueSet::getBestTechniqueIds(unsigned gen) const {

  TechniqueSet::const_iterator it;
  long double bestQuality = 0.0;
  long double bestImprovement = 0.0;
  std::vector<techIdType> bestTechs;
  std::vector<techIdType> mostImprovedTechs;

  for (it=_techniqueSet.begin(); it!=_techniqueSet.end(); it++) {
    if ((it->second->getQuality(true) == bestQuality) && (bestQuality != 0)) {
      bestTechs.push_back(it->first);
    }
    if (it->second->getQuality(true) > bestQuality) {
      bestTechs.clear();
      bestTechs.push_back(it->first);
      bestQuality = it->second->getQuality(true);
    }
  }

  for (it=_techniqueSet.begin(); it!=_techniqueSet.end(); it++) {
    if ((it->second->getRatioImproved() == bestImprovement) && (bestImprovement != 0)) {
      mostImprovedTechs.push_back(it->first);
    }
    if (it->second->getRatioImproved() > bestImprovement) {
      mostImprovedTechs.clear();
      mostImprovedTechs.push_back(it->first);
      bestImprovement = it->second->getRatioImproved();
    }
  }

  if (bestTechs.size() != mostImprovedTechs.size() && gen > 10) {
    MOSTechnique::improvement_override = true;
    return mostImprovedTechs;
  }
  else {

    bool equal = true;

    for (unsigned i = 0; i < bestTechs.size(); i++)
      equal = equal && bestTechs[i] == mostImprovedTechs[i];

    if (equal || gen <= 10) {
      MOSTechnique::improvement_override = false;
      return bestTechs;
    }
    else {
      MOSTechnique::improvement_override = true;
      return mostImprovedTechs;
    }

  }
}

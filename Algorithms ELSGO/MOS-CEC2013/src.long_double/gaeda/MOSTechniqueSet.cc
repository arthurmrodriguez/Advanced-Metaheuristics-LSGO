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
#include "MOSTechniqueES.h"
#include "MOSConfig.h"
#include "MOSParticipationFunc.h"

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
  _hasDE = 0;

}

/**
 * Destructor
 */
MOSTechniqueSet::~MOSTechniqueSet () {
  TechniqueSet::iterator it;
  for (it = _techniqueSet.begin(); it != _techniqueSet.end(); it++)
    delete it->second;
  _techniqueSet.clear();
  _esTechs.clear();
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

    if (technique->classID()==GAID::TechniqueES)
      _esTechs.push_back(technique);
    else if (technique->classID()==GAID::TechniqueDE)
      _hasDE++;

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

  if (iter->second->classID() == GAID::TechniqueDE)
    _hasDE--;

   if (iter == _techniqueSet.end()) {
    return false;
   }
   else {
    std::vector<MOSTechnique*>::iterator it = std::find(_esTechs.begin(), _esTechs.end(), iter->second);
    delete iter->second;
    _esTechs.erase(it);
    _techniqueSet.erase(iter);
    return true;
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
 * Recorrer todas las técnicas y evolucionar una generación en todas en función de los ratios
 * @param origPop Población origen
 * @param destPop Población destino
 */
void MOSTechniqueSet::offspringAll(GAPopulation* origPop, GAPopulation* destPop, unsigned maxEvals) {

  unsigned offset = 0;
  long double ratiosAcum = sumPartRatios();
  MOSTechnique* currentTechnique = NULL;

  TechniqueSet::iterator it;

  if (maxEvals < destPop->size())
    destPop->size(maxEvals);

  for (it=_techniqueSet.begin(); it!=_techniqueSet.end(); it++) {
    int size = (int) (getPartRatio(it->first) * destPop->size() / ratiosAcum);
    if (size > 0) {
      currentTechnique = it->second;
      currentTechnique->offspring (origPop, destPop, size, offset);
      offset += size;
    }
  }

  // Descendencia con la ultima tecnica de los individuos que queden descolgados
  if (destPop->size() > offset) {
    currentTechnique->offspring (origPop, destPop, destPop->size()-offset, offset);
  }

  return;

}


void MOSTechniqueSet::offspringAll(std::vector<GAPopulation*>& origPops, std::vector<GAPopulation*> destPops) {

  unsigned offset = 0;
  long double ratiosAcum = (GAEDAConfig::handle()->getPartFunction() != dynQualityMaxPartPF) ? sumPartRatios() : 1.0;
  MOSTechnique* currentTechnique = NULL;

  TechniqueSet::iterator it;
  unsigned popSize = GAEDAConfig::handle()->getPopSize();

  for (it=_techniqueSet.begin(); it!=_techniqueSet.end(); it++) {

    int size = (int) (getPartRatio(it->first) * popSize / ratiosAcum);

    if (it->first == (_techniqueSet.size() - 1))
      size += (popSize - offset - size);

    if (size > 0) {
      destPops[it->first]->size(size);
      currentTechnique = it->second;
      currentTechnique->offspring (origPops[it->first], destPops[it->first], size, 0);
      offset += size;
    }

  }

  return;
}


/**
 * Premiar una técnica a través de un vector (map) de probabilidades
 */
MOSProbVector& MOSTechniqueSet::bonusTechnique(MOSProbVector& probVector, techIdType techniqueId) {

  // Calcular bonus previsto para la tecnica ganadora
  long double incBonusTech = probVector[techniqueId] * _bonus;
  // Numero inicial de tecnicas a penalizar
  int decList = probVector.size() - 1;
  // Acumulador de penalización no aplicada (ya han llegado a 0)
  long double penNotAply = 0.0;
  // Calcular penalización inicial
  long double decTech = incBonusTech / decList;

  MOSProbVector::iterator it;
  for(it = probVector.begin(); it!=probVector.end(); it++) {
    // Ajuste de ganadora a 1.0
    if (probVector[techniqueId] + incBonusTech >= 1.0) {
      if (it->first == techniqueId) {
        // Tech ganadora
        probVector[it->first] = 1.0;
      } else {
        // Resto tech a 0
        probVector[it->first] = 0.0;
      }
    // Ajuste de ganadora
    } else {
      if ((it->second >= decTech) && (it->first != techniqueId)) {
        // Decrementar perdedoras
        probVector[it->first] -= decTech;
      } else {
        if (it->first != techniqueId) {
          // Acumular penalizacion no aplicada
          penNotAply += probVector[it->first];
          probVector[it->first] = 0.0;
        }
      }
    }
  }

  if (probVector[techniqueId] + incBonusTech < 1.0 ) {
    // Restar penalizacion no aplicada al bonus de la ganadora
    incBonusTech -= penNotAply;

    // Bonificar técnica ganadora
    probVector[techniqueId] += incBonusTech;
  }

  return probVector;
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
    (*qualityVector)[i]=(it->second->getQuality());
    i++;
  }
  return qualityVector;
}

/**
 * Devuelve el identificador de la tecnica con menos calidad
 */
techIdType MOSTechniqueSet::getBestTechniqueId() const {
  TechniqueSet::const_iterator it;
  long double bestQuality = 0.0;
  techIdType bestTech = -1;

  for (it=_techniqueSet.begin(); it!=_techniqueSet.end(); it++) {
    if (it->second->getQuality() > bestQuality) {
      bestTech = it->first;
      bestQuality = it->second->getQuality();
    }
  }
  return bestTech;
}

/**
 * Devuelve el identificador de la tecnica con menos calidad
 */
techIdType MOSTechniqueSet::getWorstTechniqueId() const {
  TechniqueSet::const_iterator it;
  long double worstQuality = 999999999999.0;
  techIdType worstTech = -1;

  for (it=_techniqueSet.begin(); it!=_techniqueSet.end(); it++) {
    if (it->second->getQuality() <= worstQuality) {
      worstTech = it->first;
      worstQuality = it->second->getQuality();
    }
  }
  return worstTech;
}

std::vector<techIdType> MOSTechniqueSet::getBestTechniqueIds(unsigned gen) const {

  TechniqueSet::const_iterator it;
  long double bestQuality = 0.0;
  long double bestImprovement = 0.0;
  std::vector<techIdType> bestTechs;
  std::vector<techIdType> mostImprovedTechs;

  for (it=_techniqueSet.begin(); it!=_techniqueSet.end(); it++) {
    if ((it->second->getQuality() == bestQuality)) {
      bestTechs.push_back(it->first);
    }
    if (it->second->getQuality() > bestQuality) {
      bestTechs.clear();
      bestTechs.push_back(it->first);
      bestQuality = it->second->getQuality();
    }
  }

  for (it=_techniqueSet.begin(); it!=_techniqueSet.end(); it++) {
    if ((it->second->getRatioImproved() == bestImprovement)) {
      mostImprovedTechs.push_back(it->first);
    }
    if (it->second->getRatioImproved() > bestImprovement) {
      mostImprovedTechs.clear();
      mostImprovedTechs.push_back(it->first);
      bestImprovement = it->second->getRatioImproved();
    }
  }

//  if (bestTechs.size() != mostImprovedTechs.size()) {
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

std::vector<techIdType> MOSTechniqueSet::getWorstTechniqueIds() const {
  TechniqueSet::const_iterator it;
  long double worstQuality = 999999999999.0;
  std::vector<techIdType> worstTechs;

  for (it=_techniqueSet.begin(); it!=_techniqueSet.end(); it++) {
    if ((it->second->getQuality() == worstQuality) && (worstQuality != 0)) {
      worstTechs.push_back(it->first);
    }
    if ((it->second->getQuality() < worstQuality)  && (worstQuality != 0)) {
      worstTechs.clear();
      worstTechs.push_back(it->first);
      worstQuality = it->second->getQuality();
    }
  }
  return worstTechs;
}

bool MOSTechniqueSet::updateStrategies(const std::vector<GAGenome*>& parents, GAGenome* child) const {

  if (_esTechs.size() > 0) {
    long double parts[_esTechs.size()];
    long double totalPart=0.0;

    for (unsigned i = 0; i < _esTechs.size(); i++) {
      parts[i]=_esTechs[i]->getPartRatio();
      totalPart+=parts[i];
    }

    for (unsigned i = 0; i < _esTechs.size(); i++)
      parts[i]=(parts[i] / totalPart) + ((i == 0) ? 0 : parts [i-1]);

    long double prob=GARandomDouble(0,1);

    for (unsigned j = 0; j < _esTechs.size(); j++)
      if (prob <= parts[j]) {

        MOSTechniqueES& es = DYN_CAST (MOSTechniqueES&, *_esTechs[j]);
        es.updateStrategies(parents, child);
        break;

      }
  }

  return true;

}

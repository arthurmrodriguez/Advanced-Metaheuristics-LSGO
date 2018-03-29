/**
 * @file
 * @brief MOSTechniqueSet class hdr.
 *
 */

#ifndef MOSTechniqueSet_H
#define MOSTechniqueSet_H

#include <map>
#include <vector>
#include <string>

#include "MOSConfig.h"

class GAGenome;
class GAPopulation;
class MOSTechnique;

/**
 * @brief Clase que representa el conjunto de técnicas
 */
class MOSTechniqueSet {

public:

  typedef std::map<techIdType, MOSTechnique*> TechniqueSet;
  typedef TechniqueSet::iterator MOSTechniqueSetIterator;
  typedef std::map<encodingType, GAGenome*> EncodingSet;

  // Registrar técnica en el conjunto
  bool registerTechnique(MOSTechnique* technique);

  // Eliminar técnica del conjunto
  bool unregisterTechnique(const techIdType techniqueId);

  // Obtener una técnica
  MOSTechnique* getTechnique(const techIdType techniqueId);

  // Número de técnicas registradas
  unsigned nTechniques() const;

  // Begin
  MOSTechniqueSetIterator begin() {return _techniqueSet.begin();}

  // End
  MOSTechniqueSetIterator end() {return _techniqueSet.end();}

  // Recorrer todas las técnicas y evolucionar una generación en todas en función de los ratios
  void offspringAll(GAPopulation* origPop, GAPopulation* destPop, unsigned maxEvals);
  void offspringAll(std::vector<GAPopulation*>& origPops, std::vector<GAPopulation*> destPops);

  // Fijar el ratio de participación de una técnica
  bool setPartRatio(const techIdType techniqueId, const double partRatio);

  // Obtener el ratio de participación de una técnica
  double getPartRatio(const techIdType techniqueId) const;

  // Iniciar ratios de participacion de las técnicas de forma equiprobable
  void initPartRatios();

  // Recupera los ratios de participación en un vector
  std::vector<double>* getPartRatiosVector(std::vector<double>* ratiosVector) const;

  // Poner a 0 la calidad de todas las técnicas
  void resetTechQuality();

  // Obtienes la calidad de la técnica
  double getTechQuality(const techIdType techniqueId);

  // Fijar la calidad de la técnica
  void setTechQuality(const techIdType techniqueId, const double quality);

  // Devolver Id de la tecnica con mejor calidad
  techIdType getBestTechniqueId() const;

  // Devolver Id de la tecnica con peor calidad
  techIdType getWorstTechniqueId() const;

  std::vector<techIdType> getBestTechniqueIds(unsigned gen) const;
  std::vector<techIdType> getWorstTechniqueIds() const;

  // Premiar una técnica a través de un vector (map) de probabilidades
  MOSProbVector& bonusTechnique(MOSProbVector& probVector, techIdType techniqueId);

  // Recupera la de las técnicas en un vector
  std::vector<double>* getQualityVector(std::vector<double>* qualityVector) const;

  // Obtener estadísticas
  unsigned long int getAllSelections();
  unsigned long int getAllCrossovers();
  unsigned long int getAllMutations();
  unsigned long int getAllReplacements();
  unsigned long int getAllEvals();

  bool updateStrategies(const std::vector<GAGenome*>& parents, GAGenome* child) const;

  // Método para obtener instancia única
  static MOSTechniqueSet* handle();

  static void destroy();

  bool hasDE() const {return _hasDE != 0;}

  // Sumatorio de los ratios de participación de las técnicas
  double sumPartRatios();

protected:

  // Constructor
  MOSTechniqueSet();
  virtual ~MOSTechniqueSet();

  TechniqueSet _techniqueSet; // Conjunto de técnicas (map)
  static MOSTechniqueSet* pSelf; // Puntero a la instancia si existe

  std::vector<MOSTechnique*> _esTechs; // ES techniques (for easy access)
  double _bonus;
  unsigned _hasDE;

};

#endif

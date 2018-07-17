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

  // Fijar el ratio de participación de una técnica
  bool setPartRatio(const techIdType techniqueId, const long double partRatio);

  // Obtener el ratio de participación de una técnica
  long double getPartRatio(const techIdType techniqueId) const;

  // Iniciar ratios de participacion de las técnicas de forma equiprobable
  void initPartRatios();

  // Recupera los ratios de participación en un vector
  std::vector<long double>* getPartRatiosVector(std::vector<long double>* ratiosVector) const;

  // Poner a 0 la calidad de todas las técnicas
  void resetTechQuality();

  // Obtienes la calidad de la técnica
  long double getTechQuality(const techIdType techniqueId);

  // Fijar la calidad de la técnica
  void setTechQuality(const techIdType techniqueId, const long double quality);

  std::vector<techIdType> getBestTechniqueIds(unsigned gen) const;

  // Recupera la de las técnicas en un vector
  std::vector<long double>* getQualityVector(std::vector<long double>* qualityVector) const;

  // Obtener estadísticas
  unsigned long int getAllSelections();
  unsigned long int getAllCrossovers();
  unsigned long int getAllMutations();
  unsigned long int getAllReplacements();
  unsigned long int getAllEvals();

  // Método para obtener instancia única
  static MOSTechniqueSet* handle();

  static void destroy();

  // Sumatorio de los ratios de participación de las técnicas
  long double sumPartRatios();

protected:

  // Constructor
  MOSTechniqueSet();
  virtual ~MOSTechniqueSet();

  TechniqueSet _techniqueSet; // Conjunto de técnicas (map)
  static MOSTechniqueSet* pSelf; // Puntero a la instancia si existe

  long double _bonus;

};

#endif

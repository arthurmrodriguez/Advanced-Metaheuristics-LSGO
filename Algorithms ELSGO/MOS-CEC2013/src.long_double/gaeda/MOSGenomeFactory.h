/**
 * @file
 * @brief MOSGenomeFactory class hdr.
 *
 * Header file for the MOS Genome Factory
 */

#ifndef MOSGenomeFactory_H
#define MOSGenomeFactory_H

#include "genomes/GAGenome.h"
#include "MOSConfig.h"

/**
 * Class: MOSGenomeFactory
 *
 * @brief Descripción de la clase MOSGenomeFactory
 */
class MOSGenomeFactory {

 public:

  GAGenome* getGenome(encodingType encoding);

  void registerGenome(encodingType encoding, GAGenome* genome);

  // Método para obtener instancia única
  static MOSGenomeFactory* handle();

  static void destroy();

 protected:

  std::map<encodingType, GAGenome*> genomeSet;
  static MOSGenomeFactory* pSelf; // Puntero a la instancia si existe
  MOSGenomeFactory();
  virtual ~MOSGenomeFactory();

};

#endif

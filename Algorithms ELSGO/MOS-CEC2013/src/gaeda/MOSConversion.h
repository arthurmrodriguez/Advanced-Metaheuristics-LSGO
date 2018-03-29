/**
 * @file
 * @brief MOSConversion class hdr.
 *
 * Header file for the MOS Conversion class
 */

#ifndef MOSConversion_H
#define MOSConversion_H

#include <map>

#include "quicksort.h"
#include "MOSConfig.h"
#include "genomes/GAGenome.h"

/**
 * Class: MOSConversion
 *
 * @brief Descripción de la clase MOSConversion
 */
class MOSConversion {

 public:

  typedef bool (*ConversionFunction) (GAGenome* genOrig, GAGenome* genDest);
  typedef std::pair<encodingType, encodingType> convertKey;

  // Registrar función de conversión
  bool registerConvFunction(const encodingType codOrig, const encodingType codDest, ConversionFunction convFunction);

  // Convertir genoma de una codificación origen a una destino
  bool convertGenome(encodingType codOrig, encodingType codDest, GAGenome* genOrig, GAGenome* genDest);

  // Método para obtener instancia única
  static MOSConversion* handle();

  static void destroy();

 protected:

  std::map<convertKey, ConversionFunction> functionSet;
  static MOSConversion* pSelf; // Puntero a la instancia si existe
  MOSConversion();
  virtual ~MOSConversion() {}

};

#endif

/**
 * @file
 * @brief MOSConversion class impl.
 *
 * Implementation of the MOSConversion class
 */

#include "MOSConversion.h"

MOSConversion* MOSConversion::pSelf = NULL;

/**
 * Registrar función de conversión
 */
bool MOSConversion::registerConvFunction(encodingType codOrig, encodingType codDest, ConversionFunction convFunction) {
  functionSet[std::make_pair(codOrig, codDest)] = convFunction;
  return true;
}


/**
 * Convertir genoma de una codificación origen a una destino
 */
bool MOSConversion::convertGenome(encodingType codOrig, encodingType codDest, GAGenome* genOrig, GAGenome* genDest) {
  return functionSet[std::make_pair(codOrig, codDest)](genOrig, genDest);
}

/**
 * Manejador de instancia única del singleton
 */
MOSConversion* MOSConversion::handle() {
  if( !pSelf ) {
    pSelf = new MOSConversion;
  }
  return pSelf;
}

void MOSConversion::destroy() {
  if (pSelf)
   delete pSelf;
  pSelf = NULL;
}

/**
 * Constructor
 */
MOSConversion::MOSConversion() {
}

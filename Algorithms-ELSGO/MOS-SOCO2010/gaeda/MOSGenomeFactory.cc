/**
 * @file
 * @brief MOSGenomeFactory class impl.
 *
 * Implementation of the MOSGenomeFactory class
 */

#include "MOSGenomeFactory.h"

MOSGenomeFactory* MOSGenomeFactory::pSelf = NULL;

/**
 * Manejador de instancia única del singleton
 */
MOSGenomeFactory* MOSGenomeFactory::handle() {
  if( !pSelf ) {
    pSelf = new MOSGenomeFactory;
  }
  return pSelf;
}

void MOSGenomeFactory::destroy() {
  if (pSelf)
   delete pSelf;
  pSelf = NULL;
}

/**
 * Constructor
 */
MOSGenomeFactory::MOSGenomeFactory() {
}

/**
 * Destructor
 */
MOSGenomeFactory::~MOSGenomeFactory() {
  std::map<encodingType, GAGenome*>::iterator it;
  for (it = genomeSet.begin(); it != genomeSet.end(); it++)
    delete it->second;
  genomeSet.clear();
}

/**
 *
 */
GAGenome* MOSGenomeFactory::getGenome(encodingType encoding) {
  return genomeSet[encoding]->clone();
}

/**
 *
 */
void MOSGenomeFactory::registerGenome(encodingType encoding, GAGenome* genome) {
  genomeSet[encoding] = genome->clone();
  return;
}

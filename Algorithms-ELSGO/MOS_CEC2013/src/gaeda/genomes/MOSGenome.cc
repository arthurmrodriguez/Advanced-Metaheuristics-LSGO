#include <algorithm>
#include <iostream>
#include <iomanip>

#include "MOSGenome.h"

#include "../std_stream.h"
#include "../GAEDAConfig.h"
#include "../MOSConversion.h"
#include "../MOSTechnique.h"
#include "../MOSTechniqueSet.h"
#include "../MOSGenomeFactory.h"

const char* MOSGenome::className () const {return "MOSGenome";}
int MOSGenome::classID () const {return GAID::GenomeMOS;}

/**
 * MOSGenome
 * @param tech
 */
MOSGenome::MOSGenome(MOSTechnique* technique)  {
  _technique = technique;
  _defaultEncoding = technique->getEncoding();
  _defaultGenome = technique->getGenome();
  _encodingSet[_defaultEncoding] = _defaultGenome;
  _fit_inc = 0;
  _mustbenulled = false;
  _mustquality = false;

  MOSTechniqueSet::MOSTechniqueSetIterator it;
  MOSTechniqueSet* techniqueSet = MOSTechniqueSet::handle();

  for (it=techniqueSet->begin(); it!=techniqueSet->end(); it++)
    _probVector[it->second->getId()] = 0.0;

  _probVector[technique->getId()] = 1.0;

  return;
}

/**
 * MOSGenome
 * @param orig
 */
MOSGenome::MOSGenome (const GAGenome& origin) {
  copy(origin);
  return;
}

/**
 * ~MOSGenome
 */
MOSGenome::~MOSGenome () {
  /* Aqui hay que recorrer el map de codificaciones y cepillarse a todos los genomas* que hay dentro*/
  MOSEncodingSet::iterator it;
  for(it=_encodingSet.begin(); it!=_encodingSet.end(); it++) {
    delete it->second;
  }
  _encodingSet.clear();
  return;
}

/**
 * clone
 */
MOSGenome* MOSGenome::clone (GAGenome::CloneMethod) const {
  return new MOSGenome((GAGenome&)*this);
}

/**
 * PRE: Tiene solo una codificacion. Previamente se ha purgado
 * El inherit sirve para heredar la información de un padre a un hijo, recuperando la información de la representación por defecto del hijo en el padre.
 * inherit
 * @param origin
 */
void MOSGenome::inherit (const MOSGenome& origin) {
  GAGenome::copy (origin);
  _encodingSet[_defaultEncoding]->copy(*(origin.getGenome(_defaultEncoding)));
  _defaultGenome = _encodingSet[_defaultEncoding];
  return;
}

/**
 * copy
 * @param origin
 * POST: Se copia la codificacion por defecto, si hay otras codificaciones no se copian.
 * Se debería ejecutar primero el purgue, y no haber ninguna codificación más que la de por defecto,
 * que debe ser la misma para mi y para el origen
 * NOTA: Se usa tambien en la construccion de la poblacion inicial
 */
void MOSGenome::copy (const GAGenome& orig) {

  const MOSGenome& origin = DYN_CAST (const MOSGenome&, orig);

  _defaultEncoding = origin._defaultEncoding;
  _technique = origin._technique;

  purgeGenome(_technique);

  GAGenome::copy (orig);

  // Copiamos la codificación por defecto
  if (!existEncoding(_defaultEncoding)) {
    _encodingSet[_defaultEncoding] = _technique->getGenome();
  }
  _encodingSet[_defaultEncoding]->copy(*(origin.getDefaultGenome()));
  _defaultGenome = _encodingSet[_defaultEncoding];

  // Copiamos las probabilidades de las técnicas
  _probVector=origin.getProbVector();

  _fit_inc = origin._fit_inc;
  _mustbenulled = origin._mustbenulled;
  _mustquality  = origin._mustquality;

  return;
}


/**
 * getTechProb
 */
double MOSGenome::getTechProb(techIdType techniqueId) const {
  MOSProbVector::const_iterator it = _probVector.find(techniqueId);
  if (it == _probVector.end()) {
    return 0.0;
  } else {
    return it->second;
  }
}


/**
 * addEncoding
 */
void MOSGenome::addEncoding(encodingType encoding, GAGenome* genome) {
  MOSEncodingSet::iterator it = _encodingSet.find(encoding);
  if (it == _encodingSet.end()) {
    _encodingSet[encoding] = genome;
  } else {
    delete _encodingSet[encoding];
    _encodingSet[encoding] = genome;
    if (encoding == _defaultEncoding) {
      _defaultGenome = genome;
    }
  }
  return;
}

/**
 * updateProbVector
 * @param newProbVector
 */
void MOSGenome::updateProbVector(MOSProbVector& probVector) {
  _probVector = probVector;
  return;
}

/**
 * getTechniqueId
 */
techIdType MOSGenome::getTechniqueId() const {
  return _technique->getId();
}

/**
 * existEncoding
 * @param encoding
 */
bool MOSGenome::existEncoding(encodingType encoding) {
  MOSEncodingSet::iterator it = _encodingSet.find(encoding);
  if (it == _encodingSet.end())
    return false;
  else
    return true;
}

/**
 * cleanEncoding
 * Elimina todos los genomas contenidos salvo el de la codificación que se pasa como argumento que pasa a ser la codificación y el genoma por defecto. Debe usarse solo en los hijos cuando hacemos copia de genomas de generaciones anteriores y queremos limpiarlos.
 * @param encoding
 */
void MOSGenome::cleanEncodings(encodingType encoding) {
  MOSEncodingSet::iterator it;
  for(it=_encodingSet.begin(); it!=_encodingSet.end(); ) {
    if (it->first != encoding) {
      delete it->second;
      _encodingSet.erase(it++);
    } else {
      _defaultGenome = _encodingSet[encoding];
      it++;
    }
  }
  return;
}

/**
 * purgeGenome
 * Se usa para inicializar los hijos antes de hacer el cruce
 */
void MOSGenome::purgeGenome(MOSTechnique* technique) {
  _evaluated=gaFalse;
  _neval=0;
  _technique=technique;
  _defaultEncoding = technique->getEncoding();
  cleanEncodings(_defaultEncoding);
  if (!existEncoding(_defaultEncoding)) {
    _defaultGenome = technique->getGenome();
    addEncoding(_defaultEncoding, _defaultGenome);
  }
  else
    _defaultGenome->resetGenomeInfo();
  _fit_inc = 0;
  _mustbenulled = false;
  _mustquality  = false;
  return;
}

/**
 * initialize
 */
void MOSGenome::initialize() {
  _evaluated=gaFalse;
  _neval=0;
  _technique->initGenome(*_defaultGenome);
  _fit_inc = 0;
  _mustbenulled = false;
  _mustquality  = false;
  return;
}

/**
 * evaluate
 * @param flag
 */
double MOSGenome::evaluate(GABoolean flag) const {
  if (_evaluated == gaFalse || flag == gaTrue) {
    MOSGenome* This = (MOSGenome*) this;
    This->_score = _technique->evalGenome(*_defaultGenome);
    This->_neval++;
    This->_evaluated = gaTrue;

    MOSEncodingSet::const_iterator it = _encodingSet.begin();
    for (;it!=_encodingSet.end();it++){
      it->second->score(This->_score);
    }
  }
  return _score;
}

/**
 * getGenome
 */
GAGenome* MOSGenome::getGenome(encodingType encoding) const {
  MOSEncodingSet::const_iterator it = _encodingSet.find(encoding);
  if (it == _encodingSet.end())
    return NULL;
  else
    return it->second;
}


/**
 * write
 * @param os
 */
int MOSGenome::write(STD_OSTREAM& os) const {
  MOSEncodingSet::const_iterator it;

  if (GAEDAConfig::handle()->getMOSGenomePrinting() == GAEDAConfig::FullPrinting) {
    for(it=_encodingSet.begin(); it!=_encodingSet.end(); it++) {
      if (_defaultEncoding == it->first) {
	//os << "default: ";
      } else {
	os << "alternative: ";
      }
      //os << it->first << ": ";
      it->second->write(os);
      os << "\n";
    }
  }
  else {
    for(it=_encodingSet.begin(); it!=_encodingSet.end(); it++) {
      if (_defaultEncoding == it->first) {
	//os << "default: ";
	//os << it->first << ": ";
	it->second->write(os);
	os << "\n";
      }
    }
  }

  if (GAEDAConfig::handle()->getEvolutiveApproach() == AutonomicEvolution) {
    MOSProbVector::const_iterator itp;
    os << "probabilities: [";
    for(itp=_probVector.begin(); itp!=_probVector.end(); itp++) {
      os << itp->second;
      if (itp != --_probVector.end())
        os << ", ";
    }
    os << "]" << std::endl;
  }

  return 1;
}

/**
 * Wrapper getValueOfNominalGene
 */
long MOSGenome::getValueOfNominalGene(int g) const {
  MOSEncodingSet::const_iterator it = _encodingSet.find(_nextEdaTech->getEncoding());
  if (it == _encodingSet.end()) {
    return -1;
  } else {
    return it->second->getValueOfNominalGene(g);
  }
}

/**
 * Wrapper setValueOfNominalGene
 */
long MOSGenome::setValueOfNominalGene(int g, long value) {
  MOSEncodingSet::const_iterator it = _encodingSet.find(_nextEdaTech->getEncoding());
  if (it == _encodingSet.end()) {
    return -1;
  } else {
    return it->second->setValueOfNominalGene(g,value);
  }
}

/**
 * Wrapper getValueOfContinuosGene
 */
double MOSGenome::getValueOfContinuousGene(int g) const {
  MOSEncodingSet::const_iterator it = _encodingSet.find(_nextEdaTech->getEncoding());
  if (it == _encodingSet.end()) {
    return -1.0;
  } else {
    return it->second->getValueOfContinuousGene(g);
  }
}

/**
 * Wrapper setValueOfContinuosGene
 */
double MOSGenome::setValueOfContinuousGene(int g, double value) {
  MOSEncodingSet::const_iterator it = _encodingSet.find(_nextEdaTech->getEncoding());
  if (it == _encodingSet.end()) {
    return -1.0;
  } else {
    return it->second->setValueOfContinuousGene(g, value);
  }
}

/**
 * Wrapper fixedSize
 */
int MOSGenome::fixedSize() const {
  MOSEncodingSet::const_iterator it = _encodingSet.find(_nextEdaTech->getEncoding());
  if (it == _encodingSet.end()) {
    return -1;
  } else {
    return it->second->fixedSize();
  }
}

/**
 * Serialización del genoma
 */
void MOSGenome::writeObject(STD_OSTREAM& os) const {
  GAGenome::writeObject(os);
  // Serializar id tecnica con la que se contruyo
  techIdType techId = this->getTechniqueId();
  os.write((char*)&techId,sizeof(techIdType));

  os.write((char*)&_fit_inc,sizeof(double));

  // Serializar numero de codificaciones
  int encodingSetSize = _encodingSet.size();
  os.write((char*)&encodingSetSize,sizeof(int));

  // Serializar conjunto de codificaciones
  MOSGenome::MOSEncodingSet::const_iterator itE;
  for(itE=_encodingSet.begin(); itE!=_encodingSet.end(); itE++) { // Recorrer set de encodings, llamar serializadores
    os.write((char*)&itE->first,sizeof(encodingType)); // Serializar identificador de la codificacion
    itE->second->writeObject(os);
  }

  // Serializar tamano vector probabilidades
  int probVectorSize = _probVector.size();
  os.write((char*)&probVectorSize,sizeof(int));

  // Serializar vector de probabilidades
  MOSProbVector::const_iterator itP;
  for(itP=_probVector.begin(); itP!=_probVector.end(); itP++) { // Recorrer el vector de probabilidades
    os.write((char*)&itP->first,sizeof(techIdType)); // Guardar id tecnica
    os.write((char*)&itP->second,sizeof(double)); // Guardar prob tecnica
  }
  os.write((char*)&_mustbenulled,sizeof(bool));
  os.write((char*)&_mustquality, sizeof(bool));
}

/**
 * Deserialización del genoma
 */
void MOSGenome::readObject(STD_ISTREAM& is) {
  GAGenome::readObject(is);

  // Deserializar identificador de tecnica
  techIdType techId;
  is.read((char*)&techId,sizeof(techIdType));
  // Inicializar estructura basica del MOS genome en funcion de la tecnica
  _technique=MOSTechniqueSet::handle()->getTechnique(techId);
  _defaultEncoding = _technique->getEncoding();
  _defaultGenome = _technique->getGenome();
  _encodingSet[_defaultEncoding] = _defaultGenome;
  _probVector[_technique->getId()] = 1.0;

  is.read((char*)&_fit_inc,sizeof(double));

  // Deserializar numero de codificaciones
  int encodingSetSize;
  is.read((char*)&encodingSetSize,sizeof(int));

  // Deserializar conjunto de codificaciones
  for(int i=0; i<encodingSetSize; i++) {
    encodingType encoding;
    is.read((char*)&encoding,sizeof(encodingType)); // Recoger id encoding
    GAGenome* genome = MOSGenomeFactory::handle()->getGenome(encoding); // Fabricar genoma adecuado a codificación
    genome->readObject(is); // deserializacion genome concreto
    this->addEncoding(encoding,genome); // añadir a conjunto representaciones
  }

  // Deserializar tamaño vector probabilidad
  int probVectorSize;
  is.read((char*)&probVectorSize,sizeof(int));

  // Deserializar vector probabilidad
  for(int i=0; i<probVectorSize; i++) {
    techIdType techId;
    is.read((char*)&techId,sizeof(techIdType));
    double prob;
    is.read((char*)&prob,sizeof(double));
    this->updateTechProb(techId,prob);
  }
  is.read((char*)&_mustbenulled,sizeof(bool));
  is.read((char*)&_mustquality,sizeof(bool));
}

double MOSGenome::compare (const GAGenome& g) const {
  MOSGenome* genome = const_cast<MOSGenome*> ( dynamic_cast<const MOSGenome*>(&g) );

  if (!genome->existEncoding(_defaultEncoding)) {
    GAGenome* newGenome = _technique->getGenome();
    genome->addEncoding(_defaultEncoding, newGenome);
    MOSConversion::handle()->convertGenome(genome->getDefaultEncoding(), _defaultEncoding, genome->getDefaultGenome(), newGenome);
  }

  return _defaultGenome->compare(* genome->getGenome(_defaultEncoding) );
}

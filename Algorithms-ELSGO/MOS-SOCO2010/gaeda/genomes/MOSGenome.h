#ifndef _GA_MOS_GENOME_H_
#define _GA_MOS_GENOME_H_

#include <map>

#include "GAGenome.h"
#include "../MOSConfig.h"

class MOSTechnique;

class MOSGenome : public GAGenome {

 public:

  GADeclareIdentity();

  // Definición tipo vector con las diferentes representaciones del genoma
  typedef std::map<encodingType, GAGenome*> MOSEncodingSet;

  MOSGenome(MOSTechnique* tech);
  MOSGenome(const GAGenome& orig);

  ~MOSGenome();

  MOSGenome* clone(GAGenome::CloneMethod flag = GAGenome::CONTENTS) const;

  void copy(const GAGenome& orig);

  void inherit(const MOSGenome& origin);

  int write(STD_OSTREAM& os) const;

  const MOSProbVector& getProbVector() const {return _probVector;}

  long double getTechProb(techIdType techniqueId) const;

  void updateProbVector(MOSProbVector& probVector);

  void updateTechProb(techIdType techniqueId, long double prob) {_probVector[techniqueId] = prob;}

  techIdType getTechniqueId() const;

  // Devuelve true si existe una representacion del genoma en la codificación dada
  bool existEncoding(encodingType encoding);

  void addEncoding(encodingType encoding, GAGenome* genome);

  int getDefaultEncoding() const {return _defaultEncoding;}

  GAGenome* getDefaultGenome() const {return _defaultGenome;}

  GAGenome* getGenome(encodingType encoding) const;

  const MOSEncodingSet& getEncodingSet() const {return _encodingSet;}

  MOSTechnique* getTechnique(                 ) const {return _technique;}
  void          setTechnique(MOSTechnique* technique) {_technique = technique;}

  void purgeGenome(MOSTechnique* technique);

  void initialize();

  long double evaluate(GABoolean flag = gaFalse) const;

  void writeObject(STD_OSTREAM& os) const;
  void readObject (STD_ISTREAM& is);

  long double compare(const GAGenome&) const;

  long double setFitnessIncrement(long double fit_inc) {return _fit_inc = fit_inc;}
  long double getFitnessIncrement(              ) {return _fit_inc;          }

  bool mustComputeQuality (bool state) {return _mustquality = state;}
  bool mustComputeQuality (          ) {return _mustquality;        }

protected:
  void cleanEncodings(encodingType encoding);

  MOSTechnique* _technique;    // Tecnica con la que se creo el individuo
  GAGenome* _defaultGenome;      // Puntero al genoma inicial
  encodingType _defaultEncoding; // Codificación inicial
  MOSProbVector _probVector;   // Vector con las probabilidades de las tecnicas
  MOSEncodingSet _encodingSet; // Map con las representaciones
  long double _fit_inc;
  bool _mustquality;
};

#endif

/**
 * @file
 * @brief MOSTechnique class hdr.
 *
 */

#ifndef MOSTechnique_H
#define MOSTechnique_H

#include "gaid.h"
#include "MOSConfig.h"
#include "GAStatistics.h"
#include "Recombinator.h"

class GAPopulation;
class MOSGenome;
class MOSQuality;

/**
 * @brief Clase que representa un técnica
 */
class MOSTechnique : public GAID {

 public:

  GADefineIdentity("MOSTechnique", GAID::Technique);

  MOSTechnique();

  /* Destructor */
  virtual ~MOSTechnique() {
   if (_genomeBase)
      delete _genomeBase;
   if (_selector)
      delete _selector;
   if (_recombinator)
      delete _recombinator;
  }

  virtual long double evolve(GAPopulation*& pop, unsigned maxEvals, MOSQuality* qualityFunction, bool& converged);
  virtual unsigned offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged) = 0;

  /* Obtiene el identificador de la técnica */
  techIdType getId() const {return _id;}

  // Get/Set the quality of the technique
  long double getQuality(bool force = false ) const {if (!MOSTechnique::improvement_override || force) return _quality;
                                                     else                                              return getRatioImproved();}
  long double setQuality(long double quality)       {return _quality = quality;}

  // Get/Set the participation ration of the technique
  long double getPartRatio(                ) const {return _partRatio;}
  long double setPartRatio(long double partRatio)       {return _partRatio = partRatio;}

  /* Devuelve el identificador de la codificación de los genomas que usa esta técnica */
  int getEncoding() {return _encoding;}

  /* Devuelve una muestra de un genoma de los que usa esta técnica (usa GenomeFactory) */
  GAGenome* getGenome() {return _genomeBase->clone();}

  // Methods to retrieve local statistics from the main algorithm
  unsigned long getSelections  () {return _stats.selections  ();}
  unsigned long getCrossovers  () {return _stats.crossovers  ();}
  unsigned long getMutations   () {return _stats.mutations   ();}
  unsigned long getReplacements() {return _stats.replacements();}
  unsigned long getEvals       () {return _stats.indEvals    ();}

  // Initializes the given genome with the individual initializer associated to this technique
  void initGenome(GAGenome& gen) {(*_initializer)(gen);}

  // Evaluates the given genome with the objective function associated to this technique
  long double evalGenome(GAGenome& gen) {return (*_evaluator)(gen);}

  // Returns the description of this technique
  const std::string& getDescription() {return _description;}

  unsigned    setTimesImproved(unsigned times)       {return _times_improved = times;}
  unsigned    getTimesImproved(              ) const {return _times_improved;        }
  long double getRatioImproved(              ) const {return (long double) _times_improved / (long double) (_evals_last_gen);}

 protected:

  techIdType                _id;            // Identificador
  std::string               _description;   // Descripción
  GAGenome::Initializer     _initializer;   // Operador de inicialización
  GAGenome::Evaluator       _evaluator;     // Operador de evaluación
  long double               _quality;       // Evaluación de la calidad de la técnica en función de lo buenos que son los genomas que genera
  long double               _partRatio;     // Ratio de participación de la técnica
  GAGenome*                 _genomeBase;    // Puntero al genoma base
  GAStatistics              _stats;         // Estadísticas
  encodingType              _encoding;      // Codificación de los genomas que usa la técnica
  GASelectionScheme*        _selector;      // Selector de padres para el offspring
  Recombinator*             _recombinator;  // Elitism
  unsigned                  _times_improved;
  unsigned                  _evals_last_gen;

 public:
  static bool               improvement_override;
};

#endif

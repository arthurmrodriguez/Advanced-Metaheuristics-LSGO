#ifndef _MOS_QUALITY_FUNC
#define _MOS_QUALITY_FUNC

#include "gaid.h"
#include "MOSConfig.h"

class GAPopulation;

// Base quality functor
class MOSQuality {
  public:
    GADefineIdentity("MOSQuality", GAID::Quality);
    virtual long double update(const GAPopulation& pop, techIdType id, unsigned newInds) = 0;
};


// Average fitness quality functor
class AverageFitnessQuality : virtual public MOSQuality {
  public:
    GADefineIdentity("QualityFitAvg", GAID::QualityFit);
    virtual long double update         (const GAPopulation& pop, techIdType id, unsigned newInds) {return averageFitness(pop, id, newInds);}
  protected:
            long double averageFitness (const GAPopulation& pop, techIdType id, unsigned newInds);
};


// Average fitness increment quality functor
class AverageFitnessIncrementQuality : virtual public MOSQuality {
  public:
    GADefineIdentity("QualityFitInc", GAID::QualityFitInc);
    virtual long double update                  (const GAPopulation& pop, techIdType id, unsigned newInds) {return averageFitnessIncrement(pop, id, newInds);}
  protected:
            long double averageFitnessIncrement (const GAPopulation& pop, techIdType id, unsigned newInds);
};


// Average diversity quality functor
class AverageDiversityQuality : virtual public MOSQuality {
  public:
    GADefineIdentity("QualityDiv", GAID::QualityDiv);
    virtual long double update           (const GAPopulation& pop, techIdType id, unsigned newInds) {return averageDiversity(pop, id, newInds);}
  protected:
            long double averageDiversity (const GAPopulation& pop, techIdType id, unsigned newInds);
};


// Compass quality functor
class CompassQuality: public AverageFitnessQuality, AverageDiversityQuality {
  public:
    GADefineIdentity("QualityCompass", GAID::QualityCompass);
    virtual long double update(const GAPopulation& pop, techIdType id, unsigned newInds) {return (averageFitness(pop, id, newInds) + averageDiversity(pop, id, newInds)) / 2.0;}
};

#endif

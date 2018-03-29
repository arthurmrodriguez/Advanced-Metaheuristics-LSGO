#ifndef MOSPARTICIPATIONFUNCTION_H
#define MOSPARTICIPATIONFUNCTION_H

#include "gaid.h"
#include <stdexcept>

class MOSEA2;
class MOSQuality;

#define MINPART    0.05
#define RESTOREGEN 0
#define PFDATA     (void*) 0
#define ADJUST     0.05

class MOSParticipation {
  public:
    GADefineIdentity("MOSParticipation", GAID::Participation);
    MOSParticipation(long double minPart = MINPART, unsigned restoreGen = RESTOREGEN, void* pfData = PFDATA) : _minPart(minPart), _restoreGen(restoreGen), _pfData(pfData) {}

    virtual long double setBaseQual (MOSEA2& alg) { 
      throw std::runtime_error("Error: setBaseQual should not be called in MOSParticipation class"); 
      return 0.0;
    }

    virtual void update(MOSEA2& algorithm) = 0;

  protected:
    long double   _minPart;
    unsigned _restoreGen;
    void*    _pfData;
};


class ConstantParticipation : public MOSParticipation {
  public:
    GADefineIdentity("ConstantParticipation", GAID::ConstantPF);
    ConstantParticipation(long double minPart = MINPART, unsigned restoreGen = RESTOREGEN, void* pfData = PFDATA);

    long double setBaseQual (MOSEA2& alg) {return 0.0;}
    virtual void update(MOSEA2& algorithm) {}
};


class DynamicParticipation : public MOSParticipation {
  public:
    GADefineIdentity("DynamicParticipation", GAID::DynamicPF);
     DynamicParticipation(long double adjust = ADJUST, long double minPart = MINPART, unsigned restoreGen = RESTOREGEN, void* pfData = PFDATA);
    ~DynamicParticipation();

    long double setBaseQual (MOSEA2& alg);

    virtual void update(MOSEA2& algorithm);

  protected:
    long double      _adjust;
    long double      _base;
};

#endif

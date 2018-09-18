#ifndef PSOGENOME_H_
#define PSOGENOME_H_

#include "GARealGenome.h"

class PSOGenome : public GARealGenome {
protected:
  GARealGenome* best_;
  GARealGenome* vel_;
public:
  PSOGenome(const PSOGenome&   );
  PSOGenome(const GARealGenome&);

  virtual ~PSOGenome();

	GARealGenome*  best      (                          ) const;
	GARealGenome*  velocity  (                          ) const;
	virtual void   updateBest(                          );
	virtual void   copy      (const PSOGenome&          );
	PSOGenome*     clone     (GAGenome::CloneMethod flag) const;
};



#endif /*PSOGENOME_H_*/

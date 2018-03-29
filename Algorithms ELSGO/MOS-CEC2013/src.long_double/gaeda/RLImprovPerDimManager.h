#ifndef RLIMPROVPERDIMMANAGER_H_
#define RLIMPROVPERDIMMANAGER_H_

#include "ImprovPerDimManager.h"

class RLImprovPerDimManager: public ImprovPerDimManager {

protected:

  long double         alpha_;
  vector<long double> values_;
  vector<bool>   explored_dims_;
  bool           all_dims_explored_;

  virtual vector<long double> dimPerfValues() const;

public:

  RLImprovPerDimManager(int num_dims, long double alpha=0.7);
  virtual ~RLImprovPerDimManager();

  virtual void updateValue  (int pos, long double new_value );
};

#endif /* RLIMPROVPERDIMMANAGER_H_ */

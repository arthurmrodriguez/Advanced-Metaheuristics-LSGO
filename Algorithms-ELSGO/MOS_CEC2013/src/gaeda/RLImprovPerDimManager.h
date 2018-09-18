#ifndef RLIMPROVPERDIMMANAGER_H_
#define RLIMPROVPERDIMMANAGER_H_

#include "ImprovPerDimManager.h"

class RLImprovPerDimManager: public ImprovPerDimManager {

protected:

  double         alpha_;
  vector<double> values_;
  vector<bool>   explored_dims_;
  bool           all_dims_explored_;

  virtual vector<double> dimPerfValues() const;

public:

  RLImprovPerDimManager(int num_dims, double alpha=0.7);
  virtual ~RLImprovPerDimManager();

  virtual void updateValue  (int pos, double new_value );
};

#endif /* RLIMPROVPERDIMMANAGER_H_ */

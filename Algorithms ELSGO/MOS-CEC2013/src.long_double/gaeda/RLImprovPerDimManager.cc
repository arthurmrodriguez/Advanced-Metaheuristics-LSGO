#include "RLImprovPerDimManager.h"
#include "logger/GALogger.h"
#include <sstream>

RLImprovPerDimManager::RLImprovPerDimManager(int num_dims, long double alpha) : ImprovPerDimManager(num_dims),
                                                                           alpha_(alpha),
                                                                           values_(num_dims),
                                                                           explored_dims_(num_dims),
                                                                           all_dims_explored_(false) {
  assert(alpha > 0.0 && alpha <= 1.0);                                                                            
  for (int i=0; i<numDims(); i++) {
    values_[i]        = 1.0;
    explored_dims_[i] = false;
  }
}

RLImprovPerDimManager::~RLImprovPerDimManager() {}

vector<long double> RLImprovPerDimManager::dimPerfValues() const {
  if (all_dims_explored_) return values_;
  else {
    vector<long double> mod_values = values_; 

    long double max_value = values_[0];
    for (int i=1; i<values_.size(); i++) {
      if (values_[i] > max_value) max_value = values_[i];
    }

    for (int i=0; i<explored_dims_.size(); i++) {
      if (explored_dims_[i] == false) mod_values[i] = max_value + 1;
    }
    
    return mod_values;
  }
}

void RLImprovPerDimManager::updateValue (int pos, long double new_value) {
  // /*LOG*/ stringstream msg; msg << "dim: " << pos << " old RL value:" << values_[pos] << " new pased value=" << new_value;
  values_[pos]        = (explored_dims_[pos] == false) ? new_value : values_[pos]*alpha_ + (1-alpha_)*new_value;
  explored_dims_[pos] = true;
  // /*LOG*/ msg << " new computed RL value:" << values_[pos] << endl;
  // /*LOG*/ GALogger::instance()->appendLogMessage("RLImprovPerDimManager", msg.str());

  if (!all_dims_explored_) {
    all_dims_explored_ = true;
    for (int i=0; i<explored_dims_.size(); i++) { 
      if (explored_dims_[i] == false) {
        all_dims_explored_ = false;
        break;
      }
    }
  }
}

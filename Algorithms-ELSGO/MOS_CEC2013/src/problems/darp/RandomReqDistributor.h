#ifndef RANDOMREQDISTRIBUTOR_H_
#define RANDOMREQDISTRIBUTOR_H_

#include "ReqDistributor.h"

class RandomReqDistributor : public ReqDistributor {
public:
  RandomReqDistributor() : ReqDistributor() {}
  virtual ~RandomReqDistributor() {}

  vector< vector<Request> > distributeReqs(vector<Request>& reqs, ReqDistMatrix& reqdistmatrix, int node);

  virtual string name() { return "Random request distributor"; }
};

#endif /* RANDOMREQDISTRIBUTOR_H_ */

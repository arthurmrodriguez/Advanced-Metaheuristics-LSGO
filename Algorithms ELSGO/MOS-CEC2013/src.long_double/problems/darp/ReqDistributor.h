#ifndef REQDISTRIBUTOR_H_
#define REQDISTRIBUTOR_H_

#include "VerticesList.h"
#include "ReqDistMatrix.h"
#include <vector>

using namespace std;

class ReqDistributor {
public:
  ReqDistributor() {}
  virtual ~ReqDistributor() {}

  virtual vector< vector<Request> > distributeReqs(vector<Request>& reqs, ReqDistMatrix& reqdistmatrix, int node) = 0;

  virtual string name() = 0;
};

#endif /* REQDISTRIBUTOR_H_ */

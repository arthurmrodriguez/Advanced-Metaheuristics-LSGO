#ifndef KMEDIODSREQDISTRIBUTOR_H_
#define KMEDIODSREQDISTRIBUTOR_H_

#include "ReqDistributor.h"
#include "VerticesList.h"
#include "ReqDistMatrix.h"
#include <vector>
#include <sstream>

#include "extras/cluster.h"

using namespace std;

class KMedoidsReqDistributor : public ReqDistributor {
protected:
  int          npass_;
  kmedoidsMode mode_;

public:

  KMedoidsReqDistributor(int npass, kmedoidsMode mode);
  virtual ~KMedoidsReqDistributor() {}

  vector< vector<Request> > distributeReqs(vector<Request>& reqs, ReqDistMatrix& reqdistmatrix, int nodes);

  virtual string name() {
    stringstream msg;
    msg << "K-Medoids request distributor npass: " << npass_ << " mode: ";
    switch (mode_) {
      case KMEDNORMAL:          msg << "normal";  break;
      case KMEDSELBESTBALANCED: msg << "best balanced"; break;
      case KMEDFORCEBALANCE:    msg << "force balanced"; break;
    }
    msg << endl;
    return msg.str();
  }
};

#endif /* KMEDIODSREQDISTRIBUTOR_H_ */

#ifndef KMEDOIDSROUTEDISTRIBUTOR_H_
#define KMEDOIDSROUTEDISTRIBUTOR_H_

#include "ReqDistMatrix.h"

#include <GAEDAlib.h>
#include <vector>

#include "extras/cluster.h"


class KMedoidsRouteDistributor : public RouteDistributor {
  // For efficiency we are going to store the same pointer as the one received in the constructor instead of performing
  // a copy.
  ReqDistMatrix* reqdistmatrix_;
  int            npass_;
  kmedoidsMode   mode_;

protected:
  virtual vector<int> createRoutesAssignment(vector<GAGenome*>& orig_routes);

public:
  KMedoidsRouteDistributor(ReqDistMatrix* reqdistmatrix, int npass, kmedoidsMode mode, RoutingGenome& gen, int nnodes);
  KMedoidsRouteDistributor(KMedoidsRouteDistributor& other);

  virtual ~KMedoidsRouteDistributor() {}

  virtual string name() { return "KMedoidsRouteDistributor"; }

  virtual RouteDistributor* clone() {return new KMedoidsRouteDistributor(*this); }

};

#endif /* KMEDOIDSROUTEDISTRIBUTOR_H_ */

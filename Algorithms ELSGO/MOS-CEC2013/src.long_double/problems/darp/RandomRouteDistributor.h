#ifndef RANDOMROUTEDISTRIBUTOR_H_
#define RANDOMROUTEDISTRIBUTOR_H_

#include <GAEDAlib.h>

class RandomRouteDistributor : public RouteDistributor {
protected:
  virtual vector<int> createRoutesAssignment(vector<GAGenome*>& orig_routes);

public:
  RandomRouteDistributor(RoutingGenome& gen, int nnodes) : RouteDistributor(gen,nnodes) {}
  RandomRouteDistributor(RandomRouteDistributor& other) : RouteDistributor(other) {}

  virtual ~RandomRouteDistributor() {}

  virtual string name() { return "RandomRouteDistributor"; }

  virtual RouteDistributor* clone() {return new RandomRouteDistributor(*this); }

};

#endif /* RANDOMROUTEDISTRIBUTOR_H_ */

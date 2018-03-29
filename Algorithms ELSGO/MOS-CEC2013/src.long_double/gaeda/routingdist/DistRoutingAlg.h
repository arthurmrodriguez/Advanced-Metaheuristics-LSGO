#ifndef DISTROUTINGALG_H
#define DISTROUTINGALG_H

#include "../islands/CommManager.h"
#include "RouteDistributor.h"
#include "../RoutingAlg.h"
#include "../genomes/GAGenome.h"
#include "../genomes/RoutingGenome.h"


/**
 * @brief Distributed Routing Algorithm Model
 *
 * This class represents the algorithm used for conducting the parallel distributed version
 * used with routing algorithms. It follows a pseudo master-slave approach where each slave
 * executes the routing algorithm over a subset of the whole set of routes. At certain points,
 * the slaves send their routes to the master node so a reassignment of the routes can be carried out.
 */

class DistRoutingAlg : public Algorithm {
public:
  enum emigrantRouteType {BEST,ACTUAL};
  // TODO refactore routingalg class so that both routing alg and routingdistclass inherit from an interface that
  // contains the methods that are common to both classes
protected:

  RoutingAlg*       alg_;
  RouteDistributor* rdistributor_;
  int               mig_period_;
  emigrantRouteType emigrant_route_;
  CommManager&      comm_manager_;

  bool haveAllIslandsConverged();

  void doMigration();

  RoutingGenome& selectIndForSendingRoutes();

  void uniteRoutes();

public:

  GADefineIdentity ("DistRoutingAlg", GAID::DistRoutingAlg);

  // This algorithm receives a genome although it should not need to receive this values but, until a refactorization
  // is conducted to transform several parts or RoutingAlg into an abstract class it is necessary to initialize
  // the RoutingAlg class
  DistRoutingAlg(RoutingAlg&       alg,
                 int               mig_period,
                 emigrantRouteType emigrant_route,
                 GAGenome&         gen,
                 CommManager&      comm_manager);

  DistRoutingAlg(const DistRoutingAlg& other);

  virtual ~DistRoutingAlg();

  DistRoutingAlg& operator= (const DistRoutingAlg& oth) {copy(oth); return *this;}
  void            copy(const DistRoutingAlg&);
  DistRoutingAlg* clone();

  void setRouteDistributor(RouteDistributor& rdist) { assert(&rdist); rdistributor_ = rdist.clone(); }

  RoutingAlg* getRoutingAlg() const {return alg_;}

  // Inherited. The algorithm redirects them to the correspondant method in alg_
  void evolve();

  void                initialize         ()       { return alg_->initialize(); }
  GAGenome*           best               () const { return alg_->best();       }
  const GAStatistics& statistics         () const { return alg_->statistics(); }
  const GAPopulation& population         () const { return alg_->population(); }
  void                step               ()       { return alg_->step(); }

};

#endif

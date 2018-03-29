/*
 * RouteDistributor.cc
 *
 *  Created on: Jun 19, 2012
 *      Author: santi
 */

#include "RouteDistributor.h"

RouteDistributor::RouteDistributor(RoutingGenome& gen, int nnodes)  : num_nodes_(nnodes) {
  // The initial set of cached values are created
  for (int i=0; i<nnodes; i++) {
    RoutingGenome* tmp = dynamic_cast<RoutingGenome*> (gen.clone()); assert(tmp);
    cached_routes_.push_back( tmp );
  }
}

RouteDistributor::RouteDistributor(const RouteDistributor& other) {
  copy(other);
}

RouteDistributor::~RouteDistributor() {
  for (int i=0; i<cached_routes_.size(); i++) delete cached_routes_[i];
}

void RouteDistributor::copy(const RouteDistributor& other) {
  num_nodes_ = other.num_nodes_;

  for (int i=0; i<cached_routes_.size(); i++) delete cached_routes_[i];
  cached_routes_.clear();
  for (int i=0; i<other.cached_routes_.size(); i++) {
    RoutingGenome* tmp = dynamic_cast<RoutingGenome*>( other.cached_routes_[i]->clone() ); assert(tmp);
    cached_routes_.push_back(tmp);
  }
}

void RouteDistributor::clearCachedRoutes() {
  for (int i=0; i<cached_routes_.size(); i++) {
    cached_routes_[i]->emptyRoutes();
  }
}

int RouteDistributor::numRoutes(vector<GAGenome*>& orig_routes) {
  int num_routes = 0;
  for (int i=0; i<orig_routes.size(); i++) {
    RoutingGenome& gen = dynamic_cast<RoutingGenome&>(* orig_routes[i]); assert(&gen);
    num_routes += gen.numRoutes();
  }

  return num_routes;
}

/*
 * Since we are using a vector of genomes for representing the whole set of routes sent/received to/from the nodes,
 * sometimes is necessary to obtain that the i-th routes corresponds to the j-th position of the k-th genome. This
 * method is used for computing this last values.
 */
void RouteDistributor::computeGenomeAndRoutePositionsOfRoutePos(vector<GAGenome*>& routes, int routepos,
                                                                /*out*/ int& genpos, /*out*/ int& genroutepos ) {
  assert(routes.size() > 0);
  assert(routepos >= 0 && routepos < numRoutes(routes));
  int sumgenroutes;
  for (genpos=0,sumgenroutes=0; genpos<routes.size(); genpos++) {
    RoutingGenome* gen = dynamic_cast<RoutingGenome*>(routes[genpos]);

    if (sumgenroutes + gen->numRoutes() > routepos) {
      genroutepos = routepos - sumgenroutes;
      break;
    }
    else {
      sumgenroutes += gen->numRoutes();
    }
  }

  assert(genpos < routes.size());
  assert(genroutepos >= 0 && genroutepos < dynamic_cast<RoutingGenome*>(routes[genpos])->numRoutes());
  assert(routes[genpos]);
}

void RouteDistributor::assignRoutesToCachedRoutes(vector<GAGenome*>& orig_routes, vector<int>& routes_assignment) {
  clearCachedRoutes();

  int pos = 0;
  for (int i=0; i<orig_routes.size(); i++) {
    RoutingGenome* orig_gen = dynamic_cast<RoutingGenome*>(orig_routes[i]);
    assert(orig_gen);

    for (int j=0; j<orig_gen->numRoutes(); j++, pos++) {
      assert(pos < routes_assignment.size());
      assert(orig_gen->routeLength(j) > 0 ); // Do not add empty routes

      cached_routes_[ routes_assignment[pos] ]->addRoute( *orig_gen, j);
    }
  }
}

/*
 * Return a collection of cached objects so there is no need to delete them
 */
vector<RoutingGenome*> RouteDistributor::assignRoutes2Nodes(vector<GAGenome*>& orig_routes) {
   assert(orig_routes.size() == num_nodes_ && cached_routes_.size() == num_nodes_);

   vector<int> routes_assignment = createRoutesAssignment(orig_routes);

   assignRoutesToCachedRoutes(orig_routes, routes_assignment);

   return cached_routes_;
}

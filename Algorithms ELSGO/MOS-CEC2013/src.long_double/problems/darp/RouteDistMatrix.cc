#include "RouteDistMatrix.h"
#include "DARPGenome.h"
#include "VerticesList.h"
#include <routingdist/RouteDistributor.h>
#include <limits>

long double scoreOfRoutes(GAGenome& g1, int route1, GAGenome& g2, int route2, ReqDistMatrix& reqdistmatrix) {
  assert(&g1);
  assert(&g2);
  DARPGenome& gen1 = dynamic_cast<DARPGenome&>(g1); assert(&gen1);
  DARPGenome& gen2 = dynamic_cast<DARPGenome&>(g2); assert(&gen2);

  if (gen1.routeLength(route1) == 0 || gen2.routeLength(route2) == 0) return 0; // Tol solve the problem when a route is empty

  long double min_score = std::numeric_limits<long double>::max();
  for (int i=0; i<gen1.routeLength(route1); i++) {
    if ( Vertex::isVertIdDelivery( gen1.gene(route1,i) ) ) continue; // Since we have a distance matrix of requests we only use
                                                                     // as input for the distance matrix the pickup ids (which have
                                                                     // stored as request ids in the req distance matrix
    for (int j=0; j<gen2.routeLength(route2); j++) {
      assert(gen1.gene(route1,i) != gen2.gene(route2,j));
      if ( Vertex::isVertIdDelivery( gen2.gene(route2,j) ) ) continue; // Same comment as before
      long double tmp_score = reqdistmatrix.getDistOfReqsPickupIds( gen1.gene(route1,i),gen2.gene(route2,j) );
      if (tmp_score < min_score) min_score = tmp_score;
    }
  }

  return min_score;
}

long double scoreOfRoutes(vector<GAGenome*>& routes, int route1pos, int route2pos, ReqDistMatrix& reqdistmatrix) {
  // -1 to check earlier if there is a mistake since they are output values and should be rewritten
  int genpos1=-1, genroutepos1=-1;
  int genpos2=-1, genroutepos2=-1;
  RouteDistributor::computeGenomeAndRoutePositionsOfRoutePos(routes,route1pos,genpos1,genroutepos1);
  RouteDistributor::computeGenomeAndRoutePositionsOfRoutePos(routes,route2pos,genpos2,genroutepos2);
  assert(genpos1 >= 0 & genpos2 >= 0);
  assert(genroutepos1 >= 0 & genroutepos1 >= 0);

  return scoreOfRoutes( * routes[genpos1], genroutepos1, * routes[genpos2], genroutepos2, reqdistmatrix);
}

RouteDistMatrix::RouteDistMatrix(vector<GAGenome*>& routes, ReqDistMatrix& reqdistmatrix) {
  // TODO: Think where to put the numRoutes method
  int numroutes = RouteDistributor::numRoutes(routes);

  resize(numroutes);

  // Initialized as specified in cluster.h
  for (int i=1; i<numroutes; i++) {
    resize(i,i);
    for (int j=0; j<i; j++) {
      setDist(i,j, scoreOfRoutes(routes,i,j,reqdistmatrix) );
    }
  }
}

RouteDistMatrix::~RouteDistMatrix() {}

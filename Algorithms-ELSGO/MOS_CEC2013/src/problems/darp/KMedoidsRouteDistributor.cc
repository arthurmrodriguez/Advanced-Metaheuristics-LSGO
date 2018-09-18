#include "KMedoidsRouteDistributor.h"
#include "RouteDistMatrix.h"

KMedoidsRouteDistributor::KMedoidsRouteDistributor(ReqDistMatrix* reqdistmatrix, int npass, kmedoidsMode mode,
                                                   RoutingGenome& gen, int nnodes) :
                                                                                        reqdistmatrix_(reqdistmatrix),
                                                                                        npass_(npass),
                                                                                        mode_(mode),
                                                                                        RouteDistributor(gen,nnodes) {
  assert(npass > 0);
}

KMedoidsRouteDistributor::KMedoidsRouteDistributor(KMedoidsRouteDistributor& other) : reqdistmatrix_(other.reqdistmatrix_),
                                                                                      npass_(other.npass_),
                                                                                      mode_(other.mode_),
                                                                                      RouteDistributor(other) {}

/*
 * Note that the received routes are in the format received from the communicator manager. Therefore, each genome
 * contains a number of routes and the total number of the routes corresponds to the sum of the routes that each genome contains
 */
vector<int> KMedoidsRouteDistributor::createRoutesAssignment(vector<GAGenome*>& routes) {
  setClusterFuncsSeed(GAGetRandomNoRankSeed()); // all the islands must use the same seed

  RouteDistMatrix routedistmatrix (routes, *reqdistmatrix_);

  int num_routes = numRoutes(routes);

  int    clusterids[num_routes]; // the ids of the cluster that each individual belong to
  double error;                  // The error of the best solution
  int    ifound;                 // Not used (the number of times the optimal solution is found

  kmedoids(num_nodes_,num_routes,routedistmatrix,npass_,clusterids,&error,&ifound,mode_);

  vector<int> dist_routes (clusterids, clusterids + sizeof(clusterids) / sizeof(int));

  assert(dist_routes.size() == num_routes);

#ifdef DEBUG
  vector<int> sum_assignments(num_nodes_,0);
  for (int i=0; i<dist_routes.size(); i++) {
    assert(dist_routes[i] >= 0 && dist_routes[i] < num_routes);
    sum_assignments[ dist_routes[i] ] += 1;
  }
  for (int i=0; i<num_routes; i++) {
    assert( sum_assignments[i] > 0 );
  }
#endif

  return dist_routes;
}

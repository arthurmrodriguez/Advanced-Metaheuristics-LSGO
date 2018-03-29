/*
 * RouteDistributor.h
 *
 *  Created on: Jun 19, 2012
 *      Author: santi
 */

#ifndef ROUTEDISTRIBUTOR_H_
#define ROUTEDISTRIBUTOR_H_

#include "../genomes/GAGenome.h"
#include "../genomes/RoutingGenome.h"
#include <vector>
#include <string>

using namespace std;

class RouteDistributor {
protected:

  int                    num_nodes_;
  vector<RoutingGenome*> cached_routes_;

  void clearCachedRoutes();
  /*
   * Each genome from orig_routes contains the information that is sent by each node to the master, i.e.,
   * the set of routes selected for migration whereas routes_assignment contains the final assignment of each
   * route to each node. It must be taken into account that routes_assignment assigns to each route a number
   * between 0..num_nodes -1 whereas its size is orig_routes[0].size() + orig_routes[1].size() + ...
   * orig_routes[num_nodes-1].size(). For example, if orig_routes contains two genomes with two and three
   * routes respectively, routes_assignment will have a size of five where the first two positions correspond to
   * the new assignments of the routes of the first genome from orig_routes whereas the las three positions
   * contain the assignments of the routes of the second genome.
   */
  void assignRoutesToCachedRoutes(vector<GAGenome*>& orig_routes, vector<int>& routes_assignment);

  virtual vector<int> createRoutesAssignment(vector<GAGenome*>& orig_routes) = 0;

public:

  RouteDistributor(RoutingGenome& gen, int nnodes);
  RouteDistributor(const RouteDistributor& other);
  virtual ~RouteDistributor();

  RouteDistributor& operator= (const RouteDistributor& oth) {copy(oth); return *this;}

  virtual void              copy(const RouteDistributor& other);
  virtual RouteDistributor* clone() = 0;

  /*
   * Returns cached genomes, no need to delete them! Not very clear but more efficient than having
   * to create and delete genomes each time is called.
   * Another confusing part is that the received genomes are the ones corresponding to the values received by
   * the communicator manager. Therefore, each genome from orig_routes contains the routes sent by each node
   * to the master node while each genome from the returned vector contains the routes that are finally
   * assigned to each node
   */
  virtual vector<RoutingGenome*> assignRoutes2Nodes(vector<GAGenome*>& orig_routes);

  virtual string name() = 0;

  static int  numRoutes(vector<GAGenome*>& orig_routes);
  static void computeGenomeAndRoutePositionsOfRoutePos(vector<GAGenome*>& routes, int routepos,
                                                       /*out*/ int& genpos, /*out*/ int& genroutepos );

};

#endif /* ROUTEDISTRIBUTOR_H_ */

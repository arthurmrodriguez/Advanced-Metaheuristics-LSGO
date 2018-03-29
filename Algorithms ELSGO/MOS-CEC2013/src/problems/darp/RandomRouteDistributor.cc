#include "RandomRouteDistributor.h"
#include <algorithm>

vector<int> RandomRouteDistributor::createRoutesAssignment(vector<GAGenome*>& orig_routes) {
  int num_routes = numRoutes(orig_routes);

  assert(num_routes >= num_nodes_);

  int routes_per_node = num_routes/num_nodes_;

  vector<int> routes_assignment;
  for (int i=0; i<num_nodes_; i++) {
    for (int j=0; j<routes_per_node; j++) {
      routes_assignment.push_back(i);
    }
  }
  int inserts2made = num_routes-routes_assignment.size();
  for (int i=0; i<inserts2made; i++) routes_assignment.push_back(GARandomInt(0,num_nodes_-1));

  assert(routes_assignment.size() == num_routes);

  random_shuffle(routes_assignment.begin(),routes_assignment.end());

  return routes_assignment;
}


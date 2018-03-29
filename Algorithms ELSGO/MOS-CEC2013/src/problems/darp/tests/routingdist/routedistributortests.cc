#include "../daropstests_aux.h"

// Hack for testing private members
#define protected public
#define private public
#include "../../../../gaeda/routingdist/RouteDistributor.cc"
#undef protected
#undef private

using namespace std;

/*
 * Note: Due to a refactoring process  assignRoutesToCachedRoutes is being tested on randomroutedistributortests
 */

struct RouteDistributorDummy : public RouteDistributor {
  RouteDistributorDummy(RoutingGenome& gen, int nnodes) : RouteDistributor(gen,nnodes) {}

  ~RouteDistributorDummy() {}

  RouteDistributor* clone() { throw runtime_error("It should not be called"); }

  virtual vector<int> createRoutesAssignment(vector<GAGenome*>& orig_routes){ throw runtime_error("It should not be called"); }

  string name() { return "RouteDistributorDummy"; }

};

class RouteDistributorTesting : public ::testing::Test {
public:
  VerticesList*    requestslist;
  DummyCostMatrix* dummymatrix;

  /*
   * Revisar si no se puede aligerar considerablemente
   */
  virtual void SetUp() {
    requestslist = new VerticesList();

    DummyDistMatrix* distmatrix = new DummyDistMatrix(MAXPOINTS);
    dummymatrix = new DummyCostMatrix(distmatrix,*requestslist);

    DARPTestsVars::initData(requestslist,dummymatrix);
  }

  virtual void TearDown() {
    delete dummymatrix;
    delete requestslist;
  }

};

TEST_F(RouteDistributorTesting,checkNumRoutes) {
  int num_routes = 3;
  int main_route = 0;

  int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );
  DARPGenome* gen        = createGen(num_routes,  main_route,gen_values);
  DARPGenome* other_gen  = createGen(num_routes*2,main_route,gen_values);
  DARPGenome* other_gen3 = createGen(num_routes*3,main_route,gen_values);

  vector<GAGenome*> all_gens; all_gens.push_back(gen); all_gens.push_back(other_gen); all_gens.push_back(other_gen3);

  RouteDistributor::numRoutes(all_gens);

  ASSERT_EQ(num_routes+num_routes*2+num_routes*3,RouteDistributor::numRoutes(all_gens));

  for(int i=0; i<all_gens.size(); i++) delete all_gens[i];
}

TEST_F(RouteDistributorTesting,checkClearCachedRoutes) {
  int num_routes = 3;
  int main_route = 0;

  int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );
  DARPGenome* gen    = createGen(num_routes,  main_route,gen_values);

  RouteDistributorDummy dumdist(*gen,5);

  for (int i=0; i<dumdist.cached_routes_.size(); i++) {
    ASSERT_NE(0,dumdist.cached_routes_[i]->numRoutes() );
  }

  dumdist.clearCachedRoutes();

  vector<GAGenome*> routes;
  for (int i=0; i<dumdist.cached_routes_.size(); i++) {
    ASSERT_EQ(0,dumdist.cached_routes_[i]->numRoutes() );
    routes.push_back( dumdist.cached_routes_[i] );
  }

  ASSERT_EQ(0,RouteDistributor::numRoutes(routes));

  delete gen;
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  GAGenome::optCriterion(GAGenome::MINIMIZATION);

  setUpTestingEnvironment();

  return RUN_ALL_TESTS();
}

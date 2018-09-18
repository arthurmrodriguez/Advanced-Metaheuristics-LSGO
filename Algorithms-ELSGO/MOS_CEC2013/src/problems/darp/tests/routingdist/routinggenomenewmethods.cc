#include "../daropstests_aux.h"

// Hack for testing private members

using namespace std;

class RoutingGenomeNewMethodsTesting : public ::testing::Test {
public:
  VerticesList*    requestslist;
  DummyCostMatrix* dummymatrix;

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

TEST_F(RoutingGenomeNewMethodsTesting,checkNumRoutes) {
  int num_routes = 3;
  int main_route = 0;

  int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );
  DARPGenome* gen    = createGen(num_routes,  main_route,gen_values);

  ASSERT_EQ(num_routes,gen->numRoutes());
  delete gen;

  num_routes = 1;
  gen = createGen(num_routes,  main_route,gen_values);
  ASSERT_EQ(num_routes,gen->numRoutes());
  delete gen;
}

TEST_F(RoutingGenomeNewMethodsTesting,checkEmptyRoutes) {
  int num_routes = 3;
  int main_route = 0;

  int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );
  DARPGenome* gen    = createGen(num_routes,  main_route,gen_values);

  gen->emptyRoutes();

  ASSERT_EQ(0,gen->numRoutes());
  delete gen;
}

TEST_F(RoutingGenomeNewMethodsTesting,checkaddRoute) {
  int num_routes = 3;
  int main_route = 0;

  int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );
  DARPGenome* gen    = createGen(num_routes,  main_route,gen_values);

  int other_routes = 2;
  int route2add    = 0;
  int a_othergen_values[] = {4,-4,5,6,-6,5}; vector<int> othergen_values (a_othergen_values, a_othergen_values + sizeof(a_othergen_values)/ sizeof(int) );

  DARPGenome* other_gen  = createGen(num_routes*2,main_route,othergen_values);

  // Add the first route
  gen->addRoute(*other_gen,route2add);
  // check sizes, number of routes and values
  ASSERT_EQ(num_routes+1,gen->size());
  ASSERT_EQ(other_gen->routeLength(route2add),gen->routeLength(gen->size()-1));
  for (int i=0; i<other_gen->routeLength(route2add); i++) {
    ASSERT_EQ(other_gen->gene(route2add,i), gen->gene(gen->size()-1,i));
  }

  // Add the last route
  route2add = other_gen->size()-1;
  gen->addRoute(*other_gen, route2add);
  ASSERT_EQ(num_routes+2,gen->size());
  // check sizes, number of routes and values
  ASSERT_EQ(other_gen->routeLength(route2add),gen->routeLength(gen->size()-1));
  for (int i=0; i<other_gen->routeLength(route2add); i++) {
    ASSERT_EQ(other_gen->gene(route2add,i), gen->gene(gen->size()-1,i));
  }

  // Try to add an inexistent route
  ASSERT_DEATH(gen->addRoute(*other_gen,100),"Assertion failed");

  delete gen;
  delete other_gen;
}

TEST_F(RoutingGenomeNewMethodsTesting,checkaddRoutes) {
  int num_routes = 3;
  int main_route = 0;

  int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );
  DARPGenome* gen    = createGen(num_routes,  main_route,gen_values);

  int other_routes = 2;
  int route2add    = 0;
  int a_othergen_values[] = {4,-4,5,6,-6,5}; vector<int> othergen_values (a_othergen_values, a_othergen_values + sizeof(a_othergen_values)/ sizeof(int) );

  DARPGenome* other_gen  = createGen(num_routes*2,main_route,othergen_values);

  int orig_size = gen->size();
  gen->addRoutes(*other_gen);

  ASSERT_EQ(orig_size+other_gen->size(),gen->size());

  for (int routepos=orig_size; routepos<gen->size(); routepos++) {
    for (int i=0; i<gen->routeLength(i); i++) {
      ASSERT_EQ(other_gen->gene(routepos-orig_size,i),gen->gene(routepos,i));
    }
  }

  delete gen;
  delete other_gen;
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  GAGenome::optCriterion(GAGenome::MINIMIZATION);

  setUpTestingEnvironment();

  return RUN_ALL_TESTS();
}

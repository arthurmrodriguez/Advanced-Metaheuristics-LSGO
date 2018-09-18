#include "../daropstests_aux.h"

// Hack for testing private members
#define protected public
#define private public
#include "../../../../gaeda/routingdist/RouteDistributor.h"
#include "../../RandomRouteDistributor.cc"
#undef protected
#undef private

using namespace std;

class RandomRouteDistributorTesting : public ::testing::Test {
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

  void checkRoutesAssignment(int num_nodes, int num_routes) {
    int gen_numroutes = 3; // default values, it is not being used actually

    int main_route = 0;
    int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );
    DARPGenome* gen    = createGen(gen_numroutes,  main_route,gen_values);

    RandomRouteDistributor a(*gen,num_nodes);

    // Both dummygen and dummy_routes are used just for testing createRoutesAssignment which only uses
    // the size of the vector
    DARPGenome* dummygen = createGen(1, 0,gen_values);
    vector<GAGenome*> dummy_routes;
    for (int i=0; i<num_routes; i++) dummy_routes.push_back(dummygen);
    vector<int> assignment  = a.createRoutesAssignment(dummy_routes);
    ASSERT_EQ(num_routes,assignment.size());
    delete dummygen;

    // We register the frequencies and we test that all nodes have at least a num_routes/num_node frequency
    vector<int> frequencies (num_nodes);
    for (int i=0; i<assignment.size(); i++) {
      ASSERT_LT(assignment[i],num_nodes);
      frequencies[ assignment[i] ] += 1;
    }
    for (int i=0; i<frequencies.size();i++) ASSERT_GE(frequencies[i],num_routes/num_nodes);

    delete gen;
  }

  virtual void TearDown() {
    delete dummymatrix;
    delete requestslist;
  }

  bool areSameRoutes(DARPGenome& gen1, int route1, DARPGenome& gen2, int route2) {
    if (gen1.routeLength(route1) != gen2.routeLength(route2) ) return false;

    bool res = true;

    for (int i=0; i<gen1.routeLength(route1); i++) {
      if (gen1.gene(route1,i) != gen2.gene(route2,i)) {
        res = false;
        break;
      }
    }

    return res;
  }

  bool isRouteInGenome(GAGenome& orig_g, int orig_routepos, GAGenome& g_tosearch ) {
    DARPGenome& orig_gen = dynamic_cast<DARPGenome&>(orig_g);
    assert(&orig_gen);
    DARPGenome& gen_tosearch = dynamic_cast<DARPGenome&>(g_tosearch);
    assert(&g_tosearch);


    bool res = false;
    for (int i=0; i<gen_tosearch.size(); i++) {
      if (areSameRoutes(orig_gen,orig_routepos,gen_tosearch,i)) {
        res = true;
        break;
      }
    }

    return res;
  }

};

TEST_F(RandomRouteDistributorTesting,checkRoutesAssignment) {
  checkRoutesAssignment(3, 3); // First test same number of routes as nodes
  checkRoutesAssignment(3, 6);
  checkRoutesAssignment(3, 5); // we test a greater number (not multiple of the number of routes)
  checkRoutesAssignment(5, 7); // we test a greater number (not multiple of the number of routes)
  checkRoutesAssignment(5, 8); // we test a greater number (not multiple of the number of routes)

#ifdef DEBUG
  ASSERT_DEATH(checkRoutesAssignment(5, 3),"Assertion failed");
#endif
}

TEST_F(RandomRouteDistributorTesting,checkAssignRoutesToCachedRoutesSimpleCase) {
  //First we create the temporal gen that is going to be used as base for the distributor
  // with a random number of routes (10 in this case) that is recommended to be different to the posterior
  // number of routes used
  int gen_numroutes = 10;
  int main_route = 0;
  vector<int> gen_values; for (int i=1; i<4; i++) {gen_values.push_back(i); gen_values.push_back(-i);}
  DARPGenome* tmpgen = createGen(gen_numroutes,  main_route,gen_values);
  // Then we create the distributor
  int num_nodes = 3;
  RandomRouteDistributor distrib(*tmpgen,num_nodes);

  //We create four genomes with five routes
  gen_values.clear(); for (int i=1; i<4; i++) {gen_values.push_back(i); gen_values.push_back(-i);}
  gen_numroutes = 1;
  DARPGenome* routetmpgen = createGen(gen_numroutes,  main_route,gen_values);

  gen_values.clear(); for (int i=5; i<10; i++) {gen_values.push_back(i); gen_values.push_back(-i);}
  gen_numroutes = 1;
  DARPGenome* routetmpgen2 = createGen(gen_numroutes,  main_route,gen_values);

  gen_values.clear(); for (int i=10; i<15; i++) {gen_values.push_back(i); gen_values.push_back(-i);}
  gen_numroutes = 1;
  DARPGenome* routetmpgen3 = createGen(gen_numroutes,  main_route,gen_values);

  gen_values.clear(); for (int i=15; i<20; i++) {gen_values.push_back(i); gen_values.push_back(-i);}
  gen_numroutes = 1;
  DARPGenome* routetmpgen4 = createGen(gen_numroutes,  main_route,gen_values);

  gen_values.clear(); for (int i=20; i<23; i++) {gen_values.push_back(i); gen_values.push_back(-i);}
  gen_numroutes = 1;
  DARPGenome* routetmpgen5 = createGen(gen_numroutes,  main_route,gen_values);


  // to be saved in routes (the vector with the genomes)
  vector<GAGenome*> routes;
  routes.push_back(routetmpgen);routes.push_back(routetmpgen2);routes.push_back(routetmpgen3);routes.push_back(routetmpgen4); routes.push_back(routetmpgen5);

  // then we create a simple assignment
  vector<int> routes_assignment;
  for (int i=0; i<routes.size(); i++) routes_assignment.push_back(i%num_nodes);

  distrib.assignRoutesToCachedRoutes(routes,routes_assignment);

  // finally we check that every route is suppose to be where it actually is
  for (int i=0; i<routes.size(); i++) {
    ASSERT_TRUE( isRouteInGenome(* routes[i], 0, *(distrib.cached_routes_[ routes_assignment[i] ]) ) );
  }

  delete tmpgen;
  delete routetmpgen;
  delete routetmpgen2;
  delete routetmpgen3;
  delete routetmpgen4;
}

TEST_F(RandomRouteDistributorTesting,checkAssignRoutesToCachedRoutesComplexCase) {
  //First we create the temporal gen that is going to be used as base for the distributor
  // with a random number of routes (10 in this case) that is recommended to be different to the posterior
  // number of routes used
  int gen_numroutes = 10;
  int main_route = 0;
  vector<int> gen_values; for (int i=1; i<4; i++) {gen_values.push_back(i); gen_values.push_back(-i);}
  DARPGenome* tmpgen = createGen(gen_numroutes,  main_route,gen_values);
  // Then we create the distributor
  int num_nodes = 3;
  RandomRouteDistributor distrib(*tmpgen,num_nodes);

  //We create five genomes with a different number of routes per genome
  gen_values.clear(); for (int i=1; i<5; i++) {gen_values.push_back(i); gen_values.push_back(-i);}
  gen_numroutes = 2;
  DARPGenome* routetmpgen = createGen(gen_numroutes,  main_route,gen_values);

  gen_values.clear(); for (int i=5; i<10; i++) {gen_values.push_back(i); gen_values.push_back(-i);}
  gen_numroutes = 1;
  DARPGenome* routetmpgen2 = createGen(gen_numroutes,  main_route,gen_values);

  gen_values.clear(); for (int i=15; i<20; i++) {gen_values.push_back(i); gen_values.push_back(-i);}
  gen_numroutes = 3;
  DARPGenome* routetmpgen3 = createGen(gen_numroutes,  main_route,gen_values);

  gen_values.clear(); for (int i=21; i<25; i++) {gen_values.push_back(i); gen_values.push_back(-i);}
  gen_numroutes = 2;
  DARPGenome* routetmpgen4 = createGen(gen_numroutes,  main_route,gen_values);

  gen_values.clear(); for (int i=26; i<30; i++) {gen_values.push_back(i); gen_values.push_back(-i);}
  gen_numroutes = 3;
  DARPGenome* routetmpgen5 = createGen(gen_numroutes,  main_route,gen_values);


  // to be saved in routes (the vector with the genomes)
  vector<GAGenome*> routes;
  routes.push_back(routetmpgen);routes.push_back(routetmpgen2);routes.push_back(routetmpgen3);routes.push_back(routetmpgen4); routes.push_back(routetmpgen5);

  int num_routes = 0;
  for (int i=0; i<routes.size(); i++) num_routes += dynamic_cast<DARPGenome*>(routes[i])->size();

  // then we create a simple assignment
  vector<int> routes_assignment;
  for (int i=0; i<num_routes; i++) routes_assignment.push_back(i%num_nodes);

  distrib.assignRoutesToCachedRoutes(routes,routes_assignment);


  // finally we check that every route is suppose to be where it actually is
  for (int i=0, routenum=0; i<routes.size(); i++) {
    DARPGenome& ind = dynamic_cast<DARPGenome&>(*routes[i]);
    for (int route=0; route<ind.size(); route++, routenum++) {
      ASSERT_TRUE( isRouteInGenome(ind, route, *(distrib.cached_routes_[ routes_assignment[routenum] ]) ) );
    }
  }

  delete tmpgen;
  delete routetmpgen;
  delete routetmpgen2;
  delete routetmpgen3;
  delete routetmpgen4;
}


int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  GAGenome::optCriterion(GAGenome::MINIMIZATION);

  setUpTestingEnvironment();

  return RUN_ALL_TESTS();
}

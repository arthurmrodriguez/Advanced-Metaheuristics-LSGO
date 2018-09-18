#include "../daropstests_aux.h"

// Hack for testing private members
#define protected public
#define private public
#include "../../ReqDistMatrix.cc"
#include "../../RouteDistMatrix.cc"

#undef protected
#undef private

#include <time.h>

using namespace std;

// Used for avoiding the call to the objective function for computing the score (easier to conduct the tests)
class ReqDistMatrixDummy : public ReqDistMatrix {
public:
  ReqDistMatrixDummy(CostMatrix& costmatrix, DARPEvaluator& eval) : ReqDistMatrix(costmatrix,eval,true,false,0,0){}

  double bestScoreFromAllCombs (Request& req1, Request& req2) {
    return req1.pickup_vert->id_ + req2.pickup_vert->id_;
  }

};

class RouteDistMatrixTesting : public ::testing::Test {
public:
  VerticesList*    requestslist;
  DummyCostMatrix* dummymatrix;
  DummyDistMatrix* dummydistmatrix;
  ReqDistMatrix*   reqdistmatrixdummy;
  vector<Request>  reqs;

  virtual void SetUp() {
    requestslist = new VerticesList();

    dummydistmatrix = new DummyDistMatrix(MAXPOINTS);

    for (int i=0; i<MAXPOINTS; i++) {
      for (int j=0; j<MAXPOINTS; j++) {
        dummydistmatrix->setCost( i,j,GAGenome::worstPossibleScore()/MAXPOINTS);
        dummydistmatrix->setCost( j,i,GAGenome::worstPossibleScore()/MAXPOINTS);
      }
    }

    int gen_numroutes = 3; // default values, it is not being used actually

    int main_route = 0;
    int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );
    DARPGenome* gen    = createGen(gen_numroutes,  main_route,gen_values);

    dummymatrix = new DummyCostMatrix(dummydistmatrix,*requestslist);

    DARPTestsVars::initData(requestslist, dummymatrix);

    assert(DARPTestsVars::darpeval);
    reqdistmatrixdummy = new ReqDistMatrixDummy(*dummymatrix,*DARPTestsVars::darpeval);


    int num_reqs = 20;

    int id,pos,load;
    bool crit;
    long fb,fe;
    Vertex::VertexType t;

    for (int i=1; i<num_reqs; i++) {
      id = i; pos = 1; crit = true; load = 1; fb = 5; fe = 10; // Only the id is meaningful
      Vertex* pv = new Vertex(id,pos,crit,load,fb,fe,Vertex::PICKUP);
      Vertex* dv = new Vertex(Vertex::getSiblingVertId(id),pos,crit,load,fb,fe,Vertex::DELIVERY);
      Request r(pv,dv);
      reqs.push_back(r);
    }

    delete gen;
  }

  virtual void TearDown() {
    if (dummymatrix)        delete dummymatrix;
    if (dummydistmatrix)    delete dummydistmatrix;
    if (requestslist)       delete requestslist;
    if (reqdistmatrixdummy) delete reqdistmatrixdummy;

    for (int i=0; i<reqs.size(); i++) {
      delete reqs[i].pickup_vert;
      delete reqs[i].delivery_vert;
    }
  }

};


// Basic test where we check the configuration of a simple genome with two routes, one with two requests and
// a second one with only one request
TEST_F(RouteDistMatrixTesting,oneGenomeThreeRoutes) {

  reqdistmatrixdummy->setUpDistMatrix(reqs);
  reqdistmatrixdummy->initializeDistMatrix(reqs);

  int num_routes = 2;
  DARPGenome gen1 (DARPGenome(num_routes,0,DummyObjFunc) );

  gen1.pushBackVertex(0,reqs[2].pickup_vert->id_);
  gen1.pushBackVertex(0,reqs[2].delivery_vert->id_);


  gen1.pushBackVertex(1,reqs[0].pickup_vert->id_);
  gen1.pushBackVertex(1,reqs[1].pickup_vert->id_);
  gen1.pushBackVertex(1,reqs[0].delivery_vert->id_);
  gen1.pushBackVertex(1,reqs[1].delivery_vert->id_);


  vector<GAGenome*> routes;
  routes.push_back(&gen1);

  RouteDistMatrix routedistmatrix(routes, *reqdistmatrixdummy);

  ASSERT_EQ( GAMin( reqdistmatrixdummy->getDistOfReqsPickupIds(gen1.gene(1,0),
                                                               gen1.gene(0,0) ),
                    reqdistmatrixdummy->getDistOfReqsPickupIds(gen1.gene(1,1),
                                                               gen1.gene(0,0) )
                  ),
              routedistmatrix.getDist(1,0)
           );

  ASSERT_NE( GAMin( reqdistmatrixdummy->getDistOfReqsPickupIds(gen1.gene(1,0),
                                                               gen1.gene(1,1) ),
                    reqdistmatrixdummy->getDistOfReqsPickupIds(gen1.gene(1,1),
                                                               gen1.gene(1,0) )
                  ),
             routedistmatrix.getDist(1,0)
           );
}

// Here we check two genomes one with two routes and one with one
TEST_F(RouteDistMatrixTesting,twoGenomesThreeRoutes) {

  reqdistmatrixdummy->setUpDistMatrix(reqs);
  reqdistmatrixdummy->initializeDistMatrix(reqs);

  int num_routes = 2;
  DARPGenome gen1 (DARPGenome(num_routes,0,DummyObjFunc) );
  gen1.pushBackVertex(0,reqs[2].pickup_vert->id_);
  gen1.pushBackVertex(0,reqs[2].delivery_vert->id_);

  gen1.pushBackVertex(1,reqs[0].pickup_vert->id_);
  gen1.pushBackVertex(1,reqs[0].delivery_vert->id_);

  num_routes = 1;
  DARPGenome gen2 (DARPGenome(num_routes,0,DummyObjFunc) );
  gen2.pushBackVertex(0,reqs[1].pickup_vert->id_);
  gen2.pushBackVertex(0,reqs[1].delivery_vert->id_);


  vector<GAGenome*> routes;
  routes.push_back(&gen1);
  routes.push_back(&gen2);

  RouteDistMatrix routedistmatrix(routes, *reqdistmatrixdummy);

  ASSERT_EQ( reqdistmatrixdummy->getDistOfReqsPickupIds(gen1.gene(0,0), gen1.gene(1,0) ) , routedistmatrix.getDist(1,0) );
  ASSERT_EQ( reqdistmatrixdummy->getDistOfReqsPickupIds(gen1.gene(0,0), gen2.gene(0,0) ) , routedistmatrix.getDist(2,0) );
  ASSERT_EQ( reqdistmatrixdummy->getDistOfReqsPickupIds(gen1.gene(1,0), gen2.gene(0,0) ) , routedistmatrix.getDist(2,1) );
}

// Test where we check three genomes, the first one with two routes , the second one with one route and the third one
// with three routes. Since each route have several vertices, in this test we make use of the fact that the dummy
// objective function is req1.pickup_vert->id_ + req2.pickup_vert->id_;
TEST_F(RouteDistMatrixTesting,threeGenomesFiveRoutes) {

  reqdistmatrixdummy->setUpDistMatrix(reqs);
  reqdistmatrixdummy->initializeDistMatrix(reqs);

  int num_routes = 2;
  DARPGenome gen1 (DARPGenome(num_routes,0,DummyObjFunc) );
  gen1.pushBackVertex(0,reqs[1].pickup_vert->id_);
  gen1.pushBackVertex(0,reqs[2].pickup_vert->id_);
  gen1.pushBackVertex(0,reqs[0].pickup_vert->id_);
  gen1.pushBackVertex(0,reqs[2].delivery_vert->id_);
  gen1.pushBackVertex(0,reqs[1].delivery_vert->id_);
  gen1.pushBackVertex(0,reqs[0].delivery_vert->id_);

  gen1.pushBackVertex(1,reqs[3].pickup_vert->id_);
  gen1.pushBackVertex(1,reqs[3].delivery_vert->id_);


  num_routes = 1;
  DARPGenome gen2 (DARPGenome(num_routes,0,DummyObjFunc) );
  gen2.pushBackVertex(0,reqs[5].pickup_vert->id_);
  gen2.pushBackVertex(0,reqs[5].delivery_vert->id_);
  gen2.pushBackVertex(0,reqs[6].pickup_vert->id_);
  gen2.pushBackVertex(0,reqs[6].delivery_vert->id_);
  gen2.pushBackVertex(0,reqs[4].pickup_vert->id_);
  gen2.pushBackVertex(0,reqs[4].delivery_vert->id_);


  num_routes = 3;
  DARPGenome gen3 (DARPGenome(num_routes,0,DummyObjFunc) );
  gen3.pushBackVertex(0,reqs[7].pickup_vert->id_);
  gen3.pushBackVertex(0,reqs[7].delivery_vert->id_);

  gen3.pushBackVertex(1,reqs[9].pickup_vert->id_);
  gen3.pushBackVertex(1,reqs[8].pickup_vert->id_);
  gen3.pushBackVertex(1,reqs[8].delivery_vert->id_);
  gen3.pushBackVertex(1,reqs[9].delivery_vert->id_);

  gen3.pushBackVertex(2,reqs[11].pickup_vert->id_);
  gen3.pushBackVertex(2,reqs[12].pickup_vert->id_);
  gen3.pushBackVertex(2,reqs[12].delivery_vert->id_);
  gen3.pushBackVertex(2,reqs[11].delivery_vert->id_);
  gen3.pushBackVertex(2,reqs[10].pickup_vert->id_);
  gen3.pushBackVertex(2,reqs[10].delivery_vert->id_);



  vector<GAGenome*> routes;
  routes.push_back(&gen1);
  routes.push_back(&gen2);
  routes.push_back(&gen3);

  RouteDistMatrix routedistmatrix(routes, *reqdistmatrixdummy);


  // Note that the minimum distance between two requests corresponds to the sum of their pickup ids (according to the
  // implemented dummy function). Based on this fact, we first compute the minimum pickup id per route to reuse
  // this values for computing the distance

  vector<double> routes_min_pick_id;
  routes_min_pick_id.push_back(gen1.gene(0,2));
  routes_min_pick_id.push_back(gen1.gene(1,0));
  routes_min_pick_id.push_back(gen2.gene(0,4));
  routes_min_pick_id.push_back(gen3.gene(0,0));
  routes_min_pick_id.push_back(gen3.gene(1,1));
  routes_min_pick_id.push_back(gen3.gene(2,4));

  for (int i=1; i<routes_min_pick_id.size(); i++) {
    for (int j=0; j<i; j++) {
      ASSERT_EQ( reqdistmatrixdummy->getDistOfReqsPickupIds(routes_min_pick_id[j], routes_min_pick_id[i] ), routedistmatrix.getDist(i,j) );
    }
  }
}


int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  GAGenome::optCriterion(GAGenome::MINIMIZATION);

  setUpTestingEnvironment();

  return RUN_ALL_TESTS();
}

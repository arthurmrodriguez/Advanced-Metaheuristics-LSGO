#include "../daropstests_aux.h"

// Hack for testing private members
#define protected public
#define private public
#include "../../ReqDistMatrix.cc"
#include "../../KMedoidsRouteDistributor.cc"
#include "../../RouteDistMatrix.cc"

#undef protected
#undef private

#include <time.h>

using namespace std;

// Used for avoiding the call to the objective function for computing the score (easier to conduct the tests)
class ReqDistMatrixDummy : public ReqDistMatrix {
public:
  ReqDistMatrixDummy(CostMatrix& costmatrix, DARPEvaluator& eval) : ReqDistMatrix(costmatrix,eval,true,false,0,0){}

  //dummy function where it returns 0 if one request is a multiple of the other or 100 otherwise
  long double bestScoreFromAllCombs (Request& req1, Request& req2) {
    // this assert makes use of how the distance matrix is constructed so that this method should only
    // when constructing the matrix, be called with req1.id >= req2.id
    assert(req1.pickup_vert->id_ >= req2.pickup_vert->id_);

    if (req1.pickup_vert->id_ % req2.pickup_vert->id_ == 0) {
//      cout << "req1:" << req1.pickup_vert->id_ << " req2:" << req2.pickup_vert->id_ << " req1%req2= " << req1.pickup_vert->id_ % req2.pickup_vert->id_ << endl;
      return 0;
    }
    else return 100;
  }

};

class RouteDistMatrixTesting : public ::testing::Test {
public:
  VerticesList*    requestslist;
  DummyCostMatrix* dummymatrix;
  DummyDistMatrix* dummydistmatrix;
  ReqDistMatrix*   reqdistmatrixdummy;
  DARPGenome*      gen;
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
    gen = createGen(gen_numroutes,  main_route,gen_values);

    dummymatrix = new DummyCostMatrix(dummydistmatrix,*requestslist);

    DARPTestsVars::initData(requestslist, dummymatrix);

    reqdistmatrixdummy = new ReqDistMatrixDummy(*dummymatrix,*DARPTestsVars::darpeval);


    int num_reqs = 300;

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
  }

  virtual void TearDown() {
    if (dummymatrix)        delete dummymatrix;
    if (dummydistmatrix)    delete dummydistmatrix;
    if (requestslist)       delete requestslist;
    if (reqdistmatrixdummy) delete reqdistmatrixdummy;
    if (gen)                delete gen;

    for (int i=0; i<reqs.size(); i++) {
      delete reqs[i].pickup_vert;
      delete reqs[i].delivery_vert;
    }
  }

};



// Test where we have six routes in three genomes: the first one with two routes , the second one with one route and the third one
// with three routes. The ids have been selected so that some requests are multiple of the others. The idea
// is that the kmedoids algorithm returns the logical clusters.
TEST_F(RouteDistMatrixTesting,threeGenomesFiveRoutes) {

  reqdistmatrixdummy->setUpDistMatrix(reqs);
  reqdistmatrixdummy->initializeDistMatrix(reqs);

  int num_routes = 2;
  DARPGenome gen1 (DARPGenome(num_routes,0,DummyObjFunc) );
  gen1.pushBackVertex(0,2);
  gen1.pushBackVertex(0,3);
  gen1.pushBackVertex(0,5);
  gen1.pushBackVertex(0,-5);
  gen1.pushBackVertex(0,-2);
  gen1.pushBackVertex(0,-3);

  gen1.pushBackVertex(1,7);
  gen1.pushBackVertex(1,-7);


  num_routes = 1;
  DARPGenome gen2 (DARPGenome(num_routes,0,DummyObjFunc) );
  gen2.pushBackVertex(0,4);
  gen2.pushBackVertex(0,-4);
  gen2.pushBackVertex(0,13);
  gen2.pushBackVertex(0,-13);
  gen2.pushBackVertex(0,17);
  gen2.pushBackVertex(0,-17);


  num_routes = 3;
  DARPGenome gen3 (DARPGenome(num_routes,0,DummyObjFunc) );
  gen3.pushBackVertex(0,259); // id is 7*37
  gen3.pushBackVertex(0,-259);

  gen3.pushBackVertex(1,52); // id is 4*13
  gen3.pushBackVertex(1,19);
  gen3.pushBackVertex(1,-52);
  gen3.pushBackVertex(1,-19);

  gen3.pushBackVertex(2,23);
  gen3.pushBackVertex(2,29);
  gen3.pushBackVertex(2,-23);
  gen3.pushBackVertex(2,-29);
  gen3.pushBackVertex(2,31);
  gen3.pushBackVertex(2,-31);


  vector<GAGenome*> routes;
  routes.push_back(&gen1);
  routes.push_back(&gen2);
  routes.push_back(&gen3);

  int          npass  = 1000;
  kmedoidsMode mode   = KMEDNORMAL;
  int          nnodes = 3; // We have labeled the ids so that three groups of two could be created
  KMedoidsRouteDistributor kmeddist(reqdistmatrixdummy, npass, mode, *gen, nnodes);

  vector<int> res = kmeddist.createRoutesAssignment(routes);

  ASSERT_EQ(res[0],res[2]);
  ASSERT_EQ(res[0],res[4]);
  ASSERT_EQ(res[1],res[3]);
  ASSERT_NE(res[0],res[5]);
  ASSERT_NE(res[1],res[5]);
}


int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  GAGenome::optCriterion(GAGenome::MINIMIZATION);

  setUpTestingEnvironment();

  return RUN_ALL_TESTS();
}

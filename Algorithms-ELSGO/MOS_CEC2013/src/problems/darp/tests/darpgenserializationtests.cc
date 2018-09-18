#include "daropstests_aux.h"
#include <iostream>
#include <sstream>

class SerializationTesting : public ::testing::Test {
public:
  VerticesList*     requestslist;
  DummyCostMatrix*  dummymatrix;
  int               num_routes;

    virtual void SetUp() {
      requestslist = new VerticesList();

      num_routes  = 3;

      vector<DummyCostPoint> costinf;
      addRequestRoute(*requestslist, 1,1,2,5,5,  costinf, 15);
      addRequestRoute(*requestslist, 2,3,4,26,35, costinf,  9);
      addRequestRoute(*requestslist, 3,5,6,45,55, costinf,  5);

      addRequestRoute(*requestslist, 10,1,4,5,6,  costinf, 15);
      addRequestRoute(*requestslist, 20,2,5,26,36, costinf,  9);
      addRequestRoute(*requestslist, 30,6,1,45,56, costinf,  5);

      addRequestRoute(*requestslist, 100,3,1,5,7,  costinf, 15);
      addRequestRoute(*requestslist, 200,2,4,26,37, costinf,  9);
      addRequestRoute(*requestslist, 300,1,5,45,57, costinf,  5);

      DummyDistMatrix* distmatrix = new DummyDistMatrix(MAXPOINTS,costinf);
      dummymatrix = new DummyCostMatrix(distmatrix,*requestslist);
      DARPTestsVars::initData(requestslist,dummymatrix);

    }

    virtual void TearDown() {
      delete dummymatrix;
      delete requestslist;
    }

    bool areEqual(DARPGenome& gen1, DARPGenome& gen2) {
#define CHECKSINGLEATTRIB(name)\
  if (gen1.name != gen2.name) return false;

  CHECKSINGLEATTRIB(_score);
  CHECKSINGLEATTRIB(age_);
  CHECKSINGLEATTRIB(gen_of_birth_);
  CHECKSINGLEATTRIB(origin_);
  CHECKSINGLEATTRIB(id);

  CHECKSINGLEATTRIB(island);
  CHECKSINGLEATTRIB(myage);

  CHECKSINGLEATTRIB(feasible_);
  CHECKSINGLEATTRIB(nx);
  CHECKSINGLEATTRIB(minX);
  CHECKSINGLEATTRIB(maxX);

  for (int routepos=0; routepos<gen1.size(); routepos++) {
    int routesize = gen1.routeLength(routepos);
    if (routesize != gen2.routeLength(routepos)) return false;

    for (int i=0; i<routesize; i++) {
      if (gen1.gene(routepos,i) != gen2.gene(routepos,i)) return false;
    }
  }

#define CHECKATTRIBVECTOR(name) \
  if ( gen1.name.size() != gen2.name.size() ) return false;\
  for (int i=0; i<gen1.name.size(); i++) {\
    if ( !(gen1.name[i] == gen2.name[i]) ) return false;\
  }\

  CHECKATTRIBVECTOR(modifiedRoutes_);

  CHECKATTRIBVECTOR(costPerRoute_);
  CHECKATTRIBVECTOR(loadVPerRoute_);
  CHECKATTRIBVECTOR(TWVPerRoute_);
  CHECKATTRIBVECTOR(rideVPerRoute_);

  CHECKATTRIBVECTOR(pickupDelayPerRoute_);
  CHECKATTRIBVECTOR(deliveryDelayPerRoute_);

  CHECKATTRIBVECTOR(localSearchedRoutes_);
  CHECKATTRIBVECTOR(feasibleRoutes_);

  CHECKSINGLEATTRIB(evalWeightsUpdateValue_);

  return true;
  }
};

TEST_F(SerializationTesting,serializationtest) {
  // Simple test, we set arbitrarily the values of the genome, we serialize it and the we check that the values
  // are still the same
  num_routes = 3;
  DARPGenome* gen = new DARPGenome(num_routes,0,DummyObjFunc);

  int route = 2;
  int seq[] = {1,-1,2,-2,3,-3};
  vector<int> seqv (seq, seq + sizeof(seq)/ sizeof(int) );
  for (int i=0; i<seqv.size(); i++) {
    gen->pushBackVertex(route  ,seqv[i]);
    gen->pushBackVertex(route-1,10*seqv[i]);
    gen->pushBackVertex(route-2,100*seqv[i]);
  }

  gen->score(); // to recompute some values

  //We set the values manually
  for (int i=0; i<gen->modifiedRoutes_.size(); i++)      gen->modifiedRoutes_[i] =      i % 2 == 0;
  for (int i=0; i<gen->localSearchedRoutes_.size(); i++) gen->localSearchedRoutes_[i] = (DARPGenome::lsstatus) (i % 3);
  for (int i=0; i<gen->feasibleRoutes_.size(); i++)      gen->feasibleRoutes_[i] =      i % 4 == 0;

  int num=0;
  for (int i=0; i<gen->costPerRoute_.size(); i++)          gen->costPerRoute_[i] = num;
  for (int i=0; i<gen->loadVPerRoute_.size(); i++)         gen->loadVPerRoute_[i] = num;
  for (int i=0; i<gen->TWVPerRoute_.size(); i++)           gen->TWVPerRoute_[i] = num;
  for (int i=0; i<gen->rideVPerRoute_.size(); i++)         gen->rideVPerRoute_[i] = num;
  for (int i=0; i<gen->pickupDelayPerRoute_.size(); i++)   gen->pickupDelayPerRoute_[i] = num;
  for (int i=0; i<gen->deliveryDelayPerRoute_.size(); i++) gen->deliveryDelayPerRoute_[i] = num;

  gen->evalWeightsUpdateValue_ = 2.123234324; // Random values


  DARPGenome* othergen = new DARPGenome(num_routes+1,0,DummyObjFunc);

  ASSERT_FALSE(areEqual(*gen,*othergen));

  ostringstream os(ios::binary);
  gen->writeObject(os);
  const string os_str = os.str();
  istringstream is(os_str, ios::binary);

  othergen->readObject(is);

  ASSERT_TRUE(areEqual(*gen,*othergen));

  gen->gene(0,0,0);

  ASSERT_FALSE(areEqual(*gen,*othergen));

  delete gen;
  delete othergen;
}


int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  GAGenome::optCriterion(GAGenome::MINIMIZATION);

  setUpTestingEnvironment();

  return RUN_ALL_TESTS();
}


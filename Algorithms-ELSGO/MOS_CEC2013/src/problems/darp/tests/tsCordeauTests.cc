
#include "daropstests_aux.h"

#define protected public
#define private public
#include "../TSDARPGenome.cc"
#undef protected
#undef private


class TSDARPGenomeTest : public ::testing::TestWithParam<int*> {
protected:
  static const int MAXPOINTS = 300; 

public:
  VerticesList*    requestslist;
  DummyCostMatrix* dummymatrix;
  int              num_routes;
  int              route     ;

    virtual void SetUp() {
      num_routes  = 3;
      route       = 0;
      int max_id  = 9;

      requestslist = new VerticesList();

      DummyDistMatrix* distmatrix = new DummyDistMatrix(MAXPOINTS);

      // We insert the possible positions. The *10 refer to the destination positions
      // Therefore the request with id 2 would have two vertices: the pickup vertec
      // with the position number id and the delivery vertex with the position number i*10
      vector<int> allposv;
      for (int i=1; i<=max_id; i++) {
        allposv.push_back(i);
        allposv.push_back(i*10);
      }

      // The costs are set so that the optimum sequence is 1 -1 2 -2 3 -3
      // The remaining combinations are set to a higher cost
      for (vector<int>::iterator it=allposv.begin(); it!=allposv.end(); it++) {
        for (vector<int>::iterator itj=it+1; itj!=allposv.end(); itj++) {
          distmatrix->setCost( *it,*itj,3);
          distmatrix->setCost( *itj,*it,3);
        }
      }

      for (int i=0; i<allposv.size()-1; i++) distmatrix->setCost(allposv[i],allposv[i+1],1);

      // wend = wbegin + MAXUSERWAITTIME (20)
      // delivery vertex wbegin = wbegin + cost (origin dest), wend = wend + MAXUSERRIDETIME (60)

      for (int i=1; i<=max_id; i++) {
        for (int j=0;j<num_routes;j++) {
                                    //id           orig_pos  dest_pos  starting time
          addRequest(*requestslist, (int) i*pow(10.0,j),    i ,      i*10,   5+(i-1)*15,  *distmatrix);
        }
      }

      dummymatrix = new DummyCostMatrix(distmatrix,*requestslist);
      DARPTestsVars::initData(requestslist,dummymatrix);

    }

    virtual void TearDown() {
      delete dummymatrix;
      delete requestslist;
    }

    TSDARPGenome* createGen(int num_routes,int main_route, vector<int>& main_route_values) {
      assert(main_route >=0 && main_route < num_routes);
      TSDARPGenome* gen = new TSDARPGenome(num_routes,DummyObjFunc);

      vector< vector<int> > values(num_routes);
      for (int i=0; i<main_route_values.size(); i++) {
        values[main_route].push_back(main_route_values[i]);

        for (int j=1; j<num_routes; j++) {
          int next_route = (main_route+j) % num_routes;
          values[next_route].push_back( (int) pow((double)10,j) * main_route_values[i] );
        }

      }

      return createGen(values);
    }


    TSDARPGenome* createGen(vector< vector<int> >& values) {
      int num_routes = values.size();

      TSDARPGenome* gen = new TSDARPGenome(num_routes,DummyObjFunc);
      for (int route=0; route<values.size(); route++) {
        vector<int>& route_values = values[route];
        for (int pos=0; pos<route_values.size(); pos++) {
          gen->pushBackVertex(route,route_values[pos]);
        }
      }

      return gen;
    }

    void checkTabuMemUpdate(TSDARPGenome& gen) {
      map< pair<int,int> , int > oldTabuMem = gen.tabuMem_;
      gen.ageMemory();

      for ( map< pair<int,int> , int >::iterator it=oldTabuMem.begin(); it!=oldTabuMem.end(); it++) {
        pair<int,int> key = it->first;
        int oldTabuValue  = it->second;

        bool found = gen.tabuMem_.find(key) != gen.tabuMem_.end();

        if (oldTabuValue == 1) ASSERT_EQ(false,found);
        else                   ASSERT_EQ(gen.tabuMem_[key],oldTabuValue-1);
      }
    }

    void checkMoveToNeighbor(vector<int> gen_values, int orig_route, int dest_route, int vertid) {
      int main_route  = 0;
      int num_routes  = 3;

      TSDARPGenome* gen  = createGen(num_routes,main_route, gen_values);

      int    orig_route_oldsize = gen->routeLength(orig_route);
      int    dest_route_oldsize = gen->routeLength(dest_route);
      double old_score          = gen->score();
      int    orig_freq          = gen->opFreq_[TSDARPGenome::routeVertIdKey(dest_route,vertid)];
      int    verticesNum_old    = 0; for (int i=0; i<gen->length(); i++) verticesNum_old += gen->routeLength(i);

      gen->moveToNeighbor(orig_route,dest_route,vertid);

      int verticesNum = 0; for (int i=0; i<gen->length(); i++) verticesNum += gen->routeLength(i);

      ASSERT_EQ(verticesNum_old,verticesNum);

      ASSERT_EQ(orig_route_oldsize-2,gen->routeLength(orig_route));
      ASSERT_EQ(dest_route_oldsize+2,gen->routeLength(dest_route));

      ASSERT_EQ(-1,gen->findPosOfVertex(orig_route,vertid));
      ASSERT_EQ(-1,gen->findPosOfVertex(orig_route,Vertex::getSiblingVertId(vertid)));
      ASSERT_NE(-1,gen->findPosOfVertex(dest_route,vertid));
      ASSERT_NE(-1,gen->findPosOfVertex(dest_route,Vertex::getSiblingVertId(vertid)));

      ASSERT_NE(old_score,gen->score()); // We are going to try moves that change the score

      pair<int,int> key = TSDARPGenome::routeVertIdKey(dest_route,vertid);

      ASSERT_EQ(0,gen->divPenalization(dest_route,vertid));

      delete gen;
    }

    void checkBestFromNeighborhood(TSDARPGenome& gen, int vertid, int origroute, int bestdestroute) {
      ASSERT_NE(-1,gen.findPosOfVertex(origroute,vertid));

      pair<RoutingGenome*,int> neighborData = gen.bestFromNeighborhood();
      TSDARPGenome*            neighbor     = dynamic_cast<TSDARPGenome*>(neighborData.first);
      int                      used_evals   = neighborData.second;

      ASSERT_TRUE(neighbor != 0);
      ASSERT_GT(used_evals,0);

      ASSERT_EQ(gen.evalWeightsUpdateValue_,neighbor->evalWeightsUpdateValue_);
      ASSERT_EQ(gen.lackOfDivPenalization_,neighbor->lackOfDivPenalization_);
      ASSERT_EQ(gen.maxTabuLife_,neighbor->maxTabuLife_);

      ASSERT_NE(-1,neighbor->findPosOfVertex(bestdestroute,vertid));
      ASSERT_EQ(-1,neighbor->findPosOfVertex(origroute,vertid));

      pair<int,int> key = TSDARPGenome::routeVertIdKey(bestdestroute,vertid);

      ASSERT_EQ(gen.opFreq_[key]+1,neighbor->opFreq_[key]);
      ASSERT_EQ(false,neighbor->allowedMove(origroute,bestdestroute,vertid));
      ASSERT_EQ(round(gen.maxTabuLife_),neighbor->tabuMem_[key]);

      delete neighbor;
    }
};

TEST_F(TSDARPGenomeTest,checkInitialLifeValue) {
  int gen_values_a[] = {1,-1,-7,7,-8,8,2,-2,3,-3};
  vector<int> gen_values (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int main_route  = 0;
  int num_routes  = 3;

  TSDARPGenome* gen = createGen(num_routes,main_route, gen_values);
  pair<int,int> key = TSDARPGenome::routeVertIdKey(0,1);

  gen->addTabuMovement(key);

  ASSERT_NE(0,gen->tabuMem_.size());
  ASSERT_EQ(round(gen->maxTabuLife_),gen->tabuMem_[key]);

  delete gen;
}

TEST_F(TSDARPGenomeTest,updateTabuMemNoValues) {
  int gen_values_a[] = {1,-1,-7,7,-8,8,2,-2,3,-3};
  vector<int> gen_values (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int main_route  = 1;
  int num_routes  = 3;

  TSDARPGenome* gen = createGen(num_routes,main_route, gen_values);

  gen->ageMemory();

  checkTabuMemUpdate(*gen);

  delete gen;
}

TEST_F(TSDARPGenomeTest,updateTabuMemValuesFromInitialTabuLife) {
  int gen_values_a[] = {1,-1,-7,7,-8,8,2,-2,3,-3};
  vector<int> gen_values (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int main_route  = 0;
  int num_routes  = 3;

  TSDARPGenome* gen = createGen(num_routes,main_route, gen_values);
  pair<int,int> key1 = TSDARPGenome::routeVertIdKey(0,1);
  pair<int,int> key2 = TSDARPGenome::routeVertIdKey(0,-7);
  pair<int,int> key3 = TSDARPGenome::routeVertIdKey(1,-10);
  pair<int,int> key4 = TSDARPGenome::routeVertIdKey(1,-30);
  pair<int,int> key5 = TSDARPGenome::routeVertIdKey(1,100);
  pair<int,int> key6 = TSDARPGenome::routeVertIdKey(1,300);

  gen->addTabuMovement(key1);
  gen->addTabuMovement(key2);
  gen->addTabuMovement(key3);
  gen->addTabuMovement(key4);
  gen->addTabuMovement(key5);
  gen->addTabuMovement(key6);

  ASSERT_NE(0,gen->tabuMem_.size());

  gen->ageMemory();

  checkTabuMemUpdate(*gen);

  delete gen;
}

TEST_F(TSDARPGenomeTest,updateTabuMemValuesWithDifferentLifeValues) {
  int gen_values_a[] = {1,-1,-7,7,-8,8,2,-2,3,-3};
  vector<int> gen_values (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int main_route  = 0;
  int num_routes  = 3;

  TSDARPGenome* gen = createGen(num_routes,main_route, gen_values);
  pair<int,int> key1 = TSDARPGenome::routeVertIdKey(0,1);
  pair<int,int> key2 = TSDARPGenome::routeVertIdKey(0,-7);
  pair<int,int> key3 = TSDARPGenome::routeVertIdKey(1,-10);
  pair<int,int> key4 = TSDARPGenome::routeVertIdKey(1,-30);
  pair<int,int> key5 = TSDARPGenome::routeVertIdKey(1,100);
  pair<int,int> key6 = TSDARPGenome::routeVertIdKey(1,300);

  gen->tabuMem_[key1] = 1;
  gen->tabuMem_[key2] = 2;
  gen->tabuMem_[key3] = 3;
  gen->tabuMem_[key4] = 1;
  gen->tabuMem_[key5] = 4;
  gen->tabuMem_[key6] = 1;

  ASSERT_NE(0,gen->tabuMem_.size());

  gen->ageMemory();

  checkTabuMemUpdate(*gen);

  delete gen;
}

TEST_F(TSDARPGenomeTest,checkFreqValues) {
  int gen_values_a[] = {1,-1,-7,7,-8,8,2,-2,3,-3};
  vector<int> gen_values (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int main_route  = 0;
  int num_routes  = 3;

  TSDARPGenome* gen = createGen(num_routes,main_route, gen_values);
  pair<int,int> key1 = TSDARPGenome::routeVertIdKey(0,1);
  pair<int,int> key2 = TSDARPGenome::routeVertIdKey(1,-10);
  pair<int,int> key3 = TSDARPGenome::routeVertIdKey(2,-300);

  ASSERT_EQ(0,gen->opFreq_[key1]); ASSERT_EQ(0,gen->opFreq_[key2]); ASSERT_EQ(0,gen->opFreq_[key3]);

  gen->increaseFreq(key1); ASSERT_EQ(1,gen->opFreq_[key1]);
  gen->increaseFreq(key2); ASSERT_EQ(1,gen->opFreq_[key2]);
  gen->increaseFreq(key3); ASSERT_EQ(1,gen->opFreq_[key3]);
}

TEST_F(TSDARPGenomeTest,checkNeighborhoodMoveBasicTests) {
  int gen_values_a[] = {1,-1,-7,7,-8,8,2,-2,3,-3};
  vector<int> gen_values (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int orig_route = 0;
  int dest_route = 1;
  int vertid     = 1;
  checkMoveToNeighbor(gen_values,orig_route,dest_route,vertid);

  orig_route = 1;
  dest_route = 2;
  vertid     = 30;
  checkMoveToNeighbor(gen_values,orig_route,dest_route,vertid);

  orig_route = 2;
  dest_route = 0;
  vertid     = 800;
  checkMoveToNeighbor(gen_values,orig_route,dest_route,vertid);

// Only when compyling when debugging information
#ifdef DEBUG 
  // Checking a delivery vertex id
  orig_route = 0;
  dest_route = 1;
  vertid     = -1;
  ASSERT_DEATH(checkMoveToNeighbor(gen_values,orig_route,dest_route,vertid),"Assertion failed");

  // Checking for a vertex that does not exist
  orig_route = 0;
  dest_route = 1;
  vertid     = 198;
  ASSERT_DEATH(checkMoveToNeighbor(gen_values,orig_route,dest_route,vertid),"Assertion failed");
#endif

  // Check to move the last vertices pair
  int gen_values_b[] = {1,-1};
  vector<int> gen_values_b2 (gen_values_b, gen_values_b + sizeof(gen_values_b)/ sizeof(int) );
  orig_route = 0;
  dest_route = 2;
  vertid     = 1;
  checkMoveToNeighbor(gen_values_b2,orig_route,dest_route,vertid);
}


TEST_F(TSDARPGenomeTest,checkBestFromNeighborhoodSimpleTest) {
  int gen_values_a[] = {30,-30,2,-2,3,-3};
  vector<int> gen_values1 (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int gen_values_b[] = {10,-10,20,-20,100,-100};
  vector<int> gen_values2 (gen_values_b, gen_values_b + sizeof(gen_values_b)/ sizeof(int) );

  int gen_values_c[] = {200,-200,300,-300};
  vector<int> gen_values3 (gen_values_c, gen_values_c + sizeof(gen_values_c)/ sizeof(int) );

  vector< vector<int> > gen_values;
  gen_values.push_back(gen_values1); gen_values.push_back(gen_values2); gen_values.push_back(gen_values3);

  TSDARPGenome* gen  = createGen(gen_values);

  int vertid    = 30;
  int origroute = 0;
  int bestroute = 2;

  checkBestFromNeighborhood(*gen,vertid,origroute,bestroute);

  delete gen;
  // hay que probar lo mismo pero con un not allowedMove
  // y otro con un allowed move pero con frecuencia de manera que se penalice el movimiento y se esocja otro
}

TEST_F(TSDARPGenomeTest,checkBestFromNeighborhoodWithTabuTest) {
  int gen_values_a[] = {1,-1,3,-3};
  vector<int> gen_values1 (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int gen_values_b[] = {5,-5,7,-7};
  vector<int> gen_values2 (gen_values_b, gen_values_b + sizeof(gen_values_b)/ sizeof(int) );

  int gen_values_c[] = {9,-9};
  vector<int> gen_values3 (gen_values_c, gen_values_c + sizeof(gen_values_c)/ sizeof(int) );

  vector< vector<int> > gen_values;
  gen_values.push_back(gen_values1); gen_values.push_back(gen_values2); gen_values.push_back(gen_values3);

  TSDARPGenome* gen  = createGen(gen_values);

  // All the movements should produce the same score
  // We forbid the first two movements and expect to have third second movement executed
  gen->tabuMem_[TSDARPGenome::routeVertIdKey(1,1)] = 1;
  gen->tabuMem_[TSDARPGenome::routeVertIdKey(2,1)] = 1;

  int vertid    = 3;
  int origroute = 0;
  int bestroute = 1;

  checkBestFromNeighborhood(*gen,vertid,origroute,bestroute);

  // We check that they still exist in the tabu Memory
  ASSERT_NE(gen->tabuMem_.end(),gen->tabuMem_.find(TSDARPGenome::routeVertIdKey(1,1)));
  ASSERT_NE(gen->tabuMem_.end(),gen->tabuMem_.find(TSDARPGenome::routeVertIdKey(2,1)));

  // We recheck that the original tabu values have not been changed
  ASSERT_EQ(1, (gen->tabuMem_.find(TSDARPGenome::routeVertIdKey(1,1)))->second);
  ASSERT_EQ(1, (gen->tabuMem_.find(TSDARPGenome::routeVertIdKey(2,1)))->second);

  delete gen;
  // hay que probar lo mismo pero con un not allowedMove
  // y otro con un allowed move pero con frecuencia de manera que se penalice el movimiento y se esocja otro
}

TEST_F(TSDARPGenomeTest,checkBestFromNeighborhoodWithTabuTestAspirationCrit) {
  int gen_values_a[] = {30,-30,2,-2,3,-3};
  vector<int> gen_values1 (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int gen_values_b[] = {10,-10,20,-20,100,-100};
  vector<int> gen_values2 (gen_values_b, gen_values_b + sizeof(gen_values_b)/ sizeof(int) );

  int gen_values_c[] = {200,-200,300,-300};
  vector<int> gen_values3 (gen_values_c, gen_values_c + sizeof(gen_values_c)/ sizeof(int) );

  vector< vector<int> > gen_values;
  gen_values.push_back(gen_values1); gen_values.push_back(gen_values2); gen_values.push_back(gen_values3);

  TSDARPGenome* gen  = createGen(gen_values);

  // We forbid the first two best movements
  // However, since moving the 30 to the second route improves the score, the movement is allowed
  gen->tabuMem_[TSDARPGenome::routeVertIdKey(2,30)] = 1;
  gen->tabuMem_[TSDARPGenome::routeVertIdKey(1,30)] = 1;

  int vertid    = 30;
  int origroute = 0;
  int bestroute = 2;

  checkBestFromNeighborhood(*gen,vertid,origroute,bestroute);

  // We check that they still exist in the tabu Memory
  ASSERT_NE(gen->tabuMem_.end(),gen->tabuMem_.find(TSDARPGenome::routeVertIdKey(2,30)));
  ASSERT_NE(gen->tabuMem_.end(),gen->tabuMem_.find(TSDARPGenome::routeVertIdKey(1,30)));

  // We recheck that the original tabu values have not been changed
  ASSERT_EQ(1, (gen->tabuMem_.find(TSDARPGenome::routeVertIdKey(1,30)))->second);
  ASSERT_EQ(1, (gen->tabuMem_.find(TSDARPGenome::routeVertIdKey(2,30)))->second);

  delete gen;
  // hay que probar lo mismo pero con un not allowedMove
  // y otro con un allowed move pero con frecuencia de manera que se penalice el movimiento y se esocja otro
}

TEST_F(TSDARPGenomeTest,checkBestFromNeighborhoodWithDivPenalizationToBestSolTest) {
  int gen_values_a[] = {30,-30,2,-2,3,-3};
  vector<int> gen_values1 (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int gen_values_b[] = {10,-10,20,-20,100,-100};
  vector<int> gen_values2 (gen_values_b, gen_values_b + sizeof(gen_values_b)/ sizeof(int) );

  int gen_values_c[] = {200,-200,300,-300};
  vector<int> gen_values3 (gen_values_c, gen_values_c + sizeof(gen_values_c)/ sizeof(int) );

  vector< vector<int> > gen_values;
  gen_values.push_back(gen_values1); gen_values.push_back(gen_values2); gen_values.push_back(gen_values3);

  TSDARPGenome* gen  = createGen(gen_values);

  // We forbid the first two best movements
  int freq = 10;
  gen->opFreq_[TSDARPGenome::routeVertIdKey(2,30)] = freq;

  int vertid    = 30;
  int origroute = 0;
  int bestroute = 2;

  checkBestFromNeighborhood(*gen,vertid,origroute,bestroute);

  delete gen;
}

TEST_F(TSDARPGenomeTest,checkBestFromNeighborhoodWithInitialOptSolAndDivPenalizationTest) {
  int gen_values_a[] = {1,-1,2,-2};
  vector<int> gen_values1 (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int gen_values_b[] = {40,-40,50,-50};
  vector<int> gen_values2 (gen_values_b, gen_values_b + sizeof(gen_values_b)/ sizeof(int) );

  int gen_values_c[] = {700,-700,800,-800};
  vector<int> gen_values3 (gen_values_c, gen_values_c + sizeof(gen_values_c)/ sizeof(int) );

  vector< vector<int> > gen_values;
  gen_values.push_back(gen_values1); gen_values.push_back(gen_values2); gen_values.push_back(gen_values3);

  TSDARPGenome* gen  = createGen(gen_values);

  // We forbid the first movement
  gen->opFreq_[TSDARPGenome::routeVertIdKey(1,1)] = 10;

  int vertid    = 1;
  int origroute = 0;
  int bestroute = 2;

  checkBestFromNeighborhood(*gen,vertid,origroute,bestroute);

  delete gen;
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  GAGenome::optCriterion(GAGenome::MINIMIZATION);

  setUpTestingEnvironment();

  return RUN_ALL_TESTS();
}


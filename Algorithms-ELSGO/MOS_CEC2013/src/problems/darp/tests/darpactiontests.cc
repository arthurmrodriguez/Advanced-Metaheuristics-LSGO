#define USEDUMMYRANDOM

#include <vector>
#include <iostream>
#include <assert.h>

using namespace std;

vector<int>    randomintvalues;
vector<double> randomdoublevalues;

int randomintpos, randomdoublepos;

void resetRandomValues() {
  randomintpos = randomdoublepos = 0;
  randomintvalues.clear();
  randomdoublevalues.clear();
}

int GARandomIntDummy(int a, int b) {
  assert(randomintpos < randomintvalues.size());
  return randomintvalues[randomintpos++];
}

double GARandomDoubleDummy(double a=0, double b=1) {
  assert(randomdoublepos < randomdoublevalues.size());
  return randomdoublevalues[randomdoublepos++];
}

float GARandomFloatDummy(float low, float high) {
  return GARandomDoubleDummy(low,high);
}

#include "daropstests_aux.h"

class DARPActionTest : public ::testing::TestWithParam<int*> {
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
          addRequest(*requestslist, i*pow(10.0,j),    i ,      i*10,   5+(i-1)*15,  *distmatrix);
        }
      }

      dummymatrix = new DummyCostMatrix(distmatrix,*requestslist);
      DARPTestsVars::initData(requestslist,dummymatrix);

      resetRandomValues();
    }

    struct RandomData{
      virtual void placeRandomData() = 0;
    };

    struct NoRandomData : public RandomData {
      virtual void placeRandomData() {}
    };

    void setRandomData(vector<RandomData*>& rdata) {
      resetRandomValues();

      for (int i=0; i<rdata.size(); i++) {
        rdata[i]->placeRandomData();
      }
    }

    /*
     * This method is flexible enough to allow the testing of the application repeatedly of a shaker to
     * the same genome in order to compare the results achieve with the execution of the associated
     * action (generated each time the shaker is executed)
     */
    template<typename T>
    void checkShakerAction(vector<int>& gen_values, int main_route, int shaker_size, vector<RandomData*>& rdata) {
      DARPGenome* orig_gen = createGen(num_routes,main_route,gen_values);

      checkShakerAction<T>(*orig_gen, main_route, shaker_size, rdata);

      delete orig_gen;
    }

    template<typename T>
    void checkShakerAction(vector< vector<int> >& gen_values, int main_route, int shaker_size, vector<RandomData*>& rdata) {
      DARPGenome* orig_gen = createGen(gen_values);

      checkShakerAction<T>(*orig_gen, main_route, shaker_size, rdata);

      delete orig_gen;
    }

    template<typename T>
    void checkShakerAction(DARPGenome& orig_gen, int main_route, int shaker_size, vector<RandomData*>& rdata) {

      DARPGenome* sk_gen   = dynamic_cast<DARPGenome*>(orig_gen.clone());

      T sk(shaker_size);

      DARPGenome* act_gen = dynamic_cast<DARPGenome*>(orig_gen.clone());

      setRandomData(rdata);

      int num_shakers = rdata.size();

      vector<VNSOpAction*> actions;
      for (int i=0; i<num_shakers; i++) {
        VNSOpAction* act = sk(*sk_gen);
        act->executeOver(*act_gen);
        actions.push_back(act);
      }

      for (int route=0; route<sk_gen->size(); route++) {
        ASSERT_EQ(sk_gen->routeLength(route),act_gen->routeLength(route));
        for (int i=0; i<sk_gen->routeLength(route); i++) {
          ASSERT_EQ(sk_gen->gene(route,i),act_gen->gene(route,i));
        }
      }

/*
      cout << "original route" << endl;
      cout << *orig_gen << endl;
      cout << "shaker gen" << endl;
      cout << *sk_gen << endl;
      cout << "action" << endl;
      cout << *act_gen << endl;
*/

      delete sk_gen;
      delete act_gen;

      for (int i=0;i<actions.size(); i++) delete actions[i];
    }

    struct SwapData {
      int route;
      int vertex;
      int length;
    };

    struct SwapShakerRandomData : RandomData {
      SwapData r1;
      SwapData r2;

      virtual void placeRandomData() {
        randomintvalues.push_back(r1.route);
        randomintvalues.push_back(r2.route);
        randomintvalues.push_back(r1.vertex);
        randomintvalues.push_back(r1.length);
        randomintvalues.push_back(r2.vertex);
        randomintvalues.push_back(r2.length);
      }
    };

    struct ChainIterStepData {
      int size;
      int route;
    };

    struct ChainDataRandomData : RandomData {
      int starting_route_orig;
      int starting_route_dest;
      int starting_vertex;
      int starting_length;

      vector<ChainIterStepData> iterations_data;

      virtual void placeRandomData() {
        randomintvalues.push_back(starting_route_orig);
        randomintvalues.push_back(starting_route_dest);
        randomintvalues.push_back(starting_vertex);
        randomintvalues.push_back(starting_length);

        for (int i=0; i<iterations_data.size(); i++) {
          randomintvalues.push_back(iterations_data[i].size);
          randomintvalues.push_back(iterations_data[i].route);
        }
      }

    };

    struct GreedyMoveRandomData : RandomData {
      double select_value1;
      double select_value2;

      virtual void placeRandomData() {
        randomdoublevalues.push_back(select_value1);
        randomdoublevalues.push_back(select_value2);
      }
    };

    struct BadAndGoodRoutesRandomValues : RandomData {
      double b_value1;
      double b_value2;
      double g_value1;
      double g_value2;

      virtual void placeRandomData() {
        randomdoublevalues.push_back(b_value1);
        randomdoublevalues.push_back(b_value2);
        randomdoublevalues.push_back(g_value1);
        randomdoublevalues.push_back(g_value2);
      }
    };


    virtual void TearDown() {
      delete dummymatrix;
      delete requestslist;
    }

};



TEST_F(DARPActionTest,simpleSwap) {
  int gen_values_a[] = {1,-1,-7,7,-8,8,2,-2,3,-3};
  vector<int> gen_values (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int main_route  = 1;
  int shaker_size = 4;

  vector<RandomData*> swap_data;
  SwapShakerRandomData s1;
  s1.r1.route  = 0;
  s1.r1.vertex = 1;
  s1.r1.length = 3;
  s1.r2.route  = 2;
  s1.r2.vertex = gen_values.size() / 2;
  s1.r2.length = gen_values.size() - s1.r2.vertex;

  swap_data.push_back(&s1);

  checkShakerAction<SwapNeighborhood>(gen_values, main_route, shaker_size, swap_data);
}

TEST_F(DARPActionTest,TwoSwapsChained) {
  int gen_values_a[] = {1,-1,-7,7,-8,8,2,-2,3,-3};
  vector<int> gen_values (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int main_route  = 1;
  int shaker_size = 4;

  vector<RandomData*> swap_data;

  SwapShakerRandomData s1;
  s1.r1.route  = 0;
  s1.r1.vertex = 1;
  s1.r1.length = 3;
  s1.r2.route  = 2;
  s1.r2.vertex = gen_values.size() / 2;
  s1.r2.length = gen_values.size() - s1.r2.vertex;

  SwapShakerRandomData s2;
  s2.r1.route  = 1;
  s2.r1.vertex = 1;
  s2.r1.length = 4;
  s2.r2.route  = 2;
  s2.r2.vertex = gen_values.size() -3;
  s2.r2.length = 3;

  swap_data.push_back(&s1);
  swap_data.push_back(&s2);


  checkShakerAction<SwapNeighborhood>(gen_values,main_route, shaker_size,swap_data);
}

TEST_F(DARPActionTest,SimpleChainNeighborhood) {
  int gen_values_a[] = {1,-1,-7,7,-8,8,2,-2,3,-3};
  vector<int> gen_values (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int main_route  = 1;
  int shaker_size = 2;

  vector<RandomData*> chain_data;

  ChainDataRandomData ch1;
  ch1.starting_route_orig = 1;
  ch1.starting_route_dest = 0;
  ch1.starting_vertex     = 0;
  ch1.starting_length     = 2;

  ChainIterStepData iter2;
  iter2.route = 2;
  iter2.size  = 2;

  ch1.iterations_data.push_back(iter2);

  chain_data.push_back(&ch1);

  checkShakerAction<ChainNeighborhood>(gen_values,main_route, shaker_size,chain_data);
}

TEST_F(DARPActionTest,ComplexChainNeighborhood) {
  int gen_values_a[] = {1,-1,-7,7,-8,8,2,-2,3,-3};
  vector<int> gen_values (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int main_route  = 1;
  int shaker_size = 3;

  vector<RandomData*> chain_data;

  ChainDataRandomData ch1;
  ch1.starting_route_orig = 1;
  ch1.starting_route_dest = 0;
  ch1.starting_vertex     = 1;
  ch1.starting_length     = 2;

  ChainIterStepData iter2;
  iter2.route = 1;
  iter2.size  = 3;
  ChainIterStepData iter3;
  iter3.route = 2;
  iter3.size  = 4;

  ch1.iterations_data.push_back(iter2);
  ch1.iterations_data.push_back(iter3);

  chain_data.push_back(&ch1);

  checkShakerAction<ChainNeighborhood>(gen_values,main_route, shaker_size,chain_data);
}

TEST_F(DARPActionTest,TwoChainNeighborhoodChained) {
  int gen_values_a[] = {1,-1,-7,7,-8,8,2,-2,3,-3};
  vector<int> gen_values (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int main_route  = 1;
  int shaker_size = 3;

  vector<RandomData*> chain_data;


  ChainDataRandomData ch1;
  ch1.starting_route_orig = 1;
  ch1.starting_route_dest = 0;
  ch1.starting_vertex     = 1;
  ch1.starting_length     = 2;

  ChainIterStepData ch1_iter2;
  ch1_iter2.route = 1;
  ch1_iter2.size  = 3;
  ChainIterStepData ch1_iter3;
  ch1_iter3.route = 2;
  ch1_iter3.size  = 4;

  ch1.iterations_data.push_back(ch1_iter2);
  ch1.iterations_data.push_back(ch1_iter3);


  ChainDataRandomData ch2;
  ch2.starting_route_orig = 2;
  ch2.starting_route_dest = 1;
  ch2.starting_vertex     = 1;
  ch2.starting_length     = 4;

  ChainIterStepData ch2_iter2;
  ch2_iter2.route = 0;
  ch2_iter2.size  = 2;
  ChainIterStepData ch2_iter3;
  ch2_iter3.route = 2;
  ch2_iter3.size  = 1;


  ch2.iterations_data.push_back(ch2_iter2);
  ch2.iterations_data.push_back(ch2_iter3);

  chain_data.push_back(&ch1);
  chain_data.push_back(&ch2);

  checkShakerAction<ChainNeighborhood>(gen_values,main_route, shaker_size,chain_data);
}

TEST_F(DARPActionTest,SimpleGreedyMove) {
  int gen_values_a[] = {1,-1,-7,7,-8,8,2,-2,3,-3};
  vector<int> gen_values (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int main_route  = 1;
  int shaker_size = 4;

  vector<RandomData*> data;

  GreedyMoveRandomData gm1;
  gm1.select_value1 = 0.2;
  gm1.select_value1 = 0.6;

  data.push_back(&gm1);

  checkShakerAction<GreedyMoveNeighborhood>(gen_values,main_route, shaker_size,data);
}

TEST_F(DARPActionTest,SimpleGreedyMoveChained) {
  int gen_values_a[] = {1,-1,-7,7,-8,8,2,-2,3,-3};
  vector<int> gen_values (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int main_route  = 1;
  int shaker_size = 4;

  vector<RandomData*> data;

  GreedyMoveRandomData gm1;
  gm1.select_value1 = 0.2;
  gm1.select_value1 = 0.6;

  GreedyMoveRandomData gm2;
  gm2.select_value1 = 0.9;
  gm2.select_value1 = 0.7;


  data.push_back(&gm1);
  data.push_back(&gm2);

  checkShakerAction<GreedyMoveNeighborhood>(gen_values,main_route, shaker_size,data);
}

TEST_F(DARPActionTest,SimpleGreedyMoveDestCenteredAndGreedySwapChained) {
  int gen_values_a[] = {1,-1,-7,7,-8,8,2,-2,3,-3};
  vector<int> gen_values (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int main_route  = 1;
  int shaker_size = 4;

  vector<RandomData*> data;

  // Since the routes are equivalent, each one is going to have the same probability
  // best route = GAMin (value1, value2)
  // worst route = GAMax (1-value1, 1-value2)
  BadAndGoodRoutesRandomValues gm1;
  gm1.b_value1 = 0.2;
  gm1.b_value2 = 0.6;  // best route 0
  gm1.g_value1 = 0.3;
  gm1.g_value2 = 0.1;  // best route 2

  BadAndGoodRoutesRandomValues gm2;
  gm2.b_value1 = 0.9;
  gm2.b_value2 = 0.9;  // best route 2
  gm2.g_value1 = 0.9;
  gm2.g_value2 = 0.9;  // worst route 0

  data.push_back(&gm1);
  data.push_back(&gm2);

  checkShakerAction<GreedyMoveNeighborhoodDestCentered>(gen_values,main_route, shaker_size,data);
  checkShakerAction<GreedySwapNeighborhood>(gen_values,main_route, shaker_size,data);
}

TEST_F(DARPActionTest,SimpleCheckAllNaturalSeqsCombsNeighborhoodTest) {
  int gen_values_a1[] = {1,10,-1,-10,7,70,-70,-7};
  vector<int> gen_values1 (gen_values_a1, gen_values_a1 + sizeof(gen_values_a1)/ sizeof(int) );

  int gen_values_a2[] = {2,20,-20,-2,6,60,-6,-60};
  vector<int> gen_values2 (gen_values_a2, gen_values_a2 + sizeof(gen_values_a2)/ sizeof(int) );

  int gen_values_a3[] = {100,-100,3,300,-300,-3};
  vector<int> gen_values3 (gen_values_a3, gen_values_a3 + sizeof(gen_values_a3)/ sizeof(int) );

  vector< vector<int> > gen_values;
  gen_values.push_back(gen_values1);
  gen_values.push_back(gen_values2);
  gen_values.push_back(gen_values3);

  int main_route  = 1;

  vector<RandomData*> data;
  NoRandomData n1;
  data.push_back(&n1);
  data.push_back(&n1);

  checkShakerAction<CheckAllNaturalSeqsCombsNeighborhood>(gen_values,main_route,0,data);
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  GAGenome::optCriterion(GAGenome::MINIMIZATION);

  setUpTestingEnvironment();

  return RUN_ALL_TESTS();
}


#include "daropstests_aux.h"
#include "../VNSDARPGenome.cc"


//NOTE: These tests need to be updated to the new defintion of the local search where no restarts are being made

extern "C" double ReturnZeroObjFunc(GAGenome& g) {

  return 0;
}

extern "C" double LSDummyObjFunc(GAGenome& g) {
  double score = DummyObjFunc(g);
  DARPGenome& gen = dynamic_cast<DARPGenome&>(g);

  VerticesList& requestslist = *DARPTestsVars::verticeslist;
 
  const int NONFEASIBLE_penalization = 10;

  // For each infeasible arc (only that a non critical vertex comes before than a critical one for this case) a penalization is added
  for (int route=0; route<gen.length(); route++){
    for (int i=0; i<gen.routeLength(route); i++) {
      Vertex& req = requestslist.getVertex( gen.gene(route,i) );
      if (req.critic_) {
        int noncrit_vert_id  = req.getSiblingVertexId();
        int noncrit_vert_pos = gen.findPosOfVertex(route,noncrit_vert_id);
        if (noncrit_vert_pos < i ) score += NONFEASIBLE_penalization;
      }
    }
  }
  


  //cout << "final score: " << score << " obj end" << endl;

  return score;
}

class LocalSearchTest : public ::testing::TestWithParam<int*> {
public:
  static const int MAXPOINTS;
  static const int MINCOST;
  static const int MAXCOST;

  VerticesList* requestslist;
  DummyCostMatrix*  dummymatrix;
  int           num_routes;

    virtual void SetUp() {
      num_routes   = 3;

      requestslist = new VerticesList();

      DummyDistMatrix* distmatrix = new DummyDistMatrix(MAXPOINTS);

      // The costs are set so that the optimum sequence is composed of the sequence of positions set in allpos
      // The remaining combinations are set to a higher cost
      int allpos[] = {1,10,2,20,3,30};
      vector<int> allposv (allpos, allpos + sizeof(allpos)/ sizeof(int) );
      for (int i=0; i<allposv.size()-1; i++) {
        for (int j=i+1; j<allposv.size(); j++) {
          int cost = (j==i+1) ? MINCOST : MAXCOST;
          distmatrix->setCost(allposv[i],allposv[j],cost);
        }
      }

      // We create some requests with the same positions but different ids
      // Note that the requests 3,30 and 300 cannot go before 1,10 and 100 because 30+cost > 5+MAXWAITTIME (which is 20)
      addRequest(*requestslist, 1, 1, 10, 5,  *distmatrix);
      addRequest(*requestslist, 2, 2, 20, 15, *distmatrix);
      addRequest(*requestslist, 3, 3, 30, 30, *distmatrix);

      addRequest(*requestslist, 10, 1, 10, 5,  *distmatrix);
      addRequest(*requestslist, 20, 2, 20, 15, *distmatrix);
      addRequest(*requestslist, 30, 3, 30, 30, *distmatrix);

      addRequest(*requestslist, 100, 1, 10, 5,  *distmatrix);
      addRequest(*requestslist, 200, 2, 20, 15, *distmatrix);
      addRequest(*requestslist, 300, 3, 30, 30, *distmatrix);

      dummymatrix = new DummyCostMatrix(distmatrix,*requestslist);
      DARPTestsVars::initData(requestslist,dummymatrix);
    }

    virtual void TearDown() {
      delete dummymatrix;
      delete requestslist;
    }
};
const int LocalSearchTest::MAXPOINTS = 300;
const int LocalSearchTest::MINCOST   = 1;
const int LocalSearchTest::MAXCOST   = 2;


TEST_F(LocalSearchTest,TestNoImprovement) {
  int main_route = 0;

  DARPGenome* gen = new DARPGenome(num_routes,0,LSDummyObjFunc);

  int seq[] = {1,-1,2,-2,3,-3};
  vector<int> seqv (seq, seq + sizeof(seq)/ sizeof(int) );
  for (int i=0; i<seqv.size(); i++) {
    gen->pushBackVertex(main_route  ,seqv[i]);
    gen->pushBackVertex((main_route+1)% num_routes,10*seqv[i]);
    gen->pushBackVertex((main_route+2)% num_routes,100*seqv[i]);
  }

  DARPGenome* copy_gen = new DARPGenome(*gen);
  DARPLocalSearch ls;
  for (int i=0; i<copy_gen->length(); i++) ls.localSearch(*copy_gen,i);
  
  for (int route=0; route<copy_gen->length(); route++) {
    for (int pos=0; pos<copy_gen->routeLength(route); pos++) {
      ASSERT_EQ(gen->gene(route,pos),copy_gen->gene(route,pos)) << "route: " << route << " pos: " << pos << " gen: " << gen->gene(route,pos) << " copy_gen: " << copy_gen->gene(route,pos);
    }
  }

  ASSERT_EQ(gen->score(),copy_gen->score());

  delete copy_gen;
  delete gen;
}

TEST_P(LocalSearchTest,CheckWithBestSeq) {
//  cout << "windows are: " << *requestslist << endl;
  int main_route = 0;

  DARPGenome* gen = new DARPGenome(num_routes,0,LSDummyObjFunc);

  int  values_size = GetParam()[0];
  int  iterations  = GetParam()[1];
  int* values      = GetParam();
  values += 2; // We need to advance the pointer to avoid the values_size and iterations parameter

  for (int i=0; i<values_size; i++) {
    gen->pushBackVertex(main_route  ,values[i]);
    gen->pushBackVertex((main_route+1)% num_routes,10*values[i]);
    gen->pushBackVertex((main_route+2)% num_routes,100*values[i]);
  }

  DARPLocalSearch ls;

  for (int i=0; i<iterations; i++) {
    if (i>0) { // We need this step if we like that the ls is applied again, otherwise the status is going to be EXPLORED
      for (int i=0; i<num_routes; i++) gen->localSearchedRoute(i,DARPGenome::LSUNEXPLORED);
    }
    ls(*gen); // as many iterations as specified for the test
  }

  for (int i=0; i<num_routes; i++) {
    ASSERT_EQ(DARPGenome::LSEXPLORED,gen->localSearchedRoute(i));

    for (int j=0; j<gen->routeLength(main_route)-1; j++) {
      int cost_consecutive = dummymatrix->getCost(gen->gene(i,j), gen->gene(i,j+1) );

      ASSERT_EQ( LocalSearchTest::MINCOST, cost_consecutive );
    }
  }

  delete gen;
}




int size          = 6;
int oneiteration  = 1;
int twoiterations = 2;
int threeiterations = 3;

int oneswap[] = {size,oneiteration,2,-2,1,-1,3,-3};
INSTANTIATE_TEST_CASE_P(OneSwap,LocalSearchTest,::testing::Values(oneswap));

int oneswap2[] = {size,oneiteration,2,-2,3,-3,1,-1};
INSTANTIATE_TEST_CASE_P(OneSwap2,LocalSearchTest,::testing::Values(oneswap2));

int inverseseq[] = {size,oneiteration,3,-3,2,-2,1,-1};
INSTANTIATE_TEST_CASE_P(InverseSeq,LocalSearchTest,::testing::Values(inverseseq));

int rotatedseq[] = {size,oneiteration,3,-3,1,-1,2,-2,};
INSTANTIATE_TEST_CASE_P(RotatedSeq,LocalSearchTest,::testing::Values(rotatedseq));


// Due to the change in the local search, now it stops as soon as it discover an improvement so
// in order to solve these problems, it needs to do several passes

int shuffledseq[] = {size,twoiterations,2,3,1,-3,-2,-1,};
INSTANTIATE_TEST_CASE_P(ShuffleSeq,LocalSearchTest,::testing::Values(shuffledseq));

int shuffledseq2[] = {size,twoiterations,2,-2,1,3,-3,-1,};
INSTANTIATE_TEST_CASE_P(ShuffleSeq2,LocalSearchTest,::testing::Values(shuffledseq2));

int shuffledseq3[] = {size,twoiterations,2,1,-2,-1,3,-3};
INSTANTIATE_TEST_CASE_P(ShuffleSeq3,LocalSearchTest,::testing::Values(shuffledseq3));


/*
 * Due to the non extreme greedyness of the ls, this permutation needs the following three iterations
1 2 3 -3 -2 -1

1 3 -3 -1 2 -2
3 -3 1 -1 2 -2
1 -1 2 -2 3 -3
*/

int shuffledseq4[] = {size,threeiterations,1,3,2,-3,-2,-1};
INSTANTIATE_TEST_CASE_P(ShuffleSeq4,LocalSearchTest,::testing::Values(shuffledseq4));

TEST_F(LocalSearchTest,LSRestartTest) {
  /*
   * In this test several scenarios are tested continuously. Not the cleanest way but for now serves
   * for checking the implemented mechanism
   */
  DARPGenome* gen = new VNSDARPGenome(num_routes,0,1,LSDummyObjFunc);

  int values_size = 6;
  int values[] = {2,-2,1,-1,3,-3};

  // A genome with three routes is created. This routes have the same values but multiplied by
  // a power of 10
  for (int i=0; i<values_size; i++) {
    gen->pushBackVertex(0,    values[i]);
    gen->pushBackVertex(1, 10*values[i]);
    gen->pushBackVertex(2,100*values[i]);
  }
  gen->evaluate(); // so that the number count is easier to track

  // No change should be conducted since we are not giving enough evals to perform any change
  int available_evals;
  int starting_gen_evals;

  /*
   * First test, not enough evals to do anything. Note that we have two positions, the first one, localSearchLastPosExplored,
   * which points to the next position where the local search should find the next critical vertex, and the second one,
   * localSearchLastCritVertexSearchPosExplored which marks the position (from localSearchLastPosExplored) that the search was
   * trying with the critical vertex
   * The state should remain the same
   */
  {
    available_evals = 1;

    starting_gen_evals = gen->nevals();

    DARPLocalSearch ls(available_evals);

    ls(*gen);

    for (int route=0; route<gen->length(); route++) {
      if (route==0) {
        ASSERT_EQ(DARPGenome::LSUNFINISHED,gen->localSearchedRoute(route));
        ASSERT_EQ(0,gen->localSearchLastPosExplored(route));
        ASSERT_EQ(0,gen->localSearchLastCritVertexSearchPosExplored(route));
      }
      else {
        ASSERT_EQ(DARPGenome::LSUNEXPLORED,gen->localSearchedRoute(route));
      }
      for (int pos=0; pos<gen->routeLength(route); pos++) {
        ASSERT_EQ(values[pos]*pow((float)10,route), gen->gene(route,pos) );
      }
    }

    ASSERT_EQ(starting_gen_evals+available_evals,gen->nevals());
  }


  /*
   * enough evals to test all the possible insertions of the vertex -2 when 2 is placed at the first position (position 0)
   * The next call (not this one) should start with the vertex 2 at position 1
   */
  {
    available_evals    = gen->routeLength(0); // We need to test all inserting positions for -2 if 2 is at pos 0
    starting_gen_evals = gen->nevals();       // That implies 5 positions + 1 to advance to next position (2 at pos 1) so that
                                              // critvertex search pos gets updated
    DARPLocalSearch ls(available_evals);

    ls(*gen);
    for (int route=0; route<gen->length(); route++) {
      if (route==0) {
        ASSERT_EQ(DARPGenome::LSUNFINISHED,gen->localSearchedRoute(route));
        ASSERT_EQ(0,gen->localSearchLastPosExplored(route));
        ASSERT_EQ(1,gen->localSearchLastCritVertexSearchPosExplored(route));
      }
      else {
        ASSERT_EQ(DARPGenome::LSUNEXPLORED,gen->localSearchedRoute(route));
      }
      for (int pos=0; pos<gen->routeLength(route); pos++) {
        ASSERT_EQ(values[pos]*pow((float)10,route), gen->gene(route,pos) );
      }
    }
  }


  /*
   * Here we give enough evaluations (the exact number) to find a better solution. The call should stop when this better
   * solution is found and set the appropriate state. In fact, due to the implementation details, the call stops with the next
   * position and with and unexplored state so that the next call continues where this execution finished.
   */
  {
    available_evals    = gen->routeLength(0)-1; // We need 4 evals to test -2 at all positions starting from position 2 since
    starting_gen_evals = gen->nevals();         // vertex 2 is placed at position 1 + 1 evaluation for finding the optimum
    double start_score = gen->score();          // vertex 2 at position 2 and -2 at position 3

    DARPLocalSearch ls(available_evals);

    ls(*gen);

    // In this case the ls has found a better solution so the end state is explored (and the values should be in the best position)
    for (int route=0; route<gen->length(); route++) {
      if (route==0) {
        ASSERT_TRUE( GAGenome::compareScores(gen->score(),start_score) == GAGenome::BETTER );
        ASSERT_EQ(DARPGenome::LSUNFINISHED,gen->localSearchedRoute(route));
        ASSERT_EQ(2,gen->localSearchLastPosExplored(route));
        ASSERT_EQ(-1,gen->localSearchLastCritVertexSearchPosExplored(route));

        // Check that we have the good values
        for (int pos=0; pos<gen->routeLength(route)-1; pos++) {
          int cost_consecutive = dummymatrix->getCost(gen->gene(route,pos), gen->gene(route,pos+1) );
          ASSERT_EQ( LocalSearchTest::MINCOST, cost_consecutive );
        }

      }
      else {
        // Check that the remaining routes have not been modified
        ASSERT_EQ(DARPGenome::LSUNEXPLORED,gen->localSearchedRoute(route));
        // Check that we have the same values that the original values
        for (int pos=0; pos<gen->routeLength(route); pos++) {
          ASSERT_EQ(values[pos]*pow((float)10,route), gen->gene(route,pos) );
        }

      }
    }
  }

  /*
   * Here we give enough evals to end the search of the first route. With this test the loop while ( pos < gen.routeLength(route)-1 )
   * from the localsearch method is going to be called several times testing its behavior.
   */
  {
      available_evals    = (5+4+3+2+1) + (4+3+2+1) + 1; // enough evals for exploring the remaing route. It starts with position
      starting_gen_evals = gen->nevals();               // 2 (i.e. vertex 2) and tries to insert first the vertex 2 at position 0
      double start_score = gen->score();                // trying all the insert position for vertex -2. Then, it repeats the same
                                                        // for position 1 (4+3+2+1) and then it finally founds it with the next call
      DARPLocalSearch ls(available_evals);
      ls(*gen);

      // In this case the ls has not modified the solution and should contain the same best results as with the previous case
      // However, in this case, the ls has finished with the first route so the state should be explored
      for (int route=0; route<gen->length(); route++) {
        if (route==0) {
          ASSERT_TRUE( GAGenome::compareScores(gen->score(),start_score) == GAGenome::EQUAL );
          ASSERT_EQ(DARPGenome::LSEXPLORED,gen->localSearchedRoute(route));

          // Check that we have the good values
          for (int pos=0; pos<gen->routeLength(route)-1; pos++) {
            int cost_consecutive = dummymatrix->getCost(gen->gene(route,pos), gen->gene(route,pos+1) );
            ASSERT_EQ( LocalSearchTest::MINCOST, cost_consecutive );
          }

        }
        else {
          // Check that the remaining routes have not been modified
          if (route == 1) ASSERT_EQ(DARPGenome::LSUNFINISHED,gen->localSearchedRoute(route));
          else            ASSERT_EQ(DARPGenome::LSUNEXPLORED,gen->localSearchedRoute(route));
          // Check that we have the same values that the original values
          for (int pos=0; pos<gen->routeLength(route); pos++) {
            ASSERT_EQ(values[pos]*pow((float)10,route), gen->gene(route,pos) );
          }

        }
      }
  }

  // For this test we check if the ls works with two routes by giving the exact evals to be applied to the remaining two routes
  {
      available_evals    = 2*(5+4+3+ 5+4+3+2+1 + 4+3+2+1 ) + 1; // First the 20(or 200) vertices are tried at positions 0, 1, 2
      starting_gen_evals = gen->nevals();                       // Then the best result is found so the search breaks and continues.
      double start_score = gen->score();                        // with the next position: 2 (since 1 corresponds to a delivery vertex)
                                                                // since no better solution is found all the possibilities are explored
                                                                // (5+4+3+2+1) the same happens with the next position 4 and
      DARPLocalSearch ls(available_evals);                      // finally, one last evaluation to finish with an evaluation remaining
      ls(*gen);                                                 // Otherwise, it would seem as unfinished since we do not save the
                                                                // non critical vertex position

      for (int route=0; route<gen->length(); route++) {

        ASSERT_TRUE( GAGenome::compareScores(gen->score(),start_score) == GAGenome::BETTER );
        ASSERT_EQ(DARPGenome::LSEXPLORED,gen->localSearchedRoute(route));

        // Check that we have the good values
        for (int pos=0; pos<gen->routeLength(route)-1; pos++) {
          int cost_consecutive = dummymatrix->getCost(gen->gene(route,pos), gen->gene(route,pos+1) );
          ASSERT_EQ( LocalSearchTest::MINCOST, cost_consecutive );
        }

      }

  }
  cout << "queda ordenar estos tests y repetirlos sin que el mejor orden sea 1-1 2 -2 sino con otro mas aleatorio" << endl;

  delete gen;
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  GAGenome::optCriterion(GAGenome::MINIMIZATION);

  setUpTestingEnvironment();

  return RUN_ALL_TESTS();
}

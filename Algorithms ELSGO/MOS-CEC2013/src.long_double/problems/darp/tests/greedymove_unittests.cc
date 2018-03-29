#include "daropstests_aux.h"

class GreedyMoveTest : public ::testing::TestWithParam<int*> {
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
          distmatrix->setCostOneSide( *it,*itj,3);
          distmatrix->setCostOneSide( *itj,*it,3);
        }
      }

      for (int i=0; i<allposv.size()-1; i++) distmatrix->setCostOneSide(allposv[i],allposv[i+1],1);

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
    }

    virtual void TearDown() {
      delete dummymatrix;
      delete requestslist;
    }

    bool areRouteCostEqual(DARPVNSShaker::RouteAndScore& routecost1, DARPVNSShaker::RouteAndScore& routecost2) {
      return routecost1.route == routecost2.route and routecost1.score == routecost2.score;
    }

    void checkRouteScores(vector< vector<int> >& all_values) {
      DARPGenome* gen = createGen(all_values);

      // The three alternatives are obtained, the default should be equaled to the maximization one
      vector<DARPVNSShaker::RouteAndScore> routecosts_def = DARPVNSShaker::getRouteAndScoresSorted(*gen);
      vector<DARPVNSShaker::RouteAndScore> routecosts_max = DARPVNSShaker::getRouteAndScoresSorted(*gen,GAGenome::MAXIMIZATION);
      vector<DARPVNSShaker::RouteAndScore> routecosts_min = DARPVNSShaker::getRouteAndScoresSorted(*gen,GAGenome::MINIMIZATION);

      ASSERT_EQ(routecosts_def.size(),all_values.size());
      ASSERT_EQ(routecosts_max.size(),all_values.size());
      ASSERT_EQ(routecosts_min.size(),all_values.size());

      // Check that the default option is maximization
      ASSERT_EQ(routecosts_def.size(),routecosts_max.size());
      for (int i=0; i<routecosts_def.size(); i++) {
        ASSERT_EQ(routecosts_def[i].score, routecosts_max[i].score);
        ASSERT_EQ(routecosts_def[i].route,routecosts_max[i].route);
      }

      //cout << "route costs are: " << endl;
      //cout << "max: " ; for (int i=0; i<routecosts_max.size(); i++) cout << routecosts_max[i].score << " "; cout << endl;
      //cout << "min: " ; for (int i=0; i<routecosts_min.size(); i++) cout << routecosts_min[i].score << " "; cout << endl;


      for (int i=0; i<gen->size(); i++) {
        for (int j=i+1; j<gen->size(); j++) {
          if (i==j) continue;
          if ( !areRoutesEqual(*gen,i,j) ) {
            ASSERT_EQ(false, areRouteCostEqual(routecosts_min[i],routecosts_min[j]));
            ASSERT_EQ(false, areRouteCostEqual(routecosts_max[i],routecosts_max[j]));
          }
        }
      }

      // We check the relationship between the elements is the expected one (according to each type of sorting)
      for (int i=0; i<routecosts_max.size()-1; i++) ASSERT_GE(routecosts_max[i].score,routecosts_max[i+1].score);
      for (int i=0; i<routecosts_max.size()-1; i++) ASSERT_LE(routecosts_min[i].score,routecosts_min[i+1].score);

      delete gen;
    }

    void checkProbsFromGen(vector <vector<int> >& all_values) {
      DARPGenome* gen = createGen(all_values);

      vector<DARPVNSShaker::RouteAndScore> routecosts_max = DARPVNSShaker::getRouteAndScoresSorted(*gen,GAGenome::MAXIMIZATION);
      vector<DARPVNSShaker::RouteAndScore> routecosts_min = DARPVNSShaker::getRouteAndScoresSorted(*gen,GAGenome::MINIMIZATION);

      vector<long double> probs_def = DARPVNSShaker::getSelectionProbs(routecosts_max);
      vector<long double> probs_max = DARPVNSShaker::getSelectionProbs(routecosts_max,GAGenome::MAXIMIZATION);
      vector<long double> probs_min = DARPVNSShaker::getSelectionProbs(routecosts_min,GAGenome::MINIMIZATION);

      //Check that the default option is maximization
      ASSERT_EQ(probs_def.size(),probs_max.size());
      for (int i=0; i<probs_def.size(); i++) ASSERT_EQ( probs_def[i], probs_max[i]);

      ASSERT_EQ(probs_max.size(),probs_min.size());

      //Chech that the values are always in ascending order
      for (int i=0; i<probs_max.size()-1; i++) {
        ASSERT_LT(probs_max[i],probs_max[i+1]);
        ASSERT_LT(probs_min[i],probs_min[i+1]);
      }

      // Check that max and min have different values
      for (int i=0; i<probs_max.size()-1; i++) ASSERT_NE(probs_max[i],probs_min[i]);

      // Check that the last positions are 1.0
      ASSERT_EQ(probs_max[probs_max.size()-1],1.0);
      ASSERT_EQ(probs_max[probs_min.size()-1],1.0);

      delete gen;
    }

    void checkBestSeqToRemove(int num_routes, int main_route, vector<int>& gen_values, list<int>& bad_seq, int insert_pos, int shaker_size) {
      DARPGenome* gen = createGen(num_routes,main_route,gen_values);

      gen->insertVertices(main_route,insert_pos,bad_seq);
      //cout << "depsues de insertar gen es " << *gen << endl;

      GreedyMoveNeighborhood shaker(shaker_size);

      list<int> del_seq = shaker.bestSeqToBeRemoved(*gen,main_route);

      //cout << "la secuencia extraida es: "; printNatSeq(&del_seq); cout << endl;

      // Each element from the extracted sequence is found in the bad sequence passed as argument
      for (list<int>::iterator delit=del_seq.begin(); delit!= del_seq.end(); delit++) {
        ASSERT_NE(bad_seq.end() , find(bad_seq.begin(),bad_seq.end(),*delit));
      }

      if (shaker_size == bad_seq.size()) ASSERT_EQ(shaker_size,del_seq.size());

      int orig_size = gen->routeLength(main_route);

      // Check that if we remove the vertices, the gen remains the same
      gen->removeVertices(main_route,del_seq);
      ASSERT_EQ(orig_size-del_seq.size(),gen->routeLength(main_route));
      for (int i=0; i<gen->routeLength(main_route); i++) {
        ASSERT_EQ(del_seq.end(),find(del_seq.begin(),del_seq.end(),gen->gene(main_route,i)));
      }

      delete gen;
    }

    void checkBestRouteForInsertion(vector< vector<int> >& gen_values, list<int>& ins_seq, int orig_route, long double score_orig, int expected_value) {
      DARPGenome* gen = createGen(gen_values);

      int best_route = GreedyMoveNeighborhood::bestRoutePosForInsertion(*gen,orig_route,score_orig,ins_seq);

      ASSERT_EQ(expected_value, best_route);

      delete gen;
    }
};


TEST_F(GreedyMoveTest,RouteScoresSortedSameValues) {
  // Three routes are created having the _1 route the best score followed by the _2 and being the _3 the worst route.

  int a_gen_values_1[] = {1,-1,2,-2,3,-3,4,-4,5,-5}; vector<int> gen_values_1 (a_gen_values_1, a_gen_values_1 + sizeof(a_gen_values_1)/ sizeof(int) );
  int a_gen_values_2[] = {1,-1,2,-2,3,-3,4,-4,5,-5}; vector<int> gen_values_2 (a_gen_values_2, a_gen_values_2 + sizeof(a_gen_values_2)/ sizeof(int) );
  int a_gen_values_3[] = {1,-1,2,-2,3,-3,4,-4,5,-5}; vector<int> gen_values_3 (a_gen_values_3, a_gen_values_3 + sizeof(a_gen_values_3)/ sizeof(int) );

  vector< vector<int> > all_values;
  all_values.push_back(gen_values_1); // The routes are inserted shuffled
  all_values.push_back(gen_values_2);
  all_values.push_back(gen_values_3);

  checkRouteScores(all_values);

}

TEST_F(GreedyMoveTest,RouteScoresSortedDiffValues) {
  // Three routes are created having the _1 route the best score followed by the _2 and being the _3 the worst route.

  int a_gen_values_1[] = {1,-1,2,-2,3,-3,4,-4,5,-5}; vector<int> gen_values_1 (a_gen_values_1, a_gen_values_1 + sizeof(a_gen_values_1)/ sizeof(int) );
  int a_gen_values_2[] = {1,-1,-2,2,3,-3,4,-4,5,-5}; vector<int> gen_values_2 (a_gen_values_2, a_gen_values_2 + sizeof(a_gen_values_2)/ sizeof(int) );
  int a_gen_values_3[] = {1,-1,4,-4,-3,3,2,-2,5,-5}; vector<int> gen_values_3 (a_gen_values_3, a_gen_values_3 + sizeof(a_gen_values_3)/ sizeof(int) );

  vector< vector<int> > all_values;
  all_values.push_back(gen_values_1); // The routes are inserted shuffled
  all_values.push_back(gen_values_2);
  all_values.push_back(gen_values_3);

  checkRouteScores(all_values);

}

TEST_F(GreedyMoveTest,RouteScoresSortedDiffValuesShuffled) {
  // Three routes are created having the _1 route the best score followed by the _2 and being the _3 the worst route.

  int a_gen_values_1[] = {1,-1,2,-2,3,-3,4,-4,5,-5}; vector<int> gen_values_1 (a_gen_values_1, a_gen_values_1 + sizeof(a_gen_values_1)/ sizeof(int) );
  int a_gen_values_2[] = {1,-1,-2,2,3,-3,4,-4,5,-5}; vector<int> gen_values_2 (a_gen_values_2, a_gen_values_2 + sizeof(a_gen_values_2)/ sizeof(int) );
  int a_gen_values_3[] = {1,-1,4,-4,-3,3,2,-2,5,-5}; vector<int> gen_values_3 (a_gen_values_3, a_gen_values_3 + sizeof(a_gen_values_3)/ sizeof(int) );

  vector< vector<int> > all_values;
  all_values.push_back(gen_values_3); // The routes are inserted shuffled
  all_values.push_back(gen_values_1);
  all_values.push_back(gen_values_2);

  checkRouteScores(all_values);

}

TEST_F(GreedyMoveTest,basicProbsTest) {
  // Three routes are created having the _1 route the best score followed by the _2 and being the _3 the worst route.

  int a_gen_values_1[] = {1,-1,2,-2,3,-3,4,-4,5,-5}; vector<int> gen_values_1 (a_gen_values_1, a_gen_values_1 + sizeof(a_gen_values_1)/ sizeof(int) );
  int a_gen_values_2[] = {1,-1,-2,2,3,-3,4,-4,5,-5}; vector<int> gen_values_2 (a_gen_values_2, a_gen_values_2 + sizeof(a_gen_values_2)/ sizeof(int) );
  int a_gen_values_3[] = {1,-1,4,-4,-3,3,2,-2,5,-5}; vector<int> gen_values_3 (a_gen_values_3, a_gen_values_3 + sizeof(a_gen_values_3)/ sizeof(int) );

  vector< vector<int> > all_values;
  all_values.push_back(gen_values_3); // The routes are inserted shuffled
  all_values.push_back(gen_values_1);
  all_values.push_back(gen_values_2);

  checkProbsFromGen(all_values);

}

TEST_F(GreedyMoveTest,specificValuesProbsTest) {

  int size = 4;
  int pos  = 0;
  vector<DARPVNSShaker::RouteAndScore> values(size);
  int sum = 0;
  for (int i=0; i<size; i++) {
    values[i].score = 3 + i; // 3,4,5,6
    sum += values[i].score;
  }

  vector<long double> probs = DARPVNSShaker::getSelectionProbs(values);
  for (int i=0; i<size; i++) {
    long double prev_value = (i==0) ? 0.0 : probs[i-1];
    ASSERT_GT(0.001,fabs( prev_value+(long double)values[i].score/(long double)sum - probs[i]) );
  }

}

TEST(SelectPos,TestingProbabilities) {
  long double probs[] = {0.10,0.20,0.30,0.48,0.78,0.93,1.0}; vector<long double> probsv (probs, probs + sizeof(probs)/ sizeof(long double) );

  vector<int> freqs(probsv.size(),0);

  const int ATTEMPTS = 1000000;

  for (int i=0; i<ATTEMPTS; i++) {
    int pos = DARPVNSShaker::selectPos(probsv);
    freqs[pos]++;
  }

  for (int i=0; i<freqs.size(); i++) {
    long double estim_prob = (i==0) ? probs[0] : probs[i]-probs[i-1];
    ASSERT_GT(0.01, fabs ( (long double) estim_prob - (long double)freqs[i]/ATTEMPTS) );
  }

}

// Testing the extraction of a sequence with the same size as the shaker
TEST_F(GreedyMoveTest,testingSelectBestSeqToRemoveSameSize) {

  int num_routes     = 3;
  int main_route     = 1;

  int gen_values_a[] = {1,-1,2,-2,3,-3};
  vector<int> gen_values (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int bad_seq_a[] = {-7,7,8,-8};
  list<int> bad_seq (bad_seq_a, bad_seq_a + sizeof(bad_seq_a)/ sizeof(int) );

  int shaker_size = 4;

  // Testing insertion at the beginning
  int insert_pos = 0;
  checkBestSeqToRemove(num_routes, main_route, gen_values, bad_seq, insert_pos, shaker_size);

  // Testing insertion at the middle
  insert_pos = gen_values.size()/2;
  checkBestSeqToRemove(num_routes, main_route, gen_values, bad_seq, insert_pos, shaker_size);

  // Testing insertion at the end
  insert_pos = gen_values.size();
  checkBestSeqToRemove(num_routes, main_route, gen_values, bad_seq, insert_pos, shaker_size);
}

// Testing the extraction of a sequence with the same size as the shaker
TEST_F(GreedyMoveTest,testingSelectBestSeqToRemoveDifSize1) {

  int num_routes     = 3;
  int main_route     = 1;

  int gen_values_a[] = {1,-1,2,-2,3,-3};
  vector<int> gen_values (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int bad_seq_a[] = {-7,7,8,-8};
  list<int> bad_seq (bad_seq_a, bad_seq_a + sizeof(bad_seq_a)/ sizeof(int) );

  int shaker_size = 2;

  // Testing insertion at the beginning
  int insert_pos = 0;
  checkBestSeqToRemove(num_routes, main_route, gen_values, bad_seq, insert_pos, shaker_size);

  // Testing insertion at the middle
  insert_pos = gen_values.size()/2;
  checkBestSeqToRemove(num_routes, main_route, gen_values, bad_seq, insert_pos, shaker_size);

  // Testing insertion at the end
  insert_pos = gen_values.size();
  checkBestSeqToRemove(num_routes, main_route, gen_values, bad_seq, insert_pos, shaker_size);
}

TEST_F(GreedyMoveTest,testingSelectBestSeqToRemoveDifSize2) {

  int num_routes     = 3;
  int main_route     = 1;

  int gen_values_a[] = {1,-1,2,-2,3,-3};
  vector<int> gen_values (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int bad_seq_a[] = {-7,7,8,-8,-9,9};
  list<int> bad_seq (bad_seq_a, bad_seq_a + sizeof(bad_seq_a)/ sizeof(int) );

  int shaker_size = 4;

  // Testing insertion at the beginning
  int insert_pos = 0;
  checkBestSeqToRemove(num_routes, main_route, gen_values, bad_seq, insert_pos, shaker_size);

  // Testing insertion at the middle
  insert_pos = gen_values.size()/2;
  checkBestSeqToRemove(num_routes, main_route, gen_values, bad_seq, insert_pos, shaker_size);

  // Testing insertion at the end
  insert_pos = gen_values.size();
  checkBestSeqToRemove(num_routes, main_route, gen_values, bad_seq, insert_pos, shaker_size);
}

TEST_F(GreedyMoveTest,bestRouteForInsertionTestBasicTests1) {
  int gen_values_a1[] = {1,-1,5,-5,6,-6};
  vector<int> gen_values1 (gen_values_a1, gen_values_a1 + sizeof(gen_values_a1)/ sizeof(int) );

  int gen_values_a2[] = {8,-8};
  vector<int> gen_values2 (gen_values_a2, gen_values_a2 + sizeof(gen_values_a2)/ sizeof(int) );

  int gen_values_a3[] = {7,-7,1,-1};
  vector<int> gen_values3 (gen_values_a3, gen_values_a3 + sizeof(gen_values_a3)/ sizeof(int) );

  vector< vector<int> > all_values;
  all_values.push_back(gen_values1);all_values.push_back(gen_values2); all_values.push_back(gen_values3);

  int ins_seq_a[] = {2,-2,3,-3,4,-4};
  list<int> ins_seq (ins_seq_a, ins_seq_a + sizeof(ins_seq_a)/ sizeof(int) );

  int    route_of_extraction = -1;
  long double orig_score          =  0;
  int    expected_result     =  0;

  checkBestRouteForInsertion(all_values, ins_seq, route_of_extraction, orig_score, expected_result);
}

TEST_F(GreedyMoveTest,bestRouteForInsertionTestBasicTests2) {
  int gen_values_a1[] = {7,-7,1,-1};
  vector<int> gen_values1 (gen_values_a1, gen_values_a1 + sizeof(gen_values_a1)/ sizeof(int) );

  int gen_values_a2[] = {8,-8};
  vector<int> gen_values2 (gen_values_a2, gen_values_a2 + sizeof(gen_values_a2)/ sizeof(int) );

  int gen_values_a3[] = {1,-1,5,-5,6,-6};
  vector<int> gen_values3 (gen_values_a3, gen_values_a3 + sizeof(gen_values_a3)/ sizeof(int) );

  vector< vector<int> > all_values;
  all_values.push_back(gen_values1);all_values.push_back(gen_values2); all_values.push_back(gen_values3);

  int ins_seq_a[] = {2,-2,3,-3,4,-4};
  list<int> ins_seq (ins_seq_a, ins_seq_a + sizeof(ins_seq_a)/ sizeof(int) );

  int    route_of_extraction = -1;
  long double orig_score          =  0;
  int    expected_result     =  2;

  checkBestRouteForInsertion(all_values, ins_seq, route_of_extraction, orig_score, expected_result);
}

TEST_F(GreedyMoveTest,checkStopConditionWithOrigScore) {
  int gen_values_a1[] = {7,-7,1,-1};
  vector<int> gen_values1 (gen_values_a1, gen_values_a1 + sizeof(gen_values_a1)/ sizeof(int) );

  int gen_values_a2[] = {8,-8};
  vector<int> gen_values2 (gen_values_a2, gen_values_a2 + sizeof(gen_values_a2)/ sizeof(int) );

  int gen_values_a3[] = {1,-1,5,-5,6,-6};
  vector<int> gen_values3 (gen_values_a3, gen_values_a3 + sizeof(gen_values_a3)/ sizeof(int) );

  vector< vector<int> > all_values;
  all_values.push_back(gen_values1);all_values.push_back(gen_values2); all_values.push_back(gen_values3);

  int ins_seq_a[] = {2,-2,3,-3,4,-4};
  list<int> ins_seq (ins_seq_a, ins_seq_a + sizeof(ins_seq_a)/ sizeof(int) );

  int    route_of_extraction = -1;
  long double orig_score          =  200;
  int    expected_result     =  0;

  checkBestRouteForInsertion(all_values, ins_seq, route_of_extraction, orig_score, expected_result);
}

TEST_F(GreedyMoveTest,checkStopConditionWithOrigScore2) {
  int gen_values_a1[] = {8,-8};
  vector<int> gen_values1 (gen_values_a1, gen_values_a1 + sizeof(gen_values_a1)/ sizeof(int) );

  int gen_values_a2[] = {7,-7,1,-1};
  vector<int> gen_values2 (gen_values_a2, gen_values_a2 + sizeof(gen_values_a2)/ sizeof(int) );

  int gen_values_a3[] = {1,-1,5,-5,6,-6};
  vector<int> gen_values3 (gen_values_a3, gen_values_a3 + sizeof(gen_values_a3)/ sizeof(int) );

  vector< vector<int> > all_values;
  all_values.push_back(gen_values1);all_values.push_back(gen_values2); all_values.push_back(gen_values3);

  int ins_seq_a[] = {2,-2,3,-3,4,-4};
  list<int> ins_seq (ins_seq_a, ins_seq_a + sizeof(ins_seq_a)/ sizeof(int) );

  int    route_of_extraction = -1;
  long double orig_score          =  20;
  int    expected_result     =  1;

  checkBestRouteForInsertion(all_values, ins_seq, route_of_extraction, orig_score, expected_result);
}

TEST_F(GreedyMoveTest,checkStopConditionWithOrigRoute1) {
  int gen_values_a1[] = {8,-8};
  vector<int> gen_values1 (gen_values_a1, gen_values_a1 + sizeof(gen_values_a1)/ sizeof(int) );

  int gen_values_a2[] = {7,-7,1,-1};
  vector<int> gen_values2 (gen_values_a2, gen_values_a2 + sizeof(gen_values_a2)/ sizeof(int) );

  int gen_values_a3[] = {1,-1,5,-5,6,-6};
  vector<int> gen_values3 (gen_values_a3, gen_values_a3 + sizeof(gen_values_a3)/ sizeof(int) );

  vector< vector<int> > all_values;
  all_values.push_back(gen_values1);all_values.push_back(gen_values2); all_values.push_back(gen_values3);

  int ins_seq_a[] = {2,-2,3,-3,4,-4};
  list<int> ins_seq (ins_seq_a, ins_seq_a + sizeof(ins_seq_a)/ sizeof(int) );

  int    route_of_extraction = 2;
  long double orig_score          = 0;
  int    expected_result     = 1;

  checkBestRouteForInsertion(all_values, ins_seq, route_of_extraction, orig_score, expected_result);
}

TEST_F(GreedyMoveTest,checkStopConditionWithOrigRoute2) {
  int gen_values_a1[] = {8,-8};
  vector<int> gen_values1 (gen_values_a1, gen_values_a1 + sizeof(gen_values_a1)/ sizeof(int) );

  int gen_values_a2[] = {7,-7,1,-1};
  vector<int> gen_values2 (gen_values_a2, gen_values_a2 + sizeof(gen_values_a2)/ sizeof(int) );

  int gen_values_a3[] = {1,-1,5,-5,6,-6};
  vector<int> gen_values3 (gen_values_a3, gen_values_a3 + sizeof(gen_values_a3)/ sizeof(int) );

  vector< vector<int> > all_values;
  all_values.push_back(gen_values1);all_values.push_back(gen_values2); all_values.push_back(gen_values3);

  int ins_seq_a[] = {2,-2,3,-3,4,-4};
  list<int> ins_seq (ins_seq_a, ins_seq_a + sizeof(ins_seq_a)/ sizeof(int) );

  int    route_of_extraction = 2;
  long double orig_score          = 200;
  int    expected_result     = 0;
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  GAGenome::optCriterion(GAGenome::MINIMIZATION);

  setUpTestingEnvironment();

  return RUN_ALL_TESTS();
}


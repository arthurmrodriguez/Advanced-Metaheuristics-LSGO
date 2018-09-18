#include "daropstests_aux.h"

class GreedySwapTest : public ::testing::TestWithParam<int*> {
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
    }

    virtual void TearDown() {
      delete dummymatrix;
      delete requestslist;
    }

    void checkBestAndWorstSeqToRemove(int num_routes, int main_route, vector<int>& values, list<int>& best_seq, list<int>& worst_seq, int ins_pos, int shaker_size) {
      // We insert the bad sequence at the specified pos
      vector<int> all_values = insertSeqAtPos(values,worst_seq,ins_pos);


      DARPGenome* gen = createGen(num_routes,main_route,all_values);


      GreedySwapNeighborhood shaker(shaker_size);

      list<int> best_seq_2del, worst_seq_2del;

      shaker.getBestAndWorstSeqsToRemove(*gen,main_route,main_route,best_seq_2del,worst_seq_2del);

      // First check that the best seq to extract is the worst seq passed as argument
      ASSERT_EQ( GAMin(shaker_size,worst_seq.size()) ,best_seq_2del.size());
      for (list<int>::iterator bestit2del=best_seq_2del.begin(); bestit2del!=best_seq_2del.end(); bestit2del++) {
        ASSERT_NE( worst_seq.end() , find(worst_seq.begin(),worst_seq.end(),*bestit2del));
      }

      // Then check the best seq
      ASSERT_EQ(GAMin(shaker_size,best_seq.size()),worst_seq_2del.size());
      for (list<int>::iterator worstit2del=worst_seq_2del.begin(); worstit2del!=worst_seq_2del.end(); worstit2del++) {

        ASSERT_NE( best_seq.end() , find(best_seq.begin(),best_seq.end(),*worstit2del));
      }

      delete gen;
    }

    vector<int> insertSeqAtPos(vector<int>& values, list<int>& seq, int inspos) {
      assert(inspos>=0);
      assert(inspos<=values.size());

      vector<int> new_values;
      int pos=0;
      for (;pos!=inspos; pos++) new_values.push_back(values[pos]);

      for (list<int>::iterator it=seq.begin(); it!=seq.end(); it++) new_values.push_back( *it );

      for(; pos!=values.size();pos++) new_values.push_back(values[pos]);

      return new_values;
    }

    list<int> getFirstNElems(vector<int>& values, int n) {
      list<int> ext_values;

      for (int i=0; i<n; i++) ext_values.push_back(values[i]);

      return ext_values;
    }

    void checkSwap(vector< vector<int> >& all_values, int route1, list<int>& route1_swap_seq, int route2,list<int>& route2_swap_seq) {

      DARPGenome* gen = createGen(all_values);

      DARPGenome tmp_gen(*gen);

      tmp_gen.swapSeqs(route1,route1_swap_seq,route2,route2_swap_seq);

      ASSERT_EQ(gen->routeLength(route1)-route1_swap_seq.size()+route2_swap_seq.size(),tmp_gen.routeLength(route1));
      ASSERT_EQ(gen->routeLength(route2)-route2_swap_seq.size()+route1_swap_seq.size(),tmp_gen.routeLength(route2));

      for (int i=0; i<tmp_gen.routeLength(route1); i++) {
        int  value                 = tmp_gen.gene(route1,i);
        bool isInOriginalGenRoute1 = gen->findPosOfVertex(route1,value) != -1;
        bool isInRoute2Swap        = find(route2_swap_seq.begin(),route2_swap_seq.end(),value) != route2_swap_seq.end();
        ASSERT_EQ(true,isInOriginalGenRoute1 or isInRoute2Swap);
      }

      for (int i=0; i<tmp_gen.routeLength(route2); i++) {
        int  value                 = tmp_gen.gene(route2,i);
        bool isInOriginalGenRoute2 = gen->findPosOfVertex(route2,value) != -1;
        bool isInRoute1Swap        = find(route1_swap_seq.begin(),route1_swap_seq.end(),value) != route1_swap_seq.end();
        ASSERT_EQ(true,isInOriginalGenRoute2 or isInRoute1Swap);
      }

      delete gen;
    }

};

TEST_F(GreedySwapTest,testingSelectWorstSeqToRemoveDiffSize) {
  int gen_values_a[] = {1,-1,2,-2,3,-3};
  vector<int> gen_values (gen_values_a, gen_values_a + sizeof(gen_values_a)/ sizeof(int) );

  int bad_seq_a[] = {-7,7,-8,8};
  list<int> bad_seq (bad_seq_a, bad_seq_a + sizeof(bad_seq_a)/ sizeof(int) );

  int num_routes  = 3;
  int main_route  = 1;
  int shaker_size = 2;

  int best_seq_a[] = {2,-2};
  list<int> best_seq (best_seq_a, best_seq_a + sizeof(best_seq_a)/ sizeof(int) );
  int ins_pos = 0;
  checkBestAndWorstSeqToRemove(num_routes, main_route, gen_values, best_seq, bad_seq, ins_pos, shaker_size);

  int best_seq2_a[] = {1,-1};
  list<int> best_seq2 (best_seq2_a, best_seq2_a + sizeof(best_seq2_a)/ sizeof(int) );
  ins_pos = 3;
  checkBestAndWorstSeqToRemove(num_routes, main_route, gen_values, best_seq2, bad_seq, ins_pos, shaker_size);

  ins_pos = gen_values.size();
  checkBestAndWorstSeqToRemove(num_routes, main_route, gen_values, best_seq, bad_seq, ins_pos, shaker_size);
}

TEST_F(GreedySwapTest,testingSwapSeqs) {
  int a_gen_values_1[] = { 1,-1,2,-2, 3,-3,4,-4}; vector<int> gen_values_1 (a_gen_values_1, a_gen_values_1 + sizeof(a_gen_values_1)/ sizeof(int) );
  int a_gen_values_2[] = {-7, 7,6,-6,-8, 8};      vector<int> gen_values_2 (a_gen_values_2, a_gen_values_2 + sizeof(a_gen_values_2)/ sizeof(int) );

  vector< vector<int> > all_values;
  all_values.push_back(gen_values_2);
  all_values.push_back(gen_values_1);

  int route1_swqp_seq_a[] = {1,2,-2};
  list<int> route1_swap_seq (route1_swqp_seq_a, route1_swqp_seq_a + sizeof(route1_swqp_seq_a)/ sizeof(int) );

  int route2_swqp_seq_a[] = {7,-8,8};
  list<int> route2_swap_seq (route2_swqp_seq_a, route2_swqp_seq_a + sizeof(route2_swqp_seq_a)/ sizeof(int) );

  int route1 = 1;
  int route2 = 0;

  checkSwap(all_values, route1, route1_swap_seq, route2,route2_swap_seq);
}

TEST_F(GreedySwapTest,testingSwapSeqsInverseRoutes) {
  int a_gen_values_1[] = { 1,-1,2,-2, 3,-3,4,-4}; vector<int> gen_values_1 (a_gen_values_1, a_gen_values_1 + sizeof(a_gen_values_1)/ sizeof(int) );
  int a_gen_values_2[] = {-7, 7,6,-6,-8, 8};      vector<int> gen_values_2 (a_gen_values_2, a_gen_values_2 + sizeof(a_gen_values_2)/ sizeof(int) );

  vector< vector<int> > all_values;
  all_values.push_back(gen_values_2);
  all_values.push_back(gen_values_1);

  int route1_swqp_seq_a[] = {7,-7};
  list<int> route1_swap_seq (route1_swqp_seq_a, route1_swqp_seq_a + sizeof(route1_swqp_seq_a)/ sizeof(int) );

  int route2_swqp_seq_a[] = {1,-1,4,-4};
  list<int> route2_swap_seq (route2_swqp_seq_a, route2_swqp_seq_a + sizeof(route2_swqp_seq_a)/ sizeof(int) );

  int route1 = 0;
  int route2 = 1;

  checkSwap(all_values, route1, route1_swap_seq, route2,route2_swap_seq);
}

TEST_F(GreedySwapTest,testingSwapSeqsFirstLast) {
  int a_gen_values_1[] = { 1,-1,2,-2, 3,-3,4,-4}; vector<int> gen_values_1 (a_gen_values_1, a_gen_values_1 + sizeof(a_gen_values_1)/ sizeof(int) );
  int a_gen_values_2[] = {-7, 7,6,-6,-8, 8};      vector<int> gen_values_2 (a_gen_values_2, a_gen_values_2 + sizeof(a_gen_values_2)/ sizeof(int) );

  vector< vector<int> > all_values;
  all_values.push_back(gen_values_2);
  all_values.push_back(gen_values_1);

  int route1_swqp_seq_a[] = {7,-7,6,-6};
  list<int> route1_swap_seq (route1_swqp_seq_a, route1_swqp_seq_a + sizeof(route1_swqp_seq_a)/ sizeof(int) );

  int route2_swqp_seq_a[] = {4,-4,3,-3};
  list<int> route2_swap_seq (route2_swqp_seq_a, route2_swqp_seq_a + sizeof(route2_swqp_seq_a)/ sizeof(int) );

  int route1 = 0;
  int route2 = 1;

  checkSwap(all_values, route1, route1_swap_seq, route2,route2_swap_seq);
}

TEST_F(GreedySwapTest,testingSwapSeqsLastMiddle) {
  int a_gen_values_1[] = { 1,-1,2,-2, 3,-3,4,-4}; vector<int> gen_values_1 (a_gen_values_1, a_gen_values_1 + sizeof(a_gen_values_1)/ sizeof(int) );
  int a_gen_values_2[] = {-7, 7,6,-6,-8, 8};      vector<int> gen_values_2 (a_gen_values_2, a_gen_values_2 + sizeof(a_gen_values_2)/ sizeof(int) );

  vector< vector<int> > all_values;
  all_values.push_back(gen_values_2);
  all_values.push_back(gen_values_1);

  int route1_swqp_seq_a[] = {7,-7,6,-6};
  list<int> route1_swap_seq (route1_swqp_seq_a, route1_swqp_seq_a + sizeof(route1_swqp_seq_a)/ sizeof(int) );

  int route2_swqp_seq_a[] = {2,-2,3,-3};
  list<int> route2_swap_seq (route2_swqp_seq_a, route2_swqp_seq_a + sizeof(route2_swqp_seq_a)/ sizeof(int) );

  int route1 = 0;
  int route2 = 1;

  checkSwap(all_values, route1, route1_swap_seq, route2,route2_swap_seq);
}


TEST_F(GreedySwapTest,testingSwapSeqsMiddleFirst) {
  int a_gen_values_1[] = { 1,-1,2,-2, 3,-3,4,-4}; vector<int> gen_values_1 (a_gen_values_1, a_gen_values_1 + sizeof(a_gen_values_1)/ sizeof(int) );
  int a_gen_values_2[] = {-7, 7,6,-6,-8, 8};      vector<int> gen_values_2 (a_gen_values_2, a_gen_values_2 + sizeof(a_gen_values_2)/ sizeof(int) );

  vector< vector<int> > all_values;
  all_values.push_back(gen_values_2);
  all_values.push_back(gen_values_1);

  int route1_swqp_seq_a[] = {6,-6};
  list<int> route1_swap_seq (route1_swqp_seq_a, route1_swqp_seq_a + sizeof(route1_swqp_seq_a)/ sizeof(int) );

  int route2_swqp_seq_a[] = {1,-1};
  list<int> route2_swap_seq (route2_swqp_seq_a, route2_swqp_seq_a + sizeof(route2_swqp_seq_a)/ sizeof(int) );

  int route1 = 0;
  int route2 = 1;

  checkSwap(all_values, route1, route1_swap_seq, route2,route2_swap_seq);
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  GAGenome::optCriterion(GAGenome::MINIMIZATION);

  setUpTestingEnvironment();

  return RUN_ALL_TESTS();
}


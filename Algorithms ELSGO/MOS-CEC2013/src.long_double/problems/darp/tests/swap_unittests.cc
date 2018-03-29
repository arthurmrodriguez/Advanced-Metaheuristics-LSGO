#include "daropstests_aux.h"

class SwapTest : public ::testing::TestWithParam<int*> {

public:
  VerticesList* requestslist;
  DummyCostMatrix*  dummymatrix;
  int           num_routes;
  int           route     ;
  DARPGenome*   gen ;

    virtual void SetUp() {
      requestslist = new VerticesList();

      num_routes  = 3;
      route       = 2;
      gen = new DARPGenome(num_routes,0,DummyObjFunc);

      int seq[] = {1,-1,2,-2,3,-3};
      vector<int> seqv (seq, seq + sizeof(seq)/ sizeof(int) );
      for (int i=0; i<seqv.size(); i++) {
        gen->pushBackVertex(route  ,seqv[i]);
        gen->pushBackVertex((route+1)%num_routes,10*seqv[i]);
        gen->pushBackVertex((route+2)%num_routes,100*seqv[i]);
      }

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
      delete gen;
      delete dummymatrix;
      delete requestslist;
    }
};

// For testing swap shaker operator manually
TEST_F(SwapTest,CheckSwapNeighborhood) {
  int neighborhood_size = 3;
//
//  GARandomSeed(5);
//
//  cout << *requestslist << endl;
//
//  cout << "before method call gen is : " << endl << *gen << endl;
//
  long double orig_num_vertices = 0;
  for (int route=0; route<gen->length(); route++) orig_num_vertices += gen->routeLength(route); 
  SwapNeighborhood shaker(neighborhood_size);
  cout << "before gen is " << *gen << endl;
  shaker(*gen);
  cout << "after gen is " << *gen << endl;

  long double after_num_vertices = 0;
  for (int route=0; route<gen->length(); route++) after_num_vertices += gen->routeLength(route); 

  ASSERT_EQ(orig_num_vertices,after_num_vertices);
//  cout << "after method call gen is : " << *gen << endl;

}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  GAGenome::optCriterion(GAGenome::MINIMIZATION);

  setUpTestingEnvironment();

  return RUN_ALL_TESTS();
}


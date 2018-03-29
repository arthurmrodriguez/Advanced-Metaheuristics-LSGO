#include "daropstests_aux.h"

class ChainTest : public ::testing::TestWithParam<int*> {

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
        gen->pushBackVertex(route-1,10*seqv[i]);
        gen->pushBackVertex(route-2,100*seqv[i]);
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

TEST_F(ChainTest,getAllSequences) {
  int neighborhood_size = 3;
//
//  GARandomSeed(11);
//
//  cout << *gen << endl;
  int orig_num_vertices = 0;
  for (int route=0; route<gen->length(); route++) orig_num_vertices += gen->routeLength(route); 

  ChainNeighborhood shaker(neighborhood_size);
  shaker(*gen);

  int after_num_vertices = 0;
  for (int route=0; route<gen->length(); route++) after_num_vertices += gen->routeLength(route); 

  ASSERT_EQ(orig_num_vertices,after_num_vertices);

//  cout << *gen << endl;


//  cout << *gen << endl;
//  cout << "score ofiginal " << gen->score() << endl;
//
//  list<int> seq; 
//  seq.push_back(-1);
//  seq.push_back(2);
//  cout << "scor eof removing seq " << scoreOfRemovingSeq(*gen,route,seq) << endl;
//
//  cout << "score ogi 2 " << gen->score() << endl;
//  gen->removeRequest(route,-1);
//  gen->removeRequest(route,2);
//  cout << "score of removing 2 " << gen->score() << endl;



  //vector<list<int>* > allseqs =  getAllSequences(*gen,route,3);

  //for (int i=0; i<allseqs.size(); i++) {
  //  list<int>* seq = allseqs[i];
  //  cout << "seq " << i << " ";
  //  for (list<int>::iterator it=seq->begin(); it!=seq->end(); it++) {
  //    cout << *it << " ";
  //  }
  //  cout << endl;
  //}


}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  GAGenome::optCriterion(GAGenome::MINIMIZATION);

  setUpTestingEnvironment();

  return RUN_ALL_TESTS();
}


#include "daropstests_aux.h"

class ZeroSplitTest : public ::testing::TestWithParam<int*> {

public:
  VerticesList* requestslist;
  DummyCostMatrix*  dummymatrix;
  int           num_routes;
  int           route     ;
  DARPGenome*   gen ;

    virtual void SetUp() {
      requestslist = new VerticesList();

      num_routes  = 5;
      route       = 2;
      gen = new DARPGenome(num_routes,0,DummyObjFunc);

      int seq[] = {2,1,-1,-2,4,-4,3,5,-5,-3};
      vector<int> seqv (seq, seq + sizeof(seq)/ sizeof(int) );
      for (int i=0; i<seqv.size(); i++) {
        gen->pushBackVertex(route  ,seqv[i]);
        gen->pushBackVertex(route+1,10*seqv[i]);
        gen->pushBackVertex(route+2,100*seqv[i]);
      }

      vector<DummyCostPoint> costinf;
      addRequestRoute(*requestslist, 1,1,2,5,5,  costinf, 15); 
      addRequestRoute(*requestslist, 2,3,4,26,35, costinf,  9);
      addRequestRoute(*requestslist, 3,5,6,45,55, costinf,  5);
      addRequestRoute(*requestslist, 4,2,1,55,65, costinf,  5);
      addRequestRoute(*requestslist, 5,4,3,65,85, costinf,  5);
      addRequestRoute(*requestslist, 6,5,6,95,105, costinf,  5);
      addRequestRoute(*requestslist, 7,4,3,115,120, costinf,  5);

      addRequestRoute(*requestslist, 10,1,4,5,6,  costinf, 15); 
      addRequestRoute(*requestslist, 20,2,5,26,36, costinf,  9);
      addRequestRoute(*requestslist, 30,6,1,45,56, costinf,  5);
      addRequestRoute(*requestslist, 40,2,1,55,65, costinf,  5);
      addRequestRoute(*requestslist, 50,4,3,65,85, costinf,  5);

      addRequestRoute(*requestslist, 100,3,1,5,7,  costinf, 15); 
      addRequestRoute(*requestslist, 200,2,4,26,37, costinf,  9);
      addRequestRoute(*requestslist, 300,1,5,45,57, costinf,  5);
      addRequestRoute(*requestslist, 400,2,1,55,65, costinf,  5);
      addRequestRoute(*requestslist, 500,4,3,65,85, costinf,  5);

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

vector<int> subVector(vector<int>& v, int startpos, int endpos) {
  assert(startpos >=0 && startpos <v.size());
  assert(endpos >0 && endpos <v.size());

  vector<int> subv;
  for (int i=startpos; i<=endpos; i++) {
    subv.push_back(v[i]);
  }

  return subv;
}

bool containTheSameElems(list<int>& l, vector<int>& v) {
  bool result = l.size() == v.size();

  if (result) {
    for (vector<int>::iterator it=v.begin(); it!=v.end(); it++) {
      bool notfound = find(l.begin(),l.end(),*it) == l.end();
      if (notfound) {
        result = false;
        break;
      }
    }
  }

  return result;
}

TEST_F(ZeroSplitTest, NoNatSeq) {
  int num_routes  = 3;
  int route       = 2;

  int seq[] = {1,2,3,-1,-2,-3};
  vector<int> seqv (seq, seq + sizeof(seq)/ sizeof(int) );
  DARPGenome* gen = createGen(num_routes,route,seqv);

  vector< list<int>* > natseqs;
  DARPVNSShaker::getAllNaturalSeqs(*gen,route,natseqs);

  ASSERT_EQ(0,natseqs.size());

  for (int i=0; i<natseqs.size(); i++) delete natseqs[i];
  delete gen;
}

TEST_F(ZeroSplitTest, TwoNatSeqs) {
  int num_routes  = 3;
  int route       = 2;

  int seq[] = {2,1,-1,-2,4,-4,3,5,-5,-3};
  vector<int> seqv (seq, seq + sizeof(seq)/ sizeof(int) );
  DARPGenome* gen = createGen(num_routes,route,seqv);

  vector< list<int>* > natseqs;
  DARPVNSShaker::getAllNaturalSeqs(*gen,route,natseqs);

  vector<int> firstnat  = subVector(seqv,0,3);
  vector<int> secondnat = subVector(seqv,4,5);
  vector<int> thirdnat  = subVector(seqv,6,9);

  ASSERT_EQ(3,natseqs.size());
  ASSERT_TRUE(containTheSameElems(* natseqs[0],firstnat));
  ASSERT_TRUE(containTheSameElems(* natseqs[1],secondnat));


  for (int i=0; i<natseqs.size(); i++) delete natseqs[i];
  delete gen;
}

TEST_F(ZeroSplitTest, TwoNatSeqs2) {
  int num_routes  = 3;
  int route       = 2;

  int seq[] = {2,1,-1,-2,4,-4,3,5,-5,-3,6,-6,7,-7};
  vector<int> seqv (seq, seq + sizeof(seq)/ sizeof(int) );
  DARPGenome* gen = createGen(num_routes,route,seqv);

  vector< list<int>* > natseqs;
  DARPVNSShaker::getAllNaturalSeqs(*gen,route,natseqs);

  vector<int> firstnat  = subVector(seqv,0,3);
  vector<int> secondnat = subVector(seqv,4,5);
  vector<int> thirdnat  = subVector(seqv,6,9);
  vector<int> fourthnat = subVector(seqv,10,11);
  vector<int> fifthnat  = subVector(seqv,12,13);

  ASSERT_EQ(5,natseqs.size());
  ASSERT_TRUE(containTheSameElems(* natseqs[0],firstnat));
  ASSERT_TRUE(containTheSameElems(* natseqs[1],secondnat));
  ASSERT_TRUE(containTheSameElems(* natseqs[2],thirdnat));
  ASSERT_TRUE(containTheSameElems(* natseqs[3],fourthnat));
  ASSERT_TRUE(containTheSameElems(* natseqs[4],fifthnat));


  for (int i=0; i<natseqs.size(); i++) delete natseqs[i];
  delete gen;
}

TEST_F(ZeroSplitTest, ThreeNatSeqs2) {
  int num_routes  = 3;
  int route       = 2;

  int seq[] = {2,1,-1,-2,4,-4,3,5,-5,-3,6,7,-7,-6};
  vector<int> seqv (seq, seq + sizeof(seq)/ sizeof(int) );
  DARPGenome* gen = createGen(num_routes,route,seqv);

  vector< list<int>* > natseqs;

  DARPVNSShaker::getAllNaturalSeqs(*gen,route,natseqs);

  vector<int> firstnat  = subVector(seqv,0,3);
  vector<int> secondnat = subVector(seqv,4,5);
  vector<int> thirdnat  = subVector(seqv,6,9);
  vector<int> fourthnat = subVector(seqv,10,13);

  ASSERT_EQ(4,natseqs.size());
  ASSERT_TRUE(containTheSameElems(* natseqs[0],firstnat));
  ASSERT_TRUE(containTheSameElems(* natseqs[1],secondnat));
  ASSERT_TRUE(containTheSameElems(* natseqs[2],thirdnat));
  ASSERT_TRUE(containTheSameElems(* natseqs[3],fourthnat));


  for (int i=0; i<natseqs.size(); i++) delete natseqs[i];
  delete gen;
}

TEST_F(ZeroSplitTest,BasicTest) {
  int neighborhood_size = 3;

  GARandomSeed(0);

  cout << "Zero split" << endl;
  cout << *gen << endl;

  long double orig_num_vertices = 0;
  for (int route=0; route<gen->length(); route++) orig_num_vertices += gen->routeLength(route); 

  ZeroSplitNeighborhood shaker(0);
  cout << "shaker is " << shaker << endl;
  shaker(*gen);
  cout << endl;
  cout << "after first zero split gen is " << endl;
  cout << *gen << endl;  
  shaker(*gen);

  long double after_num_vertices = 0;
  for (int route=0; route<gen->length(); route++) after_num_vertices += gen->routeLength(route); 

  ASSERT_EQ(orig_num_vertices,after_num_vertices);
//  vector< list<int>* > natseqs;
//  DARPVNSShaker::getAllNaturalSeqs(*gen,0, natseqs);
//
//  for(int i=0; i<natseqs.size(); i++) {
//    cout << "nat seq " << i << ": ";
//    printSeq( * natseqs[i] );
//  }
//  cout << endl;
  
  cout << endl;
  cout << "after zero split gen is" << endl;
  cout << *gen << endl;  
  cout << "fin" << endl;
 
}



int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  GAGenome::optCriterion(GAGenome::MINIMIZATION);

  setUpTestingEnvironment();

  return RUN_ALL_TESTS();
}


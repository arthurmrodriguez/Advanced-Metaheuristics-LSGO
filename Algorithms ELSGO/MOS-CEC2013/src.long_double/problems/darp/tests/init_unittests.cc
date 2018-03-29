#include "daropstests_aux.h"
#include "../darpInit.cc"
#include "../skybusInit/heuristicoInsercion.cc"
#include "../skybusInit/ordenacionPeticiones.cc"
#include "../skybusInit/peticionInsercion.cc"
#include "../skybusInit/vehiculo.cc"

extern "C" long double ReturnZeroObjFunc(GAGenome& g) {

  return 0;
}

class InitializationTest : public ::testing::TestWithParam<int*> {
  static const int MAXPOINTS = 300; 

public:
  VerticesList* requestslist;
  DummyCostMatrix*  dummymatrix;
  DummyDistMatrix*  distmatrix;
  int           num_routes;
  int           route     ;

    virtual void SetUp() {
      num_routes  = 3;
      route       = 0;

      int route_size = 3; 

      requestslist = new VerticesList();

      distmatrix = new DummyDistMatrix(MAXPOINTS);

      // The costs are set so that the optimum sequence is 1 -1 2 -2 3 -3
      // Idem with the 10 and 100 values
      // The remaining combinations are set to a higher cost
      //
      const int MIN_COST_VALUE = 1; 

      vector<int> allposv; 
      for (int i=1; i<=num_routes*route_size*2; i++) allposv.push_back(i);

      for (vector<int>::iterator it=allposv.begin(); it!=allposv.end(); it++) {
        for (vector<int>::iterator itj=it+1; itj!=allposv.end(); itj++) {
          distmatrix->setCost( *it,*itj,MIN_COST_VALUE+1);
          distmatrix->setCost( *itj,*it,MIN_COST_VALUE+1);
        }
      }

      // Optimum sequences
      vector< vector<int> > opt_seqs(3);
      for (int i= 1; i<= 6; i++) opt_seqs[0].push_back(i); 
      for (int i= 7; i<=12; i++) opt_seqs[1].push_back(i); 
      for (int i=13; i<=18; i++) opt_seqs[2].push_back(i); 
      
      for (int pos=0; pos<opt_seqs.size(); pos++) {
        vector<int> opt_seq = opt_seqs[pos];
        
        for (int i=0; i<opt_seq.size()-1; i++) {
          for (int j=i+1; j<opt_seq.size(); j++) {
            distmatrix->setCost( opt_seq[i],  opt_seq[j],MIN_COST_VALUE);
            distmatrix->setCost( opt_seq[j],opt_seq[i],  MIN_COST_VALUE);
          }
        }
      }
      
      addRequest(*requestslist, 1, 1, 2, 5,  *distmatrix);
      addRequest(*requestslist, 2, 3, 4, 15, *distmatrix);
      addRequest(*requestslist, 3, 5, 6, 30, *distmatrix);

      addRequest(*requestslist, 10,  7,  8, 5,  *distmatrix);
      addRequest(*requestslist, 20,  9, 10, 15, *distmatrix);
      addRequest(*requestslist, 30, 11, 12, 30, *distmatrix);

      addRequest(*requestslist, 100, 13, 14, 5,  *distmatrix);
      addRequest(*requestslist, 200, 15, 16, 15, *distmatrix);
      addRequest(*requestslist, 300, 17, 18, 30, *distmatrix);

      dummymatrix = new DummyCostMatrix(distmatrix,*requestslist);
      DARPTestsVars::initData(requestslist,dummymatrix);
    }

    virtual void TearDown() {
      delete dummymatrix;
      delete requestslist;
    }
};


TEST_F(InitializationTest,CheckingEachCriterion) {
  
  Vertex* req1_p = new Vertex(6 ,20,true,  1, 2, 10, Vertex::PICKUP);
  Vertex* req1_d = new Vertex(-6,21,false, 1, 7, 37, Vertex::DELIVERY);

  Vertex* req2_p = new Vertex(7 ,22,true,  1, 2, 10, Vertex::PICKUP);
  Vertex* req2_d = new Vertex(-7,23,false, 1, 7, 37, Vertex::DELIVERY);

  int ori_ori_cost   = 1;
  int ori_dest_cost  = 2;
  int dest_ori_cost  = 3;
  int dest_dest_cost = 4;

  distmatrix->setCost(req1_p->pos_,req2_p->pos_,ori_ori_cost);
  distmatrix->setCost(req1_p->pos_,req2_d->pos_,ori_dest_cost);
  distmatrix->setCost(req1_d->pos_,req2_p->pos_,dest_ori_cost);
  distmatrix->setCost(req1_d->pos_,req2_d->pos_,dest_dest_cost);

  vector< long double (*)(DistMatrix&, Vertex&, Vertex&, Vertex&, Vertex& ) > criteria = constructCriteriaList();
  
  ASSERT_EQ(ori_ori_cost,criteria[0](*distmatrix,*req1_p,*req1_d,*req2_p,*req2_d));
  ASSERT_EQ(ori_dest_cost,criteria[1](*distmatrix,*req1_p,*req1_d,*req2_p,*req2_d));
  ASSERT_EQ(dest_ori_cost,criteria[2](*distmatrix,*req1_p,*req1_d,*req2_p,*req2_d));
  ASSERT_EQ(dest_dest_cost,criteria[3](*distmatrix,*req1_p,*req1_d,*req2_p,*req2_d));

  delete req1_p;
  delete req1_d;
  delete req2_p;
  delete req2_d;
}

// General test
TEST_F(InitializationTest,TestNoImprovement) {
  route = 0;

  DARPGenome* gen = new DARPGenome(num_routes,0,DummyObjFunc);

  VNSInitConstructGenome(*gen,*distmatrix,*requestslist);

  delete gen;
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  GAGenome::optCriterion(GAGenome::MINIMIZATION);

  setUpTestingEnvironment();

  return RUN_ALL_TESTS();
}


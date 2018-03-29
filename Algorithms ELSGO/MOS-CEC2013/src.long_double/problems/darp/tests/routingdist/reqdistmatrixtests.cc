#include "../daropstests_aux.h"

// Hack for testing private members
#define protected public
#define private public
#include "../../ReqDistMatrix.cc"
#include "../../LazyReqDistMatrix.cc"

#undef protected
#undef private

#include <time.h>

using namespace std;


class ReqDistMatrixDummy : public LazyReqDistMatrix {
public:
  ReqDistMatrixDummy(CostMatrix& costmatrix, DARPEvaluator& eval) : LazyReqDistMatrix(costmatrix,eval,true,false,0,0){}

  long double bestScoreFromAllCombs (Request& req1, Request& req2) {
    return req1.delivery_vert->id_;
  }

};

class ReqDistMatrixTesting1 : public ::testing::Test {
public:
  VerticesList*    requestslist;
  DummyCostMatrix* dummymatrix;
  DummyDistMatrix* dummydistmatrix;
  ReqDistMatrix*   reqdistmatrixdummy;

  virtual void SetUp() {
    requestslist = new VerticesList();

    dummydistmatrix = new DummyDistMatrix(MAXPOINTS);

    for (int i=0; i<MAXPOINTS; i++) {
      for (int j=0; j<MAXPOINTS; j++) {
        dummydistmatrix->setCost( i,j,GAGenome::worstPossibleScore()/MAXPOINTS);
        dummydistmatrix->setCost( j,i,GAGenome::worstPossibleScore()/MAXPOINTS);
      }
    }

    int gen_numroutes = 3; // default values, it is not being used actually

    int main_route = 0;
    int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );
    DARPGenome* gen    = createGen(gen_numroutes,  main_route,gen_values);

    dummymatrix = new DummyCostMatrix(dummydistmatrix,*requestslist);

    DARPTestsVars::initData(requestslist, dummymatrix);

    reqdistmatrixdummy = new ReqDistMatrixDummy(*dummymatrix,*DARPTestsVars::darpeval);

    delete gen;
  }

  virtual void TearDown() {
    if (dummymatrix)        delete dummymatrix;
    if (dummydistmatrix)    delete dummydistmatrix;
    if (requestslist)       delete requestslist;
    if (reqdistmatrixdummy) delete reqdistmatrixdummy;
  }

};

TEST_F(ReqDistMatrixTesting1,checkcreateDistMatrix) {
  int id,pos,load;
  bool crit;
  long fb,fe;
  Vertex::VertexType t;

  id = 1; pos = 1; crit = true; load = 1; fb = 5; fe = 10; t = Vertex::PICKUP;
  Vertex v1(id,pos,crit,load,fb,fe,t);

  id = 2; pos = 2; crit = true; load = 1; fb = 5; fe = 10; t = Vertex::PICKUP;
  Vertex v2(id,pos,crit,load,fb,fe,t);

  id = 3; pos = 3; crit = true; load = 1; fb = 5; fe = 10; t = Vertex::DELIVERY;
  Vertex v3(id,pos,crit,load,fb,fe,t);

  id = -1; pos = 4; crit = true; load = 1; fb = 5; fe = 10; t = Vertex::DELIVERY;
  Vertex v4(id,pos,crit,load,fb,fe,t);

  Request r1(&v1,&v1);  // Dummy requests, only the delivery vertex is being used for the computation of the matrix (check bestScoreFromAllCombs from above)
  Request r2(&v2,&v2);
  Request r3(&v3,&v3);

  //Need to read ReqDistMatrixDummy code to understand the tests
  // Basically, the bestScoreFromAllCombs method has been rewritten so that it always returns
  // the id of the first vertex

  ASSERT_EQ(r1.delivery_vert->id_, reqdistmatrixdummy->bestScoreFromAllCombs(r1,r2) );  //Checking the dummy function

  vector<Request> reqs; reqs.push_back(r1); reqs.push_back(r2); reqs.push_back(r3);

  reqdistmatrixdummy->setUpDistMatrix(reqs);
  reqdistmatrixdummy->initializeDistMatrix(reqs);

  ASSERT_EQ(0,reqdistmatrixdummy->data_[0].size());
  for (int i=1; i<reqs.size(); i++) {
    for (int j=0; j<i; j++) {
      ASSERT_EQ( reqs[i].delivery_vert->id_, reqdistmatrixdummy->getDist(i,j) );
    }
  }

  // Then, we check that the reqids correspond to the vertices position in reqs
  ASSERT_EQ(0,reqdistmatrixdummy->reqid2matrixpos_[r1.pickup_vert->id_]);
  ASSERT_EQ(1,reqdistmatrixdummy->reqid2matrixpos_[r2.pickup_vert->id_]);
  ASSERT_EQ(2,reqdistmatrixdummy->reqid2matrixpos_[r3.pickup_vert->id_]);

  // Finally we check that cost from the req id is the same as using the position
  for (int i=1; i<reqs.size(); i++) {
    for (int j=0; j<i; j++) {
      ASSERT_EQ( reqdistmatrixdummy->getDist(i,j), reqdistmatrixdummy->getDistOfReqsPickupIds(reqs[i].pickup_vert->id_,reqs[j].pickup_vert->id_));
    }
  }

}


class ReqDistMatrixTesting2 : public ::testing::Test {
public:
  VerticesList*           verticeslist;
  DummyCostMatrix*        dummymatrix;
  DummyDistMatrix*        dummydistmatrix;
  ReqDistMatrix*          reqdistmatrix;


  virtual void setUpWithParams(vector<Vertex*>& vertices, int costv0v1, int costv1v2, int costv2v3, int costl1, int costl2,
      int vcapac, int sumtimeexp, int maxuserridetime, int maxpickuptime, int twv, int loadv, int ridev,
      DARPEvaluator::darpOptCriterionType optcrit, bool computedeltimes, bool useconstrpen=true, bool usewaitingpen=false, long waitingPenThreshold=0, long double waitingPenConstant=0) {
    dummydistmatrix = new DummyDistMatrix(MAXPOINTS);

    for (int i=0; i<dummydistmatrix->size(); i++) {
      for (int j=0; j<dummydistmatrix->size(); j++) {
        dummydistmatrix->setCost( i,j,GAGenome::worstPossibleScore()/dummydistmatrix->size());
        dummydistmatrix->setCost( j,i,GAGenome::worstPossibleScore()/dummydistmatrix->size());
      }
    }

    dummydistmatrix->setCost(vertices[0]->pos_,vertices[1]->pos_,costv0v1) ;
    dummydistmatrix->setCost(vertices[1]->pos_,vertices[2]->pos_,costv1v2) ;
    dummydistmatrix->setCost(vertices[2]->pos_,vertices[3]->pos_,costv2v3) ;

    vector<int> costs; costs.push_back(costl1); costs.push_back(costl2);
    int costpos=0;
    for (int i=0; i<vertices.size(); i++) {
      if (vertices[i]->isPickUp()) {
        for (int j=i+1;j<vertices.size(); j++) {
          if (vertices[i]->id_ == Vertex::getSiblingVertId( vertices[j]->id_ ) ) {
            // First we set the cost between the vertices of a request
            dummydistmatrix->setCost(vertices[i]->pos_,vertices[j]->pos_,costs[costpos++]);
            assert(costpos < 3);
          }
        }

      }
    }

    verticeslist = new VerticesList();
    for (int i=0; i<vertices.size(); i++) verticeslist->addVertex(vertices[i]);

    dummymatrix = new DummyCostMatrix(dummydistmatrix,*verticeslist);

    DARPTestsVars::initDataOriginalEval(verticeslist,dummymatrix, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit);

    if (computedeltimes) {

      // This code is a duplciated version of the above loop but we have to do it in different steps because here
      // we need the costmatrix and for creating the costmatrix we need to set the missing distances which were set in the
      // previous loop
      for (int i=0; i<vertices.size(); i++) {
        if (vertices[i]->isPickUp()) {
          for (int j=i+1;j<vertices.size(); j++) {
            if (vertices[i]->id_ == Vertex::getSiblingVertId( vertices[j]->id_ ) ) {
              // Then we adjust the times for the delivery vertices based on the values of the pickup vertices
              vertices[j]->fbegin_ = vertices[i]->fend_ + dummymatrix->getCost(*vertices[i],*vertices[j]);
              vertices[j]->fend_   = vertices[i]->fend_ + Vertex::maxRideTime(vertices[i]->id_,*dummymatrix);
            }
          }

        }
      }
    }

    int numroutes = 1;
    DARPGenome* gen = new DARPGenome(numroutes,0,originalobjective);

    reqdistmatrix = new LazyReqDistMatrix(*dummymatrix, *DARPTestsVars::darpeval,useconstrpen,usewaitingpen,waitingPenThreshold,waitingPenConstant);
    delete gen;
  }

  virtual long double computeNewScore(vector<Vertex*>& vertices, int costv0v1, int costv1v2, int costv2v3, int costl1, int costl2,
      int vcapac, int sumtimeexp, int maxuserridetime, int maxpickuptime, int twv, int loadv, int ridev,
      DARPEvaluator::darpOptCriterionType optcrit, bool computedeltimes=true, bool useconstrpen=true, bool usewaitingpen=false, long waitingPenThreshold=0, long double waitingPenConstant=0) {

    setUpWithParams(vertices, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime,
                    twv, loadv, ridev, optcrit,computedeltimes,useconstrpen,usewaitingpen,waitingPenThreshold,waitingPenConstant);

    vector<Request> requests;
    vector<int>     positions;
    for (int i=0; i<vertices.size(); i++) {
      if (vertices[i]->isPickUp()) {
        for (int j=i+1;j<vertices.size(); j++) {
           if (vertices[i]->id_ == Vertex::getSiblingVertId( vertices[j]->id_ ) ) {
             requests.push_back(Request(vertices[i],vertices[j]));
             positions.push_back(i);positions.push_back(j);
          }
        }
      }
    }
    assert(requests.size() == 2);
    assert(positions.size() == 4);
    assert(positions[0] == 0);

    long double new_score = (long double) reqdistmatrix->reqSequenceScore(requests[0], requests[1], positions[1], positions[2], positions[3]);

    return new_score;
  }


  virtual void checkValues(Vertex* v0, Vertex* v1, Vertex* v2, Vertex* v3, int costv0v1, int costv1v2, int costv2v3, int costl1, int costl2,
      int vcapac, int sumtimeexp, int maxuserridetime, int maxpickuptime, int twv, int loadv, int ridev,
      DARPEvaluator::darpOptCriterionType optcrit, bool computedeltimes=true) {
    vector<Vertex*> vertices;
    vertices.push_back(v0); vertices.push_back(v1); vertices.push_back(v2); vertices.push_back(v3);

    long double new_score = computeNewScore(vertices,costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit, computedeltimes);

    int numroutes = 1; // default values, it is not being used actually
    int main_route = 0;
    DARPGenome* gen = new DARPGenome(numroutes,0,originalobjective);
    for (int i=0; i<vertices.size(); i++) gen->pushBackVertex(0,vertices[i]->id_);

    long double old_score = gen->score();

    ASSERT_EQ(new_score,old_score);
  }

  virtual void checkTimes(Vertex* v0, Vertex* v1, Vertex* v2, Vertex* v3, int costv0v1, int costv1v2, int costv2v3, int costl1, int costl2,
      int vcapac, int sumtimeexp, int maxuserridetime, int maxpickuptime, int twv, int loadv, int ridev,
      DARPEvaluator::darpOptCriterionType optcrit) {
    vector<Vertex*> vertices;
    vertices.push_back(v0); vertices.push_back(v1); vertices.push_back(v2); vertices.push_back(v3);

    setUpWithParams(vertices, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit,true);

    vector<Request> requests;
    vector<int>     positions;
    for (int i=0; i<vertices.size(); i++) {
      if (vertices[i]->isPickUp()) {
        for (int j=i+1;j<vertices.size(); j++) {
           if (vertices[i]->id_ == Vertex::getSiblingVertId( vertices[j]->id_ ) ) {
             requests.push_back(Request(vertices[i],vertices[j]));
             positions.push_back(i);positions.push_back(j);
          }
        }
      }
    }
    assert(requests.size() == 2);
    assert(positions.size() == 4);
    assert(positions[0] == 0);

    int starttime = clock();
    for (int i=0; i<10000; i++) {
      long double new_score = reqdistmatrix->reqSequenceScore(requests[0], requests[1], positions[1], positions[2], positions[3]);
    }
    int new_time = clock() - starttime;

    int numroutes = 1; // default values, it is not being used actually
    int main_route = 0;
    DARPGenome* gen = new DARPGenome(numroutes,0,originalobjective);
    for (int i=0; i<vertices.size(); i++) gen->pushBackVertex(0,vertices[i]->id_);

    starttime = clock();
    for (int i=0; i<10000; i++) {
      for (int i=0; i<vertices.size(); i++) gen->gene(0,i,vertices[i]->id_);
      gen->evaluate(gaTrue);
      long double old_score = gen->score();
    }
    long double old_time = clock() - starttime;

    cout << "La diferencia es de " << old_time/new_time << endl;

  }


  virtual void TearDown() {
    if (dummymatrix)     delete dummymatrix;
    if (dummydistmatrix) delete dummydistmatrix;
    if (verticeslist)    delete verticeslist;
  }

};


TEST_F(ReqDistMatrixTesting2,R1P_R1D_R2P_R2D_simple) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=1;  pos=1; critic=true;  load=3; fe=5;  fl=15; type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=2;  pos=3; critic=true;  load=3; fe=20; fl=35; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);

  int costv0v1 = 3; int costv1v2 = 3; int costv2v3 = 3;
  int costl1 = 6 ; int costl2 = 6;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=13; int maxpickuptime=10;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  checkValues(req1p,req1d,req2p,req2d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit);
}

TEST_F(ReqDistMatrixTesting2,R1P_R1D_R2P_R2D_earlyarrivals) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=1;  pos=1; critic=true;  load=3; fe=5;  fl=10; type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=2;  pos=3; critic=true;  load=3; fe=20; fl=24; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);

  int costv0v1 = 3; int costv1v2 = 3; int costv2v3 = 3;
  int costl1 = 6 ; int costl2 = 6;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=11; int maxpickuptime=10;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  checkValues(req1p,req1d,req2p,req2d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit);
}


TEST_F(ReqDistMatrixTesting2,R1P_R1D_R2P_R2D_exactadjust) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=1;  pos=1; critic=true;  load=3; fe=5;  fl=15; type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=2;  pos=3; critic=true;  load=3; fe=20; fl=24; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);

  int costv0v1 = 3; int costv1v2 = 3; int costv2v3 = 3;
  int costl1 = 6 ; int costl2 = 6;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=11; int maxpickuptime=10;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  checkValues(req1p,req1d,req2p,req2d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit);
}

TEST_F(ReqDistMatrixTesting2,R1P_R1D_R2P_R2D_slackadjust) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=1;  pos=1; critic=true;  load=3; fe=5;  fl=15; type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=2;  pos=3; critic=true;  load=3; fe=16; fl=19; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);

  int costv0v1 = 3; int costv1v2 = 3; int costv2v3 = 3;
  int costl1 = 6 ; int costl2 = 6;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=11; int maxpickuptime=10;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  checkValues(req1p,req1d,req2p,req2d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit);
}

TEST_F(ReqDistMatrixTesting2,R1P_R1D_R2P_R2D_TWV_by_one) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=1;  pos=1; critic=true;  load=3; fe=5;  fl=15; type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=2;  pos=3; critic=true;  load=3; fe=16; fl=18; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);

  int costv0v1 = 3; int costv1v2 = 3; int costv2v3 = 3;
  int costl1 = 6 ; int costl2 = 6;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=11; int maxpickuptime=10;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  checkValues(req1p,req1d,req2p,req2d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit);
}

TEST_F(ReqDistMatrixTesting2,R1P_R1D_R2P_R2D_TWV_by_several) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=1;  pos=1; critic=true;  load=3; fe=20;  fl=25; type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=2;  pos=3; critic=true;  load=3; fe=16; fl=18; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);

  int costv0v1 = 3; int costv1v2 = 3; int costv2v3 = 3;
  int costl1 = 6 ; int costl2 = 6;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=11; int maxpickuptime=10;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  checkValues(req1p,req1d,req2p,req2d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit);
}

TEST_F(ReqDistMatrixTesting2,R1P_R2P_R1D_R2D_basic) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=1;  pos=1; critic=true;  load=3; fe=5;  fl=10; type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=2;  pos=3; critic=true;  load=3; fe=16; fl=18; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);

  int costv0v1 = 3; int costv1v2 = 3; int costv2v3 = 3;
  int costl1 = 6 ; int costl2 = 6;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=11; int maxpickuptime=10;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  checkValues(req1p,req2p,req1d,req2d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit);
}

TEST_F(ReqDistMatrixTesting2,R1P_R2P_R1D_R2D_adjusted) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=1;  pos=1; critic=true;  load=3; fe=5;  fl=10; type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=2;  pos=3; critic=true;  load=3; fe=6;  fl=13; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);

  int costv0v1 = 3; int costv1v2 = 13; int costv2v3 = 3;
  int costl1 = 6 ; int costl2 = 6;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=16; int maxpickuptime=10;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  checkValues(req1p,req2p,req1d,req2d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit);
}

TEST_F(ReqDistMatrixTesting2,R1P_R2P_R1D_R2D_TWV_2P) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=1;  pos=1; critic=true;  load=3; fe=5;  fl=10; type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=2;  pos=3; critic=true;  load=3; fe=6;  fl=13; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);

  int costv0v1 = 9; int costv1v2 = 12; int costv2v3 = 3;
  int costl1 = 6 ; int costl2 = 6;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=21; int maxpickuptime=10;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  checkValues(req1p,req2p,req1d,req2d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit);
}

TEST_F(ReqDistMatrixTesting2,R1P_R2P_R1D_R2D_TWV_1D_RTV_1D) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=1;  pos=1; critic=true;  load=3; fe=5;  fl=10; type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=2;  pos=3; critic=true;  load=3; fe=6;  fl=13; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);

  int costv0v1 = 8; int costv1v2 = 19; int costv2v3 = 2;
  int costl1 = 6 ; int costl2 = 6;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=21; int maxpickuptime=10;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  checkValues(req1p,req2p,req1d,req2d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit);
}


TEST_F(ReqDistMatrixTesting2,R1P_R2P_R1D_R2D_TWV_1D_RTV_1D_2D) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=1;  pos=1; critic=true;  load=3; fe=5;  fl=10; type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=2;  pos=3; critic=true;  load=3; fe=6;  fl=13; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);

  int costv0v1 = 8; int costv1v2 = 19; int costv2v3 = 3;
  int costl1 = 6 ; int costl2 = 6;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=21; int maxpickuptime=10;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  checkValues(req1p,req2p,req1d,req2d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit);
}

TEST_F(ReqDistMatrixTesting2,R1P_R2P_R1D_R2D_TWV_2P_1D_2D_RTV_1D_2D_by_several) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=1;  pos=1; critic=true;  load=3; fe=5;  fl=10; type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=2;  pos=3; critic=true;  load=3; fe=6;  fl=13; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);

  int costv0v1 = 50; int costv1v2 = 8; int costv2v3 = 3;
  int costl1 = 6 ; int costl2 = 6;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=21; int maxpickuptime=10;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  checkValues(req1p,req2p,req1d,req2d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit);
}

TEST_F(ReqDistMatrixTesting2,R1P_R2P_R2D_R1D_basic) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=1;  pos=1; critic=true;  load=3; fe=5;  fl=10; type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=2;  pos=3; critic=true;  load=3; fe=6;  fl=13; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);

  int costv0v1 = 3; int costv1v2 = 3; int costv2v3 = 3;
  int costl1 = 6 ; int costl2 = 6;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=21; int maxpickuptime=10;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  checkValues(req1p,req2p,req2d,req1d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit);
}

TEST_F(ReqDistMatrixTesting2,R1P_R2P_R2D_R1D_adjusted) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=1;  pos=1; critic=true;  load=3; fe=5;  fl=10; type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=2;  pos=3; critic=true;  load=3; fe=6;  fl=13; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);

  int costv0v1 = 3; int costv1v2 = 3; int costv2v3 = 3;
  int costl1 = 6 ; int costl2 = 6;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=12; int maxpickuptime=10;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  checkValues(req1p,req2p,req2d,req1d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit);
}

TEST_F(ReqDistMatrixTesting2,R1P_R2P_R2D_R1D_RTV_R1D) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=1;  pos=1; critic=true;  load=3; fe=5;  fl=10; type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=2;  pos=3; critic=true;  load=3; fe=6;  fl=13; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);

  int costv0v1 = 3; int costv1v2 = 3; int costv2v3 = 3;
  int costl1 = 6 ; int costl2 = 6;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=10; int maxpickuptime=10;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  checkValues(req1p,req2p,req2d,req1d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit);
}

TEST_F(ReqDistMatrixTesting2,R1P_R2P_R2D_R1D_TVW_R2P_R2D_R1D_RTV_R1D) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=1;  pos=1; critic=true;  load=3; fe=5;  fl=10; type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=2;  pos=3; critic=true;  load=3; fe=6;  fl=13; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);

  int costv0v1 = 30; int costv1v2 = 3; int costv2v3 = 3;
  int costl1 = 6 ; int costl2 = 6;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=10; int maxpickuptime=10;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  checkValues(req1p,req2p,req2d,req1d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit);
}

TEST_F(ReqDistMatrixTesting2,violations_turnaround_R1_R2) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=2;  pos=3; critic=true;  load=3; fe=6;  fl=13; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=1;  pos=1; critic=true;  load=3; fe=5;  fl=8; type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);

  int costv0v1 = 3; int costv1v2 = 3; int costv2v3 = 3;
  int costl1 = 6 ; int costl2 = 6;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=12; int maxpickuptime=10;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  checkValues(req1p,req2p,req2d,req1d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit);
}

// Based on an error detected on a test of the bernabeu 1000 i1 problem
TEST_F(ReqDistMatrixTesting2,detectedError1) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // General parameters
  int costl1 = 1467 ; int costl2 = 1162;
  int vcapac=18; int sumtimeexp=1; int maxuserridetime=1800; int maxpickuptime=1200;
  int twv=1; int loadv=1; int ridev=1;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;

  Vertex *req1p, *req1d, *req2p, *req2d;
  int costv0v1 ; int costv1v2 ; int costv2v3 ;

  // First case R1+ R1- R2+ R2-
  id=13433;  pos=1; critic=true;  load=1; fe=1317725100;  fl=1317726300; type=Vertex::PICKUP;   req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-13433; pos=2; critic=false; load=1; fe=1317726567;  fl=1317728100; type=Vertex::DELIVERY; req1d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=13419;  pos=3; critic=true;  load=1; fe=1317724200;  fl=1317725400; type=Vertex::PICKUP;   req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-13419; pos=4; critic=false; load=1; fe=1317725362;  fl=1317727200; type=Vertex::DELIVERY; req2d = new Vertex( id,pos,critic,load,fe,fl,type);

  costv0v1 = 1467; costv1v2 = 1495; costv2v3 = 1162;
  checkValues(req1p,req1d,req2p,req2d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit,false);

  // Second case R1+ R2+ R1- R2-
  id=13433;  pos=1; critic=true;  load=1; fe=1317725100;  fl=1317726300; type=Vertex::PICKUP;   req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=13419;  pos=3; critic=true;  load=1; fe=1317724200;  fl=1317725400; type=Vertex::PICKUP;   req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-13433; pos=2; critic=false; load=1; fe=1317726567;  fl=1317728100; type=Vertex::DELIVERY; req1d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-13419; pos=4; critic=false; load=1; fe=1317725362;  fl=1317727200; type=Vertex::DELIVERY; req2d = new Vertex( id,pos,critic,load,fe,fl,type);

  costv0v1 = 0; costv1v2 = 1467; costv2v3 = 702;
  checkValues(req1p,req2p,req1d,req2d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit,false);

  // Third case R1+ R2+ R2- R1-
  id=13433;  pos=1; critic=true;  load=1; fe=1317725100;  fl=1317726300; type=Vertex::PICKUP;   req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=13419;  pos=3; critic=true;  load=1; fe=1317724200;  fl=1317725400; type=Vertex::PICKUP;   req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-13419; pos=4; critic=false; load=1; fe=1317725362;  fl=1317727200; type=Vertex::DELIVERY; req2d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-13433; pos=2; critic=false; load=1; fe=1317726567;  fl=1317728100; type=Vertex::DELIVERY; req1d = new Vertex( id,pos,critic,load,fe,fl,type);

  costv0v1 = 0; costv1v2 = 1162; costv2v3 = 493;
  checkValues(req1p,req2p,req2d,req1d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit,false);


  // Fourth case R2+ R2- R1+ R1-
  id=13419;  pos=3; critic=true;  load=1; fe=1317724200;  fl=1317725400; type=Vertex::PICKUP;   req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-13419; pos=4; critic=false; load=1; fe=1317725362;  fl=1317727200; type=Vertex::DELIVERY; req2d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=13433;  pos=1; critic=true;  load=1; fe=1317725100;  fl=1317726300; type=Vertex::PICKUP;   req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-13433; pos=2; critic=false; load=1; fe=1317726567;  fl=1317728100; type=Vertex::DELIVERY; req1d = new Vertex( id,pos,critic,load,fe,fl,type);

  costv0v1 = 1467; costv1v2 = 1151; costv2v3 = 1467;
  checkValues(req2p,req2d,req1p,req1d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit,false);


  // Fith case R2+ R1+ R2- R1-
  id=13419;  pos=3; critic=true;  load=1; fe=1317724200;  fl=1317725400; type=Vertex::PICKUP;   req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=13433;  pos=1; critic=true;  load=1; fe=1317725100;  fl=1317726300; type=Vertex::PICKUP;   req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-13419; pos=4; critic=false; load=1; fe=1317725362;  fl=1317727200; type=Vertex::DELIVERY; req2d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-13433; pos=2; critic=false; load=1; fe=1317726567;  fl=1317728100; type=Vertex::DELIVERY; req1d = new Vertex( id,pos,critic,load,fe,fl,type);

  costv0v1 = 0; costv1v2 = 1162; costv2v3 = 493;
  int oldcostl1 = costl1; // Ugly but needed to avoid assert error due to the fact that going from req1p to req2d and then to req1d is smaller than going directly from 1p to 1d
  costl1 = 1162;
  checkValues(req2p,req1p,req2d,req1d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit,false);
  costl1 = oldcostl1;


  // Sixth case R2+ R1+ R1- R2-
  id=-13419; pos=4; critic=false; load=1; fe=1317725362;  fl=1317727200; type=Vertex::DELIVERY; req2d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=13433;  pos=1; critic=true;  load=1; fe=1317725100;  fl=1317726300; type=Vertex::PICKUP;   req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-13433; pos=2; critic=false; load=1; fe=1317726567;  fl=1317728100; type=Vertex::DELIVERY; req1d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=13419;  pos=3; critic=true;  load=1; fe=1317724200;  fl=1317725400; type=Vertex::PICKUP;   req2p = new Vertex( id,pos,critic,load,fe,fl,type);

  costv0v1 = 0; costv1v2 = 1467; costv2v3 = 702;
  checkValues(req2p,req1p,req1d,req2d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit,false);

}

TEST_F(ReqDistMatrixTesting2,testingWaitingPenalizationNoPen) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=1;  pos=1; critic=true;  load=3; fe=5;  fl=8;  type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=8;  fl=20; type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=2;  pos=3; critic=true;  load=3; fe=11; fl=16; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=18; fl=28; type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);


  int costv0v1 = 3; int costv1v2 = 3; int costv2v3 = 3;
  int costl1 = 6 ; int costl2 = 6;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=12; int maxpickuptime=20;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  long waitingPenThreshold=2;
  long double waitingPenConstant=200;

  vector<Vertex*> vertices; vertices.push_back(req1p); vertices.push_back(req1d); vertices.push_back(req2p); vertices.push_back(req2d);

  long double withoutwaitpen = computeNewScore(vertices, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp,
                          maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit, false, false);

  long double withwaitpen = computeNewScore(vertices, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp,
                       maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit, false, true, true, waitingPenThreshold, waitingPenConstant);

  ASSERT_EQ(withoutwaitpen,withwaitpen);
}

TEST_F(ReqDistMatrixTesting2,testingWaitingPenalizationMaxNoPen) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=1;  pos=1; critic=true;  load=3; fe=5;  fl=8;  type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=17; fl=20; type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=2;  pos=3; critic=true;  load=3; fe=18; fl=26; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=21; fl=38; type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);


  int costv0v1 = 3; int costv1v2 = 3; int costv2v3 = 3;
  int costl1 = 6 ; int costl2 = 6;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=12; int maxpickuptime=20;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  long waitingPenThreshold=3;
  long double waitingPenConstant=200;

  vector<Vertex*> vertices; vertices.push_back(req1p); vertices.push_back(req1d); vertices.push_back(req2p); vertices.push_back(req2d);

  long double withoutwaitpen = computeNewScore(vertices, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp,
                          maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit, false, false);

  long double withwaitpen = computeNewScore(vertices, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp,
                       maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit, false, true,true, waitingPenThreshold, waitingPenConstant);

  ASSERT_EQ(withoutwaitpen,withwaitpen);
}

TEST_F(ReqDistMatrixTesting2,testingWaitingPenalizationOnePen) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=1;  pos=1; critic=true;  load=3; fe=5;  fl=8;  type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=17; fl=20; type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=2;  pos=3; critic=true;  load=3; fe=18; fl=26; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=21; fl=38; type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);


  int costv0v1 = 3; int costv1v2 = 3; int costv2v3 = 3;
  int costl1 = 6 ; int costl2 = 6;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=12; int maxpickuptime=20;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  long waitingPenThreshold=2;
  long double waitingPenConstant=200;

  vector<Vertex*> vertices; vertices.push_back(req1p); vertices.push_back(req1d); vertices.push_back(req2p); vertices.push_back(req2d);

  long double withoutwaitpen = computeNewScore(vertices, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp,
                          maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit, false, false);

  // Adjusted so that there is only one pen value
  long double withwaitpen = computeNewScore(vertices, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp,
                       maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit, false, true, true, waitingPenThreshold, waitingPenConstant);

  ASSERT_EQ(withoutwaitpen+waitingPenConstant,withwaitpen);
}

TEST_F(ReqDistMatrixTesting2,testingWaitingPenalizationSeveralPen) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=1;  pos=1; critic=true;  load=3; fe=5;  fl=8;  type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=15; fl=20; type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=2;  pos=3; critic=true;  load=3; fe=20; fl=26; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=25; fl=38; type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);


  int costl1 = 5 ; int costl2 = 1;
  int costv0v1 = costl1; int costv1v2 = 2; int costv2v3 = costl2;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=12; int maxpickuptime=20;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  long waitingPenThreshold=2;
  long double waitingPenConstant=200;

  vector<Vertex*> vertices; vertices.push_back(req1p); vertices.push_back(req1d); vertices.push_back(req2p); vertices.push_back(req2d);

  long double withoutwaitpen = computeNewScore(vertices, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp,
                          maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit, false, false);

  // Adjusted so that there is only one pen value for R1-, three for R2+ and two for R2-
  long double withwaitpen = computeNewScore(vertices, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp,
                       maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit, false, true, true, waitingPenThreshold, waitingPenConstant);

  ASSERT_EQ(withoutwaitpen+3*waitingPenConstant,withwaitpen);
}

TEST_F(ReqDistMatrixTesting2,testTimes) {
  int id,pos,fe,fl,load;
  bool critic;
  Vertex::VertexType type;

  // The delivery vertices times are adjusted later
  id=2;  pos=3; critic=true;  load=3; fe=6;  fl=13; type=Vertex::PICKUP;   Vertex* req2p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=1;  pos=1; critic=true;  load=3; fe=5;  fl=8;  type=Vertex::PICKUP;   Vertex* req1p = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-2; pos=4; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req2d = new Vertex( id,pos,critic,load,fe,fl,type);
  id=-1; pos=2; critic=false; load=3; fe=0;  fl=0;  type=Vertex::DELIVERY; Vertex* req1d = new Vertex( id,pos,critic,load,fe,fl,type);

  int costv0v1 = 3; int costv1v2 = 3; int costv2v3 = 3;
  int costl1 = 6 ; int costl2 = 6;

  int vcapac=10; int sumtimeexp=1; int maxuserridetime=12; int maxpickuptime=10;
  int twv=100; int loadv=100; int ridev=100;
  DARPEvaluator::darpOptCriterionType optcrit = DARPEvaluator::BOTHDELAYS;
  checkTimes(req1p,req2p,req2d,req1d, costv0v1, costv1v2, costv2v3, costl1, costl2, vcapac,sumtimeexp, maxuserridetime, maxpickuptime, twv, loadv, ridev, optcrit);
}


int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  GAGenome::optCriterion(GAGenome::MINIMIZATION);

  setUpTestingEnvironment();

  return RUN_ALL_TESTS();
}

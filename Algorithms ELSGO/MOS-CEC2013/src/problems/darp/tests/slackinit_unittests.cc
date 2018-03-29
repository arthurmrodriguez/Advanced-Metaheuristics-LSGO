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

#include "daropstests_aux.h"

float GARandomFloatDummy(float low, float high) {
  return GARandomDoubleDummy(low,high);
}

#define protected public
#define private public
#include "../darpInit.cc"
#undef protected
#undef private

#include "../VNSDARPGenome.cc"

#include "../skybusInit/heuristicoInsercion.cc"
#include "../skybusInit/ordenacionPeticiones.cc"
#include "../skybusInit/peticionInsercion.cc"
#include "../skybusInit/vehiculo.cc"

class Clique : public ::testing::Test {
public:
  static const int dim = 6;

  vector< vector<bool> > matrix;

  virtual void SetUp() {
    matrix.resize(dim);

    for (int i=0; i<dim; i++) {
      matrix[i].resize(dim);
      for (int j=0; j<dim; j++) {
        matrix[i][j] = false;
      }
    }

  }

  Clique& setConnection(int i, int j) {
    assert(i<dim && j<dim);
    matrix[i][j] = matrix[j][i] = true;
    return *this;
  }

  void print() {
    for (int i=0; i<dim; i++) {
      for (int j=0; j<dim; j++) {
        cout << matrix[i][j] << " ";
      }
      cout << endl;
    }
  }

  virtual void TearDown() {
  }

};

TEST_F(Clique,basicExample) {
  for(int i=0;i<dim;i++) {
    vector<bool> marked_nodes_orig = DARPSlackInitC::computeClique(matrix,i);

    for(int j=0;j<dim;j++) {
      if (i==j) ASSERT_TRUE(marked_nodes_orig[j]);
      else      ASSERT_FALSE(marked_nodes_orig[j]);
    }
  }
}

TEST_F(Clique,fullConnectionExample) {
  for(int i=0;i<dim;i++) {
    for(int j=i+1;j<dim;j++) {
      setConnection(i,j);
    }
  }

  for(int i=0;i<60;i++) randomdoublevalues.push_back(0);

  for(int start_node=0;start_node<dim;start_node++) {
    vector<bool> marked_nodes_orig = DARPSlackInitC::computeClique(matrix,start_node);

    for(int j=0;j<dim;j++) ASSERT_TRUE(marked_nodes_orig[j]);
  }
}

TEST_F(Clique,articleExample) {
  setConnection(0,1).setConnection(0,2).setConnection(0,3).setConnection(0,5);
  setConnection(1,2).setConnection(1,3);
  setConnection(2,3).setConnection(2,4);
  setConnection(4,5);


  // This test represents the example depicted in the article
  // We test the same values with two random schemes, the results should be the same
  randomdoublevalues.push_back(0); randomdoublevalues.push_back(0); randomdoublevalues.push_back(0);

  vector<bool> marked_nodes_orig = DARPSlackInitC::computeClique(matrix,0);
  for (int i=0; i<marked_nodes_orig.size(); i++) {
    if (i<4) ASSERT_TRUE(marked_nodes_orig[i]);
    else     ASSERT_FALSE(marked_nodes_orig[i]);
  }

  randomdoublevalues.push_back(0); randomdoublevalues.push_back(0); randomdoublevalues.push_back(0);
  vector<bool> marked_nodes2 = DARPSlackInitC::computeClique(matrix,0);
  for (int i=0; i<marked_nodes_orig.size(); i++) ASSERT_EQ(marked_nodes_orig[i],marked_nodes2[i]);


  // Now we start with the last node, depending on the random value we should have 0 and 5 or 4 and 5

  randomdoublevalues.push_back(0);
  vector<bool>  marked_nodes = DARPSlackInitC::computeClique(matrix,5);
  for (int i=0; i<marked_nodes.size(); i++) {
    if (i>3) ASSERT_TRUE(marked_nodes[i]);
    else     ASSERT_FALSE(marked_nodes[i]);
  }

  randomdoublevalues.push_back(0.6);
  marked_nodes = DARPSlackInitC::computeClique(matrix,5);
  for (int i=0; i<marked_nodes.size(); i++) {
    if (i==0 or i==5) ASSERT_TRUE(marked_nodes[i]);
    else              ASSERT_FALSE(marked_nodes[i]);
  }

  // If we start with any of the most connected nodes, the result should be the same as the first one
  randomdoublevalues.push_back(0); randomdoublevalues.push_back(0); randomdoublevalues.push_back(0);
  marked_nodes2 = DARPSlackInitC::computeClique(matrix,1);
  for (int i=0; i<marked_nodes_orig.size(); i++) ASSERT_EQ(marked_nodes_orig[i],marked_nodes2[i]);
}


class RequestCompatibility : public ::testing::Test {
  static const int MAXPOINTS = 300;

public:
  DummyDistMatrix* distmatrix;
  Vertex *pick1, *pick2, *pick3, *pick4, *del1, *del2, *del3, *del4;
  Request *req1, *req2, *req3, *req4;

  virtual void SetUp() {
    distmatrix = new DummyDistMatrix(MAXPOINTS);

    int  id       = 1;
    int  position = 0;
    bool critical = true;
    int  load     = 1;
    int  fb       = 1;
    int  fe       = 1;

    critical = true;  int pos_p1 = 0; fb = 1; fe = 5;
    pick1 = new Vertex (id++, pos_p1, critical, load, fb, fe, Vertex::PICKUP);
    critical = false; int pos_d1 = 1; fb = 7; fe = 15;
    del1  = new Vertex(-1*id, pos_d1, critical, load, fb, fe, Vertex::DELIVERY);

    critical = true;  int pos_p2 = 2; fb = 9; fe = 13;
    pick2 = new Vertex(id++, pos_p2, critical, load, fb, fe, Vertex::PICKUP);
    critical = false; int pos_d2 = 3; fb = 12; fe = 14;
    del2  = new Vertex(-1*id, pos_d2, critical, load, fb, fe, Vertex::DELIVERY);

    critical = true;  int pos_p3 = 4; fb = 50; fe = 55;
    pick3 = new Vertex(id++, pos_p3, critical, load, fb, fe, Vertex::PICKUP);
    critical = false; int pos_d3 = 5; fb = 70; fe = 75;
    del3  = new Vertex(-1*id, pos_d3, critical, load, fb, fe, Vertex::DELIVERY);

    critical = true;  int pos_p4 = 4; fb = 2; fe = 3;
    pick4 = new Vertex(id++, pos_p4, critical, load, fb, fe, Vertex::PICKUP);
    critical = false; int pos_d4 = 5; fb = 5; fe = 7;
    del4  = new Vertex(-1*id, pos_d4, critical, load, fb, fe, Vertex::DELIVERY);


    req1 = new Request(pick1,del1);
    req2 = new Request(pick2,del2);
    req3 = new Request(pick3,del3);
    req4 = new Request(pick4,del4);
  }

  virtual void TearDown() {
    delete distmatrix;
    delete pick1;
    delete pick2;
    delete pick3;
    delete pick4;
    delete del1;
    delete del2;
    delete del3;
    delete del4;
    delete req1;
    delete req2;
    delete req3;
    delete req4;
  }
};

TEST_F(RequestCompatibility,firstCompCase) {

  // Only the vi vi+n vj vj+n should be valid
  distmatrix->setCostOneSide(pick1->pos_,del1->pos_,3);
  distmatrix->setCostOneSide(pick2->pos_,del2->pos_,3);

  distmatrix->setCostOneSide(del1->pos_, pick2->pos_,3);
  distmatrix->setCostOneSide(pick1->pos_,pick2->pos_,100);
  distmatrix->setCostOneSide(del1->pos_, del2->pos_,100);
  distmatrix->setCostOneSide(del2->pos_, del1->pos_,100);
  distmatrix->setCostOneSide(del2->pos_, pick1->pos_,100);
  distmatrix->setCostOneSide(pick2->pos_,pick1->pos_,100);
  distmatrix->setCostOneSide(pick2->pos_,del1->pos_ ,100);


  ASSERT_FALSE(DARPSlackInitC::areReqIncompatible(*req1,*req2,*distmatrix));
  ASSERT_FALSE(DARPSlackInitC::areReqIncompatible(*req1,*req2,*distmatrix));
}

TEST_F(RequestCompatibility,secondCompCase) {
  // Only the  vi vj vi+n vj+n should be valid
  distmatrix->setCostOneSide(pick1->pos_,del1->pos_,3);
  distmatrix->setCostOneSide(pick2->pos_,del2->pos_,3);

  distmatrix->setCostOneSide(del1->pos_, pick2->pos_,100);
  distmatrix->setCostOneSide(pick1->pos_,pick2->pos_,1);
  distmatrix->setCostOneSide(del1->pos_, del2->pos_ ,1);
  distmatrix->setCostOneSide(del2->pos_, del1->pos_ ,100);
  distmatrix->setCostOneSide(del2->pos_, pick1->pos_,100);
  distmatrix->setCostOneSide(pick2->pos_,pick1->pos_,100);
  distmatrix->setCostOneSide(pick2->pos_,pick1->pos_,100);
  distmatrix->setCostOneSide(pick2->pos_,del1->pos_ ,1);

  ASSERT_FALSE(DARPSlackInitC::areReqIncompatible(*req1,*req2,*distmatrix));
  ASSERT_FALSE(DARPSlackInitC::areReqIncompatible(*req2,*req1,*distmatrix));
}

TEST_F(RequestCompatibility,thirdCompCase) {
  // Only the vi vj vj+n vi+n should be valid
  distmatrix->setCostOneSide(pick1->pos_,del1->pos_,3);
  distmatrix->setCostOneSide(pick2->pos_,del2->pos_,3);

  distmatrix->setCostOneSide(del1->pos_, pick2->pos_,100);
  distmatrix->setCostOneSide(pick1->pos_,pick2->pos_,3);
  distmatrix->setCostOneSide(del1->pos_, del2->pos_ ,100);
  distmatrix->setCostOneSide(del2->pos_, del1->pos_ ,3);
  distmatrix->setCostOneSide(del2->pos_, pick1->pos_,100);
  distmatrix->setCostOneSide(pick2->pos_,pick1->pos_,100);
  distmatrix->setCostOneSide(pick2->pos_,pick1->pos_,100);
  distmatrix->setCostOneSide(pick2->pos_,del1->pos_ ,100);

  ASSERT_FALSE(DARPSlackInitC::areReqIncompatible(*req1,*req2,*distmatrix));
  ASSERT_FALSE(DARPSlackInitC::areReqIncompatible(*req2,*req1,*distmatrix));
}

TEST_F(RequestCompatibility,noCompatibility) {
  // Only the vi vj vj+n vi+n should be valid
  distmatrix->setCostOneSide(pick1->pos_,del1->pos_,3);
  distmatrix->setCostOneSide(pick2->pos_,del2->pos_,3);

  distmatrix->setCostOneSide(del1->pos_, pick2->pos_,100);
  distmatrix->setCostOneSide(pick1->pos_,pick2->pos_,100);
  distmatrix->setCostOneSide(del1->pos_, del2->pos_ ,100);
  distmatrix->setCostOneSide(del2->pos_, del1->pos_ ,100);
  distmatrix->setCostOneSide(del2->pos_, pick1->pos_,100);
  distmatrix->setCostOneSide(pick2->pos_,pick1->pos_,100);
  distmatrix->setCostOneSide(pick2->pos_,pick1->pos_,100);
  distmatrix->setCostOneSide(pick2->pos_,del1->pos_ ,100);

  ASSERT_TRUE(DARPSlackInitC::areReqIncompatible(*req1,*req2,*distmatrix));
}

TEST_F(RequestCompatibility,IncompatMatrixOneIncompatibility) {
  // two requests, same as the first case of incompatibility checking
  distmatrix->setCostOneSide(pick1->pos_,del1->pos_,3);
  distmatrix->setCostOneSide(pick2->pos_,del2->pos_,3);

  distmatrix->setCostOneSide(del1->pos_, pick2->pos_,100);
  distmatrix->setCostOneSide(pick1->pos_,pick2->pos_,100);
  distmatrix->setCostOneSide(del1->pos_, del2->pos_ ,100);
  distmatrix->setCostOneSide(del2->pos_, del1->pos_ ,100);
  distmatrix->setCostOneSide(del2->pos_, pick1->pos_,100);
  distmatrix->setCostOneSide(pick2->pos_,pick1->pos_,100);
  distmatrix->setCostOneSide(pick2->pos_,del1->pos_ ,100);

  vector<Request> requests; requests.push_back(*req1); requests.push_back(*req2);

  vector< vector<bool> > res = DARPSlackInitC::createReqInCompatibilityMatrix(requests,*distmatrix);

  for (int i=0; i<res.size(); i++) {
    for (int j=0; j<res.size(); j++) {
      if (i==j) ASSERT_FALSE(res[i][j]);
      else      ASSERT_TRUE(res[i][j]);
    }
  }
}

TEST_F(RequestCompatibility,IncompatMatrixFourReqAllIncompatibilities) {
  vector<Request> requests;
  requests.push_back(*req1); requests.push_back(*req2); requests.push_back(*req3); requests.push_back(*req4);

  for(int i=0; i<requests.size(); i++) {
    distmatrix->setCostOneSide(requests[i].pickup_vert->pos_,requests[i].delivery_vert->pos_,3);
    for (int j=0; j<requests.size(); j++) {
      if (i!=j) {
        distmatrix->setCostOneSide(requests[i].pickup_vert->pos_  ,requests[j].pickup_vert->pos_  ,100);
        distmatrix->setCostOneSide(requests[i].pickup_vert->pos_  ,requests[j].delivery_vert->pos_,100);
        distmatrix->setCostOneSide(requests[i].delivery_vert->pos_,requests[j].pickup_vert->pos_  ,100);
        distmatrix->setCostOneSide(requests[i].delivery_vert->pos_,requests[j].delivery_vert->pos_,100);
      }
    }
  }

  vector< vector<bool> > res = DARPSlackInitC::createReqInCompatibilityMatrix(requests,*distmatrix);

  ASSERT_EQ(requests.size(),res.size());
  if (requests.size() > 0) ASSERT_EQ(requests.size(),res[0].size());

  for (int i=0; i<res.size(); i++) {
    for (int j=0; j<res.size(); j++) {
      if (i!=j) {
        ASSERT_TRUE(DARPSlackInitC::areReqIncompatible(requests[i],requests[j],*distmatrix));
        ASSERT_TRUE(res[i][j]);
      }
      else ASSERT_FALSE(res[i][j]);
    }
  }
}

TEST_F(RequestCompatibility,IncompatMatrixFourReqNoIncompatibilities) {
  vector<Request> requests;
  requests.push_back(*req1); requests.push_back(*req2); requests.push_back(*req3); requests.push_back(*req4);

  for(int i=0; i<requests.size(); i++) {
    distmatrix->setCostOneSide(requests[i].pickup_vert->pos_,requests[i].delivery_vert->pos_,3);
    for (int j=0; j<requests.size(); j++) {
      if (i!=j) {
        distmatrix->setCostOneSide(requests[i].pickup_vert->pos_  ,requests[j].pickup_vert->pos_  ,1);
        distmatrix->setCostOneSide(requests[i].pickup_vert->pos_  ,requests[j].delivery_vert->pos_,1);
        distmatrix->setCostOneSide(requests[i].delivery_vert->pos_,requests[j].pickup_vert->pos_  ,1);
        distmatrix->setCostOneSide(requests[i].delivery_vert->pos_,requests[j].delivery_vert->pos_,1);
      }
    }
  }

  vector< vector<bool> > res = DARPSlackInitC::createReqInCompatibilityMatrix(requests,*distmatrix);

  ASSERT_EQ(requests.size(),res.size());
  if (requests.size() > 0) ASSERT_EQ(requests.size(),res[0].size());

  for (int i=0; i<res.size(); i++) {
    for (int j=0; j<res.size(); j++) {
      ASSERT_FALSE(DARPSlackInitC::areReqIncompatible(requests[i],requests[j],*distmatrix));
      ASSERT_FALSE(res[i][j]);
    }
  }
}

TEST_F(RequestCompatibility,IncompatMatrixFourReqSomeIncompatibilities) {
  // Incompatibilities between 0, 1, 2 (all connections)
  // Also between 0 and 3
  vector<Request> requests;
  requests.push_back(*req1); requests.push_back(*req2); requests.push_back(*req3); requests.push_back(*req4);

  for(int i=0; i<requests.size(); i++) {
    distmatrix->setCostOneSide(requests[i].pickup_vert->pos_,requests[i].delivery_vert->pos_,3);
    for (int j=0; j<requests.size(); j++) {
      if (i!=j) {
        if ( (i!=3 and j!=3) or (i==0 and j==3) ) {
          distmatrix->setCostOneSide(requests[i].pickup_vert->pos_  ,requests[j].pickup_vert->pos_  ,100);
          distmatrix->setCostOneSide(requests[i].pickup_vert->pos_  ,requests[j].delivery_vert->pos_,100);
          distmatrix->setCostOneSide(requests[i].delivery_vert->pos_,requests[j].pickup_vert->pos_  ,100);
          distmatrix->setCostOneSide(requests[i].delivery_vert->pos_,requests[j].delivery_vert->pos_,100);
        }
      }
    }
  }

  vector< vector<bool> > res = DARPSlackInitC::createReqInCompatibilityMatrix(requests,*distmatrix);

  ASSERT_EQ(requests.size(),res.size());
  if (requests.size() > 0) ASSERT_EQ(requests.size(),res[0].size());

  for (int i=0; i<res.size(); i++) {
    for (int j=0; j<res.size(); j++) {
      if (i!=j) {
        if ( (i!=3 and j!=3) or (i==0 and j==3) ) {
          ASSERT_TRUE( res[i][j] );
        }
      }
    }
  }
}



int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  GAGenome::optCriterion(GAGenome::MINIMIZATION);

  setUpTestingEnvironment();

  return RUN_ALL_TESTS();
}


#include "daropstests_aux.h"

#define CREATEREQUESTVARS(num) Vertex * req##num##_p, *req##num##_d;

#define SETREQUESTVARS(num) req##num##_p = & requestslist->getVertex(num); req##num##_d = & requestslist->getVertex(-1*num);

class CanReqGoBeforeTest : public ::testing::Test {
protected:
  DummyCostMatrix*  dummymatrix;
  VerticesList* requestslist;

  CREATEREQUESTVARS(1);
  CREATEREQUESTVARS(2);
  CREATEREQUESTVARS(3);
  CREATEREQUESTVARS(4);
  CREATEREQUESTVARS(5);

    virtual void SetUp() {
      requestslist = new VerticesList();
      vector<DummyCostPoint> costinf;

      addRequestRoute(*requestslist, 1,0,1,5,10,  costinf, 15); // 1 from position 0 to position 1 cost 15
      addRequestRoute(*requestslist, 2,1,3,26,35, costinf,  9); // 2 from position 1 to position 3 cost  9
      addRequestRoute(*requestslist, 3,3,5,45,55, costinf,  5); // 3 from position 4 to position 5 cost  5
      addRequestRoute(*requestslist, 4,0,2,5,10,  costinf, 22); // 4 from position 0 to position 2 cost 17
      addRequestRoute(*requestslist, 5,2,3,26,26, costinf,  9); // 5 from position 2 to position 3 cost 17

      DummyDistMatrix* distmatrix = new DummyDistMatrix(MAXPOINTS,costinf);
      dummymatrix = new DummyCostMatrix(distmatrix,*requestslist);


      SETREQUESTVARS(1)
      SETREQUESTVARS(2)
      SETREQUESTVARS(3)
      SETREQUESTVARS(4)
      SETREQUESTVARS(5)
    }

    virtual void TearDown() {
      delete dummymatrix;
      delete requestslist;
    }
};


TEST_F(CanReqGoBeforeTest,UsingSameRoute) {
  ASSERT_TRUE ( req1_p->canVertexGoBefore( *req1_d, *dummymatrix ) ) << "req1+ " << *req1_p << " req1- " << *req1_d;
  ASSERT_FALSE( req1_d->canVertexGoBefore( *req1_p, *dummymatrix ) ) << "req1+ " << *req1_p << " req1- " << *req1_d;

  ASSERT_TRUE ( req2_p->canVertexGoBefore( *req2_d, *dummymatrix ) ) << "req2+ " << *req2_p << " req2- " << *req2_d;
  ASSERT_FALSE( req2_d->canVertexGoBefore( *req2_p, *dummymatrix ) ) << "req2+ " << *req2_p << " req2- " << *req2_d;

  ASSERT_TRUE ( req3_p->canVertexGoBefore( *req3_d, *dummymatrix ) ) << "req3+ " << *req3_p << " req3- " << *req3_d;
  ASSERT_FALSE( req3_d->canVertexGoBefore( *req3_p, *dummymatrix ) ) << "req3+ " << *req3_p << " req3- " << *req3_d;
}

TEST_F(CanReqGoBeforeTest,UsingDifferentRoutes) {
  ASSERT_TRUE  ( req1_p->canVertexGoBefore( *req2_p, *dummymatrix ) ) << "req1+ " << *req1_p << " req2+ " << *req2_p;;
  ASSERT_TRUE  ( req1_p->canVertexGoBefore( *req3_p, *dummymatrix ) ) << "req1+ " << *req1_p << " req3+ " << *req3_p;;
  ASSERT_FALSE ( req4_p->canVertexGoBefore( *req5_p, *dummymatrix ) ) << "req4+ " << *req4_p << " req5+ " << *req5_p;;

  ASSERT_TRUE ( req1_d->canVertexGoBefore( *req2_d, *dummymatrix ) ) << "req1- " << *req1_d << " req2+ " << *req2_p;;
  ASSERT_TRUE ( req1_d->canVertexGoBefore( *req3_d, *dummymatrix ) ) << "req1- " << *req1_d << " req3+ " << *req3_p;;

  ASSERT_FALSE ( req2_p->canVertexGoBefore( *req1_p, *dummymatrix ) );
  ASSERT_TRUE  ( req2_p->canVertexGoBefore( *req1_d, *dummymatrix ) );

  ASSERT_FALSE ( req3_p->canVertexGoBefore( *req2_p, *dummymatrix ) );
  ASSERT_FALSE ( req3_d->canVertexGoBefore( *req1_d, *dummymatrix ) );
}

// Testing firstInsertionPos


class ExampleRequests1 : public ::testing::Test {

protected:
  DummyCostMatrix*  dummymatrix;
  VerticesList* requestslist;

    virtual void SetUp() {
      requestslist = new VerticesList();
      vector<DummyCostPoint> costinf;

      requestslist->addVertex( new Vertex( 1,0,true,  1,  2,  8,Vertex::PICKUP));
      requestslist->addVertex( new Vertex(-1,1,false, 1, 21, 23,Vertex::DELIVERY));
      requestslist->addVertex( new Vertex( 2,2,true,  1, 26, 35,Vertex::PICKUP));
      requestslist->addVertex( new Vertex(-2,3,false, 1, 58, 60,Vertex::DELIVERY));
      requestslist->addVertex( new Vertex( 3,4,true,  1, 55, 75,Vertex::PICKUP));
      requestslist->addVertex( new Vertex(-3,5,false, 1, 82, 84,Vertex::DELIVERY));
      requestslist->addVertex( new Vertex( 4,6,true,  1, 88,100,Vertex::PICKUP));
      requestslist->addVertex( new Vertex(-4,7,false, 1, 110,120,Vertex::DELIVERY));

      costinf.push_back( DummyCostPoint(0,1,14) );
      costinf.push_back( DummyCostPoint(1,2, 2) );
      costinf.push_back( DummyCostPoint(2,3,33) );
      costinf.push_back( DummyCostPoint(3,4,15) );
      costinf.push_back( DummyCostPoint(4,5,28) );
      costinf.push_back( DummyCostPoint(5,6,17) );
      costinf.push_back( DummyCostPoint(6,7,31) );

      DummyDistMatrix* distmatrix = new DummyDistMatrix(MAXPOINTS,costinf);
      dummymatrix = new DummyCostMatrix(distmatrix,*requestslist);

      DARPTestsVars::initData(requestslist,dummymatrix);
    }

    virtual void TearDown() {
      delete dummymatrix;
      delete requestslist;
      requestslist = 0;
      dummymatrix = 0;
    }
};

// Hay que refactorizar estos tests que estan un poco chapuzas
TEST_F(ExampleRequests1,InsertionPosTests) {
  int num_routes  = 3;
  int route       = 0;

  int seq[] = {1,-1,2,-2,3,-3,4,-4};
  vector<int> seqv (seq, seq + sizeof(seq)/ sizeof(int) );
  DARPGenome* gen = createGen(num_routes,route,seqv);
  //printRoute(*gen,route,*requestslist,*dummymatrix);
  //cout << *requestslist << endl;

  int start_insert_pos, end_insert_pos;

  int              newreqid   = 5;
  int              newreqdest = 20;
  bool             critic     = true;
  Vertex::VertexType type     = critic ? Vertex::PICKUP : Vertex::DELIVERY;
  int              start      = 3;
  int              end        = 6;

  for (int i=1;i<5; i++) {
    dummymatrix->setCost(i,newreqid,2);
    dummymatrix->setCost(-i,newreqid,2);
  }

  // First we test to insert it before the first request
  Vertex* req = new Vertex(newreqid,newreqdest,critic, 1, start,end, type);
  requestslist->addVertex(req);

  gen->getInsertionPos(route, req->id_,start_insert_pos,end_insert_pos);
  ASSERT_EQ ( start_insert_pos, 0 );
  ASSERT_EQ ( end_insert_pos, 1 );

  // Test that can only go at the first position
  req->fend_ = req->fbegin_ = 3;

  gen->getInsertionPos(route, req->id_,start_insert_pos,end_insert_pos);
  ASSERT_EQ ( start_insert_pos, 0 );
  ASSERT_EQ ( end_insert_pos, 0 );

  // start = 1 end = 2
  req->fbegin_ = 19;
  req->fend_   = 23;

  gen->getInsertionPos(route, req->id_,start_insert_pos,end_insert_pos);
  ASSERT_EQ ( 1,start_insert_pos );
  ASSERT_EQ ( 2,end_insert_pos );

  // start = 0
  // end = last
  req->fbegin_ = 1;
  req->fend_   = 111;

  gen->getInsertionPos(route, req->id_,start_insert_pos,end_insert_pos);
  ASSERT_EQ ( 0,start_insert_pos );
  ASSERT_EQ ( 7,end_insert_pos );

  // start = 0
  // end = next to last
  req->fbegin_ = 1;
  req->fend_   = 112;

  gen->getInsertionPos(route, req->id_,start_insert_pos,end_insert_pos);
  ASSERT_EQ ( 0,start_insert_pos );
  ASSERT_EQ ( 8,end_insert_pos );

  // start = 1
  // end = next to last
  req->fbegin_ = 19;
  req->fend_   = 112;

  gen->getInsertionPos(route, req->id_,start_insert_pos,end_insert_pos);
  ASSERT_EQ ( 1,start_insert_pos );
  ASSERT_EQ ( 8,end_insert_pos );

  // start = 4
  // end = 6
  req->fbegin_ = 73;
  req->fend_   = 89;

  gen->getInsertionPos(route, req->id_,start_insert_pos,end_insert_pos);
  ASSERT_EQ ( 4,start_insert_pos );
  ASSERT_EQ ( 6,end_insert_pos );
}

// Ugly, must be refactored
TEST_F(ExampleRequests1,CheckNonCritInsertions) {
  int num_routes  = 3;
  int route       = 0;

  int seq[] = {1,-1,2,-2,3,-3,4,-4,5};
  vector<int> seqv (seq, seq + sizeof(seq)/ sizeof(int) );
  DARPGenome* gen = createGen(num_routes,route,seqv);
  //printRoute(*gen,route,*requestslist,*dummymatrix);
  //cout << *requestslist << endl;

  int start_insert_pos, end_insert_pos;

  // critical
  int              critic_newreqid   = 5;
  int              critic_newreqdest = 7;
  bool             critic_critic     = true;
  Vertex::VertexType critic_type       = critic_critic ? Vertex::PICKUP : Vertex::DELIVERY;
  int              critic_start      = 90;
  int              critic_end        = 120;

  // Non critical
  int              noncritic_newreqid   = -5;
  int              noncritic_newreqdest = 20;
  bool             noncritic_critic     = false;
  Vertex::VertexType noncritic_type       = noncritic_critic ? Vertex::PICKUP : Vertex::DELIVERY;
  int              noncritic_start      = 92;
  int              noncritic_end        =180;

  for (int i=1;i<5; i++) {
    dummymatrix->setCost(i,noncritic_newreqid,2);
    dummymatrix->setCost(-i,noncritic_newreqid,2);
  }

  dummymatrix->setCost(critic_newreqid,noncritic_newreqid,1);


  Vertex* crit_req       = new Vertex(critic_newreqid,   critic_newreqdest,   critic_critic, 1,   critic_start,   critic_end,    critic_type);
  Vertex* noncritic_req  = new Vertex(noncritic_newreqid,noncritic_newreqdest,noncritic_critic, 1, noncritic_start,noncritic_end, noncritic_type);

  requestslist->addVertex(crit_req);
  requestslist->addVertex(noncritic_req);

  gen->getInsertionPos(route, noncritic_req->id_,start_insert_pos,end_insert_pos);
  ASSERT_EQ ( 9, start_insert_pos );
  ASSERT_EQ ( 9, end_insert_pos );

}


class scoreOfInsertingReqTest : public ::testing::Test {

protected:
  DummyCostMatrix*  dummymatrix;
  VerticesList* requestslist;

    virtual void SetUp() {
      requestslist = new VerticesList();
      vector<DummyCostPoint> costinf;

      addRequestRoute(*requestslist, 1,0,1,5,5,   costinf, 15); 
      addRequestRoute(*requestslist, 2,2,3,26,35, costinf,  9);
      addRequestRoute(*requestslist, 3,4,5,45,55, costinf,  5);
      addRequestRoute(*requestslist, 4,6,7,58,60, costinf,100);

      costinf.push_back( DummyCostPoint(1,2,2) );
      costinf.push_back( DummyCostPoint(3,4,3) );
      costinf.push_back( DummyCostPoint(5,6,4) );

      costinf.push_back( DummyCostPoint(0,20,1) );
      costinf.push_back( DummyCostPoint(1,20,2) );
      costinf.push_back( DummyCostPoint(2,20,3) );
      costinf.push_back( DummyCostPoint(3,20,4) );
      costinf.push_back( DummyCostPoint(4,20,5) );
      costinf.push_back( DummyCostPoint(5,20,6) );

      DummyDistMatrix* distmatrix = new DummyDistMatrix(MAXPOINTS,costinf);
      dummymatrix = new DummyCostMatrix(distmatrix,*requestslist);
      DARPTestsVars::initData(requestslist,dummymatrix);

    }

    virtual void TearDown() {
      delete dummymatrix;
      delete requestslist;
    }
};
  

TEST_F(scoreOfInsertingReqTest,BasicInsertions) {
  int num_routes  = 3;
  int route       = 2;

  DARPGenome* gen = new DARPGenome(num_routes,0,DummyObjFunc);

  int seq[] = {1,-1,2,-2,3,-3,4,-4};
  vector<int> seqv (seq, seq + sizeof(seq)/ sizeof(int) );
  for (int i=0; i<seqv.size(); i++) {
    gen->pushBackVertex(route  ,seqv[i]);
    gen->pushBackVertex(route-1,seqv[seqv.size()-1-i]); // So more than one route is added to the genome
  }
  recomputeWindowsNonCriticalVertex(*requestslist,*gen,route,*dummymatrix);


  printRoute(*gen,route,*requestslist,*dummymatrix);

  // Simple insertion, check that the size is bigger and that an element is added
  int newreqidins = 5;
  // The cost is set for every reqid of the gen
  for (int i=0; i<gen->routeLength(route); i++) dummymatrix->setCost(newreqidins,gen->gene(route,i),0);
  requestslist->addVertex(new Vertex(newreqidins,20,true, 1, 2, 31, Vertex::PICKUP)); // Only the id is meaningful


  int    orig_size  = gen->routeLength(route);
  long double orig_score = gen->score();

  vector<int> insert_positions; 
  insert_positions.push_back(0); insert_positions.push_back(gen->routeLength(route)/2); insert_positions.push_back(gen->routeLength(route)-1); 

  for (int k=0; k<insert_positions.size(); k++) {
    int insert_pos = insert_positions[k];
    gen->insertVertex(route,insert_pos,newreqidins);

    ASSERT_EQ( orig_size +1 , gen->routeLength(route) );                     // Size is greater
    ASSERT_EQ( newreqidins, gen->gene(route,insert_pos) );                      // first element is the corresponding one
    for (int i=0; i<gen->size(); i++) {                                      // Check that last elements are as expected
      if      (i<insert_pos) ASSERT_EQ( seqv.at(i)   , gen->gene(route,i) ); 
      else if (i>insert_pos) ASSERT_EQ( seqv.at(i-1) , gen->gene(route,i) ); 
    }

    //  Check the removal methods, the elements and the score should be the same as the original
    gen->removeVertexOfPos(route, insert_pos);
    ASSERT_EQ( orig_size , gen->routeLength(route) );                                 // Size is greater
    for (int i=0; i<orig_size; i++) ASSERT_EQ( seqv.at(i) , gen->gene(route,i) ); // Check that last elements are as expected

    gen->insertVertex(route,insert_pos,newreqidins);
    gen->removeVertex(route, newreqidins);
    ASSERT_EQ( orig_size , gen->routeLength(route) );                                 // Size is greater
    for (int i=0; i<orig_size; i++) ASSERT_EQ( seqv.at(i) , gen->gene(route,i) ); // Check that last elements are as expected

    // Check the score via the insertScoreMethod

    gen->insertVertex(route,insert_pos,newreqidins);
    long double new_score = gen->score();
    gen->removeVertex(route, newreqidins);
    ASSERT_EQ( new_score, gen->scoreOfInsertingVertex(route,insert_pos,newreqidins) );

    // Check manually the score insertions, removals
    int insert_pos_id = gen->gene(route,insert_pos);

    long double right_addition = dummymatrix->getCost( newreqidins, insert_pos_id);

    new_score = orig_score + right_addition;

    if (insert_pos > 0) {
      int prev_pos_id       = gen->gene(route,insert_pos-1);
      long double left_addition  = dummymatrix->getCost( prev_pos_id, newreqidins  );
      long double removal        = dummymatrix->getCost( prev_pos_id, insert_pos_id );

      new_score += left_addition - removal;
    }

    ASSERT_EQ( new_score  , gen->scoreOfInsertingVertex(route,insert_pos,newreqidins) );
  }


}

class InsertRequestInBestPosTest : public ::testing::Test {

protected:
  DummyCostMatrix*  dummymatrix;
  VerticesList* requestslist;
  DARPGenome*   gen;
  int           num_routes;
  int           route     ;

    virtual void SetUp() {
      requestslist = new VerticesList();

      num_routes  = 3;
      route       = 2;

      gen = new DARPGenome(num_routes,0,DummyObjFunc);

      int seq[] = {1,-1,2,-2};
      vector<int> seqv (seq, seq + sizeof(seq)/ sizeof(int) );
      for (int i=0; i<seqv.size(); i++) {
        gen->pushBackVertex(route  ,seqv[i]);
      }

    }

    virtual void TearDown() {
      delete dummymatrix;
      delete requestslist;
      delete gen;
    }
};

TEST_F(InsertRequestInBestPosTest,SingleRequest) {
  vector<DummyCostPoint> costinf;
  costinf.push_back( DummyCostPoint(0,1,2) );
  costinf.push_back( DummyCostPoint(0,2,3) );
  costinf.push_back( DummyCostPoint(0,3,4) );
  costinf.push_back( DummyCostPoint(1,2,4) );
  costinf.push_back( DummyCostPoint(1,3,4) );
  costinf.push_back( DummyCostPoint(2,3,4) );

  costinf.push_back( DummyCostPoint(0,4,1) );
  costinf.push_back( DummyCostPoint(1,4,20) );
  costinf.push_back( DummyCostPoint(2,4,20) );
  costinf.push_back( DummyCostPoint(3,4,20) );

  costinf.push_back( DummyCostPoint(0,5,20) );
  costinf.push_back( DummyCostPoint(1,5,1) );
  costinf.push_back( DummyCostPoint(2,5,20) );
  costinf.push_back( DummyCostPoint(3,5,20) );
  costinf.push_back( DummyCostPoint(4,5,20) );

  costinf.push_back( DummyCostPoint(0,6,20) );
  costinf.push_back( DummyCostPoint(1,6,20) );
  costinf.push_back( DummyCostPoint(2,6,20) );
  costinf.push_back( DummyCostPoint(3,6, 1) );
  costinf.push_back( DummyCostPoint(4,6,20) );
  costinf.push_back( DummyCostPoint(5,6,20) );


  addRequestRoute(*requestslist, 1,0,1,60,60,   costinf, 10);
  addRequestRoute(*requestslist, 2,2,3,100,100, costinf, 10);

  Vertex* req4 = new Vertex(4,4,true, 1, 2, 180, Vertex::PICKUP);
  requestslist->addVertex( req4 );

  Vertex* req5 = new Vertex(5,5,true, 1, 2, 180, Vertex::PICKUP);
  requestslist->addVertex( req5 );

  Vertex* req6 = new Vertex(6,6,true, 1, 2, 180, Vertex::PICKUP);
  requestslist->addVertex( req6 );



  DummyDistMatrix* distmatrix = new DummyDistMatrix(MAXPOINTS,costinf);
  dummymatrix = new DummyCostMatrix(distmatrix,*requestslist);
  DARPTestsVars::initData(requestslist,dummymatrix);

  recomputeWindowsNonCriticalVertex(*requestslist,*gen,route,*dummymatrix);

  //cout << "gen is " << *gen << endl;
  //cout << "windows:" << endl ;
  //cout << *requestslist << endl << endl;
  //printRoute(*gen,route,*requestslist,*dummymatrix);
  
  gen->insertVertexInBestPos(route,req4->id_);
  ASSERT_EQ(gen->gene(route,0), req4->id_);
  gen->removeVertexOfPos(route,0);

  gen->insertVertexInBestPos(route,req5->id_);
  ASSERT_EQ(req5->id_,gen->gene(route,1) );
  gen->removeVertexOfPos(route,1);
  
  gen->insertVertexInBestPos(route,req6->id_);
  ASSERT_EQ(gen->gene(route,4), req6->id_);
}

class ExtractSiblingsTest : public ::testing::Test {

protected:
  VerticesList* requestslist;
  int           num_routes;
  int           route     ;

    virtual void SetUp() {
      requestslist = new VerticesList();
      num_routes  = 3;
      route       = 2;
    }

    virtual void TearDown() {
      delete requestslist;
    }
};

bool hasReqID(DARPGenome& gen,int route, int reqid){
  bool found = false;
  for (int i=0; i<gen.routeLength(route) && !found; i++){
    if (gen.gene(route,i) == reqid) found = true;
  }
  return found;
}

bool hasReqID(list<int>& l, int reqid){
  bool found = false;
  for (list<int>::iterator it=l.begin(); it!=l.end() && !found; it++) {
    if (*it  == reqid) found = true;
  }
  return found;
}

void print(list<int>& l){
  for (list<int>::iterator it=l.begin(); it!=l.end() ; it++) {
    cout << *it << " ";
  }
  cout << endl;
}

TEST_F(ExtractSiblingsTest,NoSiblingsToExtract) {
  DARPGenome*   gen = new DARPGenome(num_routes,0,DummyObjFunc);

  int seq[] = {1,-1,2,-2,3,-3};
  vector<int> seqv (seq, seq + sizeof(seq)/ sizeof(int) );
  for (int i=0; i<seqv.size(); i++) {
    gen->pushBackVertex(route  ,seqv[i]);
  }

  list<int> requests_ids;
  requests_ids.push_back(4);
  requests_ids.push_back(-4);

  ASSERT_EQ  ( 2, requests_ids.size() );

  CheckAllNaturalSeqsCombsNeighborhood::extractMissingSiblings(*gen,route, requests_ids);

  ASSERT_EQ  ( 2, requests_ids.size() );
}

TEST_F(ExtractSiblingsTest,SimpleCase) {
  DARPGenome*   gen = new DARPGenome(num_routes,0,DummyObjFunc);

  int seq[] = {1,-1,2,-2,3,-3};
  vector<int> seqv (seq, seq + sizeof(seq)/ sizeof(int) );
  for (int i=0; i<seqv.size(); i++) {
    gen->pushBackVertex(route  ,seqv[i]);
  }

  list<int> requests_ids;
  requests_ids.push_back(1);

  ASSERT_FALSE ( hasReqID(requests_ids,-1) );
  ASSERT_TRUE  ( hasReqID(*gen,route,-1) );
  
  CheckAllNaturalSeqsCombsNeighborhood::extractMissingSiblings(*gen,route, requests_ids);

  ASSERT_TRUE  ( hasReqID(requests_ids,-1) );
  ASSERT_FALSE ( hasReqID(*gen,route,-1) );
}

TEST_F(ExtractSiblingsTest,SimpleCase2) {
  DARPGenome*   gen = new DARPGenome(num_routes,0,DummyObjFunc);

  int seq[] = {1,-1,2,-2,3,-3};
  vector<int> seqv (seq, seq + sizeof(seq)/ sizeof(int) );
  for (int i=0; i<seqv.size(); i++) {
    gen->pushBackVertex(route  ,seqv[i]);
  }

  list<int> requests_ids;
  requests_ids.push_back(-1);

  ASSERT_FALSE ( hasReqID(requests_ids,1) );
  ASSERT_TRUE  ( hasReqID(*gen,route,1) );
  
  CheckAllNaturalSeqsCombsNeighborhood::extractMissingSiblings(*gen,route, requests_ids);

  ASSERT_TRUE  ( hasReqID(requests_ids,1) );
  ASSERT_FALSE ( hasReqID(*gen,route,1) );
}

TEST_F(ExtractSiblingsTest,ComplexCase) {
  DARPGenome*   gen = new DARPGenome(num_routes,0,DummyObjFunc);

  int seq[] = {1,-1,2,-2,3,-3};
  vector<int> seqv (seq, seq + sizeof(seq)/ sizeof(int) );
  for (int i=0; i<seqv.size(); i++) {
    gen->pushBackVertex(route  ,seqv[i]);
  }

  list<int> requests_ids;
  requests_ids.push_back(-1);
  requests_ids.push_back(2);
  requests_ids.push_back(-3);

  ASSERT_FALSE ( hasReqID(requests_ids,1) );
  ASSERT_TRUE  ( hasReqID(*gen,route,1) );
  ASSERT_FALSE ( hasReqID(requests_ids,-2) );
  ASSERT_TRUE  ( hasReqID(*gen,route,-2) );
  ASSERT_FALSE ( hasReqID(requests_ids,3) );
  ASSERT_TRUE  ( hasReqID(*gen,route,3) );
  
  CheckAllNaturalSeqsCombsNeighborhood::extractMissingSiblings(*gen,route, requests_ids);

  ASSERT_TRUE ( hasReqID(requests_ids,1) );
  ASSERT_FALSE  ( hasReqID(*gen,route,1) );
  ASSERT_TRUE ( hasReqID(requests_ids,-2) );
  ASSERT_FALSE  ( hasReqID(*gen,route,-2) );
  ASSERT_TRUE ( hasReqID(requests_ids,3) );
  ASSERT_FALSE  ( hasReqID(*gen,route,3) );
}

class ExtractionFromGeneTestSingle : public ::testing::TestWithParam<int*> {

public:
  VerticesList* requestslist;
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
        gen->pushBackVertex(route-1,10*seqv[seqv.size()-1-i]);
        gen->pushBackVertex(route-2,100*seqv[seqv.size()-1-i]);
      }

    }

    virtual void TearDown() {
      delete gen;
      delete requestslist;
    }
};

TEST_P(ExtractionFromGeneTestSingle,ExtractReqSeqSingleElement) {

  int vertex = GetParam()[0];
  int routeLength = GetParam()[1];

  int old_value = gen->gene(route,vertex);
  int old_routeLength = gen->routeLength(route);

  list<int> requests_ids = gen->removeVerticesFrom(route,vertex,routeLength);
  ASSERT_EQ(routeLength,requests_ids.size());
  ASSERT_EQ(old_value,*(requests_ids.begin()));
  ASSERT_EQ(old_routeLength-1,gen->routeLength(route));
}

int extractsingle_firsttestvalues[] = {0,1};
INSTANTIATE_TEST_CASE_P(First,ExtractionFromGeneTestSingle,::testing::Values(extractsingle_firsttestvalues));

int extractsingle_lasttestvalues[] = {5,1};
INSTANTIATE_TEST_CASE_P(Last,ExtractionFromGeneTestSingle,::testing::Values(extractsingle_lasttestvalues));

int extractsingle_middletestvalues[] = {4,1};
INSTANTIATE_TEST_CASE_P(Middle,ExtractionFromGeneTestSingle,::testing::Values(extractsingle_middletestvalues));

class ExtractionFromGeneTestMultiple : public ExtractionFromGeneTestSingle {};

TEST_P(ExtractionFromGeneTestMultiple,ExtractReqSeqFirstTwoElements) {
  int vertex = GetParam()[0];
  int routeLength = GetParam()[1];

  int old_routeLength = gen->routeLength(route);

  vector<int> removed_values;
  for (int pos=0; pos<routeLength; pos++) {
    int value = gen->gene(route,(vertex + pos) % gen->routeLength(route) );
    removed_values.push_back(value);
  }

  list<int> requests_ids = gen->removeVerticesFrom(route,vertex,routeLength);

  ASSERT_EQ(routeLength,requests_ids.size());
  ASSERT_EQ(old_routeLength-routeLength,gen->routeLength(route));
  ASSERT_EQ(removed_values.size(),old_routeLength-gen->routeLength(route));

  // We find every removed value in the obtained list
  for(list<int>::iterator it=requests_ids.begin(); it!=requests_ids.end(); it++) {
    ASSERT_NE(find(removed_values.begin(),removed_values.end(),*it),removed_values.end());
  }

  // We do not find any removed value in the genome
  for(int i=0; i<gen->routeLength(route); i++) {
    ASSERT_EQ(find(removed_values.begin(),removed_values.end(),gen->gene(route,i)),removed_values.end());
  }
}

int extractmultiple_firsttestvalues[] = {0,2};
INSTANTIATE_TEST_CASE_P(First,ExtractionFromGeneTestMultiple,::testing::Values(extractmultiple_firsttestvalues));

int extractmultiple_secondtestvalues[] = {1,2};
INSTANTIATE_TEST_CASE_P(Second,ExtractionFromGeneTestMultiple,::testing::Values(extractmultiple_secondtestvalues));

int extractmultiple_lasttestvalues[] = {4,2};
INSTANTIATE_TEST_CASE_P(Last,ExtractionFromGeneTestMultiple,::testing::Values(extractmultiple_lasttestvalues));

int extractmultiple_middletestvalues[] = {2,2};
INSTANTIATE_TEST_CASE_P(Middle,ExtractionFromGeneTestMultiple,::testing::Values(extractmultiple_middletestvalues));

int extractmultiple_overlaptestvalues1[] = {4,4};
INSTANTIATE_TEST_CASE_P(Overlap1,ExtractionFromGeneTestMultiple,::testing::Values(extractmultiple_overlaptestvalues1));

int extractmultiple_overlaptestvalues2[] = {5,4};
INSTANTIATE_TEST_CASE_P(Overlap2,ExtractionFromGeneTestMultiple,::testing::Values(extractmultiple_overlaptestvalues2));

int extractmultiple_allvalues[] = {0,6};
INSTANTIATE_TEST_CASE_P(All,ExtractionFromGeneTestMultiple,::testing::Values(extractmultiple_allvalues));


int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  GAGenome::optCriterion(GAGenome::MINIMIZATION);

  setUpTestingEnvironment();

  return RUN_ALL_TESTS();
}


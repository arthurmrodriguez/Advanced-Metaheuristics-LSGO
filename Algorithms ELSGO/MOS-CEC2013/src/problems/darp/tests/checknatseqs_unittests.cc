#include "daropstests_aux.h"


class CheckNatSeqsTest : public ::testing::TestWithParam<int*> {
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

      requestslist = new VerticesList();

      DummyDistMatrix* distmatrix = new DummyDistMatrix(MAXPOINTS);

      // The costs are set so that the optimum sequence is 1 -1 2 -2 3 -3
      // The remaining combinations are set to a higher cost
      
      int allpos[] = {1,10,2,20,3,30,4,40,5,50,6,60};
      vector<int> allposv (allpos, allpos + sizeof(allpos)/ sizeof(int) );
      for (vector<int>::iterator it=allposv.begin(); it!=allposv.end(); it++) {
        for (vector<int>::iterator itj=it+1; itj!=allposv.end(); itj++) {
          distmatrix->setCostOneSide( *it,*itj,3);
          distmatrix->setCostOneSide( *itj,*it,3);
        }
      }

      distmatrix->setCostOneSide( 1,10,1);
      distmatrix->setCostOneSide(10, 2,1);
      distmatrix->setCostOneSide( 2,20,1);
      distmatrix->setCostOneSide(20, 3,1);
      distmatrix->setCostOneSide( 3,30,1);

      // wend = wbegin + MAXUSERWAITTIME (20)
      // delivery vertex wbegin = wbegin + cost (origin dest), wend = wend + MAXUSERRIDETIME (60)

      addRequest(*requestslist, 1, 1, 10,  5,  *distmatrix);
      addRequest(*requestslist, 2, 2, 20, 15, *distmatrix);
      addRequest(*requestslist, 3, 3, 30, 30, *distmatrix);
      addRequest(*requestslist, 4, 4, 40, 45, *distmatrix);
      addRequest(*requestslist, 5, 5, 50, 60, *distmatrix);

      addRequest(*requestslist, 10, 1, 10,  5, *distmatrix);
      addRequest(*requestslist, 20, 2, 20, 15, *distmatrix);
      addRequest(*requestslist, 30, 3, 30, 30, *distmatrix);
      addRequest(*requestslist, 40, 4, 40, 45, *distmatrix);
      addRequest(*requestslist, 50, 5, 50, 60, *distmatrix);

      addRequest(*requestslist, 100, 1, 10,  5, *distmatrix);
      addRequest(*requestslist, 200, 2, 20, 15, *distmatrix);
      addRequest(*requestslist, 300, 3, 30, 30, *distmatrix);
      addRequest(*requestslist, 400, 3, 30, 45, *distmatrix);
      addRequest(*requestslist, 500, 3, 30, 60, *distmatrix);

      dummymatrix = new DummyCostMatrix(distmatrix,*requestslist);
      DARPTestsVars::initData(requestslist,dummymatrix);
    }

    virtual void TearDown() {
      delete dummymatrix;
      delete requestslist;
    }

    vector< list<int>* > createNatSeqs(DARPGenome& gen, int route, vector<int>& cut_positions) {
      assert(cut_positions.size() > 0);
      assert(route >=0 && route < gen.size());

      if (cut_positions[cut_positions.size()-1] != gen.routeLength(route)) cut_positions.push_back(gen.routeLength(route));

      vector< list<int>* > natseqs;

      int genpos = 0;
      for (vector<int>::iterator it=cut_positions.begin(); it!=cut_positions.end() && genpos < gen.routeLength(route); it++) {

        int        cutpos = *it; assert(cutpos > 0 && cutpos <= gen.routeLength(route));
        list<int> *natseq = new list<int>;

        for (;genpos<cutpos; genpos++) natseq->push_back(gen.gene(route,genpos));

        natseqs.push_back(natseq);
      }
      return natseqs;
    }
};



TEST_F(CheckNatSeqsTest,getStartingPosInGenTest) {

  list<int> natseq1; natseq1.push_back(1); natseq1.push_back(2); natseq1.push_back(3);
  list<int> natseq2; natseq2.push_back(4); 
  list<int> natseq3; natseq3.push_back(4); natseq3.push_back(5); 

  vector< list<int>* > natseqs1;
  natseqs1.push_back(&natseq1); natseqs1.push_back(&natseq2); natseqs1.push_back(&natseq3);
  ASSERT_EQ( 0,CheckAllNaturalSeqsCombsNeighborhood::computeInsPosBasedOnTheNatSeqPos(natseqs1,0));
  ASSERT_EQ( 3,CheckAllNaturalSeqsCombsNeighborhood::computeInsPosBasedOnTheNatSeqPos(natseqs1,1));
  ASSERT_EQ( 4,CheckAllNaturalSeqsCombsNeighborhood::computeInsPosBasedOnTheNatSeqPos(natseqs1,2));
  ASSERT_EQ(-1,CheckAllNaturalSeqsCombsNeighborhood::computeInsPosBasedOnTheNatSeqPos(natseqs1,3));

  vector< list<int>* > natseqs2;
  natseqs2.push_back(&natseq2); natseqs2.push_back(&natseq1); natseqs2.push_back(&natseq3);
  ASSERT_EQ( 0,CheckAllNaturalSeqsCombsNeighborhood::computeInsPosBasedOnTheNatSeqPos(natseqs2,0));
  ASSERT_EQ( 1,CheckAllNaturalSeqsCombsNeighborhood::computeInsPosBasedOnTheNatSeqPos(natseqs2,1));
  ASSERT_EQ( 4,CheckAllNaturalSeqsCombsNeighborhood::computeInsPosBasedOnTheNatSeqPos(natseqs2,2));
  ASSERT_EQ(-1,CheckAllNaturalSeqsCombsNeighborhood::computeInsPosBasedOnTheNatSeqPos(natseqs2,3));

  vector< list<int>* > natseqs3;
  natseqs3.push_back(&natseq2); 
  ASSERT_EQ( 0,CheckAllNaturalSeqsCombsNeighborhood::computeInsPosBasedOnTheNatSeqPos(natseqs3,0));
  ASSERT_EQ(-1,CheckAllNaturalSeqsCombsNeighborhood::computeInsPosBasedOnTheNatSeqPos(natseqs3,-1));
  ASSERT_EQ(-1,CheckAllNaturalSeqsCombsNeighborhood::computeInsPosBasedOnTheNatSeqPos(natseqs3,1));
  ASSERT_EQ(-1,CheckAllNaturalSeqsCombsNeighborhood::computeInsPosBasedOnTheNatSeqPos(natseqs3,2));
  ASSERT_EQ(-1,CheckAllNaturalSeqsCombsNeighborhood::computeInsPosBasedOnTheNatSeqPos(natseqs3,3));
}

TEST_F(CheckNatSeqsTest,getPrevAndNextNatSeqsTests) {
  int num_routes   = 3;

  int a_gen_values[] = {1,-1,2,-2,3,-3,4,-4,5,-5}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );
  DARPGenome* gen = createGen(num_routes,route,gen_values);

  vector< vector< list<int>* > > natseqs; 
  for (int i=0; i<num_routes; i++) {
    int cut_pos[] = {gen->routeLength(i)/3, gen->routeLength(i)*2/3}; 
    vector<int> cut_positions(cut_pos, cut_pos+sizeof(cut_pos)/sizeof(int) );

    natseqs.push_back( createNatSeqs(*gen,i,cut_positions) );
  }

  int        route       = 0;
  int        route_pos   = 0; // First natural sequence
  list<int> *prev_natseq = 0;
  list<int> *next_natseq = 0;

  CheckAllNaturalSeqsCombsNeighborhood::getPrevAndNextNatSeqs(natseqs, route, route_pos, prev_natseq, next_natseq);
  ASSERT_EQ(0                          ,prev_natseq);
  ASSERT_EQ(natseqs[route][route_pos+1],next_natseq);

  route     = 2;
  route_pos = natseqs[route].size() -1; //Last natural sequence 

  CheckAllNaturalSeqsCombsNeighborhood::getPrevAndNextNatSeqs(natseqs, route, route_pos, prev_natseq, next_natseq);
  ASSERT_EQ(natseqs[route][route_pos-1],prev_natseq);
  ASSERT_EQ(0,                          next_natseq);

  route     = 1;
  route_pos = natseqs[route].size()/ 2; //Middle natural sequence 

  CheckAllNaturalSeqsCombsNeighborhood::getPrevAndNextNatSeqs(natseqs, route, route_pos, prev_natseq, next_natseq);
  ASSERT_EQ(natseqs[route][route_pos-1],prev_natseq);
  ASSERT_EQ(natseqs[route][route_pos+1],next_natseq);

  for (int i=0; i<num_routes; i++) {
    vector< list<int>* >& route_natseqs = natseqs[i];
    for (int j=0; j<route_natseqs.size(); j++) delete route_natseqs[j];
  }
}

class InsertVerticesTest : public CheckNatSeqsTest {
  public:
  void checkInsertion(int num_routes, int route, int pos, vector<int>& gen_values, list<int>& ins_values) {
  // Testing this function void DARPGenome::insertVertices (int route, int pos, list<int>& values) {
    DARPGenome* gen = createGen(num_routes,route,gen_values);
    DARPGenome  copy_gen (*gen);

    int old_size = gen->routeLength(route);

    gen->insertVertices(route,pos,ins_values);

    ASSERT_EQ( old_size + ins_values.size(), gen->routeLength(route) );

    // Check the values that are before the insertion point
    for (int i=0; i<pos; i++) ASSERT_EQ(copy_gen.gene(route,i), gen->gene(route,i) );

    // Check the values inserted
    list<int>::iterator it = ins_values.begin();
    for (int i=pos; i<pos+ins_values.size(); i++,it++) {
      ASSERT_EQ( gen->gene(route,i), *it);
    }

    // Check the values after the insertion
    for (int i=pos+ins_values.size(); i<gen->routeLength(route); i++) {
      ASSERT_EQ( copy_gen.gene(route,i-ins_values.size()), gen->gene(route,i) );
    }

    delete gen;
  }
};

TEST_F(InsertVerticesTest,firstPos) {
  int num_routes   = 3;
  int route        = 0;
  int pos          = 0;

  int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );

  int a_ins_values[] = {4,5,6}; list<int> ins_values (a_ins_values, a_ins_values + sizeof(a_ins_values)/ sizeof(int) );

  checkInsertion(num_routes,route,pos,gen_values,ins_values);
}

TEST_F(InsertVerticesTest,lastPos) {
  int num_routes   = 3;
  int route        = 1;

  int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );

  int a_ins_values[] = {4,5,6}; list<int> ins_values (a_ins_values, a_ins_values + sizeof(a_ins_values)/ sizeof(int) );

  int pos = gen_values.size();

  checkInsertion(num_routes,route,pos,gen_values,ins_values);
}

TEST_F(InsertVerticesTest,middlePos) {
  int num_routes   = 3;
  int route        = 1;

  int a_gen_values[] = {1,-1,2,-2,3,-3};
  vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );

  int a_ins_values[] = {4,5,6};
  list<int> ins_values (a_ins_values, a_ins_values + sizeof(a_ins_values)/ sizeof(int) );

  int pos = gen_values.size()/2;

  checkInsertion(num_routes,route,pos,gen_values,ins_values);
}

TEST_F(InsertVerticesTest,OneValueMiddlePos) {
  int num_routes   = 3;
  int route        = 1;

  int a_gen_values[] = {1,-1,2,-2,3,-3};
  vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );

  int a_ins_values[] = {4};
  list<int> ins_values (a_ins_values, a_ins_values + sizeof(a_ins_values)/ sizeof(int) );

  int pos = gen_values.size()/2;

  checkInsertion(num_routes,route,pos,gen_values,ins_values);
}

TEST_F(InsertVerticesTest,GenOneValueInsInFirst) {
  int num_routes   = 3;
  int route        = 1;

  int a_gen_values[] = {1};
  vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );

  int a_ins_values[] = {4,5,6,7};
  list<int> ins_values (a_ins_values, a_ins_values + sizeof(a_ins_values)/ sizeof(int) );

  int pos = 0;

  checkInsertion(num_routes,route,pos,gen_values,ins_values);
}

TEST_F(InsertVerticesTest,GenOneValueInsInLast) {
  int num_routes   = 3;
  int route        = 1;

  int a_gen_values[] = {1};
  vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );

  int a_ins_values[] = {4,5,6,7};
  list<int> ins_values (a_ins_values, a_ins_values + sizeof(a_ins_values)/ sizeof(int) );

  int pos = 1;

  checkInsertion(num_routes,route,pos,gen_values,ins_values);
}

class RemoveVerticesTest : public CheckNatSeqsTest {
  public:
  void checkRemoval(int num_routes, int route, int pos, int length, vector<int>& gen_values) {
  // Testing this function list<int> DARPGenome::removeVerticesFrom (int route, int pos, int length) {
    DARPGenome* gen = createGen(num_routes,route,gen_values);
    DARPGenome  copy_gen (*gen);

    int old_size = gen->routeLength(route);

    gen->removeVerticesFrom(route,pos,length);

    ASSERT_EQ( old_size - length, gen->routeLength(route) );

    // Check the values (if any) that are before the removal point
    for (int i=0; i<pos; i++) ASSERT_EQ(copy_gen.gene(route,i), gen->gene(route,i) );

    if (pos > 0) {
      // Check the values (if any) after the removed vertices 
      for (int i=pos; i<gen->size(); i++) {
        ASSERT_EQ( gen->gene(route,i), copy_gen.gene(route,i-length) );
      }
    }

    delete gen;
  }
};

TEST_F(RemoveVerticesTest,removeFirst2Values) {
  int num_routes   = 3;
  int route        = 0;
  int pos          = 0;
  int length       = 2;

  int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );

  checkRemoval(num_routes,route,pos,length,gen_values);
}

TEST_F(RemoveVerticesTest,removeTwoMiddleValues) {
  int num_routes   = 3;
  int route        = 0;

  int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );

  int pos          = gen_values.size()/2;;
  int length       = 2;

  checkRemoval(num_routes,route,pos,length,gen_values);
}

TEST_F(RemoveVerticesTest,removeLastTwoValues) {
  int num_routes   = 3;
  int route        = 0;

  int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );

  int length       = 2;
  int pos          = gen_values.size() - length;

  checkRemoval(num_routes,route,pos,length,gen_values);
}

TEST_F(RemoveVerticesTest,removeAllValues) {
  int num_routes   = 3;
  int route        = 0;

  int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );

  int length = gen_values.size();
  int pos    = 0;

  checkRemoval(num_routes,route,pos,length,gen_values);
}

TEST_F(RemoveVerticesTest,removeFirstValue) {
  int num_routes   = 3;
  int route        = 0;

  int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );

  int length = 1;
  int pos    = 0;

  checkRemoval(num_routes,route,pos,length,gen_values);
}

TEST_F(RemoveVerticesTest,removeMiddleValue) {
  int num_routes   = 3;
  int route        = 0;

  int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );

  int length = 1;
  int pos    = gen_values.size()/2;

  checkRemoval(num_routes,route,pos,length,gen_values);
}

TEST_F(RemoveVerticesTest,removeLastValue) {
  int num_routes   = 3;
  int route        = 0;

  int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );

  int length = 1;
  int pos    = gen_values.size()-1;

  checkRemoval(num_routes,route,pos,length,gen_values);
}

class SwapNatSeqsTest : public CheckNatSeqsTest {
  public:

    void checkSwap(DARPGenome& gen, vector< vector< list<int>* > >& natseqs, int fromroute, int from_natseq, int toroute, int to_natseq) {
      
      DARPGenome  copy_gen (gen);

      int old_size = gen.routeLength(route);

      list<int>* orig_from_natseq = natseqs[fromroute][from_natseq];
      list<int>* orig_to_natseq   = natseqs[toroute][to_natseq];

      vector< list<int>* >& fromroute_natseqs = natseqs[fromroute];
      int                   from_natseq_size  = fromroute_natseqs[from_natseq]->size();
      int                   from_startpos     = CheckAllNaturalSeqsCombsNeighborhood::computeInsPosBasedOnTheNatSeqPos(fromroute_natseqs,from_natseq);

      vector< list<int>* >& toroute_natseqs = natseqs[toroute];
      int                   to_natseq_size  = toroute_natseqs[to_natseq]->size();
      int                   to_startpos     = CheckAllNaturalSeqsCombsNeighborhood::computeInsPosBasedOnTheNatSeqPos(toroute_natseqs,to_natseq);

      CheckAllNaturalSeqsCombsNeighborhood::swapNaturalSequences(gen,natseqs,fromroute,from_natseq,toroute,to_natseq);

      //cout << "swap from route " << fromroute << " nat seq "; printNatSeq(fromroute_natseqs,from_natseq); 
      //cout << "to route " << toroute << " nat seq "; printNatSeq(toroute_natseqs,to_natseq); cout << endl << endl;
      //cout << "Original gen is " << copy_gen << endl << endl;
      //cout << "Swap gen is " << gen << endl << endl;

      for (int route=0; route<gen.size(); route++) {

        if (route == fromroute) {

          for (int i=0                             ; i<from_startpos;           i++) ASSERT_EQ(copy_gen.gene(fromroute,i),                                gen.gene(route,i));
          for (int i=from_startpos                 ; i<to_natseq_size;          i++) ASSERT_EQ(copy_gen.gene(toroute,i-from_startpos+to_startpos),        gen.gene(route,i));
          for (int i=from_startpos+to_natseq_size  ; i<gen.routeLength(route);  i++) ASSERT_EQ(copy_gen.gene(fromroute,i-to_natseq_size+from_natseq_size),gen.gene(route,i));
          
        }
        else if (route == toroute) {

          for (int i=0;                            i<to_startpos;             i++) ASSERT_EQ(copy_gen.gene(toroute,i),                                gen.gene(route,i));
          for (int i=to_startpos;                  i<from_natseq_size;        i++) ASSERT_EQ(copy_gen.gene(fromroute,i-to_startpos+from_startpos),    gen.gene(route,i));
          for (int i=to_startpos+from_natseq_size; i<gen.routeLength(route);  i++) ASSERT_EQ(copy_gen.gene(toroute,i-from_natseq_size+to_natseq_size),gen.gene(route,i));
          
        }
        else {
          for (int i=0; i<gen.routeLength(route); i++) {
            ASSERT_EQ(copy_gen.gene(route,i) , gen.gene(route,i));
          }
        }
      }

      // Finally the natural sequences swap are checked
      ASSERT_EQ(orig_from_natseq,toroute_natseqs[to_natseq]);
      ASSERT_EQ(orig_to_natseq,fromroute_natseqs[from_natseq]);
    }



};

TEST_F(SwapNatSeqsTest,swapSameSizeFirstNatSeqsOfFirstTwoRoutes) {
  int num_routes   = 3;

  int a_gen_values[] = {1,-1,2,-2,3,-3,4,-4,5,-5}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );
  DARPGenome* gen = createGen(num_routes,route,gen_values);

  int fromroute    = 0;
  int toroute      = 1;

  vector< vector< list<int>* > > natseqs; 
  for (int i=0; i<num_routes; i++) {
    int cut_pos[] = {gen->routeLength(i)/3, gen->routeLength(i)*2/3}; 
    vector<int> cut_positions(cut_pos, cut_pos+sizeof(cut_pos)/sizeof(int) );

    natseqs.push_back( createNatSeqs(*gen,i,cut_positions) );
  }

  int from_natseq = 0;
  int to_natseq   = 1;

  checkSwap( *gen, natseqs, fromroute, from_natseq, toroute, to_natseq ); 

  for (int i=0; i<num_routes; i++) {
    vector< list<int>* >& route_natseqs = natseqs[i];
    for (int j=0; j<route_natseqs.size(); j++) delete route_natseqs[j];
  }
}

TEST_F(SwapNatSeqsTest,swapDifferentSizesFirstTwoRoutes) {
  int num_routes   = 3;

  int a_gen_values[] = {1,-1,2,-2,3,-3,4,-4,5,-5}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );
  DARPGenome* gen = createGen(num_routes,route,gen_values);

  int fromroute    = 0;
  int toroute      = 1;

  vector< vector< list<int>* > > natseqs; 
  int first_cut_pos[] = {gen->routeLength(0)/3, gen->routeLength(0)*2/3}; 
  vector<int> first_cut_positions(first_cut_pos, first_cut_pos+sizeof(first_cut_pos)/sizeof(int) );

  int second_cut_pos[] = {gen->routeLength(1)-2}; 
  vector<int> second_cut_positions(second_cut_pos, second_cut_pos+sizeof(second_cut_pos)/sizeof(int) );

  int third_cut_pos[] = {gen->routeLength(2)/2}; 
  vector<int> third_cut_positions(third_cut_pos, third_cut_pos+sizeof(third_cut_pos)/sizeof(int) );
  
  natseqs.push_back( createNatSeqs(*gen,0,first_cut_positions) );
  natseqs.push_back( createNatSeqs(*gen,1,second_cut_positions) );
  natseqs.push_back( createNatSeqs(*gen,2,third_cut_positions) );

  int from_natseq = 0;
  int to_natseq   = 0;

  checkSwap( *gen, natseqs, fromroute, from_natseq, toroute, to_natseq ); 

  for (int i=0; i<num_routes; i++) {
    vector< list<int>* >& route_natseqs = natseqs[i];
    for (int j=0; j<route_natseqs.size(); j++) delete route_natseqs[j];
  }
}

TEST_F(SwapNatSeqsTest,swapLastNatSeqOfSize1FirstTwoRoutes) {
  int num_routes   = 3;

  int a_gen_values[] = {1,-1,2,-2,3,-3,4,-4,5,-5}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );
  DARPGenome* gen = createGen(num_routes,route,gen_values);

  int fromroute    = 0;
  int toroute      = 1;

  vector< vector< list<int>* > > natseqs; 
  int first_cut_pos[] = {gen->routeLength(0)/3, gen->routeLength(0)*2/3}; 
  vector<int> first_cut_positions(first_cut_pos, first_cut_pos+sizeof(first_cut_pos)/sizeof(int) );

  int second_cut_pos[] = {gen->routeLength(1)-2}; 
  vector<int> second_cut_positions(second_cut_pos, second_cut_pos+sizeof(second_cut_pos)/sizeof(int) );

  int third_cut_pos[] = {gen->routeLength(2)/2}; 
  vector<int> third_cut_positions(third_cut_pos, third_cut_pos+sizeof(third_cut_pos)/sizeof(int) );

  natseqs.push_back( createNatSeqs(*gen,0,first_cut_positions) );
  natseqs.push_back( createNatSeqs(*gen,1,second_cut_positions) );
  natseqs.push_back( createNatSeqs(*gen,2,third_cut_positions) );

  int from_natseq = 0;
  int to_natseq   = 1;

  checkSwap( *gen, natseqs, fromroute, from_natseq, toroute, to_natseq ); 

  for (int i=0; i<num_routes; i++) {
    vector< list<int>* >& route_natseqs = natseqs[i];
    for (int j=0; j<route_natseqs.size(); j++) delete route_natseqs[j];
  }
}

TEST_F(SwapNatSeqsTest,swapMiddleFirstNatSeqToLastLastNatSeq) {
  int num_routes   = 3;

  int a_gen_values[] = {1,-1,2,-2,3,-3,4,-4,5,-5}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );
  DARPGenome* gen = createGen(num_routes,route,gen_values);

  int fromroute    = 0;
  int toroute      = 2;

  vector< vector< list<int>* > > natseqs; 
  int first_cut_pos[] = {gen->routeLength(0)/3, gen->routeLength(0)*2/3}; 
  vector<int> first_cut_positions(first_cut_pos, first_cut_pos+sizeof(first_cut_pos)/sizeof(int) );

  int second_cut_pos[] = {gen->routeLength(1)-2}; 
  vector<int> second_cut_positions(second_cut_pos, second_cut_pos+sizeof(second_cut_pos)/sizeof(int) );

  int third_cut_pos[] = {3, 6, 9}; 
  vector<int> third_cut_positions(third_cut_pos, third_cut_pos+sizeof(third_cut_pos)/sizeof(int) );

  natseqs.push_back( createNatSeqs(*gen,0,first_cut_positions) );
  natseqs.push_back( createNatSeqs(*gen,1,second_cut_positions) );
  natseqs.push_back( createNatSeqs(*gen,2,third_cut_positions) );

  int from_natseq = 1;
  int to_natseq   = 2;

  checkSwap( *gen, natseqs, fromroute, from_natseq, toroute, to_natseq ); 

  for (int i=0; i<num_routes; i++) {
    vector< list<int>* >& route_natseqs = natseqs[i];
    for (int j=0; j<route_natseqs.size(); j++) delete route_natseqs[j];
  }
}

class isSwapFeasibleTest : public CheckNatSeqsTest {};

TEST_F(isSwapFeasibleTest,canNatSeqReplaceAnotherThreeNatSeqsPerRoute) {
  int num_routes   = 3;

  int a_gen_values[] = {1,-1,2,-2,3,-3,4,-4,5,-5}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );
  DARPGenome* gen = createGen(num_routes,route,gen_values);

  vector< vector< list<int>* > > natseqs; 
  for (int i=0; i<num_routes; i++) {

    int cut_pos[] = {gen->routeLength(i)/3, gen->routeLength(i)*2/3}; 
    vector<int> cut_positions(cut_pos, cut_pos+sizeof(cut_pos)/sizeof(int) );

    natseqs.push_back( createNatSeqs(*gen,i,cut_positions) );
  }

  //cout << *gen << endl;
  //printNatSeqs(natseqs);
  //printNatSeqsTimes(natseqs);

  int fromroute   = 0;
  int toroute     = 0;
  int from_natseq = 0;
  int to_natseq   = 0;

  ASSERT_EQ( true,  CheckAllNaturalSeqsCombsNeighborhood::canNatSeqReplaceNatSeq(natseqs, fromroute, from_natseq, toroute, to_natseq) );

  fromroute   = 0;
  toroute     = 2;
  from_natseq = 0;
  to_natseq   = 2;

  ASSERT_EQ( false,  CheckAllNaturalSeqsCombsNeighborhood::canNatSeqReplaceNatSeq(natseqs, fromroute, from_natseq, toroute, to_natseq) );

  fromroute   = 0;
  toroute     = 2;
  from_natseq = 0;
  to_natseq   = 1;

  ASSERT_EQ( true,  CheckAllNaturalSeqsCombsNeighborhood::canNatSeqReplaceNatSeq(natseqs, fromroute, from_natseq, toroute, to_natseq) );


  for (int i=0; i<num_routes; i++) {
    vector< list<int>* >& route_natseqs = natseqs[i];
    for (int j=0; j<route_natseqs.size(); j++) delete route_natseqs[j];
  }
}

TEST_F(isSwapFeasibleTest,canNatSeqReplaceAnotherFiveNatSeqsPerRoute) {
  int num_routes   = 3;

  int a_gen_values[] = {1,-1,2,-2,3,-3,4,-4,5,-5}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );
  DARPGenome* gen = createGen(num_routes,route,gen_values);

  vector< vector< list<int>* > > natseqs; 
  for (int i=0; i<num_routes; i++) {

    int cut_pos[] = {gen->routeLength(i)/5, gen->routeLength(i)*2/5,gen->routeLength(i)*3/5,gen->routeLength(i)*4/5 }; 
    vector<int> cut_positions(cut_pos, cut_pos+sizeof(cut_pos)/sizeof(int) );

    natseqs.push_back( createNatSeqs(*gen,i,cut_positions) );
  }

  //cout << *gen << endl;
  //printNatSeqs(natseqs);
  //printNatSeqsTimes(natseqs);

  int fromroute   = 0;
  int toroute     = 0;
  int from_natseq = 0;
  int to_natseq   = 0;

  ASSERT_EQ( true,  CheckAllNaturalSeqsCombsNeighborhood::canNatSeqReplaceNatSeq(natseqs, fromroute, from_natseq, toroute, to_natseq) );

  fromroute   = 0;
  toroute     = 2;
  from_natseq = 0;
  to_natseq   = 2;

  ASSERT_EQ( true,  CheckAllNaturalSeqsCombsNeighborhood::canNatSeqReplaceNatSeq(natseqs, fromroute, from_natseq, toroute, to_natseq) );

  fromroute   = 0;
  toroute     = 2;
  from_natseq = 0;
  to_natseq   = 1;

  ASSERT_EQ( true,  CheckAllNaturalSeqsCombsNeighborhood::canNatSeqReplaceNatSeq(natseqs, fromroute, from_natseq, toroute, to_natseq) );

  fromroute   = 0;
  toroute     = 1;
  from_natseq = 0;
  to_natseq   = 3;

  ASSERT_EQ( false,  CheckAllNaturalSeqsCombsNeighborhood::canNatSeqReplaceNatSeq(natseqs, fromroute, from_natseq, toroute, to_natseq) );

  fromroute   = 0;
  toroute     = 1;
  from_natseq = 0;
  to_natseq   = 4;

  ASSERT_EQ( false,  CheckAllNaturalSeqsCombsNeighborhood::canNatSeqReplaceNatSeq(natseqs, fromroute, from_natseq, toroute, to_natseq) );

  fromroute   = 1;
  toroute     = 0;
  from_natseq = 1;
  to_natseq   = 4;

  ASSERT_EQ( false,  CheckAllNaturalSeqsCombsNeighborhood::canNatSeqReplaceNatSeq(natseqs, fromroute, from_natseq, toroute, to_natseq) );

  fromroute   = 1;
  toroute     = 0;
  from_natseq = 1;
  to_natseq   = 2;

  ASSERT_EQ( true,  CheckAllNaturalSeqsCombsNeighborhood::canNatSeqReplaceNatSeq(natseqs, fromroute, from_natseq, toroute, to_natseq) );

  fromroute   = 1;
  toroute     = 0;
  from_natseq = 1;
  to_natseq   = 3;

  ASSERT_EQ( true,  CheckAllNaturalSeqsCombsNeighborhood::canNatSeqReplaceNatSeq(natseqs, fromroute, from_natseq, toroute, to_natseq) );

  for (int i=0; i<num_routes; i++) {
    vector< list<int>* >& route_natseqs = natseqs[i];
    for (int j=0; j<route_natseqs.size(); j++) delete route_natseqs[j];
  }
}

TEST_F(isSwapFeasibleTest,isSwapFeasibleFiveNatSeqsPerRoute) {
  int num_routes   = 3;

  int a_gen_values[] = {1,-1,2,-2,3,-3,4,-4,5,-5}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );
  DARPGenome* gen = createGen(num_routes,route,gen_values);

  vector< vector< list<int>* > > natseqs; 
  for (int i=0; i<num_routes; i++) {

    int cut_pos[] = {gen->routeLength(i)/5, gen->routeLength(i)*2/5,gen->routeLength(i)*3/5,gen->routeLength(i)*4/5 }; 
    vector<int> cut_positions(cut_pos, cut_pos+sizeof(cut_pos)/sizeof(int) );

    natseqs.push_back( createNatSeqs(*gen,i,cut_positions) );
  }

  //cout << *gen << endl;
  //printNatSeqs(natseqs);
  //printNatSeqsTimes(natseqs);

  int fromroute   = 0;
  int toroute     = 0;
  int from_natseq = 0;
  int to_natseq   = 0;

  ASSERT_EQ( true,  CheckAllNaturalSeqsCombsNeighborhood::isSwapFeasible(natseqs, fromroute, from_natseq, toroute, to_natseq) );

  fromroute   = 0;
  toroute     = 2;
  from_natseq = 0;
  to_natseq   = 2;

  ASSERT_EQ( true,  CheckAllNaturalSeqsCombsNeighborhood::isSwapFeasible(natseqs, fromroute, from_natseq, toroute, to_natseq) );

  fromroute   = 0;
  toroute     = 2;
  from_natseq = 0;
  to_natseq   = 1;

  ASSERT_EQ( true,  CheckAllNaturalSeqsCombsNeighborhood::isSwapFeasible(natseqs, fromroute, from_natseq, toroute, to_natseq) );

  fromroute   = 0;
  toroute     = 1;
  from_natseq = 0;
  to_natseq   = 3;

  ASSERT_EQ( false,  CheckAllNaturalSeqsCombsNeighborhood::isSwapFeasible(natseqs, fromroute, from_natseq, toroute, to_natseq) );

  fromroute   = 0;
  toroute     = 1;
  from_natseq = 0;
  to_natseq   = 4;

  ASSERT_EQ( false,  CheckAllNaturalSeqsCombsNeighborhood::isSwapFeasible(natseqs, fromroute, from_natseq, toroute, to_natseq) );

  fromroute   = 1;
  toroute     = 0;
  from_natseq = 1;
  to_natseq   = 4;

  ASSERT_EQ( false,  CheckAllNaturalSeqsCombsNeighborhood::isSwapFeasible(natseqs, fromroute, from_natseq, toroute, to_natseq) );

  fromroute   = 1;
  toroute     = 0;
  from_natseq = 1;
  to_natseq   = 2;

  ASSERT_EQ( true,  CheckAllNaturalSeqsCombsNeighborhood::isSwapFeasible(natseqs, fromroute, from_natseq, toroute, to_natseq) );

  fromroute   = 1;
  toroute     = 0;
  from_natseq = 1;
  to_natseq   = 3;

  ASSERT_EQ( false,  CheckAllNaturalSeqsCombsNeighborhood::isSwapFeasible(natseqs, fromroute, from_natseq, toroute, to_natseq) );

  for (int i=0; i<num_routes; i++) {
    vector< list<int>* >& route_natseqs = natseqs[i];
    for (int j=0; j<route_natseqs.size(); j++) delete route_natseqs[j];
  }
}

TEST_F(isSwapFeasibleTest,checkingScoreOfSwaps) {
  int num_routes   = 3;

  int a_gen_values[] = {1,-1,2,-2,3,-3,4,-4,5,-5}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );
  DARPGenome* gen = createGen(num_routes,route,gen_values);

  vector< vector< list<int>* > > natseqs; 
  for (int i=0; i<num_routes; i++) {

    int cut_pos[] = {gen->routeLength(i)/5, gen->routeLength(i)*2/5,gen->routeLength(i)*3/5,gen->routeLength(i)*4/5 }; 
    vector<int> cut_positions(cut_pos, cut_pos+sizeof(cut_pos)/sizeof(int) );

    natseqs.push_back( createNatSeqs(*gen,i,cut_positions) );
  }

  //cout << *gen << endl;
  //printNatSeqs(natseqs);
  //printNatSeqsTimes(natseqs);

  int fromroute   = 0;
  int toroute     = 1;
  int from_natseq = 0;
  int to_natseq   = 2;

  ASSERT_NE( GAGenome::worstPossibleScore(), CheckAllNaturalSeqsCombsNeighborhood::scoreFromSwapNaturalSequences(*gen,natseqs, fromroute, from_natseq, toroute, to_natseq) );

  fromroute   = 0;
  toroute     = 1;
  from_natseq = 0;
  to_natseq   = 3;

  ASSERT_EQ( GAGenome::worstPossibleScore(), CheckAllNaturalSeqsCombsNeighborhood::scoreFromSwapNaturalSequences(*gen,natseqs, fromroute, from_natseq, toroute, to_natseq) );

  for (int i=0; i<num_routes; i++) {
    vector< list<int>* >& route_natseqs = natseqs[i];
    for (int j=0; j<route_natseqs.size(); j++) delete route_natseqs[j];
  }
}

class BestSwapTest : public CheckNatSeqsTest {

    virtual void SetUp() {
      num_routes  = 3;
      route       = 0;

      requestslist = new VerticesList();

      DummyDistMatrix* distmatrix = new DummyDistMatrix(MAXPOINTS);

      // The costs are set so that the optimum sequence is 1 -1 2 -2 3 -3
      // The remaining combinations are set to a higher cost
      
      int allpos[] = {1,10,2,20,3,30,4,40,5,50,6,60};
      vector<int> allposv (allpos, allpos + sizeof(allpos)/ sizeof(int) );
      for (vector<int>::iterator it=allposv.begin(); it!=allposv.end(); it++) {
        for (vector<int>::iterator itj=it+1; itj!=allposv.end(); itj++) {
          distmatrix->setCostOneSide( *it,*itj,3);
          distmatrix->setCostOneSide( *itj,*it,3);
        }
      }

      distmatrix->setCostOneSide( 1,10,1);
      distmatrix->setCostOneSide(10, 2,1);
      distmatrix->setCostOneSide( 2,20,1);
      distmatrix->setCostOneSide(20, 3,1);
      distmatrix->setCostOneSide( 3,30,1);
      distmatrix->setCostOneSide( 30,4,1);
      distmatrix->setCostOneSide( 4,40,1);
      distmatrix->setCostOneSide( 40,5,1);
      distmatrix->setCostOneSide( 5,50,1);

      // wend = wbegin + MAXUSERWAITTIME (20)
      // delivery vertex wbegin = wbegin + cost (origin dest), wend = wend + MAXUSERRIDETIME (60)

      addRequest(*requestslist, 1, 1, 10,  5,  *distmatrix);
      addRequest(*requestslist, 2, 2, 20, 15, *distmatrix);
      addRequest(*requestslist, 3, 3, 30, 30, *distmatrix);
      addRequest(*requestslist, 4, 4, 40, 45, *distmatrix);
      addRequest(*requestslist, 5, 5, 50, 60, *distmatrix);

      addRequest(*requestslist, 10, 1, 10,  5, *distmatrix);
      addRequest(*requestslist, 20, 2, 20, 15, *distmatrix);
      addRequest(*requestslist, 30, 3, 30, 30, *distmatrix);
      addRequest(*requestslist, 40, 4, 40, 45, *distmatrix);
      addRequest(*requestslist, 50, 5, 50, 60, *distmatrix);

      addRequest(*requestslist, 100, 1, 10,  5, *distmatrix);
      addRequest(*requestslist, 200, 2, 20, 15, *distmatrix);
      addRequest(*requestslist, 300, 3, 30, 30, *distmatrix);
      addRequest(*requestslist, 400, 4, 40, 45, *distmatrix);
      addRequest(*requestslist, 500, 5, 50, 60, *distmatrix);

      dummymatrix = new DummyCostMatrix(distmatrix,*requestslist);
      DARPTestsVars::initData(requestslist,dummymatrix);
    }

};

TEST_F(BestSwapTest,NoSwapsNeeded) {
  int num_routes   = 3;

  int a_first_values[]  = {  1,  -1,  2,  -2,  3,  -3,  4,  -4,  5,  -5}; 
  int a_second_values[] = { 10, -10, 20, -20, 30, -30, 40, -40, 50, -50}; 
  int a_third_values[]  = {100,-100,200,-200,300,-300,400,-400,500,-500}; 

  vector<int> first_values  (a_first_values,  a_first_values  + sizeof(a_first_values)/  sizeof(int) );
  vector<int> second_values (a_second_values, a_second_values + sizeof(a_second_values)/ sizeof(int) );
  vector<int> third_values  (a_third_values,  a_third_values  + sizeof(a_third_values)/  sizeof(int) );

  DARPGenome* gen = createGen(num_routes);
  setRoute(*gen,0,first_values);
  setRoute(*gen,1,second_values);
  setRoute(*gen,2,third_values);

  vector< vector< list<int>* > > natseqs; 
  for (int i=0; i<num_routes; i++) {

    int cut_pos[] = {gen->routeLength(i)/5, gen->routeLength(i)*2/5,gen->routeLength(i)*3/5,gen->routeLength(i)*4/5 }; 
    vector<int> cut_positions(cut_pos, cut_pos+sizeof(cut_pos)/sizeof(int) );

    natseqs.push_back( createNatSeqs(*gen,i,cut_positions) );
  }

  int fromroute   = 0;
  int toroute     = 1;

  int    bestfrom_pos = -1;
  int    bestto_pos   = -1;
  double swap_score   = -1;

  CheckAllNaturalSeqsCombsNeighborhood::bestSwapOfNatSeqs(*gen, natseqs, fromroute, toroute, bestfrom_pos, bestto_pos, swap_score);

  ASSERT_EQ(gen->score(),swap_score);
  ASSERT_EQ(-1,  bestfrom_pos);
  ASSERT_EQ(-1,  bestto_pos);
  ASSERT_NE(-1,  swap_score);

  for (int i=0; i<num_routes; i++) {
    vector< list<int>* >& route_natseqs = natseqs[i];
    for (int j=0; j<route_natseqs.size(); j++) delete route_natseqs[j];
  }
}

TEST_F(BestSwapTest,SwapBetweenFirstTwoRoutesFistAndSecondNatSeqs) {
  int num_routes   = 3;

  int a_first_values[]  = {  2,  -2,  2,  -2,  3,  -3,  4,  -4,  5,  -5}; 
  int a_second_values[] = { 10, -10, 10, -10, 30, -30, 40, -40, 50, -50}; 
  int a_third_values[]  = {100,-100,200,-200,300,-300,400,-400,500,-500}; 

  vector<int> first_values  (a_first_values,  a_first_values  + sizeof(a_first_values)/  sizeof(int) );
  vector<int> second_values (a_second_values, a_second_values + sizeof(a_second_values)/ sizeof(int) );
  vector<int> third_values  (a_third_values,  a_third_values  + sizeof(a_third_values)/  sizeof(int) );

  DARPGenome* gen = createGen(num_routes);
  setRoute(*gen,0,first_values);
  setRoute(*gen,1,second_values);
  setRoute(*gen,2,third_values);

  vector< vector< list<int>* > > natseqs; 
  for (int i=0; i<num_routes; i++) {

    int cut_pos[] = {gen->routeLength(i)/5, gen->routeLength(i)*2/5,gen->routeLength(i)*3/5,gen->routeLength(i)*4/5 }; 
    vector<int> cut_positions(cut_pos, cut_pos+sizeof(cut_pos)/sizeof(int) );

    natseqs.push_back( createNatSeqs(*gen,i,cut_positions) );
  }

  //cout << *gen << endl;
  //printNatSeqs(natseqs);

  int fromroute   = 0;
  int toroute     = 1;

  int    bestfrom_pos = -1;
  int    bestto_pos   = -1;
  double swap_score   = -1;

  CheckAllNaturalSeqsCombsNeighborhood::bestSwapOfNatSeqs(*gen, natseqs, fromroute, toroute, bestfrom_pos, bestto_pos, swap_score);

  ASSERT_EQ(true, (GAGenome::compareScores(swap_score, gen->score()) == GAGenome::BETTER ));
  ASSERT_NE(-1,  bestfrom_pos);
  ASSERT_NE(-1,  bestto_pos);
  ASSERT_NE(-1,  swap_score);

  //cout << "best from route " << fromroute << " pos " << bestfrom_pos << "toroute: " << toroute << " pos: " << bestto_pos << " swap score " << swap_score << endl;
  //printRoute(*gen,fromroute,DARPGenome::getVerticesList(),DARPGenome::getCostMatrix());
  //printRoute(*gen,toroute,DARPGenome::getVerticesList(),DARPGenome::getCostMatrix());

  //if (bestfrom_pos >= 0)  {
  //  CheckAllNaturalSeqsCombsNeighborhood::swapNaturalSequences(*gen,natseqs,fromroute,bestfrom_pos,toroute,bestto_pos);

  //  printRoute(*gen,fromroute,DARPGenome::getVerticesList(),DARPGenome::getCostMatrix());
  //  printRoute(*gen,toroute,DARPGenome::getVerticesList(),DARPGenome::getCostMatrix());
  //  printRoute(*gen,2,DARPGenome::getVerticesList(),DARPGenome::getCostMatrix());
  //}


  for (int i=0; i<num_routes; i++) {
    vector< list<int>* >& route_natseqs = natseqs[i];
    for (int j=0; j<route_natseqs.size(); j++) delete route_natseqs[j];
  }
}

TEST_F(BestSwapTest,SwapBetweenFirstTwoRoutesFistAndFirstNatSeqsDiffSize) {
  int num_routes   = 3;

  int a_first_values[]  = {  2,  -2,  2,  -2,  3,  -3,  4,  -4,  5,  -5}; 
  int a_second_values[] = { 10,-10,40, -40, 50, -50}; 
  int a_third_values[]  = {300,-300,400,-400}; 

  vector<int> first_values  (a_first_values,  a_first_values  + sizeof(a_first_values)/  sizeof(int) );
  vector<int> second_values (a_second_values, a_second_values + sizeof(a_second_values)/ sizeof(int) );
  vector<int> third_values  (a_third_values,  a_third_values  + sizeof(a_third_values)/  sizeof(int) );

  DARPGenome* gen = createGen(num_routes);
  setRoute(*gen,0,first_values);
  setRoute(*gen,1,second_values);
  setRoute(*gen,2,third_values);

  vector< vector< list<int>* > > natseqs; 

  int first_cut_pos[] = {gen->routeLength(0)/5, gen->routeLength(0)*2/5,gen->routeLength(0)*3/5,gen->routeLength(0)*4/5 }; 
  vector<int> first_cut_positions(first_cut_pos, first_cut_pos+sizeof(first_cut_pos)/sizeof(int) );

  int second_cut_pos[] = {gen->routeLength(1)/3, gen->routeLength(1)*2/3 }; 
  vector<int> second_cut_positions(second_cut_pos, second_cut_pos+sizeof(second_cut_pos)/sizeof(int) );

  int third_cut_pos[] = {gen->routeLength(2)/2 }; 
  vector<int> third_cut_positions(third_cut_pos, third_cut_pos+sizeof(third_cut_pos)/sizeof(int) );

  natseqs.push_back( createNatSeqs(*gen,0,first_cut_positions) );
  natseqs.push_back( createNatSeqs(*gen,1,second_cut_positions) );
  natseqs.push_back( createNatSeqs(*gen,2,third_cut_positions) );

//  cout << *gen << endl;
//  printNatSeqs(natseqs);

  int fromroute   = 0;
  int toroute     = 1;

  int    bestfrom_pos = -1;
  int    bestto_pos   = -1;
  double swap_score   = -1;

  CheckAllNaturalSeqsCombsNeighborhood::bestSwapOfNatSeqs(*gen, natseqs, fromroute, toroute, bestfrom_pos, bestto_pos, swap_score);

  ASSERT_EQ(true, (GAGenome::compareScores(swap_score, gen->score()) == GAGenome::BETTER ));
  ASSERT_NE(-1,  swap_score);
  ASSERT_EQ(0, bestfrom_pos);
  ASSERT_EQ(0, bestto_pos);

  for (int i=0; i<num_routes; i++) {
    vector< list<int>* >& route_natseqs = natseqs[i];
    for (int j=0; j<route_natseqs.size(); j++) delete route_natseqs[j];
  }
}

TEST_F(BestSwapTest,SwapBetweenFirstTwoRoutesLastAndMiddleNatSeqsDiffSize) {
  int num_routes   = 3;

  int a_first_values[]  = {  1,  -1,  2,  -2,  3,  -3,  4,  -4,  30,  -30};
  int a_second_values[] = { 10,-10,50, -50, 40, -40}; 
  int a_third_values[]  = {300,-300,400,-400}; 

  vector<int> first_values  (a_first_values,  a_first_values  + sizeof(a_first_values)/  sizeof(int) );
  vector<int> second_values (a_second_values, a_second_values + sizeof(a_second_values)/ sizeof(int) );
  vector<int> third_values  (a_third_values,  a_third_values  + sizeof(a_third_values)/  sizeof(int) );

  DARPGenome* gen = createGen(num_routes);
  setRoute(*gen,0,first_values);
  setRoute(*gen,1,second_values);
  setRoute(*gen,2,third_values);

  vector< vector< list<int>* > > natseqs; 

  int first_cut_pos[] = {gen->routeLength(0)/5, gen->routeLength(0)*2/5,gen->routeLength(0)*3/5,gen->routeLength(0)*4/5 }; 
  vector<int> first_cut_positions(first_cut_pos, first_cut_pos+sizeof(first_cut_pos)/sizeof(int) );

  int second_cut_pos[] = {gen->routeLength(1)/3, gen->routeLength(1)*2/3 }; 
  vector<int> second_cut_positions(second_cut_pos, second_cut_pos+sizeof(second_cut_pos)/sizeof(int) );

  int third_cut_pos[] = {gen->routeLength(2)/2 }; 
  vector<int> third_cut_positions(third_cut_pos, third_cut_pos+sizeof(third_cut_pos)/sizeof(int) );

  natseqs.push_back( createNatSeqs(*gen,0,first_cut_positions) );
  natseqs.push_back( createNatSeqs(*gen,1,second_cut_positions) );
  natseqs.push_back( createNatSeqs(*gen,2,third_cut_positions) );

//  cout << *gen << endl;
//  printNatSeqs(natseqs);


  int fromroute   = 0;
  int toroute     = 1;

  int    bestfrom_pos = -1;
  int    bestto_pos   = -1;
  double swap_score   = -1;

  CheckAllNaturalSeqsCombsNeighborhood::bestSwapOfNatSeqs(*gen, natseqs, fromroute, toroute, bestfrom_pos, bestto_pos, swap_score);

  ASSERT_EQ(true, (GAGenome::compareScores(swap_score, gen->score()) == GAGenome::BETTER ));
  ASSERT_NE(-1,  swap_score);
  ASSERT_EQ(4, bestfrom_pos);
  ASSERT_EQ(1, bestto_pos);

  //cout << "best from route " << fromroute << " pos " << bestfrom_pos << "toroute: " << toroute << " pos: " << bestto_pos << " swap score " << swap_score << endl;
  //printRoute(*gen,fromroute,DARPGenome::getVerticesList(),DARPGenome::getCostMatrix());
  //printRoute(*gen,toroute,DARPGenome::getVerticesList(),DARPGenome::getCostMatrix());

//  if (bestfrom_pos >= 0)  {
//    printNatSeqsTimes(natseqs);
//
//    CheckAllNaturalSeqsCombsNeighborhood::swapNaturalSequences(*gen,natseqs,fromroute,bestfrom_pos,toroute,bestto_pos);
//
//    cout << "score of swap is " << gen->score() << endl;
//
//    printRoute(*gen,fromroute,DARPGenome::getVerticesList(),DARPGenome::getCostMatrix());
//    printRoute(*gen,toroute,DARPGenome::getVerticesList(),DARPGenome::getCostMatrix());
//    printRoute(*gen,2,DARPGenome::getVerticesList(),DARPGenome::getCostMatrix());
//
  //}


  for (int i=0; i<num_routes; i++) {
    vector< list<int>* >& route_natseqs = natseqs[i];
    for (int j=0; j<route_natseqs.size(); j++) delete route_natseqs[j];
  }
}

TEST_F(BestSwapTest,SwapBetweenFirstTwoRoutesSeveralAlternativesDiffSize) {
  int num_routes   = 3;

  int a_first_values[]  = {  2,  -2,  2,  -2,  2,  -2,  4,  -4,  5,  -5}; 
  int a_second_values[] = { 10,-10,50, -50, 40, -40}; 
  int a_third_values[]  = {100,-100,400,-400}; 

  vector<int> first_values  (a_first_values,  a_first_values  + sizeof(a_first_values)/  sizeof(int) );
  vector<int> second_values (a_second_values, a_second_values + sizeof(a_second_values)/ sizeof(int) );
  vector<int> third_values  (a_third_values,  a_third_values  + sizeof(a_third_values)/  sizeof(int) );

  DARPGenome* gen = createGen(num_routes);
  setRoute(*gen,0,first_values);
  setRoute(*gen,1,second_values);
  setRoute(*gen,2,third_values);

  vector< vector< list<int>* > > natseqs; 

  int first_cut_pos[] = {gen->routeLength(0)/5, gen->routeLength(0)*2/5,gen->routeLength(0)*3/5,gen->routeLength(0)*4/5 }; 
  vector<int> first_cut_positions(first_cut_pos, first_cut_pos+sizeof(first_cut_pos)/sizeof(int) );

  int second_cut_pos[] = {gen->routeLength(1)/3, gen->routeLength(1)*2/3 }; 
  vector<int> second_cut_positions(second_cut_pos, second_cut_pos+sizeof(second_cut_pos)/sizeof(int) );

  int third_cut_pos[] = {gen->routeLength(2)/2 }; 
  vector<int> third_cut_positions(third_cut_pos, third_cut_pos+sizeof(third_cut_pos)/sizeof(int) );

  natseqs.push_back( createNatSeqs(*gen,0,first_cut_positions) );
  natseqs.push_back( createNatSeqs(*gen,1,second_cut_positions) );
  natseqs.push_back( createNatSeqs(*gen,2,third_cut_positions) );

//  cout << *gen << endl;
//  printNatSeqs(natseqs);

  int fromroute   = 0;
  int toroute     = 2;

  int    bestfrom_pos = -1;
  int    bestto_pos   = -1;
  double swap_score   = -1;

  CheckAllNaturalSeqsCombsNeighborhood::bestSwapOfNatSeqs(*gen, natseqs, fromroute, toroute, bestfrom_pos, bestto_pos, swap_score);

  ASSERT_EQ(true, (GAGenome::compareScores(swap_score, gen->score()) == GAGenome::BETTER ));
  ASSERT_NE(-1,  swap_score);
  ASSERT_EQ(0, bestfrom_pos);
  ASSERT_EQ(0, bestto_pos);

  //cout << "best from route " << fromroute << " pos " << bestfrom_pos << "toroute: " << toroute << " pos: " << bestto_pos << " swap score " << swap_score << endl;
  //printRoute(*gen,fromroute,DARPGenome::getVerticesList(),DARPGenome::getCostMatrix());
  //printRoute(*gen,toroute,DARPGenome::getVerticesList(),DARPGenome::getCostMatrix());

//  if (bestfrom_pos >= 0)  {
//    printNatSeqsTimes(natseqs);
//
//    DARPGenome copy_gen(*gen);
//    cout << "orig gen" << endl << endl;
//    CheckAllNaturalSeqsCombsNeighborhood::swapNaturalSequences(*gen,natseqs,fromroute,2,toroute,1);
//    printRoute(*gen,fromroute,DARPGenome::getVerticesList(),DARPGenome::getCostMatrix());
//    printRoute(*gen,toroute,  DARPGenome::getVerticesList(),DARPGenome::getCostMatrix());
//    cout << "copy gen" << endl << endl;
//    CheckAllNaturalSeqsCombsNeighborhood::swapNaturalSequences(copy_gen,natseqs,fromroute,0,toroute,0);
//    printRoute(copy_gen,fromroute,DARPGenome::getVerticesList(),DARPGenome::getCostMatrix());
//    printRoute(copy_gen,toroute,  DARPGenome::getVerticesList(),DARPGenome::getCostMatrix());
//
//    cout << "score of swap is " << gen->score() << endl;
//
//    printRoute(*gen,fromroute,DARPGenome::getVerticesList(),DARPGenome::getCostMatrix());
//    printRoute(*gen,toroute,DARPGenome::getVerticesList(),DARPGenome::getCostMatrix());
//    printRoute(*gen,2,DARPGenome::getVerticesList(),DARPGenome::getCostMatrix());
//
  //}


  for (int i=0; i<num_routes; i++) {
    vector< list<int>* >& route_natseqs = natseqs[i];
    for (int j=0; j<route_natseqs.size(); j++) delete route_natseqs[j];
  }
}

TEST_F(BestSwapTest,CheckCompleteShakerOneAlternative) {
  int num_routes   = 3;

  int a_first_values[]  = {  4,  -4,  2,  -2,  3,  -3}; 
  int a_second_values[] = { 10,-10,20, -20, 40, -40}; 
  int a_third_values[]  = {300,-300,100,-100}; 

  vector<int> first_values  (a_first_values,  a_first_values  + sizeof(a_first_values)/  sizeof(int) );
  vector<int> second_values (a_second_values, a_second_values + sizeof(a_second_values)/ sizeof(int) );
  vector<int> third_values  (a_third_values,  a_third_values  + sizeof(a_third_values)/  sizeof(int) );

  DARPGenome* gen = createGen(num_routes);
  setRoute(*gen,0,first_values);
  setRoute(*gen,1,second_values);
  setRoute(*gen,2,third_values);

  vector< vector< list<int>* > > allnatseqs(gen->size());
  for (int i=0; i<gen->length(); i++)  CheckAllNaturalSeqsCombsNeighborhood::getAllNaturalSeqs(*gen,i, allnatseqs[i]);
  //printNatSeqs(allnatseqs);

  CheckAllNaturalSeqsCombsNeighborhood shaker(0);

  double old_score = gen->score();

  shaker(*gen);

  ASSERT_EQ(true, (GAGenome::compareScores(gen->score(), old_score ) == GAGenome::BETTER ));


  for (int i=0; i<num_routes; i++) {
    vector< list<int>* >& route_natseqs = allnatseqs[i];
    for (int j=0; j<route_natseqs.size(); j++) delete route_natseqs[j];
  }
}

TEST_F(BestSwapTest,CheckCompleteShakerSeveralAlternative) {
  int num_routes   = 3;

  int a_first_values[]  = {  5,  -5,  2,  -2, 4,  -4}; 
  int a_second_values[] = { 10,-10,50, -50, 40, -40}; 
  int a_third_values[]  = {300,-300,400,-400}; 

  vector<int> first_values  (a_first_values,  a_first_values  + sizeof(a_first_values)/  sizeof(int) );
  vector<int> second_values (a_second_values, a_second_values + sizeof(a_second_values)/ sizeof(int) );
  vector<int> third_values  (a_third_values,  a_third_values  + sizeof(a_third_values)/  sizeof(int) );

  DARPGenome* gen = createGen(num_routes);
  setRoute(*gen,0,first_values);
  setRoute(*gen,1,second_values);
  setRoute(*gen,2,third_values);

  double old_score = gen->score();

  CheckAllNaturalSeqsCombsNeighborhood shaker(0);
  shaker(*gen);

  ASSERT_EQ(true, (GAGenome::compareScores(gen->score(), old_score ) == GAGenome::BETTER ));

  //for (int i=0; i<num_routes; i++) {
  //  vector< list<int>* >& route_natseqs = allnatseqs[i];
  //  for (int j=0; j<route_natseqs.size(); j++) delete route_natseqs[j];
  //}
}

// checkear secuencias naturales con vertices con cargas distintas

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  GAGenome::optCriterion(GAGenome::MINIMIZATION);

  setUpTestingEnvironment();

  return RUN_ALL_TESTS();
}


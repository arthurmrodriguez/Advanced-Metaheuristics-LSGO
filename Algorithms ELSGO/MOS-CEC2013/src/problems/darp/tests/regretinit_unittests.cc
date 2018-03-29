#include "daropstests_aux.h"

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

class AuxFunctionInitTest : public ::testing::Test {
  static const int MAXPOINTS = 300; 

public:
  DummyDistMatrix* distmatrix;

  Request *req1, *req2, *req3, *req4;
  Vertex *req1pickup_vert,   *req2pickup_vert,   *req3pickup_vert,   *req4pickup_vert;
  Vertex *req1delivery_vert, *req2delivery_vert, *req3delivery_vert, *req4delivery_vert;


    virtual void SetUp() {
      distmatrix = new DummyDistMatrix(MAXPOINTS);

      int id, position, load, fb, fe;
      bool critical;

      id = 1; critical = true; position = 1; load = 2; fb = 5; fe = 10;
      req1pickup_vert = new Vertex(id++,position++,critical, load,fb,fe,Vertex::PICKUP);
      req2pickup_vert = new Vertex(id++,position++,critical, load,fb,fe,Vertex::PICKUP);
      req3pickup_vert = new Vertex(id++,position++,critical, load,fb,fe,Vertex::PICKUP);
      req4pickup_vert = new Vertex(id++,position++,critical, load,fb,fe,Vertex::PICKUP);

      id = -1; critical = false; position = 10; fb += 10; fe += 10;
      req1delivery_vert = new Vertex(id--,position++, critical, load,fb,fe,Vertex::DELIVERY);
      req2delivery_vert = new Vertex(id--,position++, critical, load,fb,fe,Vertex::DELIVERY);
      req3delivery_vert = new Vertex(id--,position++, critical, load,fb,fe,Vertex::DELIVERY);
      req4delivery_vert = new Vertex(id--,position++, critical, load,fb,fe,Vertex::DELIVERY);

      req1 = new Request( req1pickup_vert, req1delivery_vert);
      req2 = new Request( req2pickup_vert, req2delivery_vert);
      req3 = new Request( req3pickup_vert, req3delivery_vert);
      req4 = new Request( req4pickup_vert, req4delivery_vert);
    }

    int setCostsReq1Req2() {
      distmatrix->setCostOneSide( req1->pickup_vert->pos_,   req2->pickup_vert->pos_,   2.0);
      distmatrix->setCostOneSide( req1->pickup_vert->pos_,   req2->delivery_vert->pos_, 3.0);
      distmatrix->setCostOneSide( req1->delivery_vert->pos_, req2->pickup_vert->pos_,   4.0);
      distmatrix->setCostOneSide( req1->delivery_vert->pos_, req2->delivery_vert->pos_, 5.0);

      return 2+3+4+5;
    }

    int setCostsReq2Req1() {
      distmatrix->setCostOneSide( req2->pickup_vert->pos_,   req1->pickup_vert->pos_,   6.0);
      distmatrix->setCostOneSide( req2->pickup_vert->pos_,   req1->delivery_vert->pos_, 7.0);
      distmatrix->setCostOneSide( req2->delivery_vert->pos_, req1->pickup_vert->pos_,   8.0);
      distmatrix->setCostOneSide( req2->delivery_vert->pos_, req1->delivery_vert->pos_, 9.0);
      return 6+7+8+9;
    }

    int setCostsReq1Req3() {
      distmatrix->setCostOneSide( req1->pickup_vert->pos_,   req3->pickup_vert->pos_,   10.0);
      distmatrix->setCostOneSide( req1->pickup_vert->pos_,   req3->delivery_vert->pos_, 11.0);
      distmatrix->setCostOneSide( req1->delivery_vert->pos_, req3->pickup_vert->pos_,   12.0);
      distmatrix->setCostOneSide( req1->delivery_vert->pos_, req3->delivery_vert->pos_, 13.0);
      return 10+11+12+13;
    }

    int setCostsReq2Req3() {
      distmatrix->setCostOneSide( req2->pickup_vert->pos_,   req3->pickup_vert->pos_,   16.0);
      distmatrix->setCostOneSide( req2->pickup_vert->pos_,   req3->delivery_vert->pos_, 17.0);
      distmatrix->setCostOneSide( req2->delivery_vert->pos_, req3->pickup_vert->pos_,   18.0);
      distmatrix->setCostOneSide( req2->delivery_vert->pos_, req3->delivery_vert->pos_, 19.0);
      return 16+17+18+19;
    }

    int setCostsReq3Req1() {
      distmatrix->setCostOneSide( req3->pickup_vert->pos_,   req1->pickup_vert->pos_,   20.0);
      distmatrix->setCostOneSide( req3->pickup_vert->pos_,   req1->delivery_vert->pos_, 21.0);
      distmatrix->setCostOneSide( req3->delivery_vert->pos_, req1->pickup_vert->pos_,   22.0);
      distmatrix->setCostOneSide( req3->delivery_vert->pos_, req1->delivery_vert->pos_, 23.0);
      return 20+21+22+23;
    }

    int setCostsReq3Req2() {
      distmatrix->setCostOneSide( req3->pickup_vert->pos_,   req2->pickup_vert->pos_,   24.0);
      distmatrix->setCostOneSide( req3->pickup_vert->pos_,   req2->delivery_vert->pos_, 25.0);
      distmatrix->setCostOneSide( req3->delivery_vert->pos_, req2->pickup_vert->pos_,   26.0);
      distmatrix->setCostOneSide( req3->delivery_vert->pos_, req2->delivery_vert->pos_, 27.0);
      return 24+25+26+27;
    }


    virtual void TearDown() {
      delete distmatrix;
      delete req1; delete req2; delete req3;
      delete req1pickup_vert;   delete req2pickup_vert;   delete req3pickup_vert;
      delete req1delivery_vert; delete req2delivery_vert; delete req3delivery_vert;
    }
};

TEST_F(AuxFunctionInitTest,sumPickAndDelDistances) {

  int sum1_2 = setCostsReq1Req2();

  ASSERT_EQ( sum1_2, DARPRegretInsertionInitC::sumPickAndDelDistances(*distmatrix,*req1,*req2) );

  int sum2_1 = setCostsReq2Req1();

  ASSERT_EQ( sum2_1, DARPRegretInsertionInitC::sumPickAndDelDistances(*distmatrix,*req2,*req1) );


  distmatrix->setCostOneSide( req1->pickup_vert->pos_,   req1->pickup_vert->pos_,   -1.0);
  distmatrix->setCostOneSide( req1->pickup_vert->pos_,   req1->delivery_vert->pos_, -2.0);
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req1->pickup_vert->pos_,   -3.0);
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req1->delivery_vert->pos_, -4.0);

  ASSERT_EQ( -1 - 2 - 3 - 4, DARPRegretInsertionInitC::sumPickAndDelDistances(*distmatrix,*req1,*req1) );


  distmatrix->setCostOneSide( req2->pickup_vert->pos_,   req2->pickup_vert->pos_,   -5.0);
  distmatrix->setCostOneSide( req2->pickup_vert->pos_,   req2->delivery_vert->pos_, -6.0);
  distmatrix->setCostOneSide( req2->delivery_vert->pos_, req2->pickup_vert->pos_,   -7.0);
  distmatrix->setCostOneSide( req2->delivery_vert->pos_, req2->delivery_vert->pos_, -8.0);

  ASSERT_EQ( -5 - 6 - 7 - 8, DARPRegretInsertionInitC::sumPickAndDelDistances(*distmatrix,*req2,*req2) );
}


TEST_F(AuxFunctionInitTest,sumOfDistancesToOtherRequests) {

  vector<Request> requests;

#ifdef DEBUG
  // With no available positions it should failed to execute
  ASSERT_DEATH(sumOfDistancesToOtherRequests(*distmatrix, requests, 0),"Assertion failed");
#endif


  requests.push_back(*req1);

  // With just one request the result should be zero
  ASSERT_EQ(0, DARPRegretInsertionInitC::sumOfDistancesToOtherRequests(*distmatrix, requests, 0));

#ifdef DEBUG
  // Trying to execute with an invalid position should fail
  ASSERT_DEATH(sumOfDistancesToOtherRequests(*distmatrix, requests, 1),"Assertion failed");
#endif

  requests.push_back(*req2);

  // With two requests, depending on the pos passed, it should return one or the other sum of distances

  // If we do not set the distances it should return 0
  ASSERT_EQ(0, DARPRegretInsertionInitC::sumOfDistancesToOtherRequests(*distmatrix, requests, 0));

  int sum1_2 = setCostsReq1Req2();

  ASSERT_EQ(sum1_2, DARPRegretInsertionInitC::sumOfDistancesToOtherRequests(*distmatrix, requests, 0));


  int sum2_1 = setCostsReq2Req1();

  ASSERT_EQ(sum2_1, DARPRegretInsertionInitC::sumOfDistancesToOtherRequests(*distmatrix, requests, 1));


  // Finally, we test to push a third request, the value should be the sum of the other two
  requests.push_back(*req3);

  int sum1_3 = setCostsReq1Req3();

  ASSERT_EQ(sum1_2+sum1_3, DARPRegretInsertionInitC::sumOfDistancesToOtherRequests(*distmatrix, requests, 0));


  int sum2_3 = setCostsReq2Req3();

  ASSERT_EQ(sum2_1+sum2_3, DARPRegretInsertionInitC::sumOfDistancesToOtherRequests(*distmatrix, requests, 1));
}

TEST_F(AuxFunctionInitTest,decentralizationIndex) {
  vector<Request> requests;
  requests.push_back(*req1);
  requests.push_back(*req2);
  requests.push_back(*req3);

  double sum1_2 = setCostsReq1Req2();
  double sum1_3 = setCostsReq1Req3();
  double sum2_1 = setCostsReq2Req1();
  double sum2_3 = setCostsReq2Req3();
  double sum3_1 = setCostsReq3Req1();
  double sum3_2 = setCostsReq3Req2();

  ASSERT_EQ( (sum1_2+sum1_3) / (sum1_2+sum1_3+sum2_1+sum2_3+sum3_1+sum3_2) , DARPRegretInsertionInitC::decentralizationIndex(*distmatrix,requests,0) );
  ASSERT_EQ( (sum2_1+sum2_3) / (sum1_2+sum1_3+sum2_1+sum2_3+sum3_1+sum3_2) , DARPRegretInsertionInitC::decentralizationIndex(*distmatrix,requests,1) );
  ASSERT_EQ( (sum3_1+sum3_2) / (sum1_2+sum1_3+sum2_1+sum2_3+sum3_1+sum3_2) , DARPRegretInsertionInitC::decentralizationIndex(*distmatrix,requests,2) );
}

TEST_F(AuxFunctionInitTest,sortRequestsAccordingToDecentrIndex) {
  double sum1_2 = setCostsReq1Req2();
  double sum1_3 = setCostsReq1Req3();
  double sum2_1 = setCostsReq2Req1();
  double sum2_3 = setCostsReq2Req3();
  double sum3_1 = setCostsReq3Req1();
  double sum3_2 = setCostsReq3Req2();

  vector<Request> requests;

  // decentralized index 3 > 2 > 1

  // First test: no change should happened
  requests.push_back(*req3);
  requests.push_back(*req2);
  requests.push_back(*req1);

  DARPRegretInsertionInitC::processRequestsAccordingToDecentrIndex(*distmatrix,requests,1);

  ASSERT_EQ( req3->pickup_vert->id_, requests[0].pickup_vert->id_);
  ASSERT_EQ( req2->pickup_vert->id_, requests[1].pickup_vert->id_);
  ASSERT_EQ( req1->pickup_vert->id_, requests[2].pickup_vert->id_);


  // Second test, 3 should be swapped with 2, and then 3 with 1
  requests.clear();

  requests.push_back(*req1);
  requests.push_back(*req2);
  requests.push_back(*req3);

  DARPRegretInsertionInitC::processRequestsAccordingToDecentrIndex(*distmatrix,requests,1);

  ASSERT_EQ( req3->pickup_vert->id_, requests[0].pickup_vert->id_);
  ASSERT_EQ( req1->pickup_vert->id_, requests[1].pickup_vert->id_);
  ASSERT_EQ( req2->pickup_vert->id_, requests[2].pickup_vert->id_);


  // Third test: Same test as before but with alpha = 0, no change should happen
  requests.clear();

  requests.push_back(*req1);
  requests.push_back(*req2);
  requests.push_back(*req3);

  DARPRegretInsertionInitC::processRequestsAccordingToDecentrIndex(*distmatrix,requests,0);

  ASSERT_EQ( req1->pickup_vert->id_, requests[0].pickup_vert->id_);
  ASSERT_EQ( req2->pickup_vert->id_, requests[1].pickup_vert->id_);
  ASSERT_EQ( req3->pickup_vert->id_, requests[2].pickup_vert->id_);

  // Fourth test, 2 should be swapped with 1 but not with 3
  requests.clear();

  requests.push_back(*req3);
  requests.push_back(*req1);
  requests.push_back(*req2);

  DARPRegretInsertionInitC::processRequestsAccordingToDecentrIndex(*distmatrix,requests,1);

  ASSERT_EQ( req3->pickup_vert->id_, requests[0].pickup_vert->id_);
  ASSERT_EQ( req2->pickup_vert->id_, requests[1].pickup_vert->id_);
  ASSERT_EQ( req1->pickup_vert->id_, requests[2].pickup_vert->id_);
}


TEST_F(AuxFunctionInitTest,isTravelingFromLateDeliveryGoingToArriveBeforeEarlyPickup) {

  // Simple test where there is plenty of time to arrive before the early pickup
  req1->delivery_vert->fend_ = 5;
  req2->pickup_vert->fbegin_ = 100;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req2->pickup_vert->pos_,1);

  ASSERT_TRUE(DARPRegretInsertionInitC::isTravelingFromLateDeliveryGoingToArriveBeforeEarlyPickup(*req1,*req2,*distmatrix));

  // We repeat the previous test but changing the requests
  req2->delivery_vert->fend_ = 5;
  req1->pickup_vert->fbegin_ = 10;
  distmatrix->setCostOneSide( req2->delivery_vert->pos_, req1->pickup_vert->pos_,1);

  ASSERT_TRUE(DARPRegretInsertionInitC::isTravelingFromLateDeliveryGoingToArriveBeforeEarlyPickup(*req2,*req1,*distmatrix));

  // Simple test, here there is not enough time to arrive before the pickup (just one unit)
  req1->delivery_vert->fend_ = 5;
  req2->pickup_vert->fbegin_ = 10;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req2->pickup_vert->pos_,6);

  ASSERT_FALSE(DARPRegretInsertionInitC::isTravelingFromLateDeliveryGoingToArriveBeforeEarlyPickup(*req1,*req2,*distmatrix));

  // Similar to before but here there is exactly the time needed to arrive at the beginning of the pickup
  req1->delivery_vert->fend_ = 5;
  req2->pickup_vert->fbegin_ = 10;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req2->pickup_vert->pos_,5);

  ASSERT_TRUE(DARPRegretInsertionInitC::isTravelingFromLateDeliveryGoingToArriveBeforeEarlyPickup(*req1,*req2,*distmatrix));
}

TEST_F(AuxFunctionInitTest,findFirstFromPosThatDoesNotSatisfyTimeRequirementsWithReqPos) {
  std::vector<Request> requests;
#ifdef DEBUG
  // With an invalid position it should fail to execute
  ASSERT_DEATH(findFirstFromPosThatDoesNotSatisfyTimeRequirementsWithReqPos(-1,*distmatrix,requests);
  ASSERT_DEATH(findFirstFromPosThatDoesNotSatisfyTimeRequirementsWithReqPos(0,*distmatrix,requests);
  ASSERT_DEATH(findFirstFromPosThatDoesNotSatisfyTimeRequirementsWithReqPos(1,*distmatrix,requests);
#endif

  // With only one request, it should return -1
  requests.push_back(*req1);
  ASSERT_EQ(-1, DARPRegretInsertionInitC::findFirstFromPosThatDoesNotSatisfyTimeRequirementsWithReqPos(0,*distmatrix,requests) );

  // If no request is found -1 should be returned
  requests.clear();
  req1->delivery_vert->fend_ = 5;
  req2->pickup_vert->fbegin_ = 100;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req2->pickup_vert->pos_,1);
  requests.push_back(*req1); requests.push_back(*req2);

  ASSERT_EQ(-1, DARPRegretInsertionInitC::findFirstFromPosThatDoesNotSatisfyTimeRequirementsWithReqPos(0,*distmatrix,requests) );

  // Test with two requests. Request 2 should be returned
  requests.clear();
  req1->delivery_vert->fend_ = 5;
  req2->pickup_vert->fbegin_ = 6;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req2->pickup_vert->pos_,2);
  requests.push_back(*req1); requests.push_back(*req2);

  ASSERT_EQ(1, DARPRegretInsertionInitC::findFirstFromPosThatDoesNotSatisfyTimeRequirementsWithReqPos(0,*distmatrix,requests) );

  // Test with three requests. Request 2 should be returned although 3 does not satisfy the conditions
  requests.clear();
  req1->delivery_vert->fend_ = 5;
  req2->pickup_vert->fbegin_ = 6;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req2->pickup_vert->pos_,2);
  req3->pickup_vert->fbegin_ = 1;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req3->pickup_vert->pos_,1);
  requests.push_back(*req1); requests.push_back(*req2); requests.push_back(*req3);

  ASSERT_EQ(1, DARPRegretInsertionInitC::findFirstFromPosThatDoesNotSatisfyTimeRequirementsWithReqPos(0,*distmatrix,requests) );

  // Test with three requests. Request 3 should be returned since 2 satisfies the conditions
  requests.clear();
  req1->delivery_vert->fend_ = 5;
  req2->pickup_vert->fbegin_ = 7;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req2->pickup_vert->pos_,2);
  req3->pickup_vert->fbegin_ = 1;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req3->pickup_vert->pos_,1);
  requests.push_back(*req1); requests.push_back(*req2); requests.push_back(*req3);

  ASSERT_EQ(2, DARPRegretInsertionInitC::findFirstFromPosThatDoesNotSatisfyTimeRequirementsWithReqPos(0,*distmatrix,requests) );

  // Similar to before but here we start from request 2 and 3 does not satisfy the conditions from request 2
  requests.clear();
  req1->delivery_vert->fend_ = 5;
  req2->pickup_vert->fbegin_ = 7;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req2->pickup_vert->pos_,2);
  req2->delivery_vert->fend_ = 9;
  req3->pickup_vert->fbegin_ = 9;
  distmatrix->setCostOneSide( req2->delivery_vert->pos_, req3->pickup_vert->pos_,1);
  requests.push_back(*req1); requests.push_back(*req2); requests.push_back(*req3);

  ASSERT_EQ(2, DARPRegretInsertionInitC::findFirstFromPosThatDoesNotSatisfyTimeRequirementsWithReqPos(1,*distmatrix,requests) );

  // Similar to the previous test but with in this case request 3 satisfies the constraint, -1 should be returned
  requests.clear();
  req1->delivery_vert->fend_ = 5;
  req2->pickup_vert->fbegin_ = 7;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req2->pickup_vert->pos_,2);
  req2->delivery_vert->fend_ = 9;
  req3->pickup_vert->fbegin_ = 10;
  distmatrix->setCostOneSide( req2->delivery_vert->pos_, req3->pickup_vert->pos_,1);
  requests.push_back(*req1); requests.push_back(*req2); requests.push_back(*req3);

  ASSERT_EQ(-1, DARPRegretInsertionInitC::findFirstFromPosThatDoesNotSatisfyTimeRequirementsWithReqPos(1,*distmatrix,requests) );
}

TEST_F(AuxFunctionInitTest,processRequestsAccordingToBetterMargin) {
  std::vector<Request> requests;

#ifdef DEBUG
  // With an invalid nroutes an assertion should failed
  ASSERT_DEATH(processRequestsAccordingToBetterMargin(-1,*distmatrix,requests);
  ASSERT_DEATH(processRequestsAccordingToBetterMargin(1,*distmatrix,requests);
#endif


  // Nothing should happen with an empty requests vector
  DARPRegretInsertionInitC::processRequestsAccordingToBetterMargin(0,*distmatrix,requests);
  ASSERT_EQ(0,requests.size());

  // Similarly with only one request
  requests.push_back(*req1);
  DARPRegretInsertionInitC::processRequestsAccordingToBetterMargin(0,*distmatrix,requests);
  ASSERT_EQ(1,requests.size());
  DARPRegretInsertionInitC::processRequestsAccordingToBetterMargin(1,*distmatrix,requests);
  ASSERT_EQ(1,requests.size());

  // Testing with two requests
  // First, no change should happen
  requests.clear();
  req1->delivery_vert->fend_ = 5;
  req2->pickup_vert->fbegin_ = 7;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req2->pickup_vert->pos_,2);
  requests.push_back(*req1); requests.push_back(*req2);
  DARPRegretInsertionInitC::processRequestsAccordingToBetterMargin(2,*distmatrix,requests);
  ASSERT_EQ(req1->pickup_vert->id_,requests[0].pickup_vert->id_);
  ASSERT_EQ(req2->pickup_vert->id_,requests[1].pickup_vert->id_);

  // The constraints are not satisfied but since there are only two requests, no change should happen
  requests.clear();
  req1->delivery_vert->fend_ = 5;
  req2->pickup_vert->fbegin_ = 7;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req2->pickup_vert->pos_,3);
  requests.push_back(*req1); requests.push_back(*req2);
  DARPRegretInsertionInitC::processRequestsAccordingToBetterMargin(2,*distmatrix,requests);
  ASSERT_EQ(req1->pickup_vert->id_,requests[0].pickup_vert->id_);
  ASSERT_EQ(req2->pickup_vert->id_,requests[1].pickup_vert->id_);

  // Testing with three requests
  // The constraints are being satisfied, no change should happen
  requests.clear();
  req1->delivery_vert->fend_ = 5;
  req2->pickup_vert->fbegin_ = 7;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req2->pickup_vert->pos_,2);
  req2->delivery_vert->fend_ = 9;
  req3->pickup_vert->fbegin_ = 10;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req3->pickup_vert->pos_,1);
  distmatrix->setCostOneSide( req2->delivery_vert->pos_, req3->pickup_vert->pos_,1);
  requests.push_back(*req1); requests.push_back(*req2); requests.push_back(*req3);

  DARPRegretInsertionInitC::processRequestsAccordingToBetterMargin(3,*distmatrix,requests);
  ASSERT_EQ(req1->pickup_vert->id_,requests[0].pickup_vert->id_);
  ASSERT_EQ(req2->pickup_vert->id_,requests[1].pickup_vert->id_);
  ASSERT_EQ(req3->pickup_vert->id_,requests[2].pickup_vert->id_);

  // 2 and 3 should be swapped
  requests.clear();
  req1->delivery_vert->fend_ = 5;
  req2->pickup_vert->fbegin_ = 7;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req2->pickup_vert->pos_,2);
  req2->delivery_vert->fend_ = 9;
  req3->pickup_vert->fbegin_ = 5;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req3->pickup_vert->pos_,1);
  distmatrix->setCostOneSide( req2->delivery_vert->pos_, req3->pickup_vert->pos_,1);
  requests.push_back(*req1); requests.push_back(*req2); requests.push_back(*req3);

  DARPRegretInsertionInitC::processRequestsAccordingToBetterMargin(3,*distmatrix,requests);
  ASSERT_EQ(req1->pickup_vert->id_,requests[0].pickup_vert->id_);
  ASSERT_EQ(req3->pickup_vert->id_,requests[1].pickup_vert->id_);
  ASSERT_EQ(req2->pickup_vert->id_,requests[2].pickup_vert->id_);

  // Similar, but specifying a smaller number of routes, they should be swapped similarly to before
  requests.clear();
  requests.push_back(*req1); requests.push_back(*req2); requests.push_back(*req3);
  DARPRegretInsertionInitC::processRequestsAccordingToBetterMargin(2,*distmatrix,requests);
  ASSERT_EQ(req1->pickup_vert->id_,requests[0].pickup_vert->id_);
  ASSERT_EQ(req3->pickup_vert->id_,requests[1].pickup_vert->id_);
  ASSERT_EQ(req2->pickup_vert->id_,requests[2].pickup_vert->id_);

  // Four insertions
  // The constraints are being satisfied, no change should happen
  requests.clear();
  req1->delivery_vert->fend_ = 5;
  req2->pickup_vert->fbegin_ = 7;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req2->pickup_vert->pos_,2);
  req2->delivery_vert->fend_ = 9;
  req3->pickup_vert->fbegin_ = 10;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req3->pickup_vert->pos_,1);
  distmatrix->setCostOneSide( req2->delivery_vert->pos_, req3->pickup_vert->pos_,1);
  req3->delivery_vert->fend_ = 11;
  req4->pickup_vert->fbegin_ = 12;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req4->pickup_vert->pos_,1);
  distmatrix->setCostOneSide( req2->delivery_vert->pos_, req4->pickup_vert->pos_,1);
  distmatrix->setCostOneSide( req3->delivery_vert->pos_, req4->pickup_vert->pos_,1);

  requests.push_back(*req1); requests.push_back(*req2); requests.push_back(*req3); requests.push_back(*req4);

  DARPRegretInsertionInitC::processRequestsAccordingToBetterMargin(4,*distmatrix,requests);
  ASSERT_EQ(req1->pickup_vert->id_,requests[0].pickup_vert->id_);
  ASSERT_EQ(req2->pickup_vert->id_,requests[1].pickup_vert->id_);
  ASSERT_EQ(req3->pickup_vert->id_,requests[2].pickup_vert->id_);
  ASSERT_EQ(req4->pickup_vert->id_,requests[3].pickup_vert->id_);

  // Request two and three are incompatible but since they are contiguous, nothing should happened
  requests.clear();
  req1->delivery_vert->fend_ = 5;
  req2->pickup_vert->fbegin_ = 7;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req2->pickup_vert->pos_,2);
  req2->delivery_vert->fend_ = 9;
  req3->pickup_vert->fbegin_ = 10;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req3->pickup_vert->pos_,1);
  distmatrix->setCostOneSide( req2->delivery_vert->pos_, req3->pickup_vert->pos_,10);
  req3->delivery_vert->fend_ = 11;
  req4->pickup_vert->fbegin_ = 12;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req4->pickup_vert->pos_,1);
  distmatrix->setCostOneSide( req2->delivery_vert->pos_, req4->pickup_vert->pos_,1);
  distmatrix->setCostOneSide( req3->delivery_vert->pos_, req4->pickup_vert->pos_,1);

  requests.push_back(*req1); requests.push_back(*req2); requests.push_back(*req3); requests.push_back(*req4);

  DARPRegretInsertionInitC::processRequestsAccordingToBetterMargin(4,*distmatrix,requests);
  ASSERT_EQ(req1->pickup_vert->id_,requests[0].pickup_vert->id_);
  ASSERT_EQ(req2->pickup_vert->id_,requests[1].pickup_vert->id_);
  ASSERT_EQ(req3->pickup_vert->id_,requests[2].pickup_vert->id_);
  ASSERT_EQ(req4->pickup_vert->id_,requests[3].pickup_vert->id_);

  // Four insertions
  // The constraints are being satisfied, no change should happen
  requests.clear();
  req1->delivery_vert->fend_ = 5;
  req2->pickup_vert->fbegin_ = 7;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req2->pickup_vert->pos_,2);
  req2->delivery_vert->fend_  = 9;
  req3->pickup_vert->fbegin_  = 10;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req3->pickup_vert->pos_,1);
  distmatrix->setCostOneSide( req2->delivery_vert->pos_, req3->pickup_vert->pos_,1);
  req3->delivery_vert->fend_  = 11;
  req4->pickup_vert->fbegin_  = 12;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req4->pickup_vert->pos_,1);
  distmatrix->setCostOneSide( req2->delivery_vert->pos_, req4->pickup_vert->pos_,1);
  distmatrix->setCostOneSide( req3->delivery_vert->pos_, req4->pickup_vert->pos_,1);

  requests.push_back(*req1); requests.push_back(*req2); requests.push_back(*req3); requests.push_back(*req4);

  DARPRegretInsertionInitC::processRequestsAccordingToBetterMargin(4,*distmatrix,requests);
  ASSERT_EQ(req1->pickup_vert->id_,requests[0].pickup_vert->id_);
  ASSERT_EQ(req2->pickup_vert->id_,requests[1].pickup_vert->id_);
  ASSERT_EQ(req3->pickup_vert->id_,requests[2].pickup_vert->id_);
  ASSERT_EQ(req4->pickup_vert->id_,requests[3].pickup_vert->id_);


  // Request two and four are incompatible they should be swapped. Four and three are incompatible but since they are
  // contiguous they should not be swapped
  requests.clear();
  req1->delivery_vert->fend_   = 5;
  req2->pickup_vert->fbegin_   = 7;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req2->pickup_vert->pos_,2);
  req2->delivery_vert->fend_   = 9;
  req3->pickup_vert->fbegin_   = 10;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req3->pickup_vert->pos_,1);
  distmatrix->setCostOneSide( req2->delivery_vert->pos_, req3->pickup_vert->pos_,1);
  req3->delivery_vert->fend_   = 11;
  req4->pickup_vert->fbegin_   = 12;
  req4->delivery_vert->fend_   = 13;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req4->pickup_vert->pos_,1);
  distmatrix->setCostOneSide( req2->delivery_vert->pos_, req4->pickup_vert->pos_,10);
  distmatrix->setCostOneSide( req3->delivery_vert->pos_, req4->pickup_vert->pos_,1);
  distmatrix->setCostOneSide( req4->delivery_vert->pos_, req3->pickup_vert->pos_,100);

  requests.push_back(*req1); requests.push_back(*req2); requests.push_back(*req3); requests.push_back(*req4);

  DARPRegretInsertionInitC::processRequestsAccordingToBetterMargin(4,*distmatrix,requests);
  ASSERT_EQ(req1->pickup_vert->id_,requests[0].pickup_vert->id_);
  ASSERT_EQ(req2->pickup_vert->id_,requests[1].pickup_vert->id_);
  ASSERT_EQ(req4->pickup_vert->id_,requests[2].pickup_vert->id_);
  ASSERT_EQ(req3->pickup_vert->id_,requests[3].pickup_vert->id_);

  // Only using three routes. Request one and four are incompatible, they should be swapped. Besides, four and
  // two are incompatible so at the end we should have 1, 4, 2, 3
  requests.clear();
  req1->delivery_vert->fend_   = 5;
  req2->pickup_vert->fbegin_   = 7;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req2->pickup_vert->pos_,2);
  req2->delivery_vert->fend_   = 9;
  req3->pickup_vert->fbegin_   = 14; // compatible with 4
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req3->pickup_vert->pos_,1);
  distmatrix->setCostOneSide( req2->delivery_vert->pos_, req3->pickup_vert->pos_,1);
  req3->delivery_vert->fend_   = 11;
  req4->pickup_vert->fbegin_   = 12;
  req4->delivery_vert->fend_   = 13;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req4->pickup_vert->pos_,100);
  distmatrix->setCostOneSide( req2->delivery_vert->pos_, req4->pickup_vert->pos_,1);
  distmatrix->setCostOneSide( req3->delivery_vert->pos_, req4->pickup_vert->pos_,1);
  distmatrix->setCostOneSide( req4->delivery_vert->pos_, req2->pickup_vert->pos_,100);
  distmatrix->setCostOneSide( req4->delivery_vert->pos_, req3->pickup_vert->pos_,1);

  requests.push_back(*req1); requests.push_back(*req2); requests.push_back(*req3); requests.push_back(*req4);

  DARPRegretInsertionInitC::processRequestsAccordingToBetterMargin(4,*distmatrix,requests);
  ASSERT_EQ(req1->pickup_vert->id_,requests[0].pickup_vert->id_);
  ASSERT_EQ(req4->pickup_vert->id_,requests[1].pickup_vert->id_);
  ASSERT_EQ(req2->pickup_vert->id_,requests[2].pickup_vert->id_);
  ASSERT_EQ(req3->pickup_vert->id_,requests[3].pickup_vert->id_);


  // Testing with two routes from the total of four requests. All the routes are incompatible between them
  // except for the first one. Therefore, no changes should happen
  requests.clear();
  req1->delivery_vert->fend_   = 5;
  req2->pickup_vert->fbegin_   = 7;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req2->pickup_vert->pos_,2);
  req2->delivery_vert->fend_   = 9;
  req3->pickup_vert->fbegin_   = 10;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req3->pickup_vert->pos_,1);
  distmatrix->setCostOneSide( req2->delivery_vert->pos_, req3->pickup_vert->pos_,100);
  req3->delivery_vert->fend_   = 11;
  req4->pickup_vert->fbegin_   = 12;
  req4->delivery_vert->fend_   = 13;
  distmatrix->setCostOneSide( req1->delivery_vert->pos_, req4->pickup_vert->pos_,1);
  distmatrix->setCostOneSide( req2->delivery_vert->pos_, req4->pickup_vert->pos_,100);
  distmatrix->setCostOneSide( req3->delivery_vert->pos_, req4->pickup_vert->pos_,100);
  distmatrix->setCostOneSide( req4->delivery_vert->pos_, req2->pickup_vert->pos_,100);
  distmatrix->setCostOneSide( req4->delivery_vert->pos_, req3->pickup_vert->pos_,100);

  requests.push_back(*req1); requests.push_back(*req2); requests.push_back(*req3); requests.push_back(*req4);

  DARPRegretInsertionInitC::processRequestsAccordingToBetterMargin(2,*distmatrix,requests);
  ASSERT_EQ(req1->pickup_vert->id_,requests[0].pickup_vert->id_);
  ASSERT_EQ(req2->pickup_vert->id_,requests[1].pickup_vert->id_);
  ASSERT_EQ(req3->pickup_vert->id_,requests[2].pickup_vert->id_);
  ASSERT_EQ(req4->pickup_vert->id_,requests[3].pickup_vert->id_);
}




class RegretInsertionTesting : public AuxFunctionInitTest {
protected:
  static const int MAXPOINTS = 300;

public:
  VerticesList*    requestslist;
  DummyCostMatrix* dummymatrix;
  int              num_routes;
  int              route     ;

    virtual void SetUp() {
      GAGenome::optCriterion (GAGenome::MINIMIZATION);

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

    Request createRequestToGoAfterAndBeforeVertices(int origvertexid, int destvertexid, DARPGenome& gen) {
      Vertex& vert = requestslist->getVertex(origvertexid);

      // We use the origin vertex to adjust the window values of the new request
      int id = 50; int pos = 50; int load = 1;
      int fb, fe;
      fb = vert.fbegin_ + 1;
      fe = fb+1;

      Vertex* reqtmppickup_vert    = new Vertex(id,   pos,  true, load,fb,  fe,  Vertex::PICKUP);
      Vertex* reqtmpdelivery_vert  = new Vertex(-1*id,pos+1,false,load,fe+1,fe+2,Vertex::DELIVERY);
      Request reqtmp (reqtmppickup_vert,reqtmpdelivery_vert);

      requestslist->addVertex(reqtmppickup_vert);
      requestslist->addVertex(reqtmpdelivery_vert);

      // We set the values of the matrix so that the cost for all the combinations is fixed at a value
      for (int route=0; route<gen.size(); route++) {
        for (int pos=0; pos<gen.routeLength(route); pos++) {
          dummymatrix->setCost(gen.gene(route,pos),id,10);
          dummymatrix->setCost(gen.gene(route,pos),-1*id,10);
        }
      }

      // Then we set that the cost after the -20 and before 30 is 0
      dummymatrix->setCostOneSide(origvertexid,id,0);
      dummymatrix->setCostOneSide(-1*id,destvertexid,0);
      dummymatrix->setCostOneSide(id,-1*id,0);

      return reqtmp;
    }

    void checkComputeBestAndAllInsertionScore(int origvertexid, int destvertexid, int route, DARPGenome& gen ) {
      Request reqtmp = createRequestToGoAfterAndBeforeVertices(origvertexid,destvertexid,gen);

      double         best_score;
      int            best_ins_route;
      vector<double> ins_scores;

      DARPRegretInsertionInitC::computeBestAndAllInsertionScores(gen,reqtmp,best_score,best_ins_route,ins_scores);

      ASSERT_EQ(num_routes,ins_scores.size());
      for (int i=0; i<ins_scores.size(); i++) ASSERT_TRUE( best_score <= ins_scores[i] );

      bool all_equal = true;
      for (int i=0; i<ins_scores.size(); i++) if (ins_scores[i] != best_score) all_equal = false;
      ASSERT_FALSE(all_equal);


      ASSERT_EQ(route,best_ins_route);

    }

};


TEST_F(RegretInsertionTesting,scoreForInsertingRequestInRoute) {
  int num_routes = 3;
  int main_route = 0;


  // We test that after calling a scoreForInsertingRequest the score has not changed although it has returned the
  // new score (in this case and the dummy obj function, worse than before)

  int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );

  DARPGenome* gen       = createGen(num_routes,main_route,gen_values);
  DARPGenome* other_gen = createGen(num_routes,main_route,gen_values);

  // We check that the assertion for the route parameter is working
#ifdef DEBUG
  // With an invalid nroutes an assertion should failed
  ASSERT_DEATH(scoreForInsertingRequestInRoute(*gen,*req4,-1));
  ASSERT_DEATH(scoreForInsertingRequestInRoute(*gen,*req4,4));
#endif

  double orig_score = gen->score();
  int pickins, delins; //not being used
  double score = scoreForInsertingRequestInRoute(*gen,*req4,0,DARPObjFuncBasedInitC::BEST,pickins,delins);

  ASSERT_EQ(orig_score,gen->score());
  ASSERT_NE(orig_score, score);
  ASSERT_TRUE(orig_score < score);

  // We check that the genome has not changed from the operation of inserting a request
  for (int route=0; route<gen->length(); route++) {
    for (int pos=0; pos<gen->routeLength(route); pos++) {
      ASSERT_EQ(other_gen->gene(route,pos), gen->gene(route,pos));
    }
  }

  delete gen;
  delete other_gen;
}

TEST_F(RegretInsertionTesting,computeBestAndAllInsertionScoresSimpleCase) {
  int num_routes = 3;
  int main_route = 0;

  // We first test a simple case where the best score is equals to all the insertion scores and compare the results
  int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );

  DARPGenome* gen = createGen(num_routes,main_route,gen_values);

  double         best_score;
  int            best_ins_route;
  vector<double> ins_scores;

  DARPRegretInsertionInitC::computeBestAndAllInsertionScores(*gen,*req4,best_score,best_ins_route,ins_scores);
  ASSERT_EQ(num_routes,ins_scores.size());
  ASSERT_EQ(best_score,ins_scores[0]);
  for (int i=1; i<ins_scores.size(); i++) ASSERT_EQ(ins_scores[0],ins_scores[i]);
  ASSERT_EQ(0,best_ins_route);

  delete gen;
}

TEST_F(RegretInsertionTesting,computeBestAndAllInsertionScoresComplexCase) {
  int num_routes = 3;
  int main_route = 0;

  // test to insert a request that should not increase the score
  int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );
  DARPGenome* gen = createGen(num_routes,main_route,gen_values);

  // We create the new request that we are going to try to insert in several routes
  checkComputeBestAndAllInsertionScore(-1,    2,0,*gen);
  checkComputeBestAndAllInsertionScore(20,  -30,1,*gen);
  checkComputeBestAndAllInsertionScore(300,-300,2,*gen);

  delete gen;
}

TEST_F(RegretInsertionTesting,computeReqRegretValue) {
  int num_routes = 3;
  int main_route = 0;

  // test to insert a request that should not increase the score
  int a_gen_values[] = {1,-1,2,-2,3,-3}; vector<int> gen_values (a_gen_values, a_gen_values + sizeof(a_gen_values)/ sizeof(int) );
  DARPGenome* gen = createGen(num_routes,main_route,gen_values);

  Request reqtmp1 = createRequestToGoAfterAndBeforeVertices(-1,2,*gen);
  Request reqtmp2 = createRequestToGoAfterAndBeforeVertices(20,-30,*gen);

  int ins_route1=-1, ins_route2=-2;

  double regret_value1 = DARPRegretInsertionInitC::computeReqRegretValue(*gen,reqtmp1,ins_route1);
  double regret_value2 = DARPRegretInsertionInitC::computeReqRegretValue(*gen,reqtmp2,ins_route2);

  ASSERT_NE(0,regret_value1);
  ASSERT_NE(0,regret_value2);
  ASSERT_EQ(regret_value1,regret_value2);
  ASSERT_TRUE(regret_value1 > 0);
  ASSERT_NE(-1,ins_route1);
  ASSERT_NE(-1,ins_route2);

  // We create a vertex that should not be place well anywhere, its regret value should be higher
  int id = 60; int pos = 60; int load = 1;
  int fb, fe;
  fb = 100;
  fe = 200;;

  Vertex* reqtmppickup_vert    = new Vertex(id,   pos,  true, load,fb,  fe,  Vertex::PICKUP);
  Vertex* reqtmpdelivery_vert  = new Vertex(-1*id,pos+1,false,load,fe+1,fe+2,Vertex::DELIVERY);
  Request reqtmp3 (reqtmppickup_vert,reqtmpdelivery_vert);

  requestslist->addVertex(reqtmppickup_vert);
  requestslist->addVertex(reqtmpdelivery_vert);

  for (int route=0; route<gen->size(); route++) {
    for (int pos=0; pos<gen->routeLength(route); pos++) {
      dummymatrix->setCost(gen->gene(route,pos),id,10);
      dummymatrix->setCost(gen->gene(route,pos),-1*id,10);
    }
  }

  double regret_value3 = DARPRegretInsertionInitC::computeReqRegretValue(*gen,reqtmp3,ins_route2);

  ASSERT_TRUE(regret_value3 < regret_value1); // It should be equally bad in all the routes, so regret value = 0
  ASSERT_EQ(0,regret_value3);
  ASSERT_NE(-1,ins_route2);
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  GAGenome::optCriterion(GAGenome::MINIMIZATION);

  setUpTestingEnvironment();

  return RUN_ALL_TESTS();
}


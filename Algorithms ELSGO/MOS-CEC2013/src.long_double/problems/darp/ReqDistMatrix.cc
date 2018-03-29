#include "ReqDistMatrix.h"
#include <genomes/GAGenome.h>
#include <assert.h>
#include <algorithm>

ReqDistMatrix::ReqDistMatrix(CostMatrix& costmatrix, DARPEvaluator& eval,
                             long double TWVWeight, long double rideVWeight, long double loadVWeight, bool useconstrpen,
                             bool usewaitingpen, long waitingPenThreshold, long double waitingPenConstant) :
                               costmatrix_(costmatrix),
                               eval_(eval),
                               TWVWeight_(TWVWeight),
                               rideVWeight_(rideVWeight),
                               loadVWeight_(loadVWeight),
                               useconstrpen_(useconstrpen),
                               usewaitingpen_(usewaitingpen),
                               waitingPenThreshold_(waitingPenThreshold),
                               waitingPenConstant_(waitingPenConstant)  {}

ReqDistMatrix::ReqDistMatrix(vector<Request>& reqs, CostMatrix& costmatrix, DARPEvaluator& eval,
                             long double TWVWeight, long double rideVWeight, long double loadVWeight, bool useconstrpen,
                             bool usewaitingpen, long waitingPenThreshold, long double waitingPenConstant) :
                               costmatrix_(costmatrix),
                               eval_(eval),
                               TWVWeight_(TWVWeight),
                               rideVWeight_(rideVWeight),
                               loadVWeight_(loadVWeight),
                               useconstrpen_(useconstrpen),
                               usewaitingpen_(usewaitingpen),
                               waitingPenThreshold_(waitingPenThreshold),
                               waitingPenConstant_(waitingPenConstant)  {
  setUpDistMatrix(reqs);
  initializeDistMatrix(reqs);
}


ReqDistMatrix::ReqDistMatrix(CostMatrix& costmatrix, DARPEvaluator& eval, bool useconstrpen,
                             bool usewaitingpen, long waitingPenThreshold, long double waitingPenConstant) :
                               costmatrix_(costmatrix),
                               eval_(eval),
                               TWVWeight_(eval.TWVWeight()),
                               rideVWeight_(eval.rideVWeight()),
                               loadVWeight_(eval.loadVWeight()),
                               useconstrpen_(useconstrpen),
                               usewaitingpen_(usewaitingpen),
                               waitingPenThreshold_(waitingPenThreshold),
                               waitingPenConstant_(waitingPenConstant) {}


void ReqDistMatrix::setUpDistMatrix(vector<Request>& reqs) {
  clear();
  resize(reqs.size());
}

void ReqDistMatrix::initializeDistMatrix(vector<Request>& reqs) {
  // Initialized as specified in cluster.h
  for (int i=0; i<reqs.size(); i++) {
    if (i>0) resize(i,i);
    reqid2matrixpos_[reqs[i].pickup_vert->id_] = i;
    for (int j=0; j<i; j++) {
      setDist(i,j,bestScoreFromAllCombs(reqs[i],reqs[j]));
    }
  }
}

inline bool tmpCompareScores(long double s1, long double s2) {
  return GAGenome::compareScores(s1,s2) == GAGenome::BETTER;
}

long double ReqDistMatrix::bestScoreFromAllCombs(Request& req1, Request& req2) {
  vector<long double> scores;

  // The first position is not represented since it is always 0
  scores.push_back( reqSequenceScore(req1, req2, 1, 2, 3) ); // R1+ R1- R2+ R2-
  scores.push_back( reqSequenceScore(req1, req2, 2, 1, 3) ); // R1+ R2+ R1- R2-
  scores.push_back( reqSequenceScore(req1, req2, 3, 1, 2) ); // R1+ R2+ R2- R1-

  scores.push_back( reqSequenceScore(req2, req1, 1, 2, 3) ); // R2+ R2- R1+ R1-
  scores.push_back( reqSequenceScore(req2, req1, 2, 1, 3) ); // R2+ R1+ R2- R1-
  scores.push_back( reqSequenceScore(req2, req1, 3, 1, 2) ); // R2+ R1+ R1- R2-

  long double res = * min_element(scores.begin(),scores.end(),tmpCompareScores);

#ifdef DEBUG
  for (int i=0; i<scores.size(); i++) assert(res <= scores[i]); // Just for debugging with minimization problems
#endif

  return res;
}


/************** All these methods have been developed for optimizing the computation of the distance matrix between pair of  */
/************** requests. These are extremely adapted methods for the proposed problem and so they lack reusability and      */
/************** duplicate                                                                                                    */

/*
 * TODO: Duplicated code with the one in DARPEvaluator. All these code should be refactored
 */
void computeTimes(vector<Vertex*>& vertices,int startpos, CostMatrix& costmatrix,
                        /*out*/ vector<long>& arrivalTimes, /*out*/ vector<long>& departureTimes, /*out*/ vector<long>& waitingTimes) {
  assert( vertices.size() > 0);
  assert( startpos >= 0 && startpos < vertices.size() );
  assert( arrivalTimes.size()   == vertices.size() );
  assert( departureTimes.size() == vertices.size() );
  assert( waitingTimes.size()   == vertices.size() );

  if (startpos == 0) {
    arrivalTimes[0] = departureTimes[0] = vertices[0]->fbegin_;
    waitingTimes[0] = 0;
    startpos++;
  }
  assert(startpos > 0);

  for (int i=startpos; i<vertices.size(); i++) {
    arrivalTimes[i]   = departureTimes[i-1] + costmatrix.getCost(* vertices[i-1],* vertices[i] );
    departureTimes[i] = max( vertices[i]->fbegin_, arrivalTimes[i] );
    waitingTimes[i]   = departureTimes[i] - arrivalTimes[i];
  }
}

inline vector<long> computeTimeSlack(vector<Vertex*>& vertices, vector<long>& departureTimes, vector<long>& waitingTimes) {
  vector<long> slackTimes = vector<long>(vertices.size());

  int lastpos = vertices.size()-1;

  slackTimes[lastpos] = vertices[lastpos]->fend_ - departureTimes[lastpos];
  for (int i=lastpos-1;i>=0;i--) {
    slackTimes[i] = min( vertices[i]->fend_ - departureTimes[i], slackTimes[i+1]+waitingTimes[i+1]);
  }

  return slackTimes;
}

inline long computeTimeWindowViolation(vector<Vertex*>& vertices, vector<long>& departureTimes, int pos) {
  assert(pos >=0 && pos <vertices.size());

  return max(0L, departureTimes[pos]-vertices[pos]->fend_);
}

inline long computeTimeWindowViolation(vector<Vertex*>& vertices, vector<long>& departureTimes) {
  assert(vertices.size() == 4);
  long res=0;
  for (int i=1; i<vertices.size(); i++) {
    res += computeTimeWindowViolation(vertices,departureTimes,i);
  }

  return res;
}

/*
 * Note we are assuming that the difference between the times of a request contains the maximum ride time. This is faster than
 * the method used in DARPEvaluator but needs this strong hypothesis.
 */
inline long computeRideTimeViolation(vector<Vertex*>& vertices, vector<long>& departureTimes, vector<long>& arrivalTimes, int pickpos, int delpos) {
  assert(vertices[delpos]->type_ == Vertex::DELIVERY && vertices[pickpos]->type_ == Vertex::PICKUP);
  assert(delpos >=0  && delpos  < vertices.size());
  assert(pickpos >=0 && pickpos < vertices.size());

  long rideTime    = arrivalTimes[delpos]    - departureTimes[pickpos];
  long maxrideTime = vertices[delpos]->fend_ - vertices[pickpos]->fend_;
  assert(maxrideTime == (vertices[delpos]->fend_ - vertices[pickpos]->fend_) );

  return max(0L,rideTime - maxrideTime );
}

inline long computeLoadViolation(Request& req1, Request& req2, int vcapacity) {
  return max(0,req1.pickup_vert->load_+req2.pickup_vert->load_-vcapacity);
}

inline void computeCostValues(vector<Vertex*>& vertices, vector<long>& arrivalTimes, vector<long>& departureTimes,
                              CostMatrix& costmatrix, DARPEvaluator& eval,
                             /*out*/ long double& pickupDelay, /*out*/ long double& deliveryDelay, /*out*/ long& cost) {
  pickupDelay = deliveryDelay = cost = 0;

  for (int i=0; i < vertices.size(); i++) {
    if (i>0) cost += (long) costmatrix.getCost(* vertices[i-1], * vertices[i]);

    if (vertices[i]->isPickUp()) pickupDelay += eval.scalePickupValue(departureTimes[i],*vertices[i]);

    if (vertices[i]->isDelivery()) {
      for (int j=0; j<i; j++) {
        if (vertices[j]->id_ == Vertex::getSiblingVertId( vertices[i]->id_) ) {
          Vertex& del_vert  = *vertices[i]; assert(del_vert.type_  == Vertex::DELIVERY);
          Vertex& pick_vert = *vertices[j]; assert(pick_vert.type_ == Vertex::PICKUP);

          deliveryDelay += eval.scaleDeliveryDelay(pick_vert, del_vert, departureTimes[j], arrivalTimes[i] );
          break;
        }
      }
    }

  }

}

inline long double computeWaitingTimePenalization(vector<long>& waitingTimes,long threshold, long double penconstant) {
  long double res=0;
  for (int i=0; i<waitingTimes.size(); i++) {
    if (waitingTimes[i] > threshold) {
      res += (waitingTimes[i] - threshold) * penconstant;
    }
  }

  return res;
}

long double ReqDistMatrix::reqSequenceScore(Request& req1, Request& req2, int req1d, int req2p, int req2d) { // Positions of the inserted vertices.
  int req1p = 0; // req1 should be before than req2, therefore its pickup vertex should always be at the 0 position

  assert(req1p != req1d && req1p != req2p && req1p != req2d);
  assert(req1p < req2p); // Needed for later steps (adjuting ridetime). Req1 should be placed before req2

  vector<Vertex*> vertices(4);
  vertices[req1p] = req1.pickup_vert;
  vertices[req1d] = req1.delivery_vert;
  vertices[req2p] = req2.pickup_vert;
  vertices[req2d] = req2.delivery_vert;

  vector<long> arrivalTimes(4);
  vector<long> departureTimes(4);
  vector<long> waitingTimes(4);
  computeTimes(vertices,req1p, costmatrix_, arrivalTimes, departureTimes, waitingTimes);

  long pickupviolation       = computeTimeWindowViolation(vertices,departureTimes);
  long loadviolation         = computeLoadViolation(req1,req2,eval_.vehicleCapacity());
  long rideTimeViolationReq1 = computeRideTimeViolation(vertices,departureTimes,arrivalTimes,req1p,req1d);
  long rideTimeViolationReq2 = computeRideTimeViolation(vertices,departureTimes,arrivalTimes,req2p,req2d);

  if (pickupviolation == 0 and loadviolation == 0) {

    // Note1: This has been added so that the computation is the same as the old code but i think that it could be eliminated
    vector<long> timeslacks = computeTimeSlack(vertices, departureTimes, waitingTimes);
    long sumW=0; for (int i=0;i<waitingTimes.size();i++) sumW += waitingTimes[i];
    departureTimes[req1p] = arrivalTimes[req1p] = departureTimes[req1p] + min(timeslacks[req1p],sumW);
    computeTimes(vertices, req1p+1, costmatrix_, arrivalTimes, departureTimes, waitingTimes);
    rideTimeViolationReq1 = computeRideTimeViolation(vertices,departureTimes,arrivalTimes,req1p,req1d);
    rideTimeViolationReq2 = computeRideTimeViolation(vertices,departureTimes,arrivalTimes,req2p,req2d);
    //end of code that could be eliminated

    // There is no rideTimeViolation in the R1+R1-R2+R2- case
    if (  vertices[0]->id_ != Vertex::getSiblingVertId( vertices[1]->id_ ) &&
        ( rideTimeViolationReq1 > 0 || rideTimeViolationReq2 > 0) ) {

      vector<long> timeslacks = computeTimeSlack(vertices, departureTimes, waitingTimes);

      // If the previous part (below Note1) is removed in the future, this code should be activated
      //if (rideTimeViolationReq1 > 0 && timeslacks[req1p] > 0) {
      //  new line departureTimes[req1p] = arrivalTimes[req1p] = departureTimes[req1p] + timeslacks[req1p];
      //  computeTimes(vertices, req1p+1, costmatrix_, arrivalTimes, departureTimes, waitingTimes);
      //}

      // If there is a violation in the ride time of the request 1 we do not do anything
      if (rideTimeViolationReq2 > 0 && timeslacks[req2p] > 0 && rideTimeViolationReq1 == 0) {
        // The maximum available time to postpone the pickup vertex is also fixed by the maximum ride time of req1
        // which goes before
        long maxpostponed = min( timeslacks[req2p], vertices[req1d]->fend_ - departureTimes[req1p]);
        assert(maxpostponed >= 0);

        departureTimes[req2p] = arrivalTimes[req2p] + maxpostponed;
        waitingTimes[req2d]   = departureTimes[req2p] - arrivalTimes[req2p];
        computeTimes(vertices, req2p+1, costmatrix_, arrivalTimes, departureTimes, waitingTimes);
      }

      rideTimeViolationReq1 = computeRideTimeViolation(vertices,departureTimes,arrivalTimes,req1p,req1d);
      rideTimeViolationReq2 = computeRideTimeViolation(vertices,departureTimes,arrivalTimes,req2p,req2d);
    }

  }

  long rideTimeViolation = rideTimeViolationReq1 + rideTimeViolationReq2;
  assert( (vertices[0]->id_ != Vertex::getSiblingVertId( vertices[1]->id_) ) || rideTimeViolation == 0 ); // In the A+A-B+B- case there could not be any ride time violation

  long double pickupDelay, deliveryDelay;
  long cost;
  computeCostValues(vertices, arrivalTimes, departureTimes, costmatrix_, eval_, pickupDelay, deliveryDelay, cost);
  long double nonpenscore = eval_.nonPenalizedScore(cost,pickupDelay,deliveryDelay);

  long double resvalue = (useconstrpen_ == true)
                    ? eval_.score(nonpenscore,pickupviolation,rideTimeViolation,loadviolation, TWVWeight_,rideVWeight_,loadVWeight_)
                    : nonpenscore;

  // We added a penalization to the waiting times in order to avoid good values for pairs of requests where there is a
  // considerable waiting time between them

  if (usewaitingpen_) resvalue += computeWaitingTimePenalization(waitingTimes, waitingPenThreshold_, waitingPenConstant_);

  return resvalue;
}

void ReqDistMatrix::clear() {
  reqid2matrixpos_.clear();
  DistMatrix::clear();
}



///*
// * Old way of computing the scores
// */

//   DARPGenome   gen_; // Needs this attribute

// long double (KMedoidsReqDistributor::*scorefunc)(Request& req1, Request& req2); attribute to select between the old and the new opt way

//inline long double ReqDistMatrix::reqsBestScore(Request& req1, Request& req2) {
//  vector<long double> scores;
//
//  scores.push_back( seqScore(*req1.pickup_vert, *req1.delivery_vert, *req2.pickup_vert,   *req2.delivery_vert )  );
//  scores.push_back( seqScore(*req1.pickup_vert, *req2.pickup_vert,   *req1.delivery_vert, *req2.delivery_vert )  );
//  scores.push_back( seqScore(*req1.pickup_vert, *req2.pickup_vert,   *req2.delivery_vert, *req1.delivery_vert )  );
//
//  scores.push_back( seqScore(*req2.pickup_vert, *req2.delivery_vert, *req1.pickup_vert,   *req1.delivery_vert )  );
//  scores.push_back( seqScore(*req2.pickup_vert, *req1.pickup_vert,   *req2.delivery_vert, *req1.delivery_vert )  );
//  scores.push_back( seqScore(*req2.pickup_vert, *req1.pickup_vert,   *req1.delivery_vert, *req2.delivery_vert )  );
//
//
//  long double res = * min_element(scores.begin(),scores.end(),tmpCompareScores);
//
//#ifdef DEBUG
//  for (int i=0; i<scores.size(); i++) assert(res <= scores[i]); // Just for debugging with minimization problems
//#endif
//
//  return res;
//}


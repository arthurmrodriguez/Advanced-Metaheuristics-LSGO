#include "DARPEvaluator.h"

#include "CostMatrix.h"
#include "DistMatrix.h"
#include "DARPGenome.h"
#include <logger/GALogger.h>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <algorithm>

#include <time.h>


DARPEvaluator::DARPEvaluator(const long vehicleCapacity,
                             const VerticesList& vertlist, const CostMatrix& costMatrix,
                             int timeSumExp,
                             long double loadVW, long double TWVW, long double rideVW,
                             darpOptCriterionType optcrit) :
                                vehicleCapacity_(vehicleCapacity),
                                vertlist_(vertlist),
                                costMatrix_(costMatrix),
                                timeSumExp_(timeSumExp),
                                loadVWeight_(loadVW),
                                TWVWeight_(TWVW),
                                rideVWeight_(rideVW),
                                optCrit_(optcrit),
                                numRouteCalls_(0),
                                numRouteCallsStatsDispFreq_(numeric_limits<int>::max()) {}

int DARPEvaluator::requestsNum() const { return vertlist_.requestsNum(); }

bool DARPEvaluator::evalRoute(const vector<int>& route,
                              long& cost,
                              long& loadViolation, long& TWViolation, long& rideViolation,
                              long double& pickupDelay, long double& deliveryDelay,
                              bool debug,
                              bool debug2) {
  // Each of these variables stores the corresponding value for each vertex of the route
  map<int,long> arrivalTimes;
  map<int,long> waitingTimes;
  map<int,long> beginningServiceTimes;
  map<int,long> departureTimes;
  map<int,long> rideTimes;
  map<int,long> forwardTimeSlacks;
  map<int,long> loadsWhenLeavingVertex;

  return evalRoute(route,arrivalTimes,waitingTimes,beginningServiceTimes,departureTimes,rideTimes,forwardTimeSlacks,loadsWhenLeavingVertex,cost,loadViolation,TWViolation,rideViolation,pickupDelay,deliveryDelay,debug,debug2);
}

bool DARPEvaluator::evalRoute(const vector<int>& route,
                              map<int,long>&     arrivalTimes,
                              map<int,long>&     waitingTimes,
                              map<int,long>&     beginningServiceTimes,
                              map<int,long>&     departureTimes,
                              map<int,long>&     rideTimes,
                              map<int,long>&     forwardTimeSlacks,
                              map<int,long>&     loadsWhenLeavingVertex,
                              long&                   cost,
                              long&                   loadViolation,
                              long&                   TWViolation,
                              long&                   rideViolation,
                              long double&                 pickupDelay,
                              long double&                 deliveryDelay,
                              bool                    debug,
                              bool                    debug2) {
  // Just to be sure...
  cost          = 0;
  loadViolation = 0;
  TWViolation   = 0;
  rideViolation = 0;
  pickupDelay   = 0;
  deliveryDelay = 0;

  /*LOG*/ if (debug) cout << "=================================================================================================================" << endl << endl;

  // If vehicle is empty, return
  if (route.size() == 0) return true;

  /*LOG*/ if (debug) {
  /*LOG*/  cout << "L: " << Vertex::maxRideValue__  << endl;
  /*LOG*/  cout << "Path: ";
  /*LOG*/  for (unsigned i = 0; i < route.size(); i++)
  /*LOG*/    cout << route[i] << " ";
  /*LOG*/  cout << endl << endl;
  /*LOG*/}

  // 8-step evaluation scheme

  // Step 1: D0 = e0
  evalRouteStep1(route,arrivalTimes,waitingTimes,beginningServiceTimes,departureTimes,loadsWhenLeavingVertex);

  // Step 2: Compute arrivalTimes, waitingTimes, begining of Service Times and departureTimes for each vertex i on the route
  bool validSolution = evalRouteStep2(route,arrivalTimes,waitingTimes,beginningServiceTimes,departureTimes,rideTimes,loadsWhenLeavingVertex,debug);

  if (validSolution) {
    // Step 3: Compute Forward time slack at 0 (F0)
    evalRouteStep3(route,waitingTimes,beginningServiceTimes,arrivalTimes,departureTimes,forwardTimeSlacks,debug);


    // Step 4: Update Departure time time at 0 (D0)
    long sumW = evalRouteStep4(route,arrivalTimes,waitingTimes,beginningServiceTimes,departureTimes,rideTimes,forwardTimeSlacks,loadsWhenLeavingVertex,debug);


    // Step 5: Update arrival, waiting, beginning of service and departure times (A, W, B and D) for each vertex on the route
    evalRouteStep5(route,arrivalTimes,waitingTimes,beginningServiceTimes,departureTimes,rideTimes,forwardTimeSlacks,loadsWhenLeavingVertex,debug);


    // Step 6: Compute Li ride time for each request on the route
    bool lengthExceeded = evalRouteStep6(route,arrivalTimes,waitingTimes,beginningServiceTimes,departureTimes,rideTimes,forwardTimeSlacks,loadsWhenLeavingVertex,debug);

    if (lengthExceeded) {
      // Step 7: For every vertex j that is an origin (except the first one, that we have already computed)
      //   compute Fj and Wj
      //   update Ai,Wi,Bi and Di for each vertex i that comes after j in the route
      //   update Li for each request i whose destination is after j
      evalRouteStep7(route,arrivalTimes,waitingTimes,beginningServiceTimes,departureTimes,rideTimes,forwardTimeSlacks,loadsWhenLeavingVertex,sumW,debug);
    }
  }

  // Step 8: Compute changes in violations
  bool feasible = evalRouteStep8(route,arrivalTimes,waitingTimes,beginningServiceTimes,departureTimes,rideTimes, forwardTimeSlacks,loadsWhenLeavingVertex,
                                 cost, loadViolation, TWViolation, rideViolation, pickupDelay, deliveryDelay, debug);

  /*LOG*/ if (debug) {
  /*LOG*/   cout << "  The path is " << (feasible ? "feasible" : "unfeasible") << endl;
  /*LOG*/   cout << endl << endl;
  /*LOG*/ }


  /*LOG*/ if (debug2) {
  /*LOG*/    for (unsigned j = 0; j < route.size(); j++) {
  /*LOG*/      int next = (j == route.size() - 1) ? numeric_limits<int>::max() : route[j+1];
  /*LOG*/      printVertexInfo2(route[j], next, arrivalTimes, waitingTimes, beginningServiceTimes, departureTimes, forwardTimeSlacks, rideTimes, loadsWhenLeavingVertex);
  /*LOG*/    }
  /*LOG*/}

  /*LOG*/ if (debug) cout << "=================================================================================================================" << endl << endl;

  numRouteCalls_++;

  if (numRouteCalls_ % numRouteCallsStatsDispFreq_ == 0) GALogger::instance()->appendStats("DARPEvaluator numroutecalls stats");

  return feasible;
}

/*
 * Initialize values
 */
inline void DARPEvaluator::evalRouteStep1(const vector<int>& route,
                                          map<int,long>&     arrivalTimes,
                                          map<int,long>&     waitingTimes,
                                          map<int,long>&     beginningServiceTimes,
                                          map<int,long>&     departureTimes,
                                          map<int,long>&     loadsWhenLeavingVertex) {
  Vertex& first_vertex = vertlist_.getVertex(route[0]);

  departureTimes        [route[0]] = first_vertex.fbegin_;
  arrivalTimes          [route[0]] = departureTimes[route[0]];
  beginningServiceTimes [route[0]] = departureTimes[route[0]];
  waitingTimes          [route[0]] = 0;
  loadsWhenLeavingVertex[route[0]] = first_vertex.load_;
}

/*
 * Compute arrivalTimes, waitingTimes, begining of Service Times and departureTimes for each vertex i on the route
 */
inline bool DARPEvaluator::evalRouteStep2(const vector<int>& route,
                                          map<int,long>&     arrivalTimes,
                                          map<int,long>&     waitingTimes,
                                          map<int,long>&     beginningServiceTimes,
                                          map<int,long>&     departureTimes,
                                          map<int,long>&     rideTimes,
                                          map<int,long>&     loadsWhenLeavingVertex,
                                          bool                    debug) {
  /*LOG*/ if (debug) cout << endl << " ***** Step 2 *****" << endl << endl;

  for (unsigned i=1; i<route.size(); i++) {
    int vert_id      = route[i];
    int prev_vert_id = route[i-1];

    updateArrivalTime     (route,vert_id,prev_vert_id, arrivalTimes, departureTimes);
    updateBegOfServiceTime(route,vert_id,beginningServiceTimes,arrivalTimes);
    updateWaitingTime     (route,vert_id,waitingTimes,         beginningServiceTimes,arrivalTimes);
    updateDepartureTime   (route,vert_id,departureTimes,       beginningServiceTimes);

    Vertex& vert_i = vertlist_.getVertex(vert_id);
    loadsWhenLeavingVertex[vert_id] = (Vertex::isVertIdPickup(vert_id)) ? loadsWhenLeavingVertex[prev_vert_id] + vert_i.load_ : loadsWhenLeavingVertex[prev_vert_id] - vert_i.load_;
  }

  bool validSolution = true;

  for (unsigned i = 0; i < route.size(); i++) {
    rideTimes[route[i]] = rideTime(route[i],beginningServiceTimes,arrivalTimes, departureTimes);

    if (beginningServiceTimes[route[i]]  > vertlist_.getVertex(route[i]).fend_) validSolution = false;
    if (loadsWhenLeavingVertex[route[i]] > vehicleCapacity_                   ) validSolution = false;

    map<int,long> dummyforwardtimeslack;
    /*LOG*/ if (debug) printVertexInfo(route[i], arrivalTimes, waitingTimes, beginningServiceTimes, departureTimes, dummyforwardtimeslack, rideTimes, loadsWhenLeavingVertex);
  }

  return validSolution;
}

/*
 * Compute Forward time slack at 0 (F0)
 */
inline void DARPEvaluator::evalRouteStep3(const vector<int>& route,
                                          map<int,long>&     waitingTimes,
                                          map<int,long>&     beginningServiceTimes,
                                          map<int,long>&     arrivalTimes,
                                          map<int,long>&     departureTimes,
                                          map<int,long>&     forwardTimeSlacks,
                                          bool                    debug) {
  /*LOG*/ if (debug) cout << endl << " ***** Step 3 *****" << endl << endl;

  forwardTimeSlacks[route[0]] = computeForwardTimeSlack(0, route, beginningServiceTimes, arrivalTimes, departureTimes, waitingTimes, debug);
}

/*
 * Update Departure time time at 0 (D0)
 */
inline long DARPEvaluator::evalRouteStep4(const vector<int>& route,
                                          map<int,long>&     arrivalTimes,
                                          map<int,long>&     waitingTimes,
                                          map<int,long>&     beginningServiceTimes,
                                          map<int,long>&     departureTimes,
                                          map<int,long>&     rideTimes,
                                          map<int,long>&     forwardTimeSlacks,
                                          map<int,long>&     loadsWhenLeavingVertex,
                                          bool                    debug) {
  /*LOG*/ if (debug) cout << endl << " ***** Step 4 *****" << endl << endl;

  long sumW = 0;
  map<int,long>::iterator it;
  for (it = waitingTimes.begin(); it != waitingTimes.end(); it++) sumW += it->second;

  /*LOG*/ if (debug) {
  /*LOG*/  cout << "  sumW = " << sumW << endl;
  /*LOG*/  cout << "  minW = min(sumW, forwardTimeSlacks[" << route[0] << "]) = min(" << forwardTimeSlacks[route[0]] << ", " << sumW << ") = " << min(forwardTimeSlacks[route[0]], sumW) << endl;
  /*LOG*/ }

  departureTimes[route[0]] = vertlist_.getVertex(route[0]).fbegin_ + min(forwardTimeSlacks[route[0]], sumW);

  /*LOG*/ if (debug) {
  /*LOG*/  cout << "  D[" << route[0] << "] = e[" << route[0] <<"] + minW = " << departureTimes[route[0]] << endl << endl;
  /*LOG*/  printVertexInfo(route[0], arrivalTimes, waitingTimes, beginningServiceTimes, departureTimes, forwardTimeSlacks, rideTimes, loadsWhenLeavingVertex);
  /*LOG*/}

  // Just to be consistent
  arrivalTimes         [route[0]] = departureTimes[route[0]];
  beginningServiceTimes[route[0]] = departureTimes[route[0]];
  waitingTimes         [route[0]] = beginningServiceTimes[route[0]] - arrivalTimes[route[0]];

  return sumW;
}

/*
 * Update arrival, waiting, beginning of service and departure times (A, W, B and D)
 */
inline void DARPEvaluator::evalRouteStep5(const vector<int>& route,
                                          map<int,long>&     arrivalTimes,
                                          map<int,long>&     waitingTimes,
                                          map<int,long>&     beginningServiceTimes,
                                          map<int,long>&     departureTimes,
                                          map<int,long>&     rideTimes,
                                          map<int,long>&     forwardTimeSlacks,
                                          map<int,long>&     loadsWhenLeavingVertex,
                                          bool                    debug) {
  /*LOG*/ if (debug) cout << endl << " ***** Step 5 *****" << endl << endl;

  for (unsigned i=1; i<route.size(); i++) {
    int vert_id      = route[i];
    int prev_vert_id = route[i-1];

    updateArrivalTime     (route,vert_id,prev_vert_id, arrivalTimes, departureTimes);
    updateBegOfServiceTime(route,vert_id,beginningServiceTimes,arrivalTimes);
    updateWaitingTime     (route,vert_id,waitingTimes,         beginningServiceTimes,arrivalTimes);
    updateDepartureTime   (route,vert_id,departureTimes,       beginningServiceTimes);

    /*LOG*/ if (debug) printVertexInfo(route[i], arrivalTimes, waitingTimes, beginningServiceTimes, departureTimes, forwardTimeSlacks, rideTimes, loadsWhenLeavingVertex);
  }
}

/*
 * Compute user ride time (Li) for each request on the route
 */
inline bool DARPEvaluator::evalRouteStep6(const vector<int>& route,
                                          map<int,long>&     arrivalTimes,
                                          map<int,long>&     waitingTimes,
                                          map<int,long>&     beginningServiceTimes,
                                          map<int,long>&     departureTimes,
                                          map<int,long>&     rideTimes,
                                          map<int,long>&     forwardTimeSlacks,
                                          map<int,long>&     loadsWhenLeavingVertex,
                                          bool                    debug) {
  bool lengthExceeded = false;

  /*LOG*/ if (debug) cout << endl << " ***** Step 6 *****" << endl << endl;

  for (unsigned i = 0; i < route.size(); i++) {
    int vert_id = route[i];

    rideTimes[vert_id] = rideTime(vert_id,beginningServiceTimes,arrivalTimes,departureTimes);

    if ( computeRideTimeToMaxDiff(vert_id, rideTimes[vert_id]) > 0 ) lengthExceeded = true;

    /*LOG*/ if (debug) printVertexInfo(route[i], arrivalTimes, waitingTimes, beginningServiceTimes, departureTimes, forwardTimeSlacks, rideTimes, loadsWhenLeavingVertex);
  }

  return lengthExceeded;
}

/* For every vertex j that is an origin (except the first one, that we have already computed)
 * compute Fj and Wj
 * update Ai,Wi,Bi and Di for each vertex i that comes after j in the route
 * update Li for each request i whose destination is after j
 */
inline void DARPEvaluator::evalRouteStep7(const vector<int>& route,
                                          map<int,long>&     arrivalTimes,
                                          map<int,long>&     waitingTimes,
                                          map<int,long>&     beginningServiceTimes,
                                          map<int,long>&     departureTimes,
                                          map<int,long>&     rideTimes,
                                          map<int,long>&     forwardTimeSlacks,
                                          map<int,long>&     loadsWhenLeavingVertex,
                                          long                    sumW,
                                          bool                    debug) {
  /*LOG*/ if (debug) cout  << endl << " ***** Step 7 *****" << endl << endl;

  bool lengthExceeded;

  for (unsigned j = 1; j < route.size(); j++) {
    if (Vertex::isVertIdPickup(route[j])) {
      lengthExceeded = false;

      // Compute Fj
      forwardTimeSlacks[route[j]] = computeForwardTimeSlack(j, route, beginningServiceTimes, arrivalTimes, departureTimes, waitingTimes, debug);

      // Update Wj, Bj and Dj
      int vert_id = route[j];
      waitingTimes         [vert_id] = waitingTimes[vert_id] + min(forwardTimeSlacks[vert_id], sumW);
      beginningServiceTimes[vert_id] = arrivalTimes[vert_id] + waitingTimes[vert_id];
      departureTimes       [vert_id] = beginningServiceTimes[vert_id];

      /*LOG*/ if (debug) {
      /*LOG*/  cout << endl;
      /*LOG*/  cout << "route=" << j << endl;
      /*LOG*/  printVertexInfo(route[j], arrivalTimes, waitingTimes, beginningServiceTimes, departureTimes, forwardTimeSlacks, rideTimes, loadsWhenLeavingVertex);
      /*LOG*/  cout << endl;
      /*LOG*/}

      // Update Wi, Ai, Bi and Di for all i > j
      for (unsigned i = j + 1; i < route.size(); i++) {
        int vert_id      = route[i];
        int prev_vert_id = route[i-1];

        updateArrivalTime     (route,vert_id,prev_vert_id, arrivalTimes, departureTimes);
        updateBegOfServiceTime(route,vert_id,beginningServiceTimes,arrivalTimes);
        updateWaitingTime     (route,vert_id,waitingTimes,         beginningServiceTimes,arrivalTimes);
        updateDepartureTime   (route,vert_id,departureTimes,       beginningServiceTimes);

        if (Vertex::isVertIdDelivery(vert_id)) {
          int vert_id = route[i];
          rideTimes[vert_id] = rideTime(vert_id,beginningServiceTimes, arrivalTimes, departureTimes);

          if ( computeRideTimeToMaxDiff(vert_id,rideTimes[vert_id]) > 0 ) lengthExceeded = true;
       }

        /*LOG*/ if (debug) printVertexInfo(vert_id, arrivalTimes, waitingTimes, beginningServiceTimes, departureTimes, forwardTimeSlacks, rideTimes, loadsWhenLeavingVertex);
      }

      /*LOG*/ if (debug) cout << endl;

      if (!lengthExceeded) break;
    }
  }
}

/*
 * Compute changes in violations
 */
inline bool DARPEvaluator::evalRouteStep8(const vector<int>& route,
                                          map<int,long>&     arrivalTimes,
                                          map<int,long>&     waitingTimes,
                                          map<int,long>&     beginningServiceTimes,
                                          map<int,long>&     departureTimes,
                                          map<int,long>&     rideTimes,
                                          map<int,long>&     forwardTimeSlacks,
                                          map<int,long>&     loadsWhenLeavingVertex,
                                          long&                   cost,
                                          long&                   loadViolation,
                                          long&                   TWViolation,
                                          long&                   rideViolation,
                                          long double&                 pickupDelay,
                                          long double&                 deliveryDelay,
                                          bool                    debug) {
  /*LOG*/ if (debug) cout  << endl << " ***** Step 8 *****" << endl << endl;

  cost = 0;
  for (unsigned i = 0; i < route.size(); i++) {
    if (i > 0) cost += (long) costMatrix_.getCost(route[i-1], route[i]);

    /*LOG*/ if (debug) printVertexInfo(route[i], arrivalTimes, waitingTimes, beginningServiceTimes, departureTimes, forwardTimeSlacks, rideTimes, loadsWhenLeavingVertex);
  }

  loadViolation = computeLoadViolation      (loadsWhenLeavingVertex);
  TWViolation   = computeTimeWindowViolation(beginningServiceTimes);
  rideViolation = computeRideViolation      (rideTimes);

  pickupDelay   = computePickupDelay        (beginningServiceTimes); assert(!isnan(pickupDelay));
  deliveryDelay = computeDeliveryDelay      (beginningServiceTimes, arrivalTimes);

  /*LOG*/ if (debug) {
  /*LOG*/  cout << endl;
  /*LOG*/  cout << "  Total cost: " << cost << endl;
  /*LOG*/  cout << "  Load violation: " << loadViolation << endl;
  /*LOG*/  cout << "  TW violation: " << TWViolation << endl;
  /*LOG*/  cout << "  Ride violation: " << rideViolation << endl;
  /*LOG*/  cout << endl;
  /*LOG*/}
  assert(loadViolation >= 0);
  assert(TWViolation   >= 0);
  assert(rideViolation >= 0);

  bool feasible = (loadViolation == 0) && (TWViolation == 0) && (rideViolation == 0);

  return feasible;
}

inline long DARPEvaluator::rideTime(int vert_id, map<int,long>& beginningServiceTimes, map<int,long>& arrivalTimes, map<int,long>& departureTimes ) {
  // old line follows the original Cordeau and Parragh scheme but for our purposes there is no sense in waiting when you are droping someone
  // return (Vertex::isVertIdPickup(vert_id)) ? 0 : beginningServiceTimes[vert_id] - departureTimes[Vertex::getSiblingVertId(vert_id)];
  return (Vertex::isVertIdPickup(vert_id)) ? 0 : arrivalTimes[vert_id] - departureTimes[Vertex::getSiblingVertId(vert_id)];
}

inline void DARPEvaluator::updateArrivalTime(const vector<int>& route, int vert_id, int prev_vert_id, map<int,long>& arrivalTimes, map<int,long>& departureTimes) {
  arrivalTimes [vert_id] = departureTimes[prev_vert_id] + (long) costMatrix_.getTravelTime(prev_vert_id, vert_id);
}

inline void DARPEvaluator::updateBegOfServiceTime(const vector<int>& route, int vert_id, map<int,long>& beginningServiceTimes, map<int,long>& arrivalTimes) {
  Vertex& vert_i = vertlist_.getVertex(vert_id);

  beginningServiceTimes [vert_id] = max(arrivalTimes[vert_id], vert_i.fbegin_);
}

inline void DARPEvaluator::updateWaitingTime(const vector<int>& route, int vert_id, map<int,long>& waitingTimes, map<int,long>& beginningServiceTimes, map<int,long>& arrivalTimes) {
  waitingTimes [vert_id] = beginningServiceTimes[vert_id] - arrivalTimes[vert_id];
}

inline void DARPEvaluator::updateDepartureTime(const vector<int>& route, int vert_id, map<int,long>& departureTimes, map<int,long>& beginningServiceTimes) {
  departureTimes[vert_id] = beginningServiceTimes[vert_id]; // We assume service time = 0
}

long DARPEvaluator::computeForwardTimeSlack(int                vertexPos,
                                            const vector<int>& route,
                                            map<int,long>&     beginningServiceTimes,
                                            map<int,long>&     arrivalTimes,
                                            map<int,long>&     departureTimes,
                                            map<int,long>&     waitingTimes,
                                            bool                    debug) {

  vector<long> vals;

  /*LOG*/ if (debug) cout << " ++++ Computing F[" << route[vertexPos] << "] ++++" << endl;

  for (unsigned j = vertexPos; j <= route.size() - 1; j++) {
    long sumW = 0;

    /*LOG*/ if (debug) cout << "  Processing vertex: " << setw(6) << route[j];

    for (unsigned p = vertexPos + 1; p <= j; p++) sumW += waitingTimes[route[p]];

    /*LOG*/ if (debug) cout << " | SumWp = " << setw(6) << sumW;

    Vertex& reqj = vertlist_.getVertex(route[j]);

    long lj = reqj.fend_;
    long Pj = 0;

    if (Vertex::isVertIdDelivery(reqj.id_)) {
      long sibling = Vertex::getSiblingVertId(reqj.id_);
      bool found = false;

      for (unsigned k = 0; k < route.size() && route[k] != route[vertexPos]; k++) {
        if (route[k] == sibling) found = true;
      }

      if (found) Pj = rideTime(route[j],beginningServiceTimes,arrivalTimes, departureTimes);
    }

    /*LOG*/ if (debug) {
    /*LOG*/  cout << " | Pj = " << setw(6) << Pj;
    /*LOG*/  cout << " | L = " << Vertex::maxRideTime(route[j],costMatrix_);
    /*LOG*/  cout << " | lj - Bj = " << setw(6) << lj - beginningServiceTimes[route[j]];
    /*LOG*/  cout << " | L - Pj = " << setw(6) << Vertex::maxRideTime(route[j],costMatrix_) - Pj;
    /*LOG*/}

    long diff = min(lj - beginningServiceTimes[route[j]], Vertex::maxRideTime(route[j],costMatrix_)  - Pj);

    long val = sumW + diff;

    vals.push_back(val);

    /*LOG*/ if (debug) cout << "   ===>   val = " << setw(6) << val << endl;
  }

  long res = *min_element(vals.begin(), vals.end());

  res = max(0L, res);

  /*LOG*/ if (debug) {
  /*LOG*/   cout << endl;
  /*LOG*/   cout << "  F[" << route[vertexPos] << "] = " << res << endl;
  /*LOG*/   cout << " ++++ Finished Computing F[" << route[vertexPos] << "] ++++" << endl;
  /*LOG*/ }

  return res;
}


long DARPEvaluator::computeLoadViolation(map<int,long> loadsWhenLeavingVertex) {
  long res = 0;
  map<int,long>::const_iterator it;

  for (it = loadsWhenLeavingVertex.begin(); it != loadsWhenLeavingVertex.end(); it++) {
    res += max(0L, it->second - vehicleCapacity_);
  }

  return res;
}


long DARPEvaluator::computeTimeWindowViolation(map<int,long> beginningServiceTimes) {
  long res = 0;
  map<int,long>::const_iterator it;

  for (it = beginningServiceTimes.begin(); it != beginningServiceTimes.end(); it++) {
    long fend = vertlist_.getVertex(it->first).fend_;
    res += max(0L, it->second - fend);
  }
  assert(res >= 0);

  return res;
}

long DARPEvaluator::computeRideViolation(map<int,long> rideTimes) {
  long res = 0;
  map<int,long>::const_iterator it;

  for (it = rideTimes.begin(); it != rideTimes.end(); it++) {
    res += computeRideTimeToMaxDiff(it->first,it->second);
  }

  return res;
}

long DARPEvaluator::computeRideTimeToMaxDiff(int vert_id, long ridetime) {
  return max( 0L, ridetime - Vertex::maxRideTime(vert_id,costMatrix_) );
}

long double DARPEvaluator::computePickupDelay(map<int,long> beginningServiceTimes) {
  long double res = 0;
  map<int,long>::const_iterator it;
  int vert_id;

  for (it = beginningServiceTimes.begin(); it != beginningServiceTimes.end(); it++) {
    vert_id        = it->first;
    Vertex& vertex = vertlist_.getVertex(vert_id);

    if ( vertex.isPickUp() ) {
      int begin_time = it->second;
      assert(begin_time >= vertex.fbegin_);
      // The original values are scaled so that the minimum is 0 and the maximum corresponds to the
      // maximum pickup time constraint which is set earlier in the vertex.fend_ window
      res += scalePickupValue(begin_time,vertex);
    }
  }

  return res;
}

long double DARPEvaluator::computeDeliveryDelay(map<int,long> beginningServiceTimes, map<int, long> arrivalTimes) {
  long double res = 0;
  map<int,long>::const_iterator it;
  int vert_id, pick_vert_id;
  long direct_dist, begin_time, del_time, service_time;

  for (it = beginningServiceTimes.begin(); it != beginningServiceTimes.end(); it++) {
    vert_id        = it->first;
    Vertex& vertex = vertlist_.getVertex(vert_id);

    if ( vertex.isDelivery() ) {
      int    pick_vert_id  = Vertex::getSiblingVertId(vert_id);
      Vertex pick_vertex   = vertlist_.getVertex(pick_vert_id); 

      begin_time   = beginningServiceTimes[pick_vert_id];
      del_time     = arrivalTimes[vert_id];
                                             
      res += scaleDeliveryDelay(pick_vertex,vertex,begin_time,del_time);
    }
  }

  return res;
}


void DARPEvaluator::printVertexInfo(int vertex, map<int,long>& A, map<int,long>& W, map<int,long>& B, map<int,long>& D, map<int,long>& F, map<int,long>& L, map<int,long>& y) {
  cout << "  Vertex: " << setw(6) << vertex << " | A = " << setw(10) << A[vertex] << " | W = " << setw(10) << W[vertex];
  cout << " | B = " << setw(10) << B[vertex] << " | D = " << setw(10) << D[vertex] << " | F = " << setw(6) << F[vertex];
  cout << " | L = " << setw(6) << L[vertex] << " | y = " << setw(3) << y[vertex];
  cout << " | ej = " << setw(10) << vertlist_.getVertex(vertex).fbegin_ << " | lj = " << setw(10) << vertlist_.getVertex(vertex).fend_;
  cout << endl;
}


void DARPEvaluator::printVertexInfo2(int vertex, int nextvert, map<int,long>& A, map<int,long>& W, map<int,long>& B, map<int,long>& D, map<int,long>& F, map<int,long>& L, map<int,long>& y) {
  Vertex req1 = vertlist_.getVertex(vertex);
  Vertex req2 = vertlist_.getVertex(Vertex::getSiblingVertId(vertex));

  string Atime = asctime(localtime(&A[vertex]));
  Atime = Atime.substr(11, 8);

  string Dtime = asctime(localtime(&D[vertex]));
  Dtime = Dtime.substr(11, 8);

  string ejtime = asctime(localtime(&req1.fbegin_));
  ejtime = ejtime.substr(11, 8);

  string ljtime = asctime(localtime(&req1.fend_));
  ljtime = ljtime.substr(11, 8);

  string Ltime = asctime(gmtime(&L[vertex]));
  Ltime = Ltime.substr(11, 8);

  string Wtime = asctime(gmtime(&W[vertex]));
  Wtime = Wtime.substr(11, 8);

  long dist = (Vertex::isVertIdPickup(vertex)) ? (long) costMatrix_.distMatrix().getDist(req1.pos_, req2.pos_) : (long) costMatrix_.distMatrix().getDist(req2.pos_, req1.pos_);
  string dtime = asctime(gmtime(&dist));
  dtime = dtime.substr(11, 8);

  string dtimenext;
  if (nextvert != numeric_limits<int>::max()) {
    Vertex reqnext = vertlist_.getVertex(nextvert);
    long distnext = (long) costMatrix_.distMatrix().getDist(req1.pos_, reqnext.pos_);
    dtimenext = asctime(gmtime(&distnext));
    dtimenext = dtimenext.substr(11, 8);
  }
  else {
    dtimenext = "00:00:00";
  }

  cout << "  Vertex: " << setw(6) << vertex;
  cout << " | A = "  << Atime;
  cout << " | D = "  << Dtime;
  cout << " | ej = " << ejtime;
  cout << " | lj = " << ljtime;
  cout << " | L = "  << Ltime;
  cout << " | d = "  << dtime;
  if (Vertex::isRideTimeRelative() ) cout << " | L/d = " << (long double) L[vertex] / dist;
  cout << " | next = "  << dtimenext;
  cout << " | W = "  << Wtime;
  cout << " | y = "  << setw(2) << y[vertex];
  cout << endl;
}

long double DARPEvaluator::nonPenalizedScore(const DARPGenome& gen) const {
  long double result = 0.0;
  for (int route=0; route<gen.size(); route++) {
     result += nonPenalizedScore(gen,route);
  }

  return result;
}

// Developed for the computation of the optimized version of the distance matrix of the kmedoidsreq
long double DARPEvaluator::nonPenalizedScore(const DARPGenome& gen, int route) const {
  return nonPenalizedScore( gen.costPerRoute(route), gen.pickupDelayPerRoute(route), gen.deliveryDelayPerRoute(route) );
}

long double DARPEvaluator::nonPenalizedScore(long cost, long double pickupdelay, long double deliverydelay) const {
  long double value=-1;

  switch(optCrit_) {
    case TRAVELCOST:
      value = cost;
      break;
    case BOTHDELAYS:
      value = pickupdelay + deliverydelay;
      break;
    case DELIVERY:
      value = deliverydelay;
      break;
    default:
      throw runtime_error("Error: unrecognized darp opt criterion");
  }
  return value;
}

long double DARPEvaluator::score(DARPGenome& gen) {
  assert(!isnan(gen.nonPenalizedScore() ) ); assert(!isnan(loadVWeight_*gen.loadV()) ); assert(!isnan(rideVWeight_*gen.rideV()) ); assert(!isnan(TWVWeight_*gen.TWV()) );

  gen.evaluateAllRoutes();

  return score(gen.nonPenalizedScore(), gen.TWV(), gen.rideV(),gen.loadV() );
}

long double DARPEvaluator::scalePickupValue(long begintime, Vertex& vertex) {
  return pow( (long double) (begintime - vertex.fbegin_) / (long double) (vertex.fend_ - vertex.fbegin_) , timeSumExp_ );
}

long double DARPEvaluator::scaleDeliveryDelay(Vertex& pick_vert, Vertex& del_vert, long begin_time, long del_time) {
  long service_time = del_time - begin_time;  assert(del_time >= begin_time);
  long direct_dist  = (long) costMatrix_.getTravelTime(pick_vert.id_, del_vert.id_);

  assert(costMatrix_.getTravelTime(pick_vert.id_, del_vert.id_) == costMatrix_.distMatrix().getDist(pick_vert.pos_, del_vert.pos_) );

  // To avoid the problem where both the pickup and delivery positions are the same or when there is
  // a small numerical error where the service time is smaller than the direct time by at most 100 seconds
  long double scaled_value = ( (service_time < direct_dist && service_time >= direct_dist - 100) ||
                          direct_dist == 0 ) ?
                            0 :
                            (long double) (service_time                                                - direct_dist) /
                            (long double) (Vertex::maxRideTime(pick_vert.id_,del_vert.id_,costMatrix_) - direct_dist);

  assert(scaled_value >= 0);              // Having a smaller distance than the direct distance should not be possible

  return pow( scaled_value , timeSumExp_ );
}

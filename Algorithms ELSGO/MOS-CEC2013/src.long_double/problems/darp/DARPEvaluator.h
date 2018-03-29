#ifndef PATHEVALUATOR_H
#define PATHEVALUATOR_H

#include <vector>
#include <map>
#include <assert.h>
#include <math.h>

class VerticesList;
class Vertex;
class DistMatrix;
class DARPGenome;
class CostMatrix;

using namespace std;

class DARPEvaluator {

public:
  enum darpOptCriterionType {TRAVELCOST, DELIVERY, BOTHDELAYS};

protected:
  const long           vehicleCapacity_;
  const VerticesList&  vertlist_;
  const CostMatrix&    costMatrix_;
  const int            timeSumExp_; // The power of the sum of the delivery and route times
  long double               loadVWeight_;
  long double               TWVWeight_;
  long double               rideVWeight_;
  darpOptCriterionType optCrit_;
  int                  numRouteCalls_;
  int                  numRouteCallsStatsDispFreq_;

  void evalRouteStep1(const vector<int>& route,
                      map<int,long>&     arrivalTimes,
                      map<int,long>&     waitingTimes,
                      map<int,long>&     beginningServiceTimes,
                      map<int,long>&     departureTimes,
                      map<int,long>&     loadsWhenLeavingVertex);

  bool evalRouteStep2(const vector<int>& route,
                      map<int,long>&     arrivalTimes,
                      map<int,long>&     waitingTimes,
                      map<int,long>&     beginningServiceTimes,
                      map<int,long>&     departureTimes,
                      map<int,long>&     rideTimes,
                      map<int,long>&     loadsWhenLeavingVertex,
                      bool                    debug);

  void evalRouteStep3(const vector<int>& route,
                      map<int,long>&     waitingTimes,
                      map<int,long>&     arrivalTimes,
                      map<int,long>&     beginningServiceTimes,
                      map<int,long>&     departureTimes,
                      map<int,long>&     forwardTimeSlacks,
                      bool                    debug);

  long evalRouteStep4(const vector<int>& route,
                      map<int,long>&     arrivalTimes,
                      map<int,long>&     waitingTimes,
                      map<int,long>&     beginningServiceTimes,
                      map<int,long>&     departureTimes,
                      map<int,long>&     rideTimes,
                      map<int,long>&     forwardTimeSlacks,
                      map<int,long>&     loadsWhenLeavingVertex,
                      bool                    debug);

  void evalRouteStep5(const vector<int>& route,
                      map<int,long>&     arrivalTimes,
                      map<int,long>&     waitingTimes,
                      map<int,long>&     beginningServiceTimes,
                      map<int,long>&     departureTimes,
                      map<int,long>&     rideTimes,
                      map<int,long>&     forwardTimeSlacks,
                      map<int,long>&     loadsWhenLeavingVertex,
                      bool                    debug);

  bool evalRouteStep6(const vector<int>& route,
                      map<int,long>&     arrivalTimes,
                      map<int,long>&     waitingTimes,
                      map<int,long>&     beginningServiceTimes,
                      map<int,long>&     departureTimes,
                      map<int,long>&     rideTimes,
                      map<int,long>&     forwardTimeSlacks,
                      map<int,long>&     loadsWhenLeavingVertex,
                      bool                    debug);

  void evalRouteStep7(const vector<int>& route,
                      map<int,long>&     arrivalTimes,
                      map<int,long>&     waitingTimes,
                      map<int,long>&     beginningServiceTimes,
                      map<int,long>&     departureTimes,
                      map<int,long>&     rideTimes,
                      map<int,long>&     forwardTimeSlacks,
                      map<int,long>&     loadsWhenLeavingVertex,
                      long                    sumW,
                      bool                    debug);

  bool evalRouteStep8(const vector<int>& route,
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
                      bool                    debug);

  long rideTime (int vertex_id, map<int,long>& beginningServiceTimes, map<int,long>& arrivalTimes, map<int,long>& departureTimes );

  void updateArrivalTime     (const vector<int>& route, int vert_id, int prev_vert_id, map<int,long>& arrivalTimes, map<int,long>& departureTimes) ;
  void updateBegOfServiceTime(const vector<int>& route, int vert_id, map<int,long>& beginningServiceTimes, map<int,long>& arrivalTimes);
  void updateWaitingTime     (const vector<int>& route, int vert_id, map<int,long>& waitingTimes, map<int,long>& beginningServiceTimes, map<int,long>& arrivalTimes) ;
  void updateDepartureTime   (const vector<int>& route, int vert_id, map<int,long>& departureTimes, map<int,long>& beginningServiceTimes);

  long computeForwardTimeSlack(int i, const vector<int>& route, map<int,long>& B, map<int,long>& A, map<int,long>& D, map<int,long>& W, bool debug);

  long computeLoadViolation      (map<int,long> y);
  long computeTimeWindowViolation(map<int,long> B);
  long computeRideViolation      (map<int,long> L);

  long computeRideTimeToMaxDiff(int vertex_id, long ridetime);

  long double computePickupDelay  (map<int,long> B);
  long double computeDeliveryDelay(map<int,long> B, map<int, long> A);

  void printVertexInfo (int vertex, map<int,long>& A, map<int,long>& W, map<int,long>& B, map<int,long>& D, map<int,long>& F, map<int,long>& L, map<int,long>& y);
  void printVertexInfo2(int vertex, int nextvert, map<int,long>& A, map<int,long>& W, map<int,long>& B, map<int,long>& D, map<int,long>& F, map<int,long>& L, map<int,long>& y);


public:


  DARPEvaluator(const long vehicleCapacity,
                const VerticesList& vertlist, const CostMatrix& costMatrix,
                int timeSumExp,
                long double loadvW, long double TWVW, long double rideVW,
                darpOptCriterionType optcrit);

  virtual ~DARPEvaluator(){}

  //long double              maxRideValue()     const { return maxRideValue_;}
  long                vehicleCapacity() const { return vehicleCapacity_;}
  const VerticesList& vertList()        const { return vertlist_;}
  const CostMatrix&   costMatrix()      const { return costMatrix_;}

  int requestsNum() const;

  bool evalRoute(const vector<int>& route,
                 long& cost,
                 long& loadViolation, long& TWViolation, long& rideViolation,
                 long double& pickupDelay, long double& deliveryDelay,
                 bool debug = false,
                 bool debug2 = false);

  bool evalRoute(const vector<int>& route,
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
                 bool                    debug = false,
                 bool                    debug2 = false);

  long double         nonPenalizedScore(const DARPGenome& gen           ) const;
  virtual long double nonPenalizedScore(const DARPGenome& gen, int route) const;
  virtual long double nonPenalizedScore(long cost, long double pickupdelay, long double deliverydelay) const;
  virtual long double score(DARPGenome& gen);

  // These methods are virtual in order to be rewritten by a dummy class used in the tests

  // Developed for the computation of the optimimized version of the distance matrix of the kmedoidsreqdist matrix
  virtual long double score(long double nonpenscore, long TWV, long rideV, long loadV) {
    return score(nonpenscore, TWV, rideV, loadV, TWVWeight_, rideVWeight_, loadVWeight_);
  }

  // Developed for the computation of the optimimized version of the distance matrix of the kmedoidsreqdist matrix
  virtual long double score(long double nonpenscore, long TWV, long rideV, long loadV, long double TWVWeight, long double rideVWeight, long double loadVWeight) {
    return nonpenscore + TWVWeight*TWV + rideVWeight*rideV + loadVWeight*loadV;
  }

  long double loadVWeight() const { return loadVWeight_; }
  long double TWVWeight()   const { return TWVWeight_;   }
  long double rideVWeight() const { return rideVWeight_; }

  void loadVWeight(long double value) { loadVWeight_ = value; }
  void TWVWeight  (long double value) { TWVWeight_   = value; }
  void rideVWeight(long double value) { rideVWeight_ = value; }

  void setAllWeights(long double load, long double tw, long double ride) {
    loadVWeight_ = load;
    TWVWeight_   = tw;
    rideVWeight_ = ride;
  }

  int evalRouteCalls() const { return numRouteCalls_; }

  void displayStatsEachNumRouteCalls (int value) { numRouteCallsStatsDispFreq_ = value; }

  void decrementEvalRouteCalls() { numRouteCalls_--; if (numRouteCalls_<0) numRouteCalls_ = 0;}

  // Used for the pseudo duplicated code for quickly computating the distance matrix
  long double scalePickupValue(long begintime, Vertex& vertex) ;
  long double scaleDeliveryDelay(Vertex& pick_vert, Vertex& del_vert, long begin_time, long del_time);

};

#endif

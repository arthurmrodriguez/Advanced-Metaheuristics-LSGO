#ifndef DARPGENOME_H
#define DARPGENOME_H

#include "CostMatrix.h"
#include "VerticesList.h"
#include <genomes/GA1DArrayGenome.h>
#include <genomes/RoutingGenome.h>
#include <vector>
#include <limits>
#include <list>
#include <iostream>

using namespace std;

class DARPGenome : public GA1DArrayGenome< vector<int> >, public RoutingGenome {

public:
  enum lsstatus {LSEXPLORED, LSUNEXPLORED, LSUNFINISHED};

protected:

  vector<bool> modifiedRoutes_;

  vector<long> costPerRoute_;
  vector<long> loadVPerRoute_;
  vector<long> TWVPerRoute_;
  vector<long> rideVPerRoute_;

  vector<long double> pickupDelayPerRoute_;
  vector<long double> deliveryDelayPerRoute_;

  vector<lsstatus>        localSearchedRoutes_;
  vector< pair<int,int> > localSearchLastPosExplored_;

  vector<bool> feasibleRoutes_;

  long double evalWeightsUpdateValue_;    // delta  in Parragh's or Cordeau's paper

  static DARPEvaluator* darpeval__;

public:

  using GA1DArrayGenome< vector<int> >::length;

  DARPGenome (unsigned x, long double evalWeightsUpdateValue=0, GAGenome::Evaluator f = (GAGenome::Evaluator) 0, void* u = (void*) 0);

  DARPGenome (const DARPGenome& orig) : GA1DArrayGenome< vector<int> > (orig),
                                        RoutingGenome(orig),
                                        GAGenome(orig) { copy(orig); }

  DARPGenome& operator= (const GAGenome& arr) {copy(arr); return *this;}

  GAGenome* clone (GAGenome::CloneMethod flag = GAGenome::CONTENTS) const {return new DARPGenome (*this);}

  virtual void copy (const GAGenome& orig);

  virtual int resize (int x);

  int routeLength(int route) const;

  int gene  (int route, int pos) const;
  int gene  (int route, int pos, int value);

  bool modifiedRoute(int route) const {  assert(route >= 0 && route < this->length() ); return this->modifiedRoutes_[route];}

  bool modifiedRoute(int route, bool m);

  void setRouteWasModified(int route);

  long costPerRoute(int route) const;
  long costPerRoute(int route, long value);
  long cost() const;

  long loadVPerRoute(int route) const;
  long loadVPerRoute(int route, long value);
  long loadV() const;

  long TWVPerRoute(int route) const;
  long TWVPerRoute(int route, long value);
  long TWV() const;

  long rideVPerRoute(int route) const;
  long rideVPerRoute(int route, long value);
  long rideV() const;

  long double deliveryDelayPerRoute(int route) const;
  long double deliveryDelayPerRoute(int route, long double value);
  long double deliveryDelay() const;

  long double pickupDelayPerRoute(int route) const;
  long double pickupDelayPerRoute(int route, long double value);
  long double pickupDelay() const;

  lsstatus localSearchedRoute(int route) const;
  lsstatus localSearchedRoute(int route, lsstatus m);

  int  localSearchLastPosExplored                (int route) const;
  int  localSearchLastCritVertexSearchPosExplored(int route) const;
  void localSearchLastPosExplored                (int route, int pos, int crivertsearchpos);

  int randomRoutePos()             const;
  int randomNonEmptyRoutePos()     const;
  int emptyRouteOrRandomRoutePos() const;

  void pushBackVertex          (int route, int seq);
  void insertVertex            (int route, int pos, int value);
  void insertVertices          (int route, int pos, const list<int>& values);

  void getInsertionPos         (int route, int vert_id, int& start_pos, int& end_pos);
  void getInsertionPos         (int route, int vert_id, int& start_pos, int& end_pos, int start_pos_requested);

  int       insertVertexInBestPos             (int route, int vert_id, int start_pos, int end_pos,
                                               long double scoretoimprove=GAGenome::bestPossibleScore(), long maxevals= std::numeric_limits<int>::max() ) ;
  int       insertVertexInBestPos             (int route, int vert_id) ;
  int       insertVertexInBestPosFromStartPos (int route, int vert_id, int start_pos_requested,
                                               long double scoretoimprove=GAGenome::bestPossibleScore(), long maxevals= std::numeric_limits<int>::max() );
  list<int> insertVerticesInBestPos           (int route, const list<int> vertices_ids);

  int       insertVertexInFirstFeasiblePos            (int route, int vert_id, bool userandominsertpos, int start_pos, int end_pos);
  int       insertVertexInFirstFeasiblePos            (int route, int vert_id, bool userandominsertpos=false);
  int       insertVertexInFirstFeasiblePosFromStartPos(int route, int vert_id, int start_pos_requested, bool userandominsertpos=false);

  long double scoreOfInsertingVertex  (int route, int pos, int vert_id) const;
  long double scoreOfInsertingVertices(int route, const list<int> vertices_ids) const;

  long double scoreWithoutIncRouteEvals() { long double sc=score(); darpeval__->decrementEvalRouteCalls(); return sc;}

  void swapSeqs             (int route1, const list<int>& seq1, int route2, const list<int>& seq2);
  void swapWithInsertPosSeqs(int from_route, int from_ins_pos, const list<int>& from_seq,
                             int to_route,   int to_ins_pos,   const list<int>& to_seq);

  int       removeVertex           (int route, int value, int startpos=-1);
  int       removeVertexOfPos      (int route, int pos);
  list<int> removeVerticesFrom     (int route, int pos, int length);
  void      removeVertices         (int route, const list<int>& seq);

  long double    scoreOfRemovingVertices(int route, const list<int>& vertices_ids) const;

  void moveVertices (int fromroute, int toroute, const list<int>& seq);

  int findPosOfVertex (int route, int value) const;

  virtual int  write      (ostream& os) const;
  virtual void writeObject(ostream& os) const;
  virtual void readObject (istream& is);

  virtual long double nonPenalizedScore() const;

  virtual void evaluateAllRoutes();
  virtual void evaluateRoute    (int route);

  virtual int evalRouteCalls() const;

  virtual void updatePenalizations();

  int requestsNum() const {
    // Note: cannot be computed from the genome data since this value is used by the TSDARPGenome at the
    // construction phase when the genome is empty
    assert(darpeval__); return darpeval__->requestsNum();
  }

  int vehicleCapcity() const {
    assert(darpeval__); return darpeval__->vehicleCapacity();
  }

  long double loadVWeight() const { return darpeval__->loadVWeight(); }
  long double TWVWeight()   const { return darpeval__->TWVWeight();   }
  long double rideVWeight() const { return darpeval__->rideVWeight(); }

  int  numRoutes() const;
  bool moreThanOneRoute();
  void addRoutes(const GAGenome& gen);
  void addRoute (const GAGenome& gen, int routepos);
  void emptyRoutes();

  virtual RoutingGenome* cloneWithoutEmptyRoutes() const;

  // Debugging
  void printRoutes()              const; // Just for debugging
  void checkNoVertexIsMissing()   const; // Just for debugging
  void checkNoClonesInSameRoute() const; // Just for debugging
  void checkNoClonesInAllRoutes() const; // Just for debugging

  // static methods 
  //
  // TODO: Comprobar si los punteros ya tienen valor
  static void setDARPEvaluator(DARPEvaluator* patheval) { darpeval__ = patheval;}

};

#ifdef GALIB_USE_STREAMS
namespace __gnu_cxx  {
inline STD_OSTREAM& operator<< (STD_OSTREAM& os, const vector<int>& v) {
  for (unsigned int i = 0; i < v.size(); i++)
    os << v[i] << " ";
   os << endl;

  return (os);
}

inline STD_ISTREAM& operator>> (STD_ISTREAM& is, vector<int>& v) {return (is);}
}
#endif

#endif

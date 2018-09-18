#ifndef DARP_INIT__
#define DARP_INIT__

#include "DistMatrix.h"
#include "VerticesList.h"
#include "DARPGenome.h"
#include "DARPEvaluator.h"
#include <genomes/GAGenome.h>
#include <string>

using namespace std;

typedef void (*darpInitFuncT)(GAGenome& g, DistMatrix& dist, VerticesList& verticeslist);

void DARPVNSInit (GAGenome& g, DistMatrix& dist, VerticesList& verticeslist);

class DARPObjFuncBasedInitC {
  static void constructGenomeObjFuncBased(DARPGenome& gen, DistMatrix& distmatrix, VerticesList& verticeslist);

public:
  enum objfuncbasedmethod { BEST, FIRST };

  static objfuncbasedmethod method__;

  static void setMethod(objfuncbasedmethod method) { method__ = method; }

  static void DARPObjFuncBasedInit (GAGenome& g, DistMatrix& dist, VerticesList& verticeslist);
};

void DARPObjFuncBasedFastInit (GAGenome& g, DistMatrix& dist, VerticesList& verticeslist);

void DARPBoSSortedRandomInit (GAGenome& g, DistMatrix& dist, VerticesList& verticeslist);

void DARPRandomInit          (GAGenome& g, DistMatrix& dist, VerticesList& verticeslist);

class DARPRegretInsertionInitC {
private:
  // Decentralization Index
  static void processRequestsAccordingToDecentrIndex(DistMatrix& dist, std::vector<Request>& requests, double alpha);
  static double decentralizationIndex(DistMatrix& dist, std::vector<Request>& requests, int request_pos);
  static double sumOfDistancesToOtherRequests(DistMatrix& dist, std::vector<Request>& requests, int request_pos);
  static double sumPickAndDelDistances(DistMatrix& dist, Request& req_orig, Request& req_dest);

  // Better Margin
  static void processRequestsAccordingToBetterMargin(int nroutes, DistMatrix& dist, std::vector<Request>& requests);
  static int findFirstFromPosThatDoesNotSatisfyTimeRequirementsWithReqPos(int pos, DistMatrix& dist, std::vector<Request>& requests);
  static bool isTravelingFromLateDeliveryGoingToArriveBeforeEarlyPickup(Request& req, Request& other, DistMatrix& dist);

  // Regret Insertion
  static void   constructRoutesAccordingRegretInsertion(DARPGenome& gen, std::vector<Request>& requests_v, int maxreqsforregret);
  static void   findReqWithHighestRegret(DARPGenome& gen, std::list<Request*>& requests_toins, int maxreqsforregret,
                                         /*out*/ int& ins_route, /*out*/ list<Request*>::iterator& sel_reqpos);
  static double computeReqRegretValue(DARPGenome& gen, Request& req, /*out*/ int& ins_route);
  static void   computeBestAndAllInsertionScores(DARPGenome& gen, Request& req,
                                                /*out*/ double& best_score, /*out*/ int& best_ins_route, /*out*/ vector<double>& ins_scores);

  static double decentr_alpha__;
  static int    maxreqsforregret__;

public:
  static void DARPRegretInsertionInit (GAGenome& g, DistMatrix& dist, VerticesList& verticeslist);

  static void   setInitParameters(double decen, int maxreqs) { decentr_alpha__ = decen; maxreqsforregret__ = maxreqs; }
  static double decentrAlpha()        { return decentr_alpha__; }
  static int    maxReqsForRegretIns() { return maxreqsforregret__; }
};

/*
 * This initialization method is based on the article "A new insertion-based construction heuristic for solving the pickup and
 * delivery problem with time windows" Quan Lu, Maged M. Dessouky
 */
class DARPSlackInitC {
private:
  struct ComputedTimes {
    vector<long> arrivalTimes;
    vector<long> waitingTimes;
    vector<long> maxposttimeints;
  };

  static int maxreqs__;
  static DARPEvaluator* eval__;

  static void constructRoutesAccordingToSlackValue(DARPGenome& gen,std::vector<Request>& requests,
                                                   vector<Request>& init_requests, int maxreqs,
                                                   DistMatrix& dist, VerticesList& verticeslist);

  static vector<Request>        computeMaxClique(vector<Request>& requests, DistMatrix& dist);
  static vector< vector<bool> > createReqInCompatibilityMatrix(std::vector<Request>& requests, DistMatrix& dist);
  static bool                   areReqIncompatible (Request& req_i, Request& req_j, DistMatrix& dist);
  static bool                   areReqIncompatible2(Request& req_i, Request& req_j, DistMatrix& dist);
  static vector<bool>           computeClique(vector< vector<bool> > incompmatrix, int starting_node);
  static void                   markedNodeAndremoveNodesThatAreNotLinked( int pos, vector< vector<bool> >& imcompmatrix,
                                                          vector<bool>& removed_nodes, vector<bool>& marked_nodes);
  static void                   removeNode(int pos, vector< vector<bool> >& incompmatrix, vector<bool>& removed_nodes );
  static int                    findNonMarkedHighestIncidentEdgesNode( vector< vector<bool> >& incompmatrix,
                                                                       vector<bool>& removed_nodes, vector<bool>& marked_nodes);
  static vector<Request>        createReqListFromMarkedNodes(vector<Request>& requests, vector<bool>& marked_nodes);

  static void findReqWithSmallestSlackCost(DARPGenome& gen, std::list<Request*>& requests_toins, int maxreqs,
                                          int num_routes, DistMatrix& dist, VerticesList& verticeslist,
                                          /*out*/ int& ins_route, /*out*/ int& pick_ins_pos, /*out*/ int& del_ins_pos,
                                          /*out*/ bool& feasible, /*out*/ list<Request*>::iterator& sel_reqpos);

  static double computeSlackCostOfInsertingReq(DARPGenome& gen, Request& req, int num_routes,
                                               DistMatrix& dist, VerticesList& verticeslist, /*out*/ bool& feasible,
                                               /*out*/ int& best_route, /*out*/ int& best_pick_ins_pos, /*out*/ int& best_del_ins_pos);

  static double valueOfBestRouteAndPosInsertion(DARPGenome& gen, Vertex& newvert, int num_routes,
                                                DistMatrix& dist, VerticesList& verticelist,
                                               /*out*/ int& best_route, /*out*/ int& best_insert_pos);

  static double valueOfInsertingVertexAtBestPos(DARPGenome& gen, int route, Vertex& newvert, DistMatrix& dist, VerticesList& verticelist,
                                                /*out*/ int& best_insert_pos);
  static ComputedTimes computeTimes(DARPGenome& gen, int route, DistMatrix& dist, VerticesList& verticeslist);
  static void          computeArrivalAndWaitingTimes(DARPGenome& gen, int route, DistMatrix& dist, VerticesList& verticeslist, /*out*/ ComputedTimes& ctimes);
  static void          computeMaxPostTimeInterval(DARPGenome& gen, int route, VerticesList& verticeslist, /*in-out*/ ComputedTimes& ctimes);
  static int           computeArrivalDiff(DARPGenome& gen, int route, int pos, VerticesList& verticeslist, vector<long>& arrivalTimes);

  static long computeSlackCostOfInsertion  (DARPGenome& gen, int route, int pos, Vertex& newvert, DistMatrix& dist, VerticesList& verticeslist, ComputedTimes& ctimes);
  static bool isInsertionUnfeasible        (DARPGenome& gen, int route, int pos, Vertex& newvert, DistMatrix& dist, VerticesList& verticeslist, ComputedTimes& ctimes);
  static long computeNewArrival            (DARPGenome& gen, int route, int pos, Vertex& newvert, DistMatrix& dist, VerticesList& verticeslist, ComputedTimes& ctimes);
  static void computeValuesOfInsertingAtPos(DARPGenome& gen, int route, int pos, Vertex& newvert, DistMatrix& dist, VerticesList& verticeslist, ComputedTimes& ctimes,
                                           /*out*/ long& newarrival, /*out*/ long& newwait, /*out*/ long& newfts  );
  static long computeC1(DARPGenome& gen, int route, int pos, long newwait, long newfts, VerticesList& verticeslist, ComputedTimes& ctimes);
  static long computeC2(Vertex& newvert, long newarrival, long newfts);
  static long computeC3(DARPGenome& gen, int route, int pos, Vertex& newvert, VerticesList& verticeslist, DistMatrix& dist);

public:
  static void DARPSlackInit(GAGenome& g, DistMatrix& dist, VerticesList& verticeslist);

  static void setInitParameters(int maxreqs);

  static int maxReqsPerInitIt() { return maxreqs__; }

};

darpInitFuncT initMethod(string name);

/***************************** SKYBUS SOURCE CODE ***********************************************/
    
void DARPSkybusInit(GAGenome& g, DistMatrix& dist,VerticesList& verticeslist);

#endif

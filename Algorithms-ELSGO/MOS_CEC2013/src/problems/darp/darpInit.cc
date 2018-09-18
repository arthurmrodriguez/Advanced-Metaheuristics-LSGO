#include "darpInit.h"
#include "VerticesList.h"
#include "aux.h"
#include "darpOps.h"
#include <algorithm>
#include <limits>
#include <map>

/*
 * Comparator used with the DARP initializer
 * It is initialized with a list of requests that uses to creates an arbitrary list of possible
 * beggining of services (BoS) based on the values of the corresponding requests windows
 */
struct RequestsBoSComparator {
  bool userandompickup_;
  std::map<int, int> requests_bos_;

  RequestsBoSComparator(std::vector<Request>& requests, bool randompickup=true)  : userandompickup_(randompickup){
    // Creates a vector of possible beginning of services for the requests
    std::vector<Request>::iterator it;
    for (it=requests.begin(); it!=requests.end(); it++) {
      Vertex* pickup_vert = (*it).pickup_vert;
      requests_bos_[pickup_vert->id_] = (userandompickup_) ? GARandomInt(pickup_vert->fbegin_,pickup_vert->fend_)
                                                           : pickup_vert->fbegin_;
    }
  }

  bool operator() (Request reqi,Request reqj) {
    int bos_reqi = requests_bos_[reqi.pickup_vert->id_];
    int bos_reqj = requests_bos_[reqj.pickup_vert->id_];
    return (bos_reqi<bos_reqj);
  }
};

int routeMinCostCriterion(VerticesList& vertices,
                          DistMatrix&   distmatrix,
                          DARPGenome&   gen,
                          Request&      req,
                          double        (*crit)(DistMatrix&, Vertex&, Vertex&, Vertex&, Vertex& ) ) {

  Vertex& req_pick_vert = *req.pickup_vert;
  Vertex& req_del_vert  = *req.delivery_vert;
  assert(req_pick_vert.type_ == Vertex::PICKUP && req_del_vert.type_ == Vertex::DELIVERY);

  int best_route = -1;
  double best_value = numeric_limits<double>::max();
  int num_routes = gen.length();

  for (int route=0; route<num_routes; route++) {
    int route_size = gen.routeLength(route);
    assert( route_size >= 2); // A request must at least be inserted

    int route_pick_vert_id = gen.gene(route,route_size-2);
    int route_del_vert_id  = gen.gene(route,route_size-1);

    Vertex& route_pick_vert = vertices.getVertex(route_pick_vert_id);
    Vertex& route_del_vert  = vertices.getVertex(route_del_vert_id);
    assert(route_pick_vert.type_ == Vertex::PICKUP && route_del_vert.type_ == Vertex::DELIVERY);

    double route_value = crit(distmatrix, req_pick_vert, req_del_vert, route_pick_vert, route_del_vert);

    if (route_value < best_value) {
      best_route = route;
      best_value = route_value;
    }
  }
  if (best_route == -1) throw runtime_error("Error in bestRouteBestPosBased init method, maximum possible score reached");

  return best_route;
}

double OrigOrigCostCriterion(DistMatrix& distmatrix,
                          Vertex&    req_pick_vert,  Vertex& req_del_vert,
                          Vertex&    route_pick_vert, Vertex& route_del_vert) {
  return distmatrix.getDist(req_pick_vert.pos_,route_pick_vert.pos_);
}

double OrigDestCostCriterion(DistMatrix& distmatrix,
                          Vertex&    req_pick_vert,  Vertex& req_del_vert,
                          Vertex&    route_pick_vert, Vertex& route_del_vert) {
  return distmatrix.getDist(req_pick_vert.pos_,route_del_vert.pos_);
}

double DestOrigCostCriterion(DistMatrix& distmatrix,
                          Vertex&    req_pick_vert,  Vertex& req_del_vert,
                          Vertex&    route_pick_vert, Vertex& route_del_vert) {
  return distmatrix.getDist(req_del_vert.pos_,route_pick_vert.pos_);
}

double DestDestCostCriterion(DistMatrix& distmatrix,
                          Vertex&    req_pick_vert,  Vertex& req_del_vert,
                          Vertex&    route_pick_vert, Vertex& route_del_vert) {
  return distmatrix.getDist(req_del_vert.pos_,route_del_vert.pos_);
}

void pushBackRequest(DARPGenome& gen, int route, Request& req) {
  assert(req.pickup_vert);
  assert(req.delivery_vert);
  gen.pushBackVertex(route,req.pickup_vert->id_);
  gen.pushBackVertex(route,req.delivery_vert->id_);
}

vector< double (*)(DistMatrix&, Vertex&, Vertex&, Vertex&, Vertex& ) > constructCriteriaList () {
  vector< double (*)(DistMatrix&, Vertex&, Vertex&, Vertex&, Vertex& ) > criteria;
  criteria.push_back(OrigOrigCostCriterion);
  criteria.push_back(OrigDestCostCriterion);
  criteria.push_back(DestOrigCostCriterion);
  criteria.push_back(DestDestCostCriterion);

  return criteria;
}

/*
 * Returns the requests from the verticeslist sorted by the BoS criterion. They are returned by value (not very efficient)
 * but since this method is only call once at the initialization, it is clearer to use this approach
 */
std::vector<Request> getBoSSortedRequests(VerticesList& verticeslist, bool userandompickup=true) {
  std::vector<Request> requests = verticeslist.getReqsList();
  RequestsBoSComparator comp(requests,userandompickup);
  sort(requests.begin(), requests.end(), comp);

  return requests;
}


void initializeRoutesUsingFirstNRequests(DARPGenome& gen, std::vector<Request>& requests, int size) {
  assert(size <= gen.length());

  // Then all routes are initialized with one request using the first num_routes requests of the list
  for (int requests_pos=0; requests_pos<size && requests_pos<requests.size(); requests_pos++) {
    pushBackRequest(gen,requests_pos,requests[requests_pos]);
  }

}

void VNSInitConstructGenome(DARPGenome& gen, DistMatrix& distmatrix, VerticesList& verticeslist) {
  // First the requests are sorted according to an artificial beginning of service
  std::vector<Request> requests = getBoSSortedRequests(verticeslist);

  initializeRoutesUsingFirstNRequests(gen,requests, gen.length());

  vector< double (*)(DistMatrix&, Vertex&, Vertex&, Vertex&, Vertex& ) > criteria = constructCriteriaList();

  // After that, the requests are inserted at the end of each partial route in the order they appear in the list
  // The insertion route is selected according to a randomly selected criterion
  for (int requests_pos=gen.length(); requests_pos<requests.size(); requests_pos++) {

    Request& req = requests[requests_pos];

    int criterion_pos  = GARandomInt(0,criteria.size()-1);
    int selected_route = routeMinCostCriterion(verticeslist,distmatrix,gen, req, criteria[criterion_pos]);

    pushBackRequest(gen,selected_route,req);
  }
}

void insertRequestInRoute(DARPGenome& gen, Request& req, int route) {
  gen.insertVertexInBestPos(route,req.pickup_vert->id_);
  gen.insertVertexInBestPos(route,req.delivery_vert->id_);
}

void insertRequestInRoute(DARPGenome& gen, Request& req, int route, int pick_ins_pos, int del_ins_pos ) {
  gen.insertVertex(route,pick_ins_pos,req.pickup_vert->id_);
  gen.insertVertex(route,del_ins_pos, req.delivery_vert->id_);
}

double scoreForInsertingRequestInRoute(DARPGenome& gen, Request& req, int route, DARPObjFuncBasedInitC::objfuncbasedmethod method,
                                       /*out*/ int& pick_ins, /*out*/ int& del_ins) {
  assert(route >=0 && route < gen.length());

  DARPGenome tmpgen(gen);

  if (method == DARPObjFuncBasedInitC::BEST) {  // If more methods are added in a future this should be refactored to use an array of pointers
    pick_ins = tmpgen.insertVertexInBestPos(route,req.pickup_vert->id_);
    del_ins  = tmpgen.insertVertexInBestPosFromStartPos(route,req.delivery_vert->id_,pick_ins+1);
  }
  else{
    pick_ins = tmpgen.insertVertexInFirstFeasiblePos(route,req.pickup_vert->id_);
    del_ins  = tmpgen.insertVertexInFirstFeasiblePosFromStartPos(route,req.delivery_vert->id_,pick_ins+1);
  }

  gen.nevals(gen.nevals()+tmpgen.nevals());

  return tmpgen.score();
}

void bestScoreAndRouteForInsertion(DARPGenome& gen, Request& req, DARPObjFuncBasedInitC::objfuncbasedmethod method,
                                   /*out*/ double &best_score, /*out*/int &best_route, /*out*/ int& pickins, /*out*/ int& delins) {
  best_route  = -1;
  best_score  = GAGenome::worstPossibleScore();
  pickins     = -1;
  delins      = -1;

  for (int route=0; route<gen.length(); route++) {
    int tmp_pickins, tmp_delins;
    double tmp_score = scoreForInsertingRequestInRoute(gen,req,route,method,tmp_pickins,tmp_delins);

    if (GAGenome::compareScores(tmp_score, best_score) == GAGenome::BETTER) {
      best_route = route;
      best_score = tmp_score;
      pickins    = tmp_pickins;
      delins     = tmp_delins;
    }

  }

  assert(pickins != -1);
  assert(delins != -1);
  assert(pickins >= 0 && pickins <= gen.routeLength(best_route));
  assert(pickins >= 0 && pickins <= gen.routeLength(best_route)+1);
  assert(best_route >= 0);
}

/*
 * find the best score of removing the vertices from a request in a route with and substitute them with the ones of
 * the passed Request plus the following insertion of the extracted vertices in the best possible route. It return
 * both the best route position  to conduct the swap and the best route to insert the extracted vertices
 * TODO: Needs to be refactored
 */
//void bestScoreOfSwappingWithRequestAndSwapValues(const DARPGenome& gen, const Request& req, VerticesList& verticeslist,
//                                                /*out*/ double &best_score, /*out*/ int &route_swap, /*out*/ int &route_ins ) {
//  route_swap = route_ins = -1;
//  best_score = GAGenome::worstPossibleScore();
//
//  double tmp_score;
//  int    tmp_route_ins;
//
//  for (int route=0; route<gen.length(); route++) {
//    if ( gen.routeLength(route) == 2 ) {
//      // Carry out the swap
//      DARPGenome tmp_gen(gen);
//
//      int tmp_pickup_vert_id = tmp_gen.removeVertexOfPos(route,0);
//      int tmp_del_vert_id    = tmp_gen.removeVertexOfPos(route,0);
//
//      tmp_gen.pushBackVertex(route,req.pickup_vert->id_);
//      tmp_gen.pushBackVertex(route,req.delivery_vert->id_);
//
//      Vertex &tmp_pickup_vert = verticeslist.getVertex(tmp_pickup_vert_id);
//      Vertex &tmp_del_vert    = verticeslist.getVertex(tmp_del_vert_id);
//
//      Request tmp_req(&tmp_pickup_vert,&tmp_del_vert);
//
//      bestScoreAndRouteForInsertion(tmp_gen, tmp_req, tmp_score, tmp_route_ins);
//
//      if (GAGenome::compareScores(tmp_score, best_score) == GAGenome::BETTER) {
//        best_score = tmp_score;
//        route_swap = route;
//        route_ins  = tmp_route_ins;
//      }
//    }
//  }
//  assert(best_score == GAGenome::worstPossibleScore() or route_ins != -1);
//}

/*
 * First it inserts the vertices of the passed request in the best place inside the genome. Then it checks if performing a
 * swap with a route with only one request and inserting the extracted request in the best possible way inside the genome
 * would not get a better score. If this happens, this swap is conducted. Otherwise, the request remains in the first
 * insertion point.
 */
void insertInBestRouteObjFuncBased(DARPGenome& gen, Request& req, VerticesList& verticeslist, DARPObjFuncBasedInitC::objfuncbasedmethod method) {
  gen.evaluate(); // We first make sure that the gen has been evaluated so that the calls to the route evaluation are not wasted

  double best_score;
  int    best_route, pickinspos, delinspos;

  bestScoreAndRouteForInsertion(gen,req,method,best_score,best_route, pickinspos, delinspos);

  insertRequestInRoute(gen,req,best_route,pickinspos,delinspos);


  /* has been disabled temporarely. It is extremely slow with high dimensions
  double swap_best_score;
  int    swap_route;
  int    swap_route_ins;
  bestScoreOfSwappingWithRequestAndSwapValues(gen, req, verticeslist, swap_best_score, swap_route, swap_route_ins );

  if (GAGenome::compareScores(swap_best_score, best_score) == GAGenome::BETTER) {
    int tmp_pickup_vert_id = gen.removeVertexOfPos(swap_route,0);
    int tmp_del_vert_id    = gen.removeVertexOfPos(swap_route,0);

    gen.pushBackVertex(swap_route,req.pickup_vert->id_);
    gen.pushBackVertex(swap_route,req.delivery_vert->id_);

    gen.insertVertexInBestPos(swap_route_ins,tmp_pickup_vert_id);
    gen.insertVertexInBestPos(swap_route_ins,tmp_del_vert_id);
  }
  else {
    gen.insertVertexInBestPos(best_route,req.pickup_vert->id_);
    gen.insertVertexInBestPos(best_route,req.delivery_vert->id_);
  }
  */

}

void DARPObjFuncBasedInitC::constructGenomeObjFuncBased(DARPGenome& gen, DistMatrix& distmatrix, VerticesList& verticeslist) {
  // First the requests are sorted according to an artificial beginning of service
  std::vector<Request> requests = getBoSSortedRequests(verticeslist);

  initializeRoutesUsingFirstNRequests(gen,requests, gen.length());

  // After that, the requests are inserted at the end of each partial route in the order they appear in the list
  // The insertion route is selected according to the best value obtained by the objective function
  for (int requests_pos=gen.length(); requests_pos<requests.size(); requests_pos++) {
    insertInBestRouteObjFuncBased(gen,requests[requests_pos],verticeslist,method__);
  }

}


DARPObjFuncBasedInitC::objfuncbasedmethod DARPObjFuncBasedInitC::method__ = DARPObjFuncBasedInitC::BEST;

void DARPObjFuncBasedInitC::DARPObjFuncBasedInit(GAGenome& g, DistMatrix& dist, VerticesList& verticeslist) {
  DARPGenome& gen = dynamic_cast<DARPGenome&>(g);
  assert(&gen);

  constructGenomeObjFuncBased(gen,dist,verticeslist);
}



void insertReqInFirstFeasiblePos(DARPGenome& gen, Request& req, int route) {
  int pick_ins = gen.insertVertexInFirstFeasiblePos(route,req.pickup_vert->id_);
  gen.insertVertexInFirstFeasiblePosFromStartPos(route,req.delivery_vert->id_,pick_ins+1);
}

void insertRequestInFirstFeasibleRoute(DARPGenome& gen, Request& req) {
  gen.evaluate(); // We first make sure that the gen has been evaluated so that the calls to the route evaluation are not wasted

  DARPGenome* bestgen = dynamic_cast<DARPGenome*>(gen.clone());
  insertReqInFirstFeasiblePos(*bestgen,req,0);

  if (!bestgen->feasible()) {

    // We start at pos 1 because we have already tried the first pos;
    vector<int> routes_sorted;
    for (int route=1; route<gen.length(); route++) routes_sorted.push_back(route);
    random_shuffle(routes_sorted.begin(),routes_sorted.end());

    for (int i=0; i<routes_sorted.size(); i++) {
      int route = routes_sorted[i];

      DARPGenome* tmpgen = dynamic_cast<DARPGenome*>(gen.clone());
      insertReqInFirstFeasiblePos(*tmpgen,req,route);

      if (tmpgen->feasible() ||
          GAGenome::compareScores( tmpgen->score(), bestgen->score() ) == GAGenome::BETTER ) {
        std::swap(bestgen,tmpgen);

        if (bestgen->feasible()) {
          delete tmpgen;
          break;
        }
      }

      delete tmpgen;
    }

  }
  gen.copy(*bestgen);

  delete bestgen;
}

void constructGenomeObjFuncBasedFastInit(DARPGenome& gen, DistMatrix& distmatrix, VerticesList& verticeslist) {
  // First the requests are sorted according to an artificial beginning of service
  std::vector<Request> requests = getBoSSortedRequests(verticeslist);

  initializeRoutesUsingFirstNRequests(gen,requests, gen.length());

  // After that, the requests are inserted at the end of each partial route in the order they appear in the list
  // The insertion route is selected according to the best value obtained by the objective function
  for (int requests_pos=gen.length(); requests_pos<requests.size(); requests_pos++) {
    insertRequestInFirstFeasibleRoute(gen,requests[requests_pos]);
  }
}

void DARPObjFuncBasedFastInit (GAGenome& g, DistMatrix& dist, VerticesList& verticeslist) {
  DARPGenome& gen = dynamic_cast<DARPGenome&>(g);
  assert(&gen);

  constructGenomeObjFuncBasedFastInit(gen,dist,verticeslist);
}

// DARP Init
void DARPVNSInit(GAGenome& g, DistMatrix& dist, VerticesList& verticeslist) {
  DARPGenome& gen = dynamic_cast<DARPGenome&>(g);
  assert(&gen);

  VNSInitConstructGenome(gen,dist,verticeslist);
}

void randomlyConstructGenome(DARPGenome& gen, std::vector<Request>& requests) {
  int num_routes = gen.length();

  // After that, the requests are inserted at the end of each partial route in the order they appear in the list
  // The insertion route is selected according to a randomly selected criterion
  for (int requests_pos=0; requests_pos<requests.size(); requests_pos++) {

    Request& req = requests[requests_pos];

    int selected_route = GARandomInt(0,num_routes-1);

    pushBackRequest(gen,selected_route,req);
  }
}

void DARPBoSSortedRandomInit(GAGenome& g, DistMatrix& dist, VerticesList& verticeslist) {
  DARPGenome& gen = dynamic_cast<DARPGenome&>(g);
  assert(&gen);

  // First the requests are sorted according to an artificial beginning of service
  std::vector<Request> requests = getBoSSortedRequests(verticeslist);

  randomlyConstructGenome(gen, requests);
}

void DARPRandomInit(GAGenome& g, DistMatrix& dist, VerticesList& verticeslist) {
  DARPGenome& gen = dynamic_cast<DARPGenome&>(g);
  assert(&gen);

  std::vector<Request> requests = verticeslist.getReqsList();

  randomlyConstructGenome(gen, requests);
}

double DARPRegretInsertionInitC::sumPickAndDelDistances(DistMatrix& dist, Request& req_orig, Request& req_dest) {
  Vertex* orig_pickup   = req_orig.pickup_vert;
  Vertex* orig_delivery = req_orig.delivery_vert;
  Vertex* dest_pickup   = req_dest.pickup_vert;
  Vertex* dest_delivery = req_dest.delivery_vert;

  return dist.getDist(orig_pickup->pos_,  dest_pickup->pos_) + dist.getDist(orig_pickup->pos_,  dest_delivery->pos_) +
         dist.getDist(orig_delivery->pos_,dest_pickup->pos_) + dist.getDist(orig_delivery->pos_,dest_delivery->pos_);
}

double DARPRegretInsertionInitC::sumOfDistancesToOtherRequests(DistMatrix& dist, std::vector<Request>& requests, int request_pos) {
  assert(request_pos >=0 && request_pos < requests.size());

  double value = 0;

  for (int i=0; i<requests.size(); i++) {
    if (i==request_pos) continue;

    value += sumPickAndDelDistances(dist,requests[request_pos],requests[i]);
  }

  return value;
}

double DARPRegretInsertionInitC::decentralizationIndex(DistMatrix& dist, std::vector<Request>& requests, int request_pos) {
  double numerator = sumOfDistancesToOtherRequests(dist,requests,request_pos);

  double denominator = numerator;
  for (int i=0; i<requests.size(); i++) {
    if (i==request_pos) continue;

    denominator += sumOfDistancesToOtherRequests(dist,requests,i);
  }

  return numerator/denominator;
}

/*
 * Not optimized. Could be optimized to avoid computing each time the decentralizationIndex if the same
 * request is being swap in each iteration
 */
void DARPRegretInsertionInitC::processRequestsAccordingToDecentrIndex(DistMatrix& dist, std::vector<Request>& requests, double alpha) {
  assert(alpha >= 0 && alpha <= 1); // by increasing alpha we put a greater emphasis on the spatial aspect of the problem

  for (int i=requests.size()-1;i>0; i--) {

    if ( decentralizationIndex(dist,requests,i) - decentralizationIndex(dist,requests,i-1) >=  1.0 - alpha ) {
      swap( requests[i-1], requests[i] );
    }
  }
}

bool DARPRegretInsertionInitC::isTravelingFromLateDeliveryGoingToArriveBeforeEarlyPickup(Request& req, Request& other, DistMatrix& dist) {
  double distance = dist.getDist(req.delivery_vert->pos_,  other.pickup_vert->pos_);

  return req.delivery_vert->fend_ + distance <= other.pickup_vert->fbegin_;
}

int DARPRegretInsertionInitC::findFirstFromPosThatDoesNotSatisfyTimeRequirementsWithReqPos(int pos, DistMatrix& dist, std::vector<Request>& requests) {
  assert(pos >=0 && pos < requests.size());

  int res = -1;

  Request& req = requests[pos];

  for (int i=pos+1; i<requests.size(); i++) {
    Request& req_i = requests[i];

    if ( ! isTravelingFromLateDeliveryGoingToArriveBeforeEarlyPickup(req,req_i,dist) ) {
      res = i;
      break;
    }
  }

  return res;
}

/*
 * We try to sort the routes so that the first ones (the ones that will be used as seed for the routes) have
 * incompatibilities for being in the same route (such as not satisfying the aforementioned constraint)
 */
void DARPRegretInsertionInitC::processRequestsAccordingToBetterMargin(int nroutes, DistMatrix& dist, std::vector<Request>& requests) {
  assert(nroutes >= 0 && nroutes <= requests.size());

  for (int req_pos=0; req_pos<nroutes-1; req_pos++) {
    int pos = findFirstFromPosThatDoesNotSatisfyTimeRequirementsWithReqPos(req_pos,dist,requests);

    if (pos != -1 && pos != req_pos+1 ) {         // Only when we find a position (which is not the next one) that satisfies the
      swap( requests[req_pos+1], requests[pos] ); // constraint we conduct the swap
    }
  }
}

void DARPRegretInsertionInitC::computeBestAndAllInsertionScores(DARPGenome& gen, Request& req,
                                     /*out*/ double& best_score, /*out*/ int& best_ins_route, /*out*/ vector<double>& ins_scores) {
  ins_scores.clear();
  ins_scores.resize(gen.length());

  best_ins_route = -1;
  best_score     = GAGenome::worstPossibleScore();

  int pickins,delins; // Not being used right now, in a future they should be returned

  for (int route=0; route<gen.length(); route++) {
    double ins_score  = scoreForInsertingRequestInRoute(gen, req, route, DARPObjFuncBasedInitC::BEST, pickins, delins);
    ins_scores[route] = ins_score;

    if (GAGenome::compareScores(ins_score,best_score) == GAGenome::BETTER) {
      best_ins_route = route;
      best_score     = ins_score;
    }
  }
}

double DARPRegretInsertionInitC::computeReqRegretValue(DARPGenome& gen, Request& req, /*out*/ int& ins_route) {
  vector<double> ins_scores;
  double         best_score;

  computeBestAndAllInsertionScores(gen, req, best_score, ins_route, ins_scores);

  double regret_value = 0;
  for (int i=0; i<ins_scores.size(); i++) regret_value += fabs( ins_scores[i] - best_score );
  assert(regret_value >= 0);

  return regret_value;
}

void DARPRegretInsertionInitC::findReqWithHighestRegret(DARPGenome& gen, std::list<Request*>& requests_toins, int maxreqsforregret,
                             /*out*/ int& ins_route, /*out*/ list<Request*>::iterator& sel_reqpos) {
  double worst_regret_value = -1;

  int i=0;
  for (list<Request*>::iterator it=requests_toins.begin(); it!=requests_toins.end() && i<maxreqsforregret; it++, i++) {
    Request& req_tmp      = **it;
    int      ins_route_tmp;
    double   regret_value = computeReqRegretValue(gen,req_tmp,ins_route_tmp);
    assert(regret_value >= 0);

    if ( regret_value > worst_regret_value ) {
      worst_regret_value  = regret_value;
      sel_reqpos          = it;
      ins_route           = ins_route_tmp;
    }

  }

  assert(worst_regret_value != -1);
  assert(ins_route != -1);
}

void DARPRegretInsertionInitC::constructRoutesAccordingRegretInsertion(DARPGenome& gen, std::vector<Request>& requests_v, int maxreqsforregret) {
  assert(requests_v.size() > gen.length());

  initializeRoutesUsingFirstNRequests(gen,requests_v, gen.length());

  std::list<Request *> requests_toins;
  // Since the gen already contains the first gen.length requests, we insert the remaining requests
  for (int i=gen.length(); i<requests_v.size(); i++) requests_toins.push_back( & requests_v[i] );
  assert(requests_toins.size() > 0);


  while( ! requests_toins.empty() ) {
    int                      ins_route = -1;
    list<Request*>::iterator sel_reqpos;

    gen.evaluate(); // Before we find the request with the highest regret we make sure that the score has been evaluated
                    // otherwise, we could be wasting a lot of evaluations (think in what happens the first time, when the
                    // genome has been initializated with the first m requests but has not been evaluated

    findReqWithHighestRegret(gen, requests_toins, maxreqsforregret,ins_route, sel_reqpos);
    assert(ins_route != -1);

    insertRequestInRoute(gen, **sel_reqpos , ins_route);

    requests_toins.erase(sel_reqpos);
  }

}

double DARPRegretInsertionInitC::decentr_alpha__ = 0.5;
int DARPRegretInsertionInitC::maxreqsforregret__ = 40;

/*
 * Follows the method described in the article  A new regret insertion heuristic for solving large-scale dial-a-ride problems
 * with time windows M. Diana
 */
void DARPRegretInsertionInitC::DARPRegretInsertionInit (GAGenome& g, DistMatrix& dist, VerticesList& verticeslist) {
  DARPGenome& gen = dynamic_cast<DARPGenome&>(g);
  assert(&gen);

  // First the requests are sorted according to the earliest pickup time (false makes the method to use the actual pickup time)
  std::vector<Request> requests = getBoSSortedRequests(verticeslist,false);

  processRequestsAccordingToBetterMargin(gen.length(),dist,requests);

  processRequestsAccordingToDecentrIndex(dist, requests, decentr_alpha__);

  constructRoutesAccordingRegretInsertion(gen,requests, maxreqsforregret__); // TODO: Needs to be parametrized
}



int DARPSlackInitC::maxreqs__ = 0;

void DARPSlackInitC::setInitParameters(int maxreqs) {
  maxreqs__ = maxreqs;
}

void DARPSlackInitC::DARPSlackInit (GAGenome& g, DistMatrix& dist, VerticesList& verticeslist) {
  DARPGenome& gen = dynamic_cast<DARPGenome&>(g);
  assert(&gen);

  // First the requests are sorted according to the earliest pickup time (false makes the method to use the actual pickup time)
  std::vector<Request> requests = getBoSSortedRequests(verticeslist,false);

  // Then, the max clique is computed for being used as the initial set of requests
  // If no clique is found, we always going to have a single request
  vector<Request> init_requests = computeMaxClique(requests,dist);

  // Finally, the remaining set of requests are constructed according to the slack measure
  constructRoutesAccordingToSlackValue(gen,requests,init_requests, maxreqs__, dist, verticeslist);
}

/*
 * TODO: Needs to be refactored with the very similar method from the Regret initialization
 */
void DARPSlackInitC::constructRoutesAccordingToSlackValue(DARPGenome& gen,std::vector<Request>& requests,
                                                          std::vector<Request>& init_requests, int maxreqs,
                                                          DistMatrix& dist, VerticesList& verticeslist) {

  assert(requests.size() > gen.length());
  assert(init_requests.size() > 0);
  assert(maxreqs__ > 0);

  int routes_used = init_requests.size();
  if (routes_used > gen.length() ) {                           // If we have more init requests than routes, we sort them
    routes_used = gen.length();                                // randomly so that at each execution the prune is different
    random_shuffle(init_requests.begin(),init_requests.end());
  }

  initializeRoutesUsingFirstNRequests(gen,init_requests, routes_used);

  std::list<Request *> requests_toins;

  for (int i=0; i<requests.size(); i++) {
    if ( find(init_requests.begin(), init_requests.end(), requests[i]) == init_requests.end() ) {
      requests_toins.push_back( & requests[i] );  // we dont insert requests that were inserted at the initialization
    }
  }
  assert(requests_toins.size() > 0);

  int ins_route = -1, pick_ins_pos = -1, del_ins_pos = -1;
  bool feasible;
  int num_req_toins = requests_toins.size();
  int oldinscalls   = 0;

  while( ! requests_toins.empty() ) {
    int                      ins_route = -1;
    list<Request*>::iterator sel_reqpos;

    gen.evaluate(); // hack: to count the number of evaluations that are conducted with the other method. It should be
                    // done properly in the findReqWithSmallestSlackCost but for now we are leaving this hack


    findReqWithSmallestSlackCost(gen, requests_toins, maxreqs, routes_used, dist, verticeslist,
                                ins_route, pick_ins_pos, del_ins_pos, feasible, sel_reqpos);

    if (!feasible) {                                         // If no feasible insertion location is found, we add another route
      sel_reqpos = requests_toins.begin();                   // and insert the next request in the next empty route
      if (routes_used < gen.length()) {
        insertRequestInRoute(gen, **sel_reqpos, routes_used);
        routes_used++;
      }
      else {                                                // Otherwise we insert it in the best possible route
        oldinscalls++;
        insertInBestRouteObjFuncBased(gen,**sel_reqpos,verticeslist, DARPObjFuncBasedInitC::BEST);
      }
    }
    else {                                                  // feasible route case
      insertRequestInRoute(gen, **sel_reqpos, ins_route, pick_ins_pos, del_ins_pos );
    }
    requests_toins.erase(sel_reqpos);
  }
  cout << "Objective based insertion method calls = " << oldinscalls << " (" << (double)oldinscalls/num_req_toins << " %)" << endl;
  cout << gen << endl;
}

void displayIncomMatrix(vector< vector<bool> >& incomp_matrix) {
  cout << "Incompatibility Matrix: " << endl;
  for (int i=0; i< incomp_matrix.size(); i++) {
    int row_shown = false;

    for (int j=0; j< incomp_matrix.size(); j++) {
      if (i!=j and incomp_matrix[i][j]) {
        if (!row_shown) { cout << i << " : "; row_shown = true; }
          cout << j << " ";
        }
      }
    if (row_shown) cout << endl;
  }
  cout << endl;
}

vector<Request> DARPSlackInitC::computeMaxClique(std::vector<Request>& requests, DistMatrix& dist) {
  vector< vector<bool> > incomp_matrix = createReqInCompatibilityMatrix(requests, dist);
#ifdef DEBUG
  displayIncomMatrix(incomp_matrix);
#endif

  vector<Request> res;

  for (int i=0; i<requests.size(); i++) {
    vector<bool>    marked_nodes = computeClique(incomp_matrix,i);
    vector<Request> tmp          = createReqListFromMarkedNodes(requests,marked_nodes);
    if (tmp.size() > res.size() ) res = tmp;
  }

  return res;
}

vector< vector<bool> > DARPSlackInitC::createReqInCompatibilityMatrix(std::vector<Request>& requests, DistMatrix& dist) {
  vector< vector<bool> > incompmatrix (requests.size());

  for (int i=0; i<requests.size(); i++) incompmatrix[i].resize(requests.size(),false);

  for (int i=0; i<requests.size()-1; i++) {
    for (int j=i+1; j<requests.size(); j++) {
      if (areReqIncompatible(requests[i], requests[j], dist)) {
        assert(areReqIncompatible(requests[j],requests[i],dist));

        incompmatrix[i][j] = incompmatrix[j][i] = true;
      }
    }
  }

  return incompmatrix;
}

bool DARPSlackInitC::areReqIncompatible(Request& req_i, Request& req_j, DistMatrix& dist) {
  // only six cases need to be analyzed to determine the compatibility
  return areReqIncompatible2(req_i,req_j,dist) and areReqIncompatible2(req_j,req_i,dist);
}

inline bool DARPSlackInitC::areReqIncompatible2(Request& req_i, Request& req_j, DistMatrix& dist) {
  return ! (
/* vi vi+n vj vj+n */   req_i.delivery_vert->canVertexGoBefore(*req_j.pickup_vert,dist)        /*i+n  - j   */
                      or
/* vi vj vi+n vj+n */ ( req_i.pickup_vert->canVertexGoBefore(*req_j.pickup_vert,dist)     and  /* i   - j   */
                        req_j.pickup_vert->canVertexGoBefore(*req_i.delivery_vert,dist)   and  /* j   - i+n */
                        req_i.delivery_vert->canVertexGoBefore(*req_j.delivery_vert,dist)      /* i+n - j+n */
                      )
                      or
/* vi vj vj+n vi+n */ ( req_i.pickup_vert->canVertexGoBefore(*req_j.pickup_vert,dist)     and  /* i   - j   */
                        req_j.delivery_vert->canVertexGoBefore(*req_i.delivery_vert,dist)      /* j+n - i+n */
                      )
           );
}

/*
 * Computes the max clique according to what is defined in the paper: a new insertion-based construction heuristic for solving
 * the pickup and delivery problem with time windows - Quan Lu, Maged M. Dessouky
 * Note that that matrix is not passed by reference since it is being modified
 */
vector<bool> DARPSlackInitC::computeClique(vector< vector<bool> > incompmatrix, int starting_node) {
  vector<bool> removed_nodes(incompmatrix.size(),false);
  vector<bool> marked_nodes(incompmatrix.size(),false);

  int pos = starting_node;

  do {
    markedNodeAndremoveNodesThatAreNotLinked(pos,incompmatrix,removed_nodes,marked_nodes);

    pos = findNonMarkedHighestIncidentEdgesNode(incompmatrix,removed_nodes, marked_nodes);
  } while(pos != -1);

  return marked_nodes;
}

void DARPSlackInitC::markedNodeAndremoveNodesThatAreNotLinked( int pos, vector< vector<bool> >& incompmatrix,
                                                               vector<bool>& removed_nodes, vector<bool>& marked_nodes) {
  marked_nodes[pos] = true;
  for (int i=0; i<incompmatrix.size(); i++) {
    // Remove all the nodes where there is no edge linking node i and these nodes directly
    if (i!=pos and !removed_nodes[i] && !marked_nodes[i] && !incompmatrix[pos][i]) {
      assert(marked_nodes[i] == false);
      removeNode(i,incompmatrix,removed_nodes);
    }
  }

}

void DARPSlackInitC::removeNode(int pos, vector< vector<bool> >& incompmatrix, vector<bool>& removed_nodes ) {
  removed_nodes[pos] = true;
  for (int i=0; i<incompmatrix.size(); i++) incompmatrix[i][pos] = incompmatrix[pos][i] = false;
}

int DARPSlackInitC::findNonMarkedHighestIncidentEdgesNode( vector< vector<bool> >& incompmatrix,
                                                           vector<bool>& removed_nodes, vector<bool>& marked_nodes) {
  int max_inc = 0;
  int pos     = -1;

  for (int i=0; i<incompmatrix.size(); i++) {

    if (!removed_nodes[i] && !marked_nodes[i]) {

      int inc_num = 0;
      for (int j=0; j<incompmatrix.size(); j++) if ( incompmatrix[i][j] ) inc_num++;

      if ( inc_num > max_inc or
          (inc_num > 0 and inc_num == max_inc and GARandomDouble() <= 0.5) ) { // In the max clique algorithm of the paper, mentions
        max_inc = inc_num;                                                     // that in case of equal values, a position is
        pos     = i;                                                           // selected arbitrarily
      }

    }

  }

  return pos;
}

vector<Request> DARPSlackInitC::createReqListFromMarkedNodes(vector<Request>& requests, vector<bool>& marked_nodes) {
  vector<Request> res;
  for (int i=0; i<marked_nodes.size(); i++) if (marked_nodes[i]) res.push_back(requests[i]);

  return res;
}

void DARPSlackInitC::findReqWithSmallestSlackCost(DARPGenome& gen, std::list<Request*>& requests_toins, int maxreqs,
                                                 int num_routes, DistMatrix& dist, VerticesList& verticeslist,
                                                 /*out*/ int& ins_route, /*out*/ int& pick_ins_pos, /*out*/ int& del_ins_pos,
                                                 /*out*/ bool& feasible, /*out*/ list<Request*>::iterator& sel_reqpos) {
  double best_value   = numeric_limits<double>::max();
  bool   feasible_tmp = false;
  int    route_tmp    = -1, pick_ins_pos_tmp = -1, del_ins_pos_tmp = -1;

  feasible = false;
  ins_route = pick_ins_pos = del_ins_pos = -1;

  int i=0;
  for (list<Request*>::iterator it=requests_toins.begin(); it!=requests_toins.end() && i<maxreqs; it++, i++) {
    Request& req_tmp = **it;

    double value = computeSlackCostOfInsertingReq(gen, req_tmp, num_routes, dist, verticeslist, feasible_tmp,
                                                  route_tmp, pick_ins_pos_tmp, del_ins_pos_tmp);
    assert(value >= -1);

    if ( feasible_tmp and value < best_value ) {
      best_value   = value;
      ins_route    = route_tmp;
      feasible     = feasible_tmp;
      pick_ins_pos = pick_ins_pos_tmp;
      del_ins_pos  = del_ins_pos_tmp;
      sel_reqpos   = it;
    }

  }

  assert(feasible == false or ( ins_route != -1 && pick_ins_pos != -1 && del_ins_pos != -1) );
}

/*
 * Returns -1 if no feasible solution is found
 */
double DARPSlackInitC::computeSlackCostOfInsertingReq(DARPGenome& gen, Request& req, int num_routes,
                                                      DistMatrix& dist, VerticesList& verticeslist, /*out*/ bool& feasible,
                                                      /*out*/ int& best_route, /*out*/ int& best_pick_ins_pos, /*out*/ int& best_del_ins_pos) {

  // We assume that the insertion in feasible and if not we update the state
  feasible = true;

  // First, we find the best place (according to the slack measure) to inserting the pickup vertex
  best_route = -1; best_pick_ins_pos = -1;
  double pick_ins_value =
      valueOfBestRouteAndPosInsertion(gen, * req.pickup_vert, num_routes, dist, verticeslist, best_route, best_pick_ins_pos);
  assert(pick_ins_value == -1 or ( best_route != -1 and best_pick_ins_pos != -1) );

  // If unfeasible we return with the unfeasible value (-1)
  if (pick_ins_value == -1) {feasible = false; return -1;}

  // Then we insert the vertex in a copy of the genome in the best insertion found
  DARPGenome tmp_gen(gen);
  tmp_gen.insertVertex(best_route,best_pick_ins_pos,req.pickup_vert->id_);

  // Then we conduct a similar process for the delivery vertex
  best_del_ins_pos = -1;
  double pick_del_value =
      valueOfInsertingVertexAtBestPos(tmp_gen,best_route,*req.delivery_vert,dist,verticeslist,best_del_ins_pos);
  assert(pick_del_value == -1 or best_del_ins_pos != -1);

  // If unfeasible we return with the unfeasible value (-1)
  if (pick_del_value == -1) {feasible = false; return -1; }

  tmp_gen.insertVertex(best_route,best_del_ins_pos,req.delivery_vert->id_);

  feasible = tmp_gen.feasible();
  assert(feasible || tmp_gen.TWV() == 0); // if a TWV has occurred, we should have detected it earlier

  return pick_ins_value + pick_del_value;
}

double DARPSlackInitC::valueOfBestRouteAndPosInsertion(DARPGenome& gen, Vertex& newvert, int num_routes,
                                                   DistMatrix& dist, VerticesList& verticeslist,
                                                  /*out*/ int& best_route, /*out*/ int& best_insert_pos) {
  best_route        = -1;
  best_insert_pos   = -1;
  double best_value = -1;

  for (int route = 0; route<num_routes; route++) {
    int    tmp_pos   = -1;
    double tmp_value = valueOfInsertingVertexAtBestPos(gen,route,newvert,dist,verticeslist,tmp_pos);
    assert(tmp_value >= -1);

    if (tmp_value != -1 and (best_value == -1 or tmp_value < best_value) ) {
      best_value      = tmp_value;
      best_route      = route;
      best_insert_pos = tmp_pos;
    }
  }

  return best_value;
}

/*
 * TODO: This could should be refactored with the one found in DARPGenome::scoreOfInsertingVertex to avoid duplicating
 * code (right now this function is very similar to the one found in the aforementioned method.
 */
double DARPSlackInitC::valueOfInsertingVertexAtBestPos(DARPGenome& gen, int route, Vertex& newvert, DistMatrix& dist,
                                                       VerticesList& verticeslist, /*out*/ int& best_insert_pos) {
  int start_pos, end_pos;
  gen.getInsertionPos (route, newvert.id_, start_pos, end_pos);
  assert(start_pos >=0 && start_pos <= gen.routeLength(route) );
  assert(end_pos >= start_pos && end_pos <= gen.routeLength(route) ); // Can be inserted at the next pos to the last one

  ComputedTimes ctimes = computeTimes(gen,route,dist,verticeslist);

  best_insert_pos          = -1;
  double best_insert_value = -1;
  // From the starting pos, the remaining positions are evaluated for finding the best one
  // We also check of inserting the vertex after the last one
  for (int i=start_pos; i<=end_pos; i++) {
    double tmp_insert_value = computeSlackCostOfInsertion(gen,route,i,newvert,dist,verticeslist,ctimes);
    assert(tmp_insert_value >= -1);

    if (tmp_insert_value != -1 and (best_insert_value == -1 or tmp_insert_value < best_insert_value) ){
      best_insert_pos   = i;
      best_insert_value = tmp_insert_value;
    }
  }
  assert(best_insert_pos >= -1);

  return best_insert_value;
}

DARPSlackInitC::ComputedTimes DARPSlackInitC::computeTimes(DARPGenome& gen, int route, DistMatrix& dist, VerticesList& verticeslist) {
  ComputedTimes ctimes;

  computeArrivalAndWaitingTimes(gen,route,dist,verticeslist,ctimes);
  computeMaxPostTimeInterval   (gen,route, verticeslist, ctimes);

  return ctimes;
}

/*
 * We recompute these values (also computed in DARPEvaluator) for two reasons:
 *   1) The ones computed in DARPEvaluator have been adjusted using the forward time slack so they are not exactly the same
 *   2) In order to use the values computed in DARPEvaluator, we need to save them in DARPGenome and that implies and
 *      increment of the execution time of almost 2x
 */
void DARPSlackInitC::computeArrivalAndWaitingTimes(DARPGenome& gen, int route, DistMatrix& dist, VerticesList& verticeslist,
                                                 /*out*/ ComputedTimes& ctimes) {
  int route_size = gen.routeLength(route);
  ctimes.arrivalTimes.resize(route_size);
  ctimes.waitingTimes.resize(route_size);

  ctimes.arrivalTimes[0] = verticeslist.getVertex( gen.gene(route,0) ).fbegin_;
  ctimes.waitingTimes[0] = 0;

  for (int i=1; i<route_size; i++) {
    Vertex& vertexim1 = verticeslist.getVertex( gen.gene(route,i-1) );
    Vertex& vertex    = verticeslist.getVertex( gen.gene(route,i)   );

    int departureim1 = GAMax( ctimes.arrivalTimes[i-1], vertexim1.fbegin_ );
    ctimes.arrivalTimes[i]  = departureim1 + dist.getDist( vertexim1.pos_, vertex.pos_ );
    ctimes.waitingTimes[i]  = GAMax(0, vertex.fbegin_ - ctimes.arrivalTimes[i] );

    assert(ctimes.arrivalTimes[i] > 0);
    assert(ctimes.waitingTimes[i] >= 0);
  }
}

/*
 * Although almost similar to the forwardTimeSlack computed with the DarpEvaluator class, we need to recompute these
 * values since the forwardTimeSlack of the DARPEvaluator version are not computed for the last vertex and whenever there is
 * window violation.
 */
void DARPSlackInitC::computeMaxPostTimeInterval(DARPGenome& gen, int route, VerticesList& verticeslist, /*in-out*/ ComputedTimes& ctimes) {
  ctimes.maxposttimeints.resize(gen.routeLength(route));

  int lastpos = ctimes.maxposttimeints.size() - 1;
  ctimes.maxposttimeints[lastpos] = computeArrivalDiff(gen,route,lastpos,verticeslist,ctimes.arrivalTimes);

  for (int i=ctimes.maxposttimeints.size()-2; i>=0; i--) {
    // It could happen that a vertex in a previous iteration is inserted creating an infeasible route so that the computeArrival
    // Diff return a value <0. In that case we return 0 which represents the slack that is availabl
    ctimes.maxposttimeints[i] = GAMin ( GAMax(computeArrivalDiff(gen, route, i, verticeslist, ctimes.arrivalTimes),0) ,
                                        ctimes.maxposttimeints[i+1] + ctimes.waitingTimes[i+1] );
    assert(ctimes.maxposttimeints[i]>=0);
  }
}

inline int DARPSlackInitC::computeArrivalDiff(DARPGenome& gen, int route, int pos, VerticesList& verticeslist, vector<long>& arrivalTimes) {
  Vertex& vert = verticeslist.getVertex(gen.gene(route,pos));

  return vert.fend_ - GAMax(arrivalTimes[pos] ,vert.fbegin_);
}

long DARPSlackInitC::computeSlackCostOfInsertion(DARPGenome& gen, int route, int pos, Vertex& newvert, DistMatrix& dist,
                                                 VerticesList& verticeslist, ComputedTimes& ctimes) {
  if (isInsertionUnfeasible(gen,route,pos,newvert,dist,verticeslist,ctimes)) return -1;

  long newarrival, newwait, newfts;
  computeValuesOfInsertingAtPos(gen,route,pos,newvert,dist,verticeslist,ctimes,newarrival,newwait,newfts);

  long c1 = (pos > 0 ) ? computeC1(gen,route,pos,newwait,newfts,verticeslist,ctimes) : 0;

  long c2 = computeC2(newvert,newarrival,newfts);

  long c3 = computeC3(gen, route, pos, newvert, verticeslist, dist);

  assert(c1 >= 0);
  assert(c2 >= 0);
  // c3 represents the difference in distance. It can value less than 0

  return c1 + c2 + c3;
}

bool DARPSlackInitC::isInsertionUnfeasible (DARPGenome& gen, int route, int pos, Vertex& newvert, DistMatrix& dist,
                                            VerticesList& verticeslist, ComputedTimes& ctimes) {
  assert(pos <= gen.routeLength(route));

  long newarrival = computeNewArrival(gen,route,pos,newvert,dist,verticeslist,ctimes);

  if (newarrival > newvert.fend_) return true;

  if (pos != gen.routeLength(route)) {
    Vertex& vertexatpos = verticeslist.getVertex(gen.gene(route,pos));

    return ( GAMax(newarrival, newvert.fbegin_) + /*newvert.servicetime_*/ + dist.getDist(newvert.pos_,vertexatpos.pos_) >
             ctimes.arrivalTimes[pos] + ctimes.waitingTimes[pos] + ctimes.maxposttimeints[pos] );
  }

  return false;
}

long DARPSlackInitC::computeNewArrival(DARPGenome& gen, int route, int pos, Vertex& newvert, DistMatrix& dist, VerticesList& verticeslist, ComputedTimes& ctimes) {
  long newarrival = -1;
  if (pos == 0) newarrival = newvert.fbegin_;
  else {
    int vertexidatposm1   = gen.gene(route,pos-1);
    Vertex& vertexatposm1 = verticeslist.getVertex(vertexidatposm1);

    assert( ctimes.arrivalTimes[pos-1] != 0);
    long departureatposm1 = GAMax(ctimes.arrivalTimes[pos-1] , vertexatposm1.fbegin_);

    newarrival = departureatposm1 + /*vertexatposm1.service_ +*/
                 dist.getDist(vertexatposm1.pos_,newvert.pos_);
    assert(newarrival >= 0);
  }

  return newarrival;
}

void DARPSlackInitC::computeValuesOfInsertingAtPos(DARPGenome& gen, int route, int pos, Vertex& newvert, DistMatrix& dist,
                                                   VerticesList& verticeslist, ComputedTimes& ctimes,
                                                   /*out*/ long& newarrival, /*out*/ long& newwait, /*out*/ long& newfts  ) {
  newarrival = computeNewArrival(gen,route,pos,newvert,dist,verticeslist,ctimes);

  newwait = GAMax( 0, newvert.fbegin_ - newarrival);

  if (pos == gen.routeLength(route))  newfts =  newvert.fend_ - GAMax(newarrival,newvert.fbegin_);
  else {
    Vertex& vertexatpos = verticeslist.getVertex( gen.gene(route,pos) );

    long delta = dist.getDist(newvert.pos_,vertexatpos.pos_) + /*newvert.service_ +*/ + newwait;

    if (pos != 0) {
      Vertex& vertexatposm1 = verticeslist.getVertex( gen.gene(route,pos-1) );

      delta += dist.getDist(vertexatposm1.pos_,newvert.pos_) - dist.getDist(vertexatposm1.pos_,vertexatpos.pos_);
    }

    newfts = GAMin (
                    newvert.fend_ - GAMax(newarrival,newvert.fbegin_),
                    ctimes.waitingTimes[pos] + ctimes.maxposttimeints[pos] - delta
                   );
  }

  assert(newarrival >= 0 && newwait >=0 && newfts >= 0);
}


long DARPSlackInitC::computeC1(DARPGenome& gen, int route, int pos, long newwait, long newfts, VerticesList& verticeslist, ComputedTimes& ctimes) {
  assert(pos > 0);

  long beta = newwait + newfts; // Step 0
  long c1   = 0;
  int  k    = pos - 1;

  do {
    if (beta >= ctimes.maxposttimeints[k] or k == 0) break; // Step1

    if (k == pos -1) { // Step 2
      c1 = ctimes.maxposttimeints[k] - beta;
    }
    else { // Step 3
      long waitkp1 = ctimes.waitingTimes[k+1];
      long ftsk    = ctimes.maxposttimeints[k];
      if ( waitkp1 > 0 ) {
        c1 += GAMin( ftsk - beta, waitkp1 );
      }
    }

    beta += ctimes.waitingTimes[k]; // Step 4
    k--;
  } while (true);

  return c1;
}

long DARPSlackInitC::computeC2(Vertex& newvert, long newarrival, long newfts) {
  return newvert.fend_ - GAMax(newarrival,newvert.fbegin_) - newfts;
}

long DARPSlackInitC::computeC3(DARPGenome& gen, int route, int pos, Vertex& newvert, VerticesList& verticeslist, DistMatrix& dist) {
  if (pos == gen.routeLength(route)) return 0; // If we are inserting the vertex at the end, we do not need to compute the
                                               // the decrease of the time window slack of the verteices that come after

  Vertex vertpos = verticeslist.getVertex( gen.gene(route, pos) );

  long c3 = dist.getDist(newvert.pos_, vertpos.pos_);

  if (pos > 0) {
    Vertex vertposm1 = verticeslist.getVertex( gen.gene(route, pos - 1) );
    c3 += dist.getDist(vertposm1.pos_,newvert.pos_) - dist.getDist(vertposm1.pos_,vertpos.pos_);
  }

  return c3;
}




darpInitFuncT initMethod(string name) {

  darpInitFuncT initf = 0;

#define INITCOMPCASE(init)         if (name.compare(#init) == 0) initf = init;
#define INITCOMPCASEC(classn,init) if (name.compare(#init) == 0) initf = classn::init;

  INITCOMPCASE(DARPVNSInit);
  INITCOMPCASE(DARPBoSSortedRandomInit);
  INITCOMPCASE(DARPRandomInit);
  INITCOMPCASE(DARPSkybusInit);
  INITCOMPCASE(DARPObjFuncBasedFastInit);
  INITCOMPCASEC(DARPObjFuncBasedInitC,DARPObjFuncBasedInit);
  INITCOMPCASEC(DARPRegretInsertionInitC,DARPRegretInsertionInit);
  INITCOMPCASEC(DARPSlackInitC,DARPSlackInit);

  if (initf == 0) throw runtime_error("Unrecognized init name");

  return initf;
}




/***************************** SKYBUS SOURCE CODE ***********************************************/

#include "skybusInit/heuristicoInsercion.h"

void constructGenomeSkybus(DARPGenome& gen, DistMatrix& distmatrix, VerticesList& verticeslist,
              long vehicleCapacity,long maxDelay, double DARP_ALPHA, double DARP_BETA, long nVehicles) {

  std::vector<Request> requests=verticeslist.getReqsList();

  std::vector<Vehiculo> vehiculos;
  heuristicoInicializacion(vehiculos,vehicleCapacity,maxDelay,DARP_ALPHA,DARP_BETA,verticeslist, nVehicles,distmatrix);

  int vertices_inserted = 0;
  for(int i=0;i<vehiculos.size();i++){
    Vehiculo v=vehiculos[i];
    for(int j=0; j<v.getSolucion().size();j++){
      gen.pushBackVertex(i,v.getSolucion()[j].getIdPeticion());
      vertices_inserted++;
    }
  }
  if (vertices_inserted  != requests.size() * 2) {
    stringstream msg; msg << "different number of vertices, expected=" << requests.size() * 2 << " inserted=" << vertices_inserted;
    throw runtime_error(msg.str());
  }

}
void DARPSkybusInit(GAGenome& g, DistMatrix& dist,VerticesList& verticeslist) {
  DARPGenome& gen = dynamic_cast<DARPGenome&>(g);
  assert(&gen);

  long   vehicleCapacity = gen.vehicleCapcity();
  long   maxDelay = Vertex::maxDelay();
  double DARP_ALPHA = gen.loadVWeight();
  double DARP_BETA = gen.TWVWeight();
  long   nVehicles = gen.length();

  constructGenomeSkybus(gen,dist,verticeslist,vehicleCapacity,maxDelay,DARP_ALPHA, DARP_BETA, nVehicles);
}

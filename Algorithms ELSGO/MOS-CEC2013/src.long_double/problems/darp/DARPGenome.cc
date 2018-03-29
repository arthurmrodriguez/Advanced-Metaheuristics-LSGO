#include "DARPGenome.h"
#include "aux.h"
#include "DARPEvaluator.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <algorithm>

DARPEvaluator* DARPGenome::darpeval__ = 0;

DARPGenome::DARPGenome(unsigned x, long double evalWeightsUpdateValue, GAGenome::Evaluator f, void* u)
  : GA1DArrayGenome< vector<int> > (x, f, u), RoutingGenome(), GAGenome(),
    modifiedRoutes_             (x, false),
    costPerRoute_               (x, 0    ),
    loadVPerRoute_              (x, 0    ),
    TWVPerRoute_                (x, 0    ),
    rideVPerRoute_              (x, 0    ),
    pickupDelayPerRoute_        (x, 0.0  ),
    deliveryDelayPerRoute_      (x, 0.0  ),
    localSearchedRoutes_        (x, LSUNEXPLORED),
    localSearchLastPosExplored_ (x, make_pair(-1,-1) ),
    feasibleRoutes_             (x, true),
    evalWeightsUpdateValue_     (evalWeightsUpdateValue) {
}

void DARPGenome::copy (const GAGenome& orig) {
  const DARPGenome* other = dynamic_cast<const DARPGenome*>(&orig); assert(other);
  if (other->length() != length()) resize(other->length());

  GA1DArrayGenome< vector<int> >::copy (orig);
  RoutingGenome::onlyCopyAttributes(orig);

  modifiedRoutes_             = other->modifiedRoutes_;
  costPerRoute_               = other->costPerRoute_;
  loadVPerRoute_              = other->loadVPerRoute_;
  TWVPerRoute_                = other->TWVPerRoute_;
  rideVPerRoute_              = other->rideVPerRoute_;
  pickupDelayPerRoute_        = other->pickupDelayPerRoute_;
  deliveryDelayPerRoute_      = other->deliveryDelayPerRoute_;
  localSearchedRoutes_        = other->localSearchedRoutes_;
  localSearchLastPosExplored_ = other->localSearchLastPosExplored_;
  feasibleRoutes_             = other->feasibleRoutes_;
  evalWeightsUpdateValue_     = other->evalWeightsUpdateValue_;
}

int DARPGenome::resize (int x) {
  int orig_length = length();

  GA1DArrayGenome< vector<int> >::resize(x);

  modifiedRoutes_.resize(x);
  costPerRoute_.resize(x);
  loadVPerRoute_.resize(x);
  TWVPerRoute_.resize(x);
  rideVPerRoute_.resize(x);

  pickupDelayPerRoute_.resize(x);;
  deliveryDelayPerRoute_.resize(x);;

  localSearchedRoutes_.resize(x);
  localSearchLastPosExplored_.resize(x);
  feasibleRoutes_.resize(x);

  for (int i=orig_length; i<length(); i++) setRouteWasModified(i);

  return this->sz;
}

int DARPGenome::routeLength(int route) const {
  assert(route >= 0 && route < this->length() );
  return this->a[route].size();
}

int DARPGenome::gene(int route, int pos) const {
  assert(route >= 0 && route < this->length() );
  assert(pos >= 0 && pos < this->a[route].size() );
  return this->a[route][pos];
}

int DARPGenome::gene(int route, int pos, int value) {
  assert(route >= 0 && route < this->length() );
  assert(pos >= 0 && pos < this->a[route].size() );

  setRouteWasModified(route);

  return this->a[route][pos] = value;
}

bool DARPGenome::modifiedRoute(int route, bool m) {
  assert(route >= 0 && route < modifiedRoutes_.size() );
  return this->modifiedRoutes_[route] = m;
}

void DARPGenome::setRouteWasModified(int route) {
  assert(route >= 0 && route < modifiedRoutes_.size() && route < localSearchedRoutes_.size() && route < feasibleRoutes_.size());
  _evaluated                  = gaFalse;
  modifiedRoutes_[route]      = true;
  localSearchedRoutes_[route] = LSUNEXPLORED;
  feasibleRoutes_[route]      = false;
}

long DARPGenome::costPerRoute(int route) const {
  assert(route >= 0 && route < costPerRoute_.size() );
  const_cast<DARPGenome&>(*this).evaluateRoute(route);
  return this->costPerRoute_[route];
}

long DARPGenome::costPerRoute(int route, long value) {
  assert(route >= 0 && route < costPerRoute_.size() );
  assert(value >= 0);
  return this->costPerRoute_[route] = value;
}

long DARPGenome::cost() const {
  return accumulate(costPerRoute_.begin(), costPerRoute_.end(), 0);
}

long DARPGenome::loadVPerRoute(int route) const {
  assert(route >= 0 && route < loadVPerRoute_.size() );
  const_cast<DARPGenome&>(*this).evaluateRoute(route);
  return this->loadVPerRoute_[route];
}

long DARPGenome::loadVPerRoute(int route, long value) {
  assert(route >= 0 && route < loadVPerRoute_.size() );
  assert(value >= 0);
  return this->loadVPerRoute_[route] = value;
}

long DARPGenome::loadV() const {
  return accumulate(loadVPerRoute_.begin(), loadVPerRoute_.end(), 0);
}

long DARPGenome::TWVPerRoute(int route) const {
  assert(route >= 0 && route < TWVPerRoute_.size() );
  const_cast<DARPGenome&>(*this).evaluateRoute(route);
  return this->TWVPerRoute_[route];
}

long DARPGenome::TWVPerRoute(int route, long value) {
  assert(route >= 0 && route < TWVPerRoute_.size() );
  assert(value >= 0);
  return this->TWVPerRoute_[route] = value;
}

long DARPGenome::TWV() const {
  return accumulate(TWVPerRoute_.begin(), TWVPerRoute_.end(), 0);
}

long DARPGenome::rideVPerRoute(int route) const {
  assert(route >= 0 && route < rideVPerRoute_.size() );
  const_cast<DARPGenome&>(*this).evaluateRoute(route);
  return this->rideVPerRoute_[route];
}

long DARPGenome::rideVPerRoute(int route, long value) {
  assert(route >= 0 && route < rideVPerRoute_.size() );
  assert(value >= 0);
  return this->rideVPerRoute_[route] = value;
}

long DARPGenome::rideV() const {
  return accumulate(rideVPerRoute_.begin(), rideVPerRoute_.end(), 0);
}

long double DARPGenome::deliveryDelayPerRoute(int route) const {
  assert(route >= 0 && route < deliveryDelayPerRoute_.size() );
  const_cast<DARPGenome&>(*this).evaluateRoute(route);
  return this->deliveryDelayPerRoute_[route];
}

long double DARPGenome::deliveryDelayPerRoute(int route, long double value) {
  assert(route >= 0 && route < deliveryDelayPerRoute_.size() );
  assert(value >= 0);
  return this->deliveryDelayPerRoute_[route] = value;
}

long double DARPGenome::deliveryDelay() const {
  return accumulate(deliveryDelayPerRoute_.begin(), deliveryDelayPerRoute_.end(), 0.0);
}

long double DARPGenome::pickupDelayPerRoute(int route) const {
  assert(route >= 0 && route < pickupDelayPerRoute_.size() );
  const_cast<DARPGenome&>(*this).evaluateRoute(route);
  return this->pickupDelayPerRoute_[route];
}

long double DARPGenome::pickupDelayPerRoute(int route, long double value) {
  assert(route >= 0 && route < pickupDelayPerRoute_.size() );
  assert(value >= 0);
  return this->pickupDelayPerRoute_[route] = value;
}

long double DARPGenome::pickupDelay() const {
  return accumulate(pickupDelayPerRoute_.begin(), pickupDelayPerRoute_.end(), 0.0);
}

DARPGenome::lsstatus DARPGenome::localSearchedRoute(int route) const {
  assert(route >= 0 && route < localSearchedRoutes_.size() );
  const_cast<DARPGenome&>(*this).evaluateRoute(route);
  return this->localSearchedRoutes_[route];
}

DARPGenome::lsstatus DARPGenome::localSearchedRoute(int route, lsstatus m) {
  assert(route >= 0 && route < localSearchedRoutes_.size() );
  return this->localSearchedRoutes_[route] = m;
}

int DARPGenome::localSearchLastPosExplored(int route) const {
  assert(route >= 0 && route < localSearchedRoutes_.size() );
  assert(this->localSearchedRoutes_[route] == LSUNFINISHED);
  assert(this->localSearchLastPosExplored_[route].first >= 0 &&
         this->localSearchLastPosExplored_[route].first < routeLength(route));

  return this->localSearchLastPosExplored_[route].first;
}

int DARPGenome::localSearchLastCritVertexSearchPosExplored(int route) const {
  assert(route >= 0 && route < localSearchedRoutes_.size() );
  assert(this->localSearchedRoutes_[route] == LSUNFINISHED);
  // -1 forces to recompute the starting position to be explored. It is neccessary in certain scenarios (can be seen in the  ls
  // implementation
  assert(this->localSearchLastPosExplored_[route].second >= -1 &&
         this->localSearchLastPosExplored_[route].second < routeLength(route));

  return this->localSearchLastPosExplored_[route].second;
}

void DARPGenome::localSearchLastPosExplored(int route, int pos, int crivertsearchpos) {
  assert(route >= 0 && route < localSearchedRoutes_.size() );
  assert(pos >= 0 && pos <= routeLength(route) );                            // Could be equsl to the routelength if we are exploring
  assert(crivertsearchpos >= -1 && crivertsearchpos <= routeLength(route) );  // Could be equsl to the routelength if we are exploring
  assert(this->localSearchedRoutes_[route] == LSUNFINISHED);                 // to insert it at the next position of the last vertex

  this->localSearchLastPosExplored_[route] = make_pair(pos,crivertsearchpos);
}

int DARPGenome::randomRoutePos() const {
  int size = length();
  if ( size == 1) return 0;
  else            return GARandomInt(0,length()-1);
}

int DARPGenome::randomNonEmptyRoutePos() const {
  int route;
  do route = randomRoutePos();
  while (routeLength(route) == 0);

  return route;
}

int DARPGenome::emptyRouteOrRandomRoutePos() const {
  int route = -1;
  for (int i=0; i<length(); i++) {
    if (routeLength(i) == 0) {
      route = i;
      break;
    }
  }
  if (route == -1) route = randomRoutePos();

  return route;
}

void DARPGenome::pushBackVertex( int route, int seq) {
  assert(route >= 0 && route < length() );
  this->a[route].push_back(seq);

  setRouteWasModified(route);
}

void DARPGenome::insertVertex (int route, int pos, int value) {
  assert(route >= 0 && route < this->length() );
  assert(pos >=0 && pos <= this->a[route].size());

  this->a[route].resize(this->a[route].size()+1);

  for (int i=this->a[route].size()-1; i>pos; i--) {
    this->a[route][i] = this->a[route][i-1];
  }
  this->a[route][pos] = value;

  setRouteWasModified(route);
}

void DARPGenome::insertVertices (int route, int pos, const list<int>& values) {
  assert(route >= 0 && route < this->length() );
  assert(pos >=0 && pos <= this->a[route].size());

  int old_size = this->a[route].size();

  this->a[route].resize(this->a[route].size()+values.size());

  // First we move the old values that are going to be at the right to the rightmost positions
  int num_tocopy = old_size - pos;
  for (int i=this->a[route].size()-1; i>this->a[route].size()-1-num_tocopy; i--) {
    this->a[route][i] = this->a[route][i-values.size()];
  }

  // Then we copy the new values
  list<int>::const_iterator it=values.begin();
  for (int i=pos; i<pos+values.size(); i++, it++) {
    this->a[route][i] = *it;
  }

  setRouteWasModified(route);
}

void DARPGenome::getInsertionPos (int route, int vertex_id, int& start_pos, int& end_pos) {
  const CostMatrix&   costmatrix   = darpeval__->costMatrix();
  const VerticesList& verticeslist = darpeval__->vertList();

  assert(routeLength(route) <= verticeslist.size());

  Vertex& vertex = verticeslist.getVertex( vertex_id );

  start_pos = -1;
  end_pos   = 0;

  for (int i=0; i<routeLength(route); i++) {

    int vertextmp_id   = gene(route,i);
    Vertex& vertex_tmp = verticeslist.getVertex( vertextmp_id );

    if (vertex.canVertexGoBefore(vertex_tmp,costmatrix)) {
      if (start_pos == -1) start_pos = end_pos = i;
    }
    if (vertex.isDelivery() && vertex.getSiblingVertexId() == vertextmp_id) { // If the vertex is a delivery v, the starting position
      start_pos = i+1;                                                        // must be the next one to the pickup vertex
    }
    if (vertex.canVertexGoAfter(vertex_tmp,costmatrix)) {
      end_pos = i+1;
    }

  }

  if (start_pos == -1) { // it must be inserted at the end
    start_pos = end_pos = routeLength(route);
    assert( routeLength(route) == 0 ||
            vertex.canVertexGoBefore(
                              verticeslist.getVertex( gene(route,routeLength(route)-1) ),
                              costmatrix
                                    ) == false );
  }
  assert(start_pos != -1 && start_pos <= end_pos);
}

void DARPGenome::getInsertionPos (int route, int vertex_id, int& start_pos, int& end_pos, int start_pos_requested) {
  // The feasible positions are computed
  getInsertionPos(route,vertex_id,start_pos,end_pos);
  start_pos = GAMax(start_pos,start_pos_requested);
}

int DARPGenome::insertVertexInBestPos(int route, int vert_id, int start_pos, int end_pos, long double scoretoimprove, long maxevals) {
  assert(start_pos >=0 && start_pos <= routeLength(route) );
  assert(end_pos >= start_pos && end_pos <= routeLength(route) ); // Can be inserted at the next pos to the last one
  assert(maxevals > 0);

  int orig_evals = nevals();
  int used_evals = 0;
  int insert_pos = start_pos;

  DARPGenome* bestinsgen = dynamic_cast<DARPGenome*>( clone() );
  bestinsgen->insertVertex(route,insert_pos,vert_id);
  bestinsgen->evaluate();
  used_evals++;

  // From the starting pos, the remaining positions are evaluated for finding the best one
  // We also check of inserting the vertex after the last one
  for (int i=insert_pos+1; i<=end_pos && used_evals<maxevals ; i++) {
    DARPGenome* tmpgen = dynamic_cast<DARPGenome*>( clone() );
    tmpgen->insertVertex(route,i,vert_id);
    used_evals++;

    if (GAGenome::compareScores(tmpgen->score(), bestinsgen->score() ) == GAGenome::BETTER) { // tmp is better
      insert_pos = i;

      std::swap(bestinsgen,tmpgen);

      if (GAGenome::compareScores(bestinsgen->score(), scoretoimprove) == GAGenome::BETTER) { // the best score has been reached
        delete tmpgen;                                                                        // so we break. Note that we need
        break;                                                                                // to delete tmpgen before breaking
      }
    }

    delete tmpgen;
  }

  assert(insert_pos >= 0);

  copy(*bestinsgen);
  nevals(orig_evals + used_evals);

  delete bestinsgen;

  return insert_pos;
}

int DARPGenome::insertVertexInFirstFeasiblePos(int route, int vert_id, bool userandominsertpos, int start_pos, int end_pos) {
  int orig_evals = nevals();
  int used_evals = 0;

  vector<int> insert_positions; for (int i=start_pos; i<=end_pos; i++) insert_positions.push_back(i);
  if (userandominsertpos) random_shuffle(insert_positions.begin(),insert_positions.end());

  DARPGenome* bestinsgen    = dynamic_cast<DARPGenome*>( clone() );
  int         bestinsertpos = insert_positions[0];
  bestinsgen->insertVertex(route,bestinsertpos,vert_id);
  bestinsgen->evaluate();
  used_evals++;

  if (!bestinsgen->feasible()) {

    for (int i=1; i<insert_positions.size() ; i++) {
      DARPGenome* tmpgen = dynamic_cast<DARPGenome*>( clone() );
      tmpgen->insertVertex(route,insert_positions[i],vert_id);
      tmpgen->evaluate();
      used_evals++;

      if (tmpgen->feasible() ||
          GAGenome::compareScores(tmpgen->score(), bestinsgen->score() ) == GAGenome::BETTER) {
        std::swap(bestinsgen,tmpgen);
        bestinsertpos = insert_positions[i];

        if (bestinsgen->feasible() ) {
          delete tmpgen;
          break;
        }
      }

      delete tmpgen;
    }

  }

  assert(bestinsertpos >= 0);

  copy(*bestinsgen);
  nevals(orig_evals + used_evals);

  delete bestinsgen;

  return bestinsertpos;

}

int DARPGenome::insertVertexInFirstFeasiblePos(int route, int vert_id, bool userandominsertpos) {
  int start_pos=-1, end_pos=-1;
  getInsertionPos(route,vert_id,start_pos,end_pos);

  return insertVertexInFirstFeasiblePos(route,vert_id,userandominsertpos,start_pos,end_pos);
}

int DARPGenome::insertVertexInBestPos(int route, int vert_id) {
  int start_pos=-1, end_pos=-1;
  getInsertionPos(route,vert_id,start_pos,end_pos);

  return insertVertexInBestPos(route, vert_id, start_pos, end_pos);
}

int DARPGenome::insertVertexInBestPosFromStartPos(int route, int vert_id, int start_pos_requested,
                                                  long double scoretoimprove, long maxevals) {
  int start_pos, end_pos;
  getInsertionPos(route,vert_id,start_pos,end_pos,start_pos_requested);

  return insertVertexInBestPos(route, vert_id,start_pos,end_pos, scoretoimprove, maxevals);
}

int DARPGenome::insertVertexInFirstFeasiblePosFromStartPos(int route, int vert_id, int start_pos_requested, bool userandominsertpos) {
  int start_pos, end_pos;
  getInsertionPos(route,vert_id,start_pos,end_pos,start_pos_requested);

  return insertVertexInFirstFeasiblePos(route,vert_id,userandominsertpos,start_pos,end_pos);
}

list<int> DARPGenome::insertVerticesInBestPos (int route, list<int> vertices_ids) {
  // Note that a copy of vertices_ids is needed since we are removing or the first vertex or the critical one
  int vert_id;
  int vert_pos;
  list<int> ins_positions;

  while ( ! vertices_ids.empty() ) {
    // All vertices are inserted one-by-one.
    vert_id      = vertices_ids.front();
    Vertex& vert = darpeval__->vertList().getVertex(vert_id);

    // If the request is not critical, the corresponding critical request is found to be inserted first
    if ( ! vert.critic_ ) {
      int vert_id_sibling  = vert.getSiblingVertexId();
      if ( find(vertices_ids.begin(), vertices_ids.end(), vert_id_sibling) != vertices_ids.end() ) { // found
        vert_id = vert_id_sibling;
      }
    }
    assert(find(vertices_ids.begin(),vertices_ids.end(),vert_id) != vertices_ids.end());

    ins_positions.push_back( insertVertexInBestPos(route,vert_id) );

    // Once the vertex is inserted it is removed from the list of vertices
    vertices_ids.remove(vert_id);
  }
  // repeat until the list of vertices ids is empty

  // Note that at the end we have called evalRoute as many times as insertVertexInBestPos has needed for each vertex
  // but, at the end, there is still one evaluation missing that should be counted the next time score is called since
  // the state of the best insertion is not saved and needs to be recomputed


  return ins_positions;
}

//TODO: Analyze this method for a possible loss in performance due to the multiple DARPGenome copies
long double DARPGenome::scoreOfInsertingVertices(int route, const list<int> vertices_ids) const {
  assert(route >= 0 && route < length() );
  
  DARPGenome tmpgen(*this); // By using a new genome we make sure that the cache flags remain in the original genome

  tmpgen.insertVerticesInBestPos(route,vertices_ids);

  long double score = tmpgen.score(); // To force the computation of the score and the increase in the nevals

  DARPGenome& nevalsgen = const_cast<DARPGenome&>(*this);

  nevalsgen.nevals(nevals()+tmpgen.nevals());

  return score;
}

void DARPGenome::swapSeqs(int route1, const list<int>& seq1, int route2, const list<int>& seq2) {
  removeVertices(route1,seq1);
  removeVertices(route2,seq2);

  insertVerticesInBestPos(route1,seq2);
  insertVerticesInBestPos(route2,seq1);
}

void DARPGenome::swapWithInsertPosSeqs(int from_route, int from_ins_pos, const list<int>& from_seq,
                                       int to_route,   int to_ins_pos,   const list<int>& to_seq) {

  removeVertices(to_route,to_seq);
  insertVertices(to_route,to_ins_pos,from_seq);

  removeVertices(from_route,from_seq);
  insertVertices(from_route,from_ins_pos, to_seq);
}



//TODO: Analyze this method for a possible loss in performance due to the multiple DARPGenome copies
long double DARPGenome::scoreOfInsertingVertex(int route, int pos, int vert_id) const {
  assert(route >= 0 && route < length() );
  assert(pos >= 0 && pos <= routeLength(route));

  DARPGenome tmpgen(*this); // By using a new genome we make sure that the cache flags remain in the original genome

  tmpgen.insertVertex(route,pos,vert_id);

  DARPGenome& nevalsgen = const_cast<DARPGenome&>(*this);
  long double score = tmpgen.score(); // To force the computation of the score and the increase in the nevals

  nevalsgen.nevals(nevals()+tmpgen.nevals());

  return score;
}

/*
 * It is assumed that the value should be included in the genome
 */
int DARPGenome::removeVertex (int route, int value, int startpos) {
  assert(route >= 0 && route < this->length() );
  int pos = findPosOfValue(value,this->a[route], startpos);

  assert(pos >= 0);
  removeVertexOfPos(route,pos);

  return pos;
}

int DARPGenome::removeVertexOfPos (int route, int pos) {
  assert(route >= 0 && route < this->length() );
  assert(pos >=0 && pos <= this->a[route].size());

  int vertex_removed = this->a[route][pos];

  for (int i=pos; i<this->a[route].size()-1; i++) {
    this->a[route][i] = this->a[route][i+1];
  }

  this->a[route].resize( this->a[route].size() - 1 );

  setRouteWasModified(route);

  return vertex_removed;
}

list<int> DARPGenome::removeVerticesFrom (int route, int pos, int length) {
  assert(route >= 0 && route < this->length() );
  assert(pos >=0 && pos <= this->routeLength(route));
  assert(length <= a[route].size() );

  vector<int> new_values;
  list<int>   removed_values;

  int last_pos = (pos + length-1) % this->a[route].size();

  for (int i=0; i<this->a[route].size(); i++) {
    if ( (i>=pos && i<pos+length) || (last_pos < pos && i <= last_pos) ) {
      removed_values.push_back(this->a[route][i]);
      continue;
    }

    new_values.push_back(this->a[route][i]);
  }
  assert(new_values.size() == this->a[route].size() - length);
  assert(removed_values.size() == length);

  this->a[route].resize(new_values.size());
  for (int i=0; i<new_values.size(); i++) this->a[route][i] = new_values[i];

  setRouteWasModified(route);

  return removed_values;
}

void DARPGenome::removeVertices(int route, const list<int>& seq) {
  for (list<int>::const_iterator it=seq.begin(); it!=seq.end(); it++) {
    removeVertex(route,*it);
  }
}

long double DARPGenome::scoreOfRemovingVertices(int route, const list<int>& vertices_ids) const {
  assert(route >= 0 && route < length() );

  DARPGenome tmpgen(*this);

  for(list<int>::const_iterator it=vertices_ids.begin(); it!=vertices_ids.end(); it++) {
    tmpgen.removeVertex(route,*it);
  }

  long double score = tmpgen.score(); 

  DARPGenome& evals_gen = const_cast<DARPGenome&>(*this); //Neccesary for improving the nevals of gen
  evals_gen.nevals(nevals()+1);

  return score;
}

void DARPGenome::moveVertices (int fromroute, int toroute, const list<int>& seq) {
  removeVertices(fromroute,seq);
  insertVerticesInBestPos(toroute,seq);
}

int DARPGenome::findPosOfVertex(int route, int value) const {
  int pos=-1;
  for (int i=0; i<this->routeLength(route); i++) {
    if (this->gene(route,i) == value) {
      pos = i;
      break;
    }
  }
  return pos;
}

int DARPGenome::write (ostream& os) const {

  os << "evalWeightsUpdateValue_: " << evalWeightsUpdateValue_ << endl;
  os << "score:                   " << score() << endl;
  os << endl;
  for (int i=0; i<length(); i++) {
    os << "route " << setw(2) << i << ": ";
    os << "mod = " << modifiedRoutes_[i] << ", ls status = " << localSearchedRoutes_[i] << ", ";
    os << " ls last pos explored=" << localSearchLastPosExplored_[i].first << ", ";
    os << " ls last crit vert search pos explored=" << localSearchLastPosExplored_[i].second << ", ";
    os << "cost = " << setw(9) << costPerRoute_[i] << ", loadV = " << setw(6) << loadVPerRoute_[i] << ", ";
    os << "TWV = " << setw(6) << TWVPerRoute_[i] << ", rideV = " << setw(6) << rideVPerRoute_[i] << ", ";
    os << "pickup delay = " << setprecision(2) << setw(10) << pickupDelayPerRoute_[i] << ", delivery delay = " << setprecision(2) << setw(10) << deliveryDelayPerRoute_[i];
    os << " ===> ";
    for (int j=0; j<routeLength(i); j++) {
      os << setw(6) << gene(i,j) << " ";
    }
    os << endl;
  }

  return 0;
}

void DARPGenome::printRoutes () const {

  cout << "evalWeightsUpdateValue_: " << evalWeightsUpdateValue_ << endl;
  cout << endl;
  for (int i=0; i<length(); i++) {
    cout << "route " << setw(2) << i << ": ";
    for (int j=0; j<routeLength(i); j++) {
      cout << setw(6) << gene(i,j) << " ";
    }
    cout << endl;
  }
}


void DARPGenome::writeObject(ostream& os) const{
  GAGenome::writeObject(os);

  os.write ( (char*) (&feasible_), sizeof(feasible_) );

  os.write ( (char*) (&nx),  sizeof (nx) );
  os.write ( (char*) (&minX),sizeof (minX) );
  os.write ( (char*) (&maxX),sizeof (maxX) );

  for (int routepos=0; routepos<nx; routepos++) {
    int routesize = routeLength(routepos);
    os.write( (char*) (&routesize), sizeof(routesize));
    for (int i=0; i<routesize; i++) {
      int vertex = gene(routepos,i);
      os.write ( (char*) (&vertex), sizeof (vertex) );
    }
  }

#define WRITEVECTOR(name,type)  \
  {\
    int tmpsize = name.size();\
    os.write( (char*) ( &tmpsize ), sizeof( tmpsize ) ); \
    for (int i=0; i<name.size(); i++){ \
      type tmp = name[i];\
      os.write( (char*) (&tmp), sizeof(tmp)  );\
    }\
  }

  WRITEVECTOR(modifiedRoutes_       ,bool);

  WRITEVECTOR(costPerRoute_         ,long);
  WRITEVECTOR(loadVPerRoute_        ,long);
  WRITEVECTOR(TWVPerRoute_          ,long);
  WRITEVECTOR(rideVPerRoute_        ,long);

  WRITEVECTOR(pickupDelayPerRoute_  ,long double);
  WRITEVECTOR(deliveryDelayPerRoute_,long double);

  typedef pair<int,int> lastpostype;

  WRITEVECTOR(localSearchedRoutes_  ,     lsstatus);
  WRITEVECTOR(localSearchLastPosExplored_,lastpostype);
  WRITEVECTOR(feasibleRoutes_       ,     bool);

  os.write ( (char*) (&evalWeightsUpdateValue_), sizeof(evalWeightsUpdateValue_) );

}

void DARPGenome::readObject (istream& is) {
  GAGenome::readObject(is);

  is.read( (char*) (&feasible_), sizeof(feasible_) );

  int tmpnx;
  is.read ( (char*) (&tmpnx), sizeof (tmpnx) );
  is.read ( (char*) (&minX),  sizeof (minX) );
  is.read ( (char*) (&maxX),  sizeof (maxX) );

  resize(tmpnx);
  assert(nx==tmpnx);

  for (int routepos=0; routepos<nx; routepos++) {

    int routesize;
    is.read( (char*) (&routesize), sizeof(routesize));
    a[routepos].resize(routesize);

    for (int i=0; i<routesize; i++) {
      int vertex;
      is.read ( (char*) (&vertex), sizeof (vertex) );
      a[routepos][i] = vertex;
    }
  }

#define READVECTOR(name,type) \
  { \
    int size; is.read( (char*) (&size), sizeof(size) );\
    name.resize(size);\
    for (int i=0; i<size; i++) { \
      type tmp; \
      is.read( (char*) (&tmp), sizeof(tmp) );\
      name[i] = tmp;\
    }\
  }


  READVECTOR(modifiedRoutes_       ,bool);

  READVECTOR(costPerRoute_         ,long);
  READVECTOR(loadVPerRoute_        ,long);
  READVECTOR(TWVPerRoute_          ,long);
  READVECTOR(rideVPerRoute_        ,long);

  READVECTOR(pickupDelayPerRoute_  ,long double);
  READVECTOR(deliveryDelayPerRoute_,long double);

  typedef pair<int,int> lastpostype;

  READVECTOR(localSearchedRoutes_       ,lsstatus);
  READVECTOR(localSearchLastPosExplored_,lastpostype);
  READVECTOR(feasibleRoutes_            ,bool);

  is.read ( (char*) (&evalWeightsUpdateValue_), sizeof (evalWeightsUpdateValue_) );
}


long double DARPGenome::nonPenalizedScore() const {
  return darpeval__->nonPenalizedScore(*this);
}

void DARPGenome::evaluateAllRoutes() {
  bool isfeasible=true;

  for (unsigned route=0; route<length(); route++) {
    evaluateRoute(route);
    isfeasible &= feasibleRoutes_[route];
  }
  assert(isfeasible || loadV() > 0 || rideV() > 0 || TWV() > 0 );

  feasible(isfeasible);
}

void DARPGenome::evaluateRoute(int route) {
  if (modifiedRoute(route)) {
    vector<int> path (routeLength(route), 0);
    for (unsigned j = 0; j < routeLength(route); j++) path[j] = gene(route, j);

    map<int,long> arTimes;
    map<int,long> waitTimes;
    map<int,long> beginningServiceTimes;
    map<int,long> depTimes;
    map<int,long> rideTimes;
    map<int,long> forwTimeSlacks;
    map<int,long> loadsWhenLeavingVertex;

    long   cost;
    long   loadViolation, TWViolation, rideViolation;
    long double pickupDelay, deliveryDelay;

    darpeval__->evalRoute(path,arTimes,waitTimes,beginningServiceTimes,depTimes,rideTimes,forwTimeSlacks,loadsWhenLeavingVertex,cost,loadViolation,TWViolation,rideViolation,pickupDelay,deliveryDelay);

    costPerRoute         (route, cost         );
    loadVPerRoute        (route, loadViolation);
    TWVPerRoute          (route, TWViolation  );
    rideVPerRoute        (route, rideViolation);
    pickupDelayPerRoute  (route, pickupDelay  );
    deliveryDelayPerRoute(route, deliveryDelay);
    modifiedRoute        (route, false        );

    feasibleRoutes_[route] = ( loadViolation==0 && TWViolation==0 && rideViolation == 0 );

    assert( (loadViolation==0 && loadVPerRoute(route) == 0) || (loadViolation!=0 && loadVPerRoute(route) > 0) );
    assert( ( (TWViolation==0 && TWVPerRoute(route)   == 0) || (TWViolation  !=0 &&   TWVPerRoute(route) > 0) ));
    assert( (rideViolation==0 && rideVPerRoute(route) == 0) || (rideViolation!=0 && rideVPerRoute(route) > 0) );
  }
}

int DARPGenome::evalRouteCalls() const {
  return darpeval__->evalRouteCalls();
}

void DARPGenome::updatePenalizations() {
#define UPDATEEVALWEIGHT(measure) \
    if (measure() > 0.0 ) darpeval__->measure##Weight( darpeval__->measure##Weight() * ( 1.0 + evalWeightsUpdateValue_) ); \
    else                  darpeval__->measure##Weight( darpeval__->measure##Weight() / ( 1.0 + evalWeightsUpdateValue_) );

  UPDATEEVALWEIGHT(loadV);
  UPDATEEVALWEIGHT(TWV);
  UPDATEEVALWEIGHT(rideV);
}

int DARPGenome::numRoutes() const {
  return length();
}

bool DARPGenome::moreThanOneRoute() {
  int num_full_routes=0;
  for (int i=0; i<length(); i++) {
    if (routeLength(i) > 0) {
      num_full_routes++;
      if (num_full_routes > 1) return true;
    }
  }
  return false;
}

void DARPGenome::addRoute (const GAGenome& g, int routepos) {
  const DARPGenome& gen = dynamic_cast<const DARPGenome&>(g); assert(&gen);
  assert(routepos >=0 && routepos < gen.length() );

  int orig_length = length();
  int new_length  = orig_length + 1;

  resize( new_length ); // Sets all the cached attributes to the corresponding size and calls setRouteWasModified
                        // for the new routes (in case the size is larger)

#define COPYADDROUTE(name) this->name[orig_length] = gen.name[routepos];

  COPYADDROUTE(a);
  COPYADDROUTE(modifiedRoutes_);

  COPYADDROUTE(costPerRoute_);
  COPYADDROUTE(loadVPerRoute_);
  COPYADDROUTE(TWVPerRoute_);
  COPYADDROUTE(rideVPerRoute_);

  COPYADDROUTE(pickupDelayPerRoute_);
  COPYADDROUTE(deliveryDelayPerRoute_);

  COPYADDROUTE(localSearchedRoutes_);
  COPYADDROUTE(localSearchLastPosExplored_);

  COPYADDROUTE(feasibleRoutes_);

  bool isfeasible = (new_length == 1) ? true : feasible_;
  isfeasible &= feasibleRoutes_[orig_length];
  feasible(isfeasible);

  // If the evalWeightsUpdateValue_ is bigger than the current one, maybe it was due to
  // added route so, by default, we follow this approach that favors higher penalizations
  if (gen.evalWeightsUpdateValue_ > evalWeightsUpdateValue_) evalWeightsUpdateValue_ = gen.evalWeightsUpdateValue_;
}

/*
 * Note: does not call the previous method for efficiency
 */
void DARPGenome::addRoutes(const GAGenome& g) {
  const DARPGenome& gen = dynamic_cast<const DARPGenome&> (g); assert(&gen);

  for (int routepos=0; routepos < gen.length(); routepos++) {
    // previous line when the copy was conducted also here
    //assert(orig_length+routepos < length());
    //a[orig_length+routepos] = gen.a[routepos];
    addRoute(g,routepos);
  }
}

void DARPGenome::emptyRoutes() {
  resize(0);
}

RoutingGenome* DARPGenome::cloneWithoutEmptyRoutes() const {
  // Not the most efficient way but the most secure and quicker method to implement
  // TODO: In a future this should be implemented more efficiently
  DARPGenome* newgen = dynamic_cast<DARPGenome*>(clone());
  newgen->emptyRoutes();
  const GAGenome& thisgen = *this;

  for (int i=0; i<numRoutes(); i++) {
    if (routeLength(i) > 0) {
      newgen->addRoute(thisgen,i);
    }
  }

  return newgen;
}

void DARPGenome::checkNoVertexIsMissing() const {
  for (int route=0; route<length(); route++) {
    for (int pos=0; pos<routeLength(route); pos++) {
      int vert_id = gene(route,pos);
      int vert_id_sibling  = Vertex::getSiblingVertId(vert_id);

      bool found = false;
      int start_pos=-1;
      int end_pos=-1;

      if (Vertex::isVertIdPickup(vert_id)) {
        start_pos = pos+1;
        end_pos   = routeLength(route);
      }
      else { // delivery
        start_pos = 0;
        end_pos   = pos;
      }

      for (int i=start_pos; i<end_pos; i++) {
        if ( gene(route,i) == vert_id_sibling ) {
          found = true;
          break;
        }
      }
      if (!found) {
        cout << "Sibling Vertex " << vert_id_sibling << " not found in route " << route << endl;
        printRoutes();
        throw runtime_error("Sibling vertex not found");
      }
    }
  }
}

void DARPGenome::checkNoClonesInSameRoute() const {

  for (int route=0; route<length(); route++) {
    for (int pos=0; pos<routeLength(route); pos++) {
      int vert_id = gene(route,pos);
       for (int i=pos+1; i<routeLength(route);i++) {
         if (vert_id == gene(route,i)) {
           cout << "Vertex " << vert_id << " clone found in route=" << route << endl;
           printRoutes();
           assert(false);
         }
       }
    }
  }

}

void DARPGenome::checkNoClonesInAllRoutes()  const {

  for (int route=0; route<length(); route++) {
    for (int pos=0; pos<routeLength(route); pos++) {
      int vert_id = gene(route,pos);

       for (int route2=route+1; route2<length(); route2++) {
         for (int pos2=0; pos2<routeLength(route2); pos2++) {
           if (gene(route2,pos2) == vert_id) {
             cout << "Vertex: " << vert_id << " found on routes: " << route << " and " << route2 << endl;
             printRoutes();
             assert(false);
           }
         }
       }
    }
  }


}

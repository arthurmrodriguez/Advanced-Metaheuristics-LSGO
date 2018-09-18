#include "DARPGenome.h"
#include <genomes/GAGenome.h>
#include "darpOps.h"
#include "aux.h"
#include <list>
#include <set>
#include <algorithm>
#include <limits>

const double MAXSEQSIZERATIO = 0.7;

DARPEvaluator* DARPVNSOp::darpeval__ = 0;

void DARPVNSShaker::getRandomNonEmptyRoutes(const DARPGenome& gen, int& route1, int& route2) {
  route1 = gen.randomNonEmptyRoutePos(); // First, two routes are chosen randomly (they must be different)
  do {
    route2 = gen.randomNonEmptyRoutePos();
  } while (route2 == route1);
  assert(route2 != route1);
}

void DARPVNSShaker::getRandomBadGoodRoutes(const DARPGenome& gen, /*out*/ int& badroute, /*out*/ int& goodroute) {
  vector<RouteAndScore> routescores = getRouteAndScoresSorted(gen);
  vector<double>        probs       = getSelectionProbs(routescores);

  int b_val1 = selectPos(probs);
  int b_val2 = selectPos(probs);
  // Note that since GAMin is a macro GAMin( selectPos(...), selectPos(...) ) would call selectPos 3 times
  badroute = routescores[GAMin(b_val1,b_val2)].route;

  do {
    int g_val1 = gen.size()-1-selectPos(probs);
    int g_val2 = gen.size()-1-selectPos(probs);
    goodroute = routescores[ GAMax(g_val1,g_val2)].route;
  } while (goodroute==badroute);
}

void DARPVNSShaker::getAllNaturalSeqs(const DARPGenome& gen, int route, vector< list<int>* >& natseqs) {
  set<int> processed_values;
  int      load = 0;
  int      vert_id;
  int      start_pos = 0;

  for (int i=0; i<gen.routeLength(route); i++) {
    vert_id        = gen.gene(route,i);
    Vertex& vertex = darpeval__->vertList().getVertex( vert_id );

    if ( vertex.isPickUp() ) {
      assert(processed_values.find(vert_id) == processed_values.end());
      processed_values.insert(vert_id);

      load += vertex.load_;
    }
    else { //delivery
      assert(processed_values.find( Vertex::getSiblingVertId(vert_id) ) != processed_values.end() );
      load -= vertex.load_;

      assert(load >=0);
      if (load == 0) {
        // int prev_vert_id = gen.gene(route,i-1);
        // Previously, we did not accept nat sequences of two vertices !vertex.isSiblingOf(prev_vert_id) 
        if ( !(start_pos == 0 && i==gen.routeLength(route)-1) ) {
          natseqs.push_back( getVerticesIdsList(gen,route,start_pos,i) );
        }
        start_pos = i+1;
      }
    }

  }
}

list<int>* DARPVNSShaker::getVerticesIdsList(const DARPGenome& gen, int route, int start_pos, int end_pos) {
  assert(start_pos >= 0      && start_pos < gen.routeLength(route));
  assert(end_pos > start_pos && end_pos   < gen.routeLength(route));

  list<int>* seq = new list<int>();
  for (int i=start_pos; i<=end_pos; i++ ) {
    seq->push_back(gen.gene(route,i));
  }

  return seq;
}

/*
 * For every member in the vertices ids check that the sibling is also in the list and if not,
 * extract it from the genome and add it to the list
 */
void DARPVNSShaker::extractMissingSiblings(DARPGenome& gen, int route, list<int>& vertices_ids) {
  // this should be optimized in a future
  list<int> siblings;

  for (list<int>::iterator it=vertices_ids.begin(); it!=vertices_ids.end(); it++) {
    int vert_id         = *it;
    int vert_id_sibling = Vertex::getSiblingVertId(vert_id);

    list<int>::iterator it_found = find( vertices_ids.begin(),vertices_ids.end(), vert_id_sibling );
    if (it_found == vertices_ids.end() ) { // Sibling not found
      siblings.push_back(vert_id_sibling);
      gen.removeVertex(route,vert_id_sibling);
    }

  }

  for (list<int>::iterator it=siblings.begin(); it!=siblings.end(); it++) {
    vertices_ids.push_back(*it);
  }
}

/*
 * Randomly extracts a sequence of vertices ids from the genome
 * Not optimized right now and a list<int> is returned (an extra copy)
 */
list<int> DARPVNSShaker::extractRandomRouteSeqVertices(DARPGenome& gen, int route, int neighborhood_size) {
  assert(neighborhood_size>0);
  int route_length = gen.routeLength(route);
  int vertex       = GARandomInt(0,route_length-1);
  int max_length   = GAMin(neighborhood_size,(int)(route_length*MAXSEQSIZERATIO));
  int length       = GARandomInt(1,max_length);

  list<int> vertices_ids = gen.removeVerticesFrom(route,vertex,length);

  // If the vertices ids list does not contain both the pickup and delivery vertices, the
  // sibling vertex is also inserted.
  extractMissingSiblings(gen,route,vertices_ids);

  return vertices_ids;
}

#define BESTSEQTOBE(option)\
int DARPVNSShaker::bestSeqPosToBe##option##ed(const DARPGenome& gen, int route, const vector< list<int>* >& allseqs) { \
  int best_pos, worst_pos; \
  bestAndWorstSeqToBe##option##ed(gen,route,allseqs,best_pos,worst_pos); \
\
  return best_pos;\
}

BESTSEQTOBE(Remov)
BESTSEQTOBE(Insert)

#define WORSTSEQTOBE(option)\
int DARPVNSShaker::worstSeqPosToBe##option##ed(const DARPGenome& gen, int route, const vector< list<int>* >& allseqs) { \
  int best_pos, worst_pos;\
  bestAndWorstSeqToBe##option##ed(gen,route,allseqs,best_pos,worst_pos);\
\
  return worst_pos;\
}

WORSTSEQTOBE(Remov)
WORSTSEQTOBE(Insert)


#define BESTANDWORSTSEQTOBE(option) \
void DARPVNSShaker::bestAndWorstSeqToBe##option##ed(const DARPGenome& gen, int route, const vector< list<int>* >& allseqs, /*out*/ int &best_pos, /*out*/ int &worst_pos) { \
  best_pos   = 0; \
  worst_pos  = 0; \
\
  if (allseqs.size() == 1) return;\
\
  double best_score  = gen.scoreOf##option##ingVertices(route,*(allseqs[best_pos])); \
  double worst_score = best_score; \
  double tmp_score; \
\
  for (int i=1; i<allseqs.size(); i++) { \
    tmp_score = gen.scoreOf##option##ingVertices(route,*(allseqs[i])); \
\
    if (GAGenome::compareScores(tmp_score, best_score) == GAGenome::BETTER) { \
      best_score = tmp_score; \
      best_pos   = i; \
    } \
    if (GAGenome::compareScores(tmp_score, worst_score) == GAGenome::WORSE) { \
      worst_score = tmp_score; \
      worst_pos   = i; \
    } \
  } \
  while (worst_pos == best_pos ) worst_pos = GARandomInt(0,allseqs.size()-1);\
  assert(best_pos != worst_pos);\
\
}

BESTANDWORSTSEQTOBE(Remov)
BESTANDWORSTSEQTOBE(Insert)

/*
 * Generates all the combinations of the vertices of a specific size obtained randomly between 1 and max_size
 */
vector< list<int>* >* DARPVNSShaker::getAllSequencesOfSize(const DARPGenome& gen, int route, int max_size, bool userandomsize) {
  vector<int> allpickupvertices;
  for (int i=0; i<gen.routeLength(route); i++) {
    int vert_id = gen.gene(route,i);
    // the delivery vertices ids are not
    // inserted since they are generated at the end
    if ( Vertex::isVertIdPickup(vert_id) ) {
      allpickupvertices.push_back(vert_id);
    }
  }

  max_size = GAMin(allpickupvertices.size(),(int) round((double)max_size/2.0));  // /2.0 since later we complete it with the siblings
  assert( max_size<=gen.routeLength(route) );

  vector< list<int>* >* allseqs = new vector < list<int>* >();

  list<int> tmp_sol;
  int       size;

  if (userandomsize) size = GARandomInt(1,max_size);
  else               size = max_size;

  getAllCombinations(allpickupvertices,size,tmp_sol,0, *allseqs);

  // Finally, for each sequence we add the corresponding delivery id
  for (int i=0; i<allseqs->size(); i++) {
    list<int>* seq = (*allseqs)[i];
    list<int> tmpseq(*seq);
    for(list<int>::iterator it=tmpseq.begin(); it!=tmpseq.end(); it++) {
      seq->push_back( Vertex::getSiblingVertId(*it)); //Add the delivery vertex id;
    }
  }

  return allseqs;
}

int DARPVNSShaker::selectPos(vector<double>& probs) {
  double value = GARandomDouble();

  int    i;
  int    lower = 0;
  int    upper = probs.size()-1;

  while(upper >= lower) {
    i = (lower + upper)/2;
    assert(i >= 0 && i < probs.size());

    if  (probs[i] > value) upper = i-1;
    else                   lower = i+1;
  }
  lower = GAMin((int)probs.size() - 1, lower);
  lower = GAMax(0,lower);
  assert(lower >=0 && lower < probs.size());

  return lower;
}

string DARPVNSShaker::stringOfSeq(list<int>& seq) {
  stringstream msg;
  for (list<int>::iterator  it=seq.begin(); it!=seq.end(); it++) {
    msg << *it << " ";
  }
  msg << endl;

  return msg.str();
}

/*
 * Returns all the possible sequences of ids of maximum size max_size that could be created with the vertices ids
 * contained in the genome.
 * Note: only sequences with both pickup and delivery ids are allowed!!
 * Therefore, we first generate all the sequences with pickup ids and then we complete them with the delivery ids
 *
 */
//vector< list<int>* > getAllSequencesUpToMaxSize(DARPGenome& gen, int route, int max_size) {
//  vector< list<int>* > allseqs;
//
//  // The allseqs vector is initialized with all the possible sequences of size 1
//  for (int i=0; i<gen.routeLength(route); i++) {
//    int vert_id = gen.gene(route,i);
//    // As previously mentioned, the delivery vertices ids are not
//    // inserted since they are generated at the end
//    if ( Vertex::isVertIdPickup(vert_id) ) {
//      list<int> *seq = new list<int>();
//      seq->push_back(gen.gene(route,i));
//      allseqs.push_back(seq);
//    }
//  }
//
//  max_size = GAMin(allseqs.size(),round((double)max_size/2.0));  // /2.0 since later we complete it with the siblings
//
//  // Then the sequences of greater size are generated by adding each vertex id to the previously created sequences
//  int start_pos = 0;
//
//  while ( allseqs[allseqs.size()-1]->size() < max_size) { // While we have not reached the maximum sequence size
//
//    int last_pos = allseqs.size()-1;                      // Needs to be computed before so that it does not change
//                                                          // during the calls of the following for
//    for(int seqpos=start_pos; seqpos<=last_pos; seqpos++) {
//
//      list<int>* seq    = allseqs[seqpos];
//
//      for (int genpos=0; genpos<gen.routeLength(route); genpos++) {
//        int vertex_id = gen.gene(route,genpos);
//
//        if (Vertex::isVertIdPickup(vertex_id) &&
//            find(seq->begin(),seq->end(),vertex_id) == seq->end() ) { // Not found
//          list<int>* newseq = new list<int>(*seq);
//          newseq->push_back(vertex_id);
//          allseqs.push_back(newseq);
//        }
//      }
//
//    }
//    start_pos = last_pos+1;
//  }
//
//  // Finally, for each sequence we add the corresponding delivery id
//  for (int i=0; i<allseqs.size(); i++) {
//    list<int>* seq = allseqs[i];
//    list<int> tmpseq(*seq);
//    for(list<int>::iterator it=tmpseq.begin(); it!=tmpseq.end(); it++) {
//      seq->push_back( Vertex::getSiblingVertId(*it)); //Add the delivery vertex id;
//    }
//  }
//
//
//  /*
//   * Don't understand why this code does not work
//   */
////  for (int i=0; i<allseqs.size(); i++) {
////    list<int>* seq = allseqs[i];
////    cout << endl;
////    for(list<int>::reverse_iterator rit=seq->rbegin(); rit!=seq->rend(); rit++) {
////      seq->push_back( Request::getSiblingReqId(*rit)); //Add the delivery vertex id;
////    }
////    cout << endl;
////  }
//
//  return allseqs;
//}


// Just for testing purposes
bool eachRouteHasACompleteRequest(DARPGenome& gen) {
  for (int route=0; route<gen.length(); route++) {
    for (int pos=0; pos<gen.routeLength(route); pos++) {
      int vert_id = gen.gene(route,pos);
      if (gen.findPosOfVertex(route,Vertex::getSiblingVertId(vert_id)) == -1) return false;
    }
  }
  return true;
}

// Just for testing purposes
bool eachCritcVertexComesBefore(DARPGenome& gen) {
  for (int route=0; route<gen.length(); route++) {
    for (int pos=0; pos<gen.routeLength(route); pos++) {
      int vert_id = gen.gene(route,pos);
      if (Vertex::isVertIdPickup(vert_id)) {
        if (gen.findPosOfVertex(route,Vertex::getSiblingVertId(vert_id)) <= pos) return false;
      }
    }
  }
  return true;

}


vector<DARPVNSShaker::RouteAndScore> DARPVNSShaker::getRouteAndScoresSorted(const DARPGenome& gen, GAGenome::OptCriterion optcrit) {
  gen.evaluate(); // Compute the score (if not already)

  vector<DARPVNSShaker::RouteAndScore> routescores (gen.size());
  for (int route=0; route<gen.size(); route++) {
    routescores[route].route = route;
    routescores[route].score  = darpeval__->nonPenalizedScore(gen,route);
  }

  // Sort so that the best routes are at the first positions
  std::sort(routescores.begin(), routescores.end(), descendingRouteAndScoreOrderFunc );

  if (optcrit==GAGenome::MINIMIZATION) std::reverse(routescores.begin(), routescores.end());

  return routescores;
}

// Descending order: The ones with the higher score will be at the first positions
bool DARPVNSShaker::descendingRouteAndScoreOrderFunc(const GreedyMoveNeighborhood::RouteAndScore& a, const GreedyMoveNeighborhood::RouteAndScore& b) {
  return a.score > b.score;
}

vector<double> DARPVNSShaker::getSelectionProbs (const vector<RouteAndScore>& routescores, GAGenome::OptCriterion optcrit) {
  vector<double> probs_sum (routescores.size());

  // TODO: This code is almost the same as the one from GASelectionScheme. Both codes should be refactored
  if (optcrit == GAGenome::MAXIMIZATION) {
    probs_sum[0] = routescores[0].score;
    for (int i=1; i<probs_sum.size(); i++) probs_sum[i] = routescores[i].score + probs_sum[i-1];
    for (int i=0; i<probs_sum.size(); i++) probs_sum[i] /= probs_sum[ probs_sum.size() - 1 ];
  }
  else { // Minimization
    double minscore = routescores[0].score;
    double maxscore = minscore;

    for (vector<RouteAndScore>::const_iterator it=routescores.begin(); it!=routescores.end(); it++) {
      if      (it->score < minscore) minscore = it->score;
      else if (it->score > maxscore) maxscore = it->score;
    }

    probs_sum[0] = maxscore + minscore - routescores[0].score;
    for (int i=1; i<probs_sum.size(); i++) probs_sum[i] = maxscore + minscore - routescores[i].score + probs_sum[i-1];
    for (int i=0; i<probs_sum.size(); i++) probs_sum[i] /= probs_sum[ probs_sum.size() - 1 ];

  }

  return probs_sum;
}


list<int> DARPVNSShaker::bestSeqToBeRemoved (const DARPGenome& gen, int route) {
  vector< list<int>* >* allseqs      = getAllSequencesOfSize(gen,route,size_, false);
  int                   best_seq_pos = bestSeqPosToBeRemoved(gen,route,*allseqs);

  list<int> best_seq = *( (*allseqs)[best_seq_pos] );

  for(int i=0; i<allseqs->size(); i++) delete (*allseqs)[i];
  delete allseqs;

  return best_seq;
}

/*
 * Factory method for creating the shakers given a name and its size
 */
DARPVNSShaker* DARPVNSShaker::createNewShaker(string name, int size) {
  DARPVNSShaker* shaker = 0;

#define SHAKERCOMPCASE(shakername) if (name.compare(#shakername) == 0) shaker = new shakername(size);

  SHAKERCOMPCASE(SwapNeighborhood);
  SHAKERCOMPCASE(ChainNeighborhood);
  SHAKERCOMPCASE(GreedyMoveNeighborhood);
  SHAKERCOMPCASE(GreedyMoveNeighborhoodDestCentered);
  SHAKERCOMPCASE(GreedySwapNeighborhood);
  SHAKERCOMPCASE(GreedySwapNeighborhoodDestCentered);
  SHAKERCOMPCASE(ZeroSplitNeighborhood);
  SHAKERCOMPCASE(CheckAllNaturalSeqsCombsNeighborhood);

  if (shaker == 0) throw runtime_error("Unrecognized shaker name");

  return shaker;
}



// Shakers

VNSOpAction* SwapNeighborhood::operator() (GAGenome& g) {
  DARPGenome& gen = dynamic_cast<DARPGenome&>(g);

  if (! gen.moreThanOneRoute() ) return  new NoVNSOpAction();

  int route1, route2;
  getRandomNonEmptyRoutes(gen, route1, route2);

  // Then, on each route, a sequence to be swapped is extracted according to the option
  list<int> route1_seq = extractRandomRouteSeqVertices(gen,route1,size_);
  list<int> route2_seq = extractRandomRouteSeqVertices(gen,route2,size_);

  // And inserted one-by-one into the other route
  gen.insertVerticesInBestPos(route2, route1_seq);

  gen.insertVerticesInBestPos(route1, route2_seq);

  //assert(eachRouteHasACompleteRequest(gen));
  //assert(eachCritcVertexComesBefore(gen));
  return new SwapAction(route1,route1_seq,route2,route2_seq);
}

VNSOpAction* ChainNeighborhood::operator() (GAGenome& g) {
  DARPGenome& gen = dynamic_cast<DARPGenome&>(g);

  if (! gen.moreThanOneRoute() ) return  new NoVNSOpAction();

  int route1 = gen.randomNonEmptyRoutePos();
  int route2;
  do {
    route2 = gen.randomRoutePos();
  } while (route2 == route1);
  assert(route2 != route1);

  // Same as with the swap shaker: with extract a random seq of route1 and add it to route2
  list<int> route1_seq = extractRandomRouteSeqVertices(gen,route1,size_);

  gen.insertVerticesInBestPos(route2, route1_seq);

  MoveAction* act = new MoveAction(route1,route1_seq,route2);

  int       orig_route  = route1;
  int       selec_route = route2;
  list<int> moved_seq   = route1_seq;

  for (int mov=1; mov<size_; mov++) { // starts at 1 since a move has been already carried out

    vector< list<int>* >* allseqs = getAllSequencesOfSize(gen,selec_route,size_);

    int        best_pos  = bestSeqPosToBeRemoved(gen,selec_route,*allseqs);
    list<int>& seq_2move = *((*allseqs)[best_pos]);

    int new_route;
    do {
      new_route = gen.randomRoutePos();
    } while (  new_route == selec_route or                                    // are the same routes?
             ( new_route == orig_route and areReqsEqual(moved_seq,seq_2move)  // are we undoing the previous movement?
               and gen.size() > 2) );                                         // to avoid endless loops with two routes
    assert(selec_route != new_route);

    gen.moveVertices (selec_route, new_route, seq_2move);

    act->addDatum(selec_route,seq_2move,new_route);

    orig_route  = selec_route;
    selec_route = new_route;
    moved_seq   = seq_2move;

    for(int i=0; i<allseqs->size(); i++) delete (*allseqs)[i];

    delete allseqs;

  }
  //assert(eachRouteHasACompleteRequest(gen));
  //assert(eachCritcVertexComesBefore(gen));

  return act;
}

VNSOpAction* GreedyMoveNeighborhood::operator() (GAGenome& g) {

  DARPGenome& gen = dynamic_cast<DARPGenome&>(g);

  int orig_route = routeForExtraction(gen);

  list<int> best_seq = bestSeqToBeRemoved(gen,orig_route);

  double orig_score = gen.score();

  gen.removeVertices(orig_route,best_seq); 

  int best_ins_route = bestRoutePosForInsertion(gen,orig_route,orig_score,best_seq);
  assert(best_ins_route >= 0);

  gen.insertVerticesInBestPos(best_ins_route,best_seq);

  return new MoveAction(orig_route,best_seq,best_ins_route);
}

int GreedyMoveNeighborhood::routeForExtraction(const DARPGenome& gen) {
  vector<RouteAndScore> routescores = getRouteAndScoresSorted(gen);
  vector<double>        probs       = getSelectionProbs(routescores);

  int route1 = selectPos(probs);
  int route2 = selectPos(probs);

  int route_sel = GAMin(route1, route2);

  return routescores[route_sel].route;
}

/*
 * Orig route is the route from where the vertices sequence was extracted
 */
int GreedyMoveNeighborhood::bestRoutePosForInsertion(const DARPGenome& gen, int orig_route, double orig_score, const list<int>& vertices) {
  // Extra greedy approach: we try all the possible insertions
  int    best_ins_route = -1;
  double best_score     = GAGenome::worstPossibleScore();
  double tmp_score;

  for (int route=0; route<gen.size(); route++) {
    if (route==orig_route) continue;

    tmp_score = gen.scoreOfInsertingVertices( route, vertices );
    if ( GAGenome::compareScores(tmp_score,best_score) == GAGenome::BETTER ) {

      best_score     = tmp_score;
      best_ins_route = route;

      // We stop when we found an insertion point that improves the result
      if ( GAGenome::compareScores(tmp_score,orig_score) == GAGenome::BETTER) break;
    }
  }

  return best_ins_route;
}

VNSOpAction* GreedyMoveNeighborhoodDestCentered::operator() (GAGenome& g) {

  DARPGenome& gen = dynamic_cast<DARPGenome&>(g);

  int badroute=-1, goodroute=-1;
  getRandomBadGoodRoutes(gen,badroute,goodroute);

  list<int> bestseq = bestSeqToBeRemovedAndInserted(gen,badroute,goodroute);

  double orig_score = gen.score();

  gen.removeVertices(badroute,bestseq);

  gen.insertVerticesInBestPos(goodroute,bestseq);

  return new MoveAction(badroute,bestseq,goodroute);
}

list<int> GreedyMoveNeighborhoodDestCentered::bestSeqToBeRemovedAndInserted(const DARPGenome& gen, int fromroute, int toroute) {
  vector< list<int>* >* allseqs      = getAllSequencesOfSize(gen,fromroute,size_, false);
  int                   best_seq_pos = bestSeqPosToBeInserted(gen,toroute,*allseqs);

  list<int> best_seq = *( (*allseqs)[best_seq_pos] );

  for(int i=0; i<allseqs->size(); i++) delete (*allseqs)[i];
  delete allseqs;

  return best_seq;
}

string GreedySwapNeighborhood::getStatsResults() {
  stringstream msg;
  msg << VNSOp::getStatsResults();
  msg << " BestOpt results: bestbset=" << bestswaps_num_[bestbest] << " bestworst=" << bestswaps_num_[bestworst];
  msg << " worstbest=" << bestswaps_num_[worstbest] << " worstworst=" << bestswaps_num_[worstworst] << endl;

  return msg.str();
}

VNSOpAction* GreedySwapNeighborhood::operator() (GAGenome& g) {
  DARPGenome& gen = dynamic_cast<DARPGenome&>(g);

  if (! gen.moreThanOneRoute() ) return  new NoVNSOpAction();

  int route1=-1, route2=-1;
  getRandomBadGoodRoutes(gen,route1,route2);
  //getRandomNonEmptyRoutes(gen, route1, route2); Old way for selecting the routes

  // We need at least two sequences for each route in order to execute the shaker
  if (gen.routeLength(route1) / size_ < 2 and gen.routeLength(route2) < 2 ) return new NoVNSOpAction();

  // Then, on each route, the best and worst sequence to be swapped is extracted
  list<int> best_seq1; list<int> worst_seq1;
  list<int> best_seq2; list<int> worst_seq2;

  getBestAndWorstSeqsToRemove(gen,route1,route2,best_seq1,worst_seq1);
  getBestAndWorstSeqsToRemove(gen,route2,route1,best_seq2,worst_seq2);

  list<int> sel_seq1, sel_seq2;
  getBestSwap(gen,route1,best_seq1,worst_seq1,route2,best_seq2,worst_seq2,sel_seq1,sel_seq2);

  gen.swapSeqs(route1,sel_seq1,route2,sel_seq2);

  return new SwapAction(route1,sel_seq1,route2,sel_seq2);
}

void GreedySwapNeighborhood::getBestSwap(const DARPGenome& gen,
                                         int               route1, list<int>& best_seq1, list<int>& worst_seq1,
                                         int               route2, list<int>& best_seq2, list<int>& worst_seq2,
                                        /*out*/ list<int>& sel_seq1, /*out*/ list<int>& sel_seq2 ) {

  double best_score = GAGenome::worstPossibleScore();
  double tmp_score;
  int    bestopt = -1;

#define CHECKSCORE(gen,route1,seq1,route2, seq2, opt1, opt2) \
  tmp_score = scoreOfSwappingSeqs(gen,route1,opt1##_seq1,route2, opt2##_seq2); \
  if ( GAGenome::compareScores(tmp_score, best_score ) == GAGenome::BETTER)  { \
    best_score = tmp_score; \
    sel_seq1   = opt1##_seq1; \
    sel_seq2   = opt2##_seq2; \
    bestopt    = opt1##opt2;\
  }

  CHECKSCORE(gen,route1,seq1,route2,seq2,best, best)
  CHECKSCORE(gen,route1,seq1,route2,seq2,best, worst)
  CHECKSCORE(gen,route1,seq1,route2,seq2,worst,best)
  CHECKSCORE(gen,route1,seq1,route2,seq2,worst,worst)

  assert(bestopt >= 0);
  bestswaps_num_[bestopt]++;
}

void GreedySwapNeighborhood::getBestAndWorstSeqsToRemove(const DARPGenome&  gen,
                                                         int                fromroute,
                                                         int                toroute,
                                                         /*out*/ list<int>& best_seq,
                                                         /*out*/ list<int>& worst_seq) {
  vector< list<int>* >* allseqs = getAllSequencesOfSize(gen,fromroute,size_, false);

  int best_seq_pos = -1, worst_seq_pos = -1;
  bestAndWorstSeqToBeRemoved(gen,fromroute,*allseqs,best_seq_pos,worst_seq_pos);

  best_seq  = * (*allseqs)[best_seq_pos];
  worst_seq = * (*allseqs)[worst_seq_pos];

  for(int i=0; i<allseqs->size(); i++) delete (*allseqs)[i];
  delete allseqs;
}

double GreedySwapNeighborhood::scoreOfSwappingSeqs(const DARPGenome& gen, int route1, list<int>& seq1, int route2, list<int>& seq2) {
  DARPGenome tmp_gen(gen);

  tmp_gen.swapSeqs(route1,seq1,route2,seq2);

  double score = tmp_gen.score();
  const_cast<DARPGenome&>(gen).nevals( gen.nevals() + tmp_gen.nevals() ); // Update the nevals

  return score;
}


void GreedySwapNeighborhoodDestCentered::getBestAndWorstSeqsToRemove(const DARPGenome&  gen,
                                                                     int                fromroute,
                                                                     int                toroute,
                                                                     /*out*/ list<int>& best_seq,
                                                                     /*out*/ list<int>& worst_seq) {
  vector< list<int>* >* allseqs = getAllSequencesOfSize(gen,fromroute,size_, false);

  int best_seq_pos = -1, worst_seq_pos = -1;
  bestAndWorstSeqToBeInserted(gen,toroute,*allseqs,best_seq_pos,worst_seq_pos);
  assert(best_seq_pos !=-1 && worst_seq_pos !=-1);
  assert(best_seq_pos != worst_seq_pos);

  best_seq  = * (*allseqs)[best_seq_pos];
  worst_seq = * (*allseqs)[worst_seq_pos];

  for(int i=0; i<allseqs->size(); i++) delete (*allseqs)[i];
  delete allseqs;
}

VNSOpAction* ZeroSplitNeighborhood::operator() (GAGenome& g) {
  DARPGenome& gen = dynamic_cast<DARPGenome&>(g);
  bool anynatseq = false;
  vector< vector< list<int>* > > allnatseqs(gen.size());

  for (int i=0; i<gen.length(); i++) {
    getAllNaturalSeqs(gen,i, allnatseqs[i]);
    if (!allnatseqs[i].empty()) anynatseq = true;
  }

  VNSOpAction* act;

  if (anynatseq) { // If no natural sequences are found we don't do anything

    int selecroute;
    do {           // We select first a route that contains natural sequences
      selecroute = GARandomInt(0,allnatseqs.size()-1);
    } while (allnatseqs[selecroute].empty());

    int toroute;
    do {
      toroute = gen.emptyRouteOrRandomRoutePos();
    } while (toroute == selecroute);
    assert(selecroute != toroute);

    int natseqpos = GARandomInt(0,allnatseqs[selecroute].size()-1);

    list<int>& sel_seq = *(allnatseqs[selecroute][natseqpos]);
    gen.moveVertices (selecroute, toroute, sel_seq);

    act = new MoveAction(selecroute,sel_seq,toroute);
  }
  else act = new NoVNSOpAction();

  for (int i=0; i<gen.length(); i++) {
    for (int j=0; j<allnatseqs[i].size(); j++) {
      delete allnatseqs[i][j];
    }
  }

  //assert(eachRouteHasACompleteRequest(gen));
  //assert(eachCritcVertexComesBefore(gen));

  return act;
}

int CheckAllNaturalSeqsCombsNeighborhood::computeInsPosBasedOnTheNatSeqPos(const vector< list<int>* >& route_natseqs, int natseq_pos) {
  if (natseq_pos <0 || natseq_pos >= route_natseqs.size()) return -1;

  int natseq_starpos = 0;
  for (int i=0; i<natseq_pos; i++) natseq_starpos += route_natseqs[i]->size();

  return natseq_starpos;
}

SwapWithInsertionPosDatum
CheckAllNaturalSeqsCombsNeighborhood::swapNaturalSequences(DARPGenome&                     gen,
                                                           vector< vector< list<int>* > >& allnatseqs,
                                                           int                             fromroute,
                                                           int                             from_natseq_pos,
                                                           int                             toroute,
                                                           int                             to_natseq_pos) {
  assert(fromroute != toroute);
  assert(fromroute >= 0 && fromroute < allnatseqs.size());
  assert(toroute >= 0 && toroute < allnatseqs.size());
  assert(from_natseq_pos >= 0 && to_natseq_pos >= 0);
  
  vector< list<int>* > fromroute_natseqs = allnatseqs[fromroute];
  vector< list<int>* > toroute_natseqs   = allnatseqs[toroute];

  list<int>& fromroute_natseq = * fromroute_natseqs[from_natseq_pos];
  list<int>& toroute_natseq   = * toroute_natseqs[to_natseq_pos];

  int fromroute_ins_pos = computeInsPosBasedOnTheNatSeqPos(fromroute_natseqs, from_natseq_pos);
  int toroute_ins_pos   = computeInsPosBasedOnTheNatSeqPos(toroute_natseqs,   to_natseq_pos);
  assert(fromroute_ins_pos != -1 && toroute_ins_pos != -1);

  gen.swapWithInsertPosSeqs(fromroute, fromroute_ins_pos,fromroute_natseq,
                            toroute,   toroute_ins_pos,  toroute_natseq);

  // Finally, we swap the natural sequences
  list<int>* tmp_natseqs                 = allnatseqs[fromroute][from_natseq_pos];
  allnatseqs[fromroute][from_natseq_pos] = allnatseqs[toroute][to_natseq_pos];
  allnatseqs[toroute][to_natseq_pos]     = tmp_natseqs;

  return SwapWithInsertionPosDatum(fromroute,fromroute_ins_pos,fromroute_natseq,
                                   toroute,  toroute_ins_pos,  toroute_natseq);
}



void CheckAllNaturalSeqsCombsNeighborhood::getPrevAndNextNatSeqs(vector< vector< list<int>* > >& allnatseqs, 
                                                                 int                             route, 
                                                                 int                             pos,
                                                                 list<int>*&                     prev_natseq,
                                                                 list<int>*&                     next_natseq) {

  const vector< list<int>* >& route_natseqs = allnatseqs[route];
  assert(pos >=0 && pos < route_natseqs.size() );

  prev_natseq = next_natseq = 0;

  if (pos > 0)                      prev_natseq = route_natseqs[pos-1];
  if (pos < route_natseqs.size()-1) next_natseq = route_natseqs[pos+1];
}

bool CheckAllNaturalSeqsCombsNeighborhood::canNatSeqReplaceNatSeq(vector< vector< list<int>* > >& allnatseqs,
                                                                  int                             fromroute,
                                                                  int                             from_natseq_pos,
                                                                  int                             toroute,
                                                                  int                             to_natseq_pos) {

  assert(fromroute >=0 && fromroute < allnatseqs.size());
  assert(toroute >=0 && toroute < allnatseqs.size());

  const CostMatrix& costmatrix = darpeval__->costMatrix();

  vector< list<int>* >& fromroute_natseqs = allnatseqs[fromroute];
  list<int>*            from_natseq       = fromroute_natseqs[from_natseq_pos];

  list<int>* to_prev_natseq;
  list<int>* to_next_natseq;

  getPrevAndNextNatSeqs(allnatseqs,toroute,to_natseq_pos, to_prev_natseq, to_next_natseq);

  if (to_prev_natseq) {
    int     from_first_vertex_id = from_natseq->front();
    Vertex& from_first_vertex    = darpeval__->vertList().getVertex( from_first_vertex_id );

    int     to_prev_last_vertex_id = to_prev_natseq->back();
    Vertex& to_prev_last_vertex    = darpeval__->vertList().getVertex( to_prev_last_vertex_id );

    if ( !to_prev_last_vertex.canVertexGoBefore(from_first_vertex,costmatrix) ) return false;
  }
  if (to_next_natseq) {
    int     from_last_vertex_id = from_natseq->back();
    Vertex& from_last_vertex    = darpeval__->vertList().getVertex( from_last_vertex_id );

    int     to_next_first_vertex_id = to_next_natseq->front();
    Vertex& to_next_first_vertex    = darpeval__->vertList().getVertex( to_next_first_vertex_id );

    if ( !from_last_vertex.canVertexGoBefore(to_next_first_vertex,costmatrix) ) return false;
  }

  return true;
}

bool CheckAllNaturalSeqsCombsNeighborhood::isSwapFeasible(vector< vector< list<int>* > >& allnatseqs,
                                                          int                             fromroute,
                                                          int                             from_natseq_pos,
                                                          int                             toroute,
                                                          int                             to_natseq_pos) {

  return canNatSeqReplaceNatSeq(allnatseqs,fromroute,from_natseq_pos,toroute,  to_natseq_pos) && 
         canNatSeqReplaceNatSeq(allnatseqs,toroute,  to_natseq_pos,  fromroute,from_natseq_pos);
}


double CheckAllNaturalSeqsCombsNeighborhood::scoreFromSwapNaturalSequences(const DARPGenome&               gen,
                                                                           vector< vector< list<int>* > >& allnatseqs,
                                                                           int                             fromroute,
                                                                           int                             from_natseq_pos,
                                                                           int                             toroute,
                                                                           int                             to_natseq_pos) {
  double score;

  if ( !isSwapFeasible(allnatseqs,fromroute,from_natseq_pos,toroute,to_natseq_pos) ) {
    score = GAGenome::worstPossibleScore();
  }
  else {
    DARPGenome tmp_gen (gen);
    swapNaturalSequences(tmp_gen,allnatseqs,fromroute,from_natseq_pos,toroute,to_natseq_pos);

    score = tmp_gen.score();
    int old_evals = gen.nevals();
    const_cast<DARPGenome&>(gen).nevals(gen.nevals()+tmp_gen.nevals());  // Need to update the nevals
    assert(gen.nevals() == old_evals+ 1);
  }

  return score;
}

void CheckAllNaturalSeqsCombsNeighborhood::bestSwapOfNatSeqs(const DARPGenome&                     gen,
                                                             const vector< vector< list<int>* > >& allnatseqs,
                                                                   int                             fromroute,
                                                                   int                             toroute,
                                                                   int&    /*out*/                 bestfrom_pos,
                                                                   int&    /*out*/                 bestto_pos,
                                                                   double& /*out*/                 swap_score) {
  
  bestfrom_pos = bestto_pos = -1;
  double tmp_score;
  swap_score = gen.score();

  // habria que controlar que no se haga el swap si no es factible
  for (int from_natseq_pos=0; from_natseq_pos<allnatseqs[fromroute].size(); from_natseq_pos++) {
    for (int to_natseq_pos=0; to_natseq_pos<allnatseqs[toroute].size(); to_natseq_pos++) {

      vector< vector< list<int>* > > tmp_natseqs (allnatseqs);
      tmp_score = scoreFromSwapNaturalSequences(gen,tmp_natseqs,fromroute,from_natseq_pos,toroute,to_natseq_pos);

      if (GAGenome::compareScores(tmp_score, swap_score) == GAGenome::BETTER ) {
        bestfrom_pos = from_natseq_pos;
        bestto_pos   = to_natseq_pos;
        swap_score   = tmp_score;
      }
    }
  }
}

VNSOpAction* CheckAllNaturalSeqsCombsNeighborhood::operator() (GAGenome& g) {
  DARPGenome& gen = dynamic_cast<DARPGenome&>(g);

  // For avoiding useless computation time we only apply this shaker if the score has changed
  if (gen.score() == prev_score_) return new NoVNSOpAction();

  bool anynatseq = false;
  vector< vector< list<int>* > > allnatseqs(gen.size());

  for (int i=0; i<gen.length(); i++) {
    getAllNaturalSeqs(gen,i, allnatseqs[i]);
    if (!allnatseqs[i].empty()) anynatseq = true;
  }

  // TODO: Se puede aplicar la idea de que mientras haya una mejora siga aplicando el shaker
  // y quizas hacerlo en funcion del size

  // If no natural sequence is found we return with no action
  if (!anynatseq) return new NoVNSOpAction();

  SwapWithInsertionPosAction* act = new SwapWithInsertionPosAction();

  for (int fromroute=0; fromroute<allnatseqs.size()-1; fromroute++) {
    if (allnatseqs[fromroute].empty()) continue;

    int from_natseq_pos,  from_bestnatseq_pos;
    int to_natseq_pos,    to_bestnatseq_pos, best_toroute;

    double best_swap_score = gen.score();
    double swap_score;
    for (int toroute=fromroute+1; toroute<allnatseqs.size(); toroute++) {
      if (allnatseqs[toroute].empty()) continue;
      // Give me the best swap available for route "fromroute" along with the values of the swap
      bestSwapOfNatSeqs(gen,allnatseqs,fromroute,toroute,from_natseq_pos,to_natseq_pos,swap_score);

      if (GAGenome::compareScores(swap_score, best_swap_score) == GAGenome::BETTER ) {
        assert(from_natseq_pos >= 0 && to_natseq_pos >= 0);
        best_swap_score     = swap_score;
        from_bestnatseq_pos = from_natseq_pos;
        to_bestnatseq_pos   = to_natseq_pos;
        best_toroute        = toroute;
      }
    }
    assert( (GAGenome::compareScores(best_swap_score, gen.score()) != GAGenome::BETTER) ||
            (from_bestnatseq_pos >= 0 && to_bestnatseq_pos >= 0 && best_toroute >=0) );

    //Conduct swap
    if (GAGenome::compareScores(best_swap_score, gen.score()) == GAGenome::BETTER ) {
      SwapWithInsertionPosDatum opdatum =
          swapNaturalSequences(gen,allnatseqs,fromroute,from_bestnatseq_pos,best_toroute,to_bestnatseq_pos);
      act->addDatum(opdatum);
    }
  }

  prev_score_ = gen.score();

  for (int i=0; i<allnatseqs.size(); i++) {
    for (int j=0; j<allnatseqs[i].size(); j++) {
      delete allnatseqs[i][j];
    }
  }

  return act;
}

// Local Search Methods

void DARPLocalSearch::insertNonCriticalVertexAtBestPos(DARPGenome& gen, int route, int critvert_pos) {
  int critvert_id     = gen.gene(route,critvert_pos);
  int noncrit_vert_id = Vertex::getSiblingVertId(critvert_id);

  gen.insertVertexInBestPos(route,noncrit_vert_id);
}

/*
 * Traverses the genome from the start_pos until it find the first vertex that is critical
 * Returns -1 if no critical vertex has been found
 */
int DARPLocalSearch::firstCriticalVertexPos(DARPGenome& gen, int route, int start_pos) {
  assert(start_pos >= 0 && start_pos < gen.routeLength(route));

  for (int i=start_pos; i<gen.routeLength(route); i++) {
    Vertex& vert = darpeval__->vertList().getVertex( gen.gene(route,i) );

    if (vert.critic_) return i;
  }

  return -1;
}

void DARPLocalSearch::localSearchToCritVertexPos(DARPGenome& gen, int route, int critvert_origpos, int critvert_searchpos,
                                                 long maxevals,
                                                 /*out*/int& critvert_newpos, /*out*/bool& unfinished) {

  if (maxevals <= 0) {
    critvert_newpos = -1; // so that it get recomputed the next time is called
    unfinished      = true;
    return;
  }

  unfinished = false;

  int       orig_length = gen.routeLength(route);
  double    orig_score  = gen.score();
  GAGenome* origGenome  = gen.clone(); // so that we saved the cache with the score values computed

  int critvert_id    = gen.gene(route,critvert_origpos);
  int noncritvert_id = Vertex::getSiblingVertId(critvert_id);

  gen.removeVertexOfPos(route,critvert_origpos);                                // First, we remove the critical and non critical
  int noncritvert_origpos = gen.removeVertex(route,noncritvert_id,critvert_origpos); // vertices

  int startinsert_pos = -1, endinsert_pos = -1;
  gen.getInsertionPos(route, critvert_id, startinsert_pos, endinsert_pos); // Then we find the possible insertion positions

  if (critvert_searchpos != -1) {
    assert(critvert_searchpos >= startinsert_pos);
    startinsert_pos = critvert_searchpos;
  }

  int noncritvert_newpos;

  for (critvert_newpos=startinsert_pos; critvert_newpos<=endinsert_pos; critvert_newpos++) {
    int orig_evals = gen.nevals();

    gen.insertVertex(route,critvert_newpos,critvert_id); // First, the critical vertex

    // Then, the non critical vertex.
    Vertex& critvert = darpeval__->vertList().getVertex(critvert_id);

    // If it happens that the non critical vertex is a delivery vertex
    // we can assume that the starting pos is going to be higher than the crtical pos
    int noncrit_starting_pos = critvert.isPickUp() ? critvert_newpos+1 : 0;
    noncritvert_newpos       = gen.insertVertexInBestPosFromStartPos(route, noncritvert_id, noncrit_starting_pos, orig_score, maxevals);

    maxevals -= gen.nevals() - orig_evals;

    if (GAGenome::compareScores(gen.score(), orig_score) == GAGenome::BETTER) {
      break;
    }
    else {                                             // the new score is worse, so we remove the vertices
      gen.removeVertexOfPos(route,noncritvert_newpos); // The order of removal is crucial
      gen.removeVertexOfPos(route,critvert_newpos);
    }

    assert(maxevals >= 0);
    if (maxevals == 0) {
      unfinished = true;
      break; // Need to be here so that critvert_newpos does not get incremented in case the maxevals <= 0
    }
  }

  assert(gen.routeLength(route) == orig_length -2 || GAGenome::compareScores(gen.score(), orig_score) == GAGenome::BETTER);

  if (gen.routeLength(route) == orig_length - 2) {  // If the last position has been reached and
    int last_nevals = gen.nevals();                 // no improvement has been found, the vertices
    gen.copy(*origGenome);                          // are left at their original positions.
    gen.nevals(last_nevals);                        // Again, the insertion order is crucial

    if (!unfinished) {
      critvert_newpos = critvert_origpos; // If we stop not because of the lack of evaluation calls but because we have reached
                                          // the end, we rest the critvert_newpos so that it does the insert at the original
                                          // position (we leave it as it was)
    }
  }

  assert(gen.routeLength(route) == orig_length);

  delete origGenome;
}

// Local Search for a specific route
int DARPLocalSearch::localSearch(DARPGenome& gen, int route, int maxevals) {
  int orig_evals     = gen.nevals();
  bool ls_unfinished = false;

  // We do nothing if the route has already been explored or if we do not have any available
  // evals
  if (gen.localSearchedRoute(route) != DARPGenome::LSEXPLORED && maxevals > 0) {

    int critvert_id;
    int critvert_newpos;

    int pos = 0;
    int startcritvert_pos = -1;
    if (gen.localSearchedRoute(route) == DARPGenome::LSUNFINISHED) {
      pos               = gen.localSearchLastPosExplored(route);
      startcritvert_pos = gen.localSearchLastCritVertexSearchPosExplored(route);
      assert( Vertex::isVertIdPickup( gen.gene(route,pos)) ); // An stored pos should always point to a critical vertex
    }

    // The last pos should contain a non critical vertex
    while ( pos < gen.routeLength(route)-1 ) {
      pos = firstCriticalVertexPos(gen,route,pos);

      if (pos == -1) break; // No critical vertex has been found

      int start_nevals = gen.nevals();

      localSearchToCritVertexPos(gen, route, pos, startcritvert_pos, maxevals, critvert_newpos, ls_unfinished);
      maxevals -= gen.nevals() - start_nevals;

      if (ls_unfinished) break;

      startcritvert_pos = -1; // So that the next call does not use the cached value used for the first iteration

      pos++;

      assert(maxevals >= 0 || ls_unfinished);
      // -1 because the localSearchToCritVertex consumes one more evaluation than the ones requested

      /* In order to deal with high dimensional problems, the following greedy approach has been commented
      if (critvert_newpos < pos ) {      // If the vertex is inserted before the current pos, the process
        pos = critvert_newpos;           // restarts the process at the next position (pos will be incremented
      }                                  // at the next statement)
      else if (critvert_newpos == pos) { // If no change has happened in the position, the next vertex is analyzed
        pos = pos +1;                    // Otherwise, the same position is again analyzed since a change in the critical vertex
      }                                  // has happened
      */

    }

    if (ls_unfinished) {
      gen.localSearchedRoute(route, DARPGenome::LSUNFINISHED);
      gen.localSearchLastPosExplored(route,pos,critvert_newpos);
    }
    else {
      gen.localSearchedRoute(route, DARPGenome::LSEXPLORED);
    }

    //assert(eachRouteHasACompleteRequest(gen));
    //assert(eachCritcVertexComesBefore(gen));
  }

  return gen.nevals() - orig_evals;
}

VNSOpAction* DARPLocalSearch::operator() (GAGenome& g) {
  DARPGenome& gen = dynamic_cast<DARPGenome&>(g);

  int maxevals = maxevalspercall_;

  vector<int> routes;
  for (int i=0; i<gen.size(); i++) routes.push_back(i);

  if (randomroutes_) random_shuffle(routes.begin(),routes.end());

  for (vector<int>::iterator it=routes.begin(); it!= routes.end(); it++) {
    int consumed_evals = localSearch(gen,*it,maxevals);
    maxevals -= consumed_evals;
  }

  return 0;
}

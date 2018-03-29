#include "TSDARPGenome.h"
#include "DARPEvaluator.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <math.h>

/*
 * Cordeau's initial values
 * evalWeightsUpdateValue_ (delta)  = 0.5
 * lackOfDivPenalization_  (lambda) = 0.015
 * maxTabuLife             (teta)   = 7.5 * log10( requests_num )
 */

TSDARPGenome::TSDARPGenome(unsigned x, GAGenome::Evaluator f, void* u)
  : DARPGenome(x,0.5,f,u),
    lackOfDivPenalization_(0.015),
    maxTabuLife_( 7.5 * log10(requestsNum()) ) {}

void TSDARPGenome::copy (const GAGenome& orig) {
  DARPGenome::copy(orig);

  const TSDARPGenome* other = dynamic_cast<const TSDARPGenome*>(&orig); assert(other);

  evalWeightsUpdateValue_ = other->evalWeightsUpdateValue_;
  lackOfDivPenalization_  = other->lackOfDivPenalization_;
  maxTabuLife_            = other->maxTabuLife_;

  tabuMem_ = other->tabuMem_;
  opFreq_  = other->opFreq_;
}

int TSDARPGenome::write (STD_OSTREAM& os) const {
  os << endl;
  os << "lackOfDivPenalization: " << lackOfDivPenalization_ << endl;
  os << "maxTabuLife: " << maxTabuLife_ << endl;
  os << "Tabu Memory: " << endl;
  for (map< pair<int,int>, int>::const_iterator it=tabuMem_.begin(); it!=tabuMem_.end(); it++) {
    os << "  Route: " << it->first.first << " Vertid: " << it->first.second << " tabu life: " << it->second << endl;
  }
  os << "Op Freq: " << endl;
  for (map< pair<int,int>, int>::const_iterator it=opFreq_.begin(); it!=opFreq_.end(); it++) {
    os << "  Route: " << it->first.first << " Vertid: " << it->first.second << " op freq: " << it->second << endl;
  }


  DARPGenome::write(os);

  return 0;
}

void TSDARPGenome::writeObject(ostream& os) const{
  DARPGenome::writeObject(os);

  os.write ( (char*) (&evalWeightsUpdateValue_), sizeof(evalWeightsUpdateValue_) );
  os.write ( (char*) (&lackOfDivPenalization_),  sizeof(lackOfDivPenalization_) );
  os.write ( (char*) (&maxTabuLife_),            sizeof(maxTabuLife_) );

#define WRITETSMAP(name) \
  {\
    int size = name.size();\
    os.write( (char*) (& size), sizeof(size) ); \
    for(  map< pair<int,int>, int >::const_iterator it = name.begin(); it != name.end(); it++ ) { \
      pair<int,int> key   = it->first; \
      int           value = it->second;\
      os.write( (char*) (& key),   sizeof(key) ); \
      os.write( (char*) (& value), sizeof(value) ); \
    }\
  }

  WRITETSMAP(tabuMem_);
  WRITETSMAP(opFreq_);
}

void TSDARPGenome::readObject (istream& is) {
  DARPGenome::readObject(is);

#define READTSMAP(name)\
  {\
    int size; is.read( (char*) &size, sizeof(size) );\
    for (int i=0; i<size; i++) {\
      pair<int,int> key; \
      int           value; \
      is.read( (char*) &key,   sizeof(key) );\
      is.read( (char*) &value, sizeof(value) );\
      name[key] = value;\
    }\
  }\

  is.read ( (char*) (&evalWeightsUpdateValue_), sizeof (evalWeightsUpdateValue_) );
  is.read ( (char*) (&lackOfDivPenalization_),  sizeof (lackOfDivPenalization_) );
  is.read ( (char*) (&maxTabuLife_),            sizeof (maxTabuLife_) );

  READTSMAP(tabuMem_);
  READTSMAP(opFreq_);
}

pair<RoutingGenome*,int> TSDARPGenome::bestFromNeighborhood() {
  TSDARPGenome* bestNeighbor = 0;
  long double        bestScore    = GAGenome::worstPossibleScore(); // Needed for moving always in case the score has the same value
  int           usedEvals    = 0;

  for (int route=0; route<length(); route++) {
    for (int pos=0; pos<routeLength(route); pos++) {

      int vertid = gene(route,pos);
      if (Vertex::isVertIdDelivery(vertid)) continue; // Because they are moved in pairs

      for (int newroute=0; newroute<length(); newroute++) {

        if (newroute == route) continue;

        TSDARPGenome* neighbor = dynamic_cast<TSDARPGenome*>( clone() );
        neighbor->moveToNeighbor(route,newroute,vertid);

        long double neighborScore = neighbor->score();
        if ( GAGenome::compareScores(neighborScore,score()) == GAGenome::WORSE ) {
          neighborScore += neighbor->divPenalization(newroute,vertid);
        }
        usedEvals += neighbor->nevals(); assert(neighbor->nevals() > 0);

        if ( allowedMove(route,newroute,vertid) || GAGenome::compareScores(neighborScore,score()) == GAGenome::BETTER)  {
          neighbor->updateTSMemsWithMov(newroute,vertid); // Note: Cannot be executed before the computation of the div. penalization

          if ( GAGenome::compareScores(neighborScore,bestScore) == GAGenome::BETTER) {
            if (bestNeighbor) delete bestNeighbor;
            bestScore    = neighborScore;
            bestNeighbor = neighbor;
          }

        }

        // If the neighbor is not stored as the new bestneighbor we delete it
        if (neighbor != bestNeighbor) delete neighbor;
      }

    }
  }
  assert(usedEvals > 0);

  return pair<RoutingGenome*,int> (bestNeighbor,usedEvals);
}

bool TSDARPGenome::allowedMove(int oldroute, int newroute, int vertid) {
  assert(oldroute >= 0 and oldroute < length()); assert(newroute >= 0 and newroute < length());

  return tabuMem_.find( routeVertIdKey(newroute,vertid) ) == tabuMem_.end();
}

void TSDARPGenome::moveToNeighbor(int oldroute, int newroute, int vertid) {
  assert(oldroute >= 0 and oldroute < length()); assert(newroute >= 0 and newroute < length());
  assert(Vertex::isVertIdPickup(vertid));

  int sibling_vertid = Vertex::getSiblingVertId(vertid);
  list<int> vertices;

  vertices.push_back(vertid);
  vertices.push_back(sibling_vertid);

  moveVertices(oldroute,newroute,vertices);
  pair<int,int> key = routeVertIdKey(newroute,vertid);
}

long double TSDARPGenome::divPenalization(int newroute, int vertid) {
  int routesNum = length();
  pair<int,int> key = routeVertIdKey(newroute,vertid);

  int freq_value = opFreq_[key];

  return lackOfDivPenalization_ * nonPenalizedScore() * sqrt( (long double) requestsNum() * routesNum ) * freq_value;
}

void TSDARPGenome::updateDynParameters() {
  lackOfDivPenalization_  = GARandomDouble(0,0.015);
  evalWeightsUpdateValue_ = GARandomDouble(0.0,0.5);
  maxTabuLife_            = GARandomDouble(0.0,7.5 * log10(requestsNum()) );
}

void TSDARPGenome::updatePenalizations() {
  ageMemory();
  DARPGenome::updatePenalizations();
}

void TSDARPGenome::ageMemory() {
  for (map< pair<int,int> ,int >::iterator it = tabuMem_.begin(); it!=tabuMem_.end(); it++) {
    it->second -= 1;
    if (it->second <= 0) tabuMem_.erase(it);
  }
}

void TSDARPGenome::updateTSMemsWithMov(int newroute,int vertid) {
  // The update of these values cannot be donde before the computation of the diversity penalization
  pair<int,int> key = routeVertIdKey(newroute,vertid);
  increaseFreq(key);
  addTabuMovement(key);
  assert(opFreq_[key] > 0); assert(tabuMem_[key] > 0);
}

string TSDARPGenome::memoryDataInformation() {
  stringstream msg;

  for (map< pair<int,int>, int>::const_iterator it=opFreq_.begin(); it!=opFreq_.end(); it++) {
    if (it->second > 1) {
      msg << "  Route: " << it->first.first << " VertId: " << it->first.second << " op freq: " << it->second << endl;
    }
  }

  return msg.str();
}

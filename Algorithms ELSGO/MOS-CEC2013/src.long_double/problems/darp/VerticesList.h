#ifndef VERTICESLIST_H
#define VERTICESLIST_H

#include <assert.h>
#include <map>
#include <iostream>
#include <vector>

class DistMatrix;
class CostMatrix;

using namespace std;

struct Vertex {
  enum VertexType {PICKUP, DELIVERY};

  int  id_;
  int  pos_;
  bool critic_;
  int  load_;

  /* int service_;  If in a future is used it should be uncommented from other regions of code (right now are commented so
   * it should be easy to find them by searching for service*/

  // TODO: Tener cuidado por si se modifican las ventanas en otro sitio que no sea la funcion de fitness
  long fbegin_;
  long fend_;

  VertexType type_;

  Vertex(int id, int pos, bool crit, int load, long fb, long fe, VertexType t) : id_(id), pos_(pos), critic_(crit), load_(load), fbegin_(fb), fend_(fe), type_(t) {
    assert(id_ != 0);
  }

  Vertex(const Vertex& other) {
    copy(other);
  }

  void copy(const Vertex& other) {
    id_     = other.id_;
    pos_    = other.pos_;
    critic_ = other.critic_;
    load_   = other.load_;
    fbegin_ = other.fbegin_;
    fend_   = other.fend_;
    type_   = other.type_;

  }

  Vertex* clone() { return new Vertex(*this); }

  Vertex& operator=(const Vertex& other) { copy(other); return *this; }

  int getSiblingVertexId() { return getSiblingVertId(id_);}

  bool canVertexGoBefore(Vertex& vert, const CostMatrix& costmatrix);
  bool canVertexGoAfter (Vertex& vert, const CostMatrix& costmatrix);

  bool canVertexGoBefore(Vertex& vert, const DistMatrix& costmatrix);
  bool canVertexGoAfter (Vertex& vert, const DistMatrix& costmatrix);

  bool isPickUp   (           ) { return isVertIdPickup(id_);   }
  bool isDelivery (           ) { return isVertIdDelivery(id_); }
  bool isSiblingOf(int vert_id) { return areSiblings(id_,vert_id); }

  static int  getSiblingVertId(int vert_id)          { return vert_id * -1;}
  static bool areSiblings     (int vert1, int vert2) { return getSiblingVertId(vert1) == vert2; }
  static bool isVertIdPickup  (int vert_id)          { return vert_id > 0;}
  static bool isVertIdDelivery(int vert_id)          { return !isVertIdPickup(vert_id); }

  static bool   relativeTime__;
  static long   relativeDelayConstant__;
  static long double maxRideValue__;
  static long   maxDelay__;

  static bool isRideTimeRelative() { return relativeTime__; }

  static long maxDelay()  { return maxDelay__; }

  static void setTimeConstants(bool relativeTime, long relativeDelayConstant, long double maxRideValue, long delay);
  static long maxRideTime     (int vert_id, const CostMatrix& distmatrix);
  static long maxRideTime     (int vert_or, int vert_dest, const CostMatrix& distmatrix);
  static long maxRideTime     (int id_or, int id_dest, const DistMatrix& distmatrix);
  static long maxRideTime     (long double distance);

};

inline std::ostream& operator<< (std::ostream& os, const Vertex& vert) {
  os << "(id: " << vert.id_ << " position: "  << vert.pos_;
  if (vert.critic_) os << " is critical";
  else             os << " not critical";
  os << " e_i: " << vert.fbegin_ << " l_i: " << vert.fend_;
  if (vert.type_ == Vertex::PICKUP) os << " pickup";
  else                              os << " delivery" ;
  os << ")";


  return os;
}

struct Request {
  Vertex* pickup_vert;
  Vertex* delivery_vert;

  Request(Vertex* p_vert, Vertex* d_vert) : pickup_vert(p_vert), delivery_vert(d_vert) {}

  Request(const Request& other) { copy(other); }

  void copy(const Request& other) {
    pickup_vert   = other.pickup_vert;
    delivery_vert = other.delivery_vert;
  }

  Request& operator=(const Request& other) { copy(other); return *this; }

  bool operator==(const Request& other) const { return pickup_vert == other.pickup_vert and delivery_vert == other.delivery_vert;}
};

inline std::ostream& operator<< (std::ostream& os, const Request& vert) {
  os << "Request: " << endl;
  os << " pickup   vertex: " << * vert.pickup_vert   << endl;
  os << " delivery vertex: " << * vert.delivery_vert << endl;

  return os;
}

class VerticesList {

protected:
  std::map<int, Vertex*> vertices_;

public:
  VerticesList () {}                                    // Just for the tests
  VerticesList (std::string reqFile, DistMatrix& dist); // TODO: No serviria con CostMatrix?
  VerticesList (vector<Request>& requests);             // For the request distributor option, it clones the vertices

  ~VerticesList () {
    std::map<int, Vertex*>::iterator it;
    for (it = vertices_.begin(); it != vertices_.end(); it++)
      delete it->second;
  }

  Vertex& getVertex(int id) const {
    assert(*vertices_.find(id) != *vertices_.end());
    Vertex& tmp = *(const_cast<VerticesList&>(*this).vertices_[id]);
    return tmp;
  }

  // Gets ownership of the vertex, i.e., it deletes in its destructor
  void addVertex(Vertex* d) {vertices_[d->id_] = d;}

  std::vector<int> getPickUpIds  () const;
  std::vector<int> getDeliveryIds() const;
  std::vector<int> getVerticesds () const;

  int size() const {return vertices_.size();}

  std::ostream& write (std::ostream& os) const {
    for (std::map<int, Vertex*>::const_iterator it = vertices_.begin(); it != vertices_.end(); it++) {
      os << *( (*it).second ) << std::endl;
    }

    return os;
  }

  std::vector<Request> getReqsList();

  std::vector<Request> getReqsFromVertIds(vector<int>& ids);

  int requestsNum() const { return vertices_.size() / 2; }

};

inline std::ostream& operator<< (std::ostream& os, const VerticesList& vertlist) {
  return vertlist.write(os);
}

#endif

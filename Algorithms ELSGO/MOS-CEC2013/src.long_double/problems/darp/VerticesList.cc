#include "VerticesList.h"
#include "DistMatrix.h"
#include "CostMatrix.h"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>

#include <sys/time.h>
#include <stdlib.h>

bool   Vertex::relativeTime__          = false;
long   Vertex::relativeDelayConstant__ = 0;
long double Vertex::maxRideValue__          = 0;
long   Vertex::maxDelay__              = 0;

bool Vertex::canVertexGoBefore(Vertex& vert, const CostMatrix& costmatrix) {
  return (fbegin_ + costmatrix.getCost(*this,vert) <= vert.fend_);
}

bool Vertex::canVertexGoBefore(Vertex& vert, const DistMatrix& distmatrix) {
  return (fbegin_ + distmatrix.getDist(this->pos_,vert.pos_) <= vert.fend_);
}


// TODO: Refactor the following methods so that only one version is implemented. For example,
// Dist and CostMatrix could implement the same interface
bool Vertex::canVertexGoAfter (Vertex& vert, const CostMatrix& costmatrix) {
  return vert.canVertexGoBefore(*this,costmatrix);
}

bool Vertex::canVertexGoAfter (Vertex& vert, const DistMatrix& distmatrix) {
  return vert.canVertexGoBefore(*this,distmatrix);
}


void Vertex::setTimeConstants(bool relativeTime, long relativeDelayConstant, long double maxRideValue, long delay) {
  relativeTime__          = relativeTime;
  relativeDelayConstant__ = relativeDelayConstant;
  maxRideValue__          = maxRideValue;
  maxDelay__              = delay;
}

using namespace std;

long Vertex::maxRideTime(int vert_id, const CostMatrix& costmatrix) {

  int pick_vertid,del_vertid;
  pick_vertid = del_vertid = vert_id;

  if (Vertex::isVertIdPickup(vert_id)) del_vertid  = Vertex::getSiblingVertId(vert_id);
  else                                 pick_vertid = Vertex::getSiblingVertId(vert_id);

  return maxRideTime(pick_vertid,del_vertid,costmatrix);
}

long Vertex::maxRideTime(int vert_or, int vert_dest, const CostMatrix& costmatrix) {
  return maxRideTime(costmatrix.getTravelTime(vert_or,vert_dest));
}

long Vertex::maxRideTime(int id_or, int id_dest, const DistMatrix& distmatrix) {
  return maxRideTime(distmatrix.getDist(id_or,id_dest));
}

long Vertex::maxRideTime(long double distance) {
  return (relativeTime__) ? (long) (maxRideValue__ * distance + relativeDelayConstant__ )
                          : (long) maxRideValue__;
}

VerticesList::VerticesList(std::string reqFile, DistMatrix& dist) {
  assert(Vertex::maxDelay__ >= 0 && Vertex::maxDelay__ <= 24*60*60); // The delay cannot be longer than a day
  // Open reqFile file as an input stream
  std::ifstream fin (reqFile.c_str());

  // Check if file can be correctly open
  if (fin.fail() or fin.bad()) {
    std::stringstream msg;
    msg << "[RequestsList] Error: file '" << reqFile << "' not found or could not be open" << std::endl;
    throw std::runtime_error(msg.str());
  }

  // Process the whole file
  while (!fin.eof()) {

    // Read one line
    std::string buffer;
    getline(fin, buffer);

    // Check if this line was properly read
    if (!fin.fail()) {

      std::istringstream iss(buffer);
      std::string token;
      std::vector<std::string> tokens;

      while (getline(iss, token, '|'))
        tokens.push_back(token);

      assert(tokens.size() >= 19);

      int id       = atoi(tokens[ 0].c_str());
      int idop     = Vertex::getSiblingVertId(id);
      int load     = atoi(tokens[11].c_str());
      int idOrigen = atoi(tokens[18].c_str());
      int idDest   = atoi(tokens[17].c_str());

      std::istringstream issfechahora (tokens[9]);
      std::vector<std::string> tokensfechahora;

      while (getline(issfechahora, token, '_')) tokensfechahora.push_back(token);

      assert(tokensfechahora.size() == 2);

      std::istringstream issfecha (tokensfechahora[0]);
      std::vector<std::string> tokensfecha;

      while (getline(issfecha, token, '-')) tokensfecha.push_back(token);

      tm date;

      date.tm_year = atoi(tokensfecha[0].c_str()) - 1900; // Years start in 1900
      date.tm_mon  = atoi(tokensfecha[1].c_str()) - 1;    // Months start in 0
      date.tm_mday = atoi(tokensfecha[2].c_str());

      assert(tokensfecha.size() == 3);

      std::istringstream isshora (tokensfechahora[1]);
      std::vector<std::string> tokenshora;

      while (getline(isshora, token, ':')) tokenshora.push_back(token);

      date.tm_hour  = atoi(tokenshora[0].c_str());
      date.tm_min   = atoi(tokenshora[1].c_str());
      date.tm_sec   = atoi(tokenshora[2].c_str());
      date.tm_isdst = 1; //Avoid Daylight saving time change

      assert(tokenshora.size() == 3);

      long fbeginPickup = mktime(&date);
      long fendPickup   = fbeginPickup + Vertex::maxDelay__;

      long fbeginDelivery = fbeginPickup + (long int) dist.getDist(idOrigen, idDest);
      long fendDelivery = fendPickup + Vertex::maxRideTime(idOrigen, idDest, dist);

//       std::cout << "id: " << id << std::endl;
//       std::cout << "fbeginPickup: " << fbeginPickup << std::endl;
//       std::cout << "fendPickup: " << fendPickup << std::endl;
//       std::cout << "fbeginDelivery: " << fbeginDelivery << std::endl;
//       std::cout << "fendDelivery: " << fendDelivery << std::endl;
//       std::cout << "dist: " << dist.getDist(idOrigen, idDest) << std::endl;
//       std::cout << std::endl;

      Vertex* req = NULL;

      req = new Vertex(id, idOrigen, true, load, fbeginPickup, fendPickup, Vertex::PICKUP);
      addVertex(req);

      req = new Vertex(idop, idDest, false, load, fbeginDelivery, fendDelivery, Vertex::DELIVERY);
      addVertex(req);
    }

  }

  // Close input file
  fin.close();

}

/*
 * Note: it clones the vertices contained in the vector of requests that are passed as argument and
 * constructs the list with these vertices
 */
VerticesList::VerticesList (vector<Request>& requests) {
  for(vector<Request>::iterator it=requests.begin(); it!=requests.end(); it++) {
    Vertex* vp = it->pickup_vert->clone();
    Vertex* vd = it->delivery_vert->clone();
    assert(vp && vd);

    addVertex(vp);
    addVertex(vd);
  }

}

std::vector<int> VerticesList::getVerticesds() const{

  std::vector<int> res;
  std::map<int, Vertex*>::const_iterator it;

  res.reserve(vertices_.size());

  for (it = vertices_.begin(); it != vertices_.end(); it++)
    res.push_back(it->first);

  return res;

}


std::vector<int> VerticesList::getPickUpIds() const {

  std::vector<int> res;
  std::map<int, Vertex*>::const_iterator it;

  res.reserve(vertices_.size());

  for (it = vertices_.begin(); it != vertices_.end(); it++)
    if (Vertex::isVertIdPickup(it->first))
      res.push_back(it->first);

  return res;

}


std::vector<int> VerticesList::getDeliveryIds() const{

  std::vector<int> res;
  std::map<int, Vertex*>::const_iterator it;

  res.reserve(vertices_.size());

  for (it = vertices_.begin(); it != vertices_.end(); it++)
    if (!Vertex::isVertIdPickup(it->first))
      res.push_back(it->first);

  return res;
}

std::vector<Request> VerticesList::getReqsList() {
  std::vector<Request> requestslist;

  std::map<int, Vertex*>::iterator it;
  for (it = vertices_.begin(); it != vertices_.end(); it++) {
    Vertex* vert = (*it).second;
    if ( vert->isPickUp() ) {
      Vertex* delvert = vertices_[vert->getSiblingVertexId()];
      assert(vert && delvert);
      Request req (vert,delvert);
      requestslist.push_back(req);
    }
  }

  return requestslist;
}

std::vector<Request> VerticesList::getReqsFromVertIds(vector<int>& ids) {
  assert(ids.size() % 2 == 0);

  std::vector<Request> requestslist;

  for (int i=0; i<ids.size(); i+=2) {
    Vertex& pickvertex = getVertex(ids[i]);
    Vertex& delvertex  = getVertex(ids[i+1]);
    Request req(&pickvertex,&delvertex);

    requestslist.push_back(req);
  }

  return requestslist;
}

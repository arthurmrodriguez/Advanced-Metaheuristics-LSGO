#ifndef ROUTINGGENOME_H
#define ROUTINGGENOME_H

#include "GAGenome.h"
#include <stdexcept>
#include <utility>

/*
 * Used for specifyng a common interface of the genome that all the routing algorithms are going to need
 * The problem-dependant genome should at least inherit from this genome
 */

class RoutingGenome : virtual public GAGenome {
protected:

  bool feasible_; // Is this solution feasible? (for those algorithms that need it)

public:

  RoutingGenome () : GAGenome(), feasible_(true) {}
  RoutingGenome (const RoutingGenome& orig) : GAGenome(orig) {copy(orig);}

  virtual void copy (const GAGenome & orig) {
    GAGenome::copy(orig);
    onlyCopyAttributes(orig);
  }

  /*
   * If this changes, DARPGenome (addRoute and copy) should also change
   */
  virtual void onlyCopyAttributes (const GAGenome & orig) {
    const RoutingGenome* other = dynamic_cast<const RoutingGenome*>(&orig); assert(other);
    feasible_ = other->feasible_;
  }

  bool feasible () const { evaluate(); return feasible_; }
  bool feasible (bool f) {return feasible_ = f;}

  virtual long double nonPenalizedScore() const = 0;

  virtual int evalRouteCalls() const = 0;

  // These methods are used by the TS Algorithm

  virtual pair<RoutingGenome*,int> bestFromNeighborhood() { throw runtime_error("bestFromNeighborhood not implemented"); }

  virtual void updateDynParameters() { throw runtime_error("bestFromNeighborhood not implemented"); }
  virtual void updatePenalizations() { throw runtime_error("bestFromNeighborhood not implemented"); }

  virtual string memoryDataInformation() { throw runtime_error("bestFromNeighborhood not implemented"); }

  virtual long double deliveryDelay() const = 0;
  virtual long double pickupDelay  () const = 0;

  // These methods are neccessary for the distributed model
  virtual int  numRoutes() const = 0;
  virtual void addRoutes(const GAGenome& gen) = 0;
  virtual void addRoute (const GAGenome& gen, int routepos) = 0;
  virtual void emptyRoutes() = 0;
  virtual int  routeLength(int route) const = 0;

  virtual RoutingGenome* cloneWithoutEmptyRoutes() const = 0;

  virtual long TWV()   const = 0;
  virtual long loadV() const = 0;
  virtual long rideV() const = 0;

  virtual void printRoutes()              const = 0; // Just for debugging
  virtual void checkNoVertexIsMissing()   const = 0; // Just for debugging
  virtual void checkNoClonesInSameRoute() const = 0; // Just for debugging
  virtual void checkNoClonesInAllRoutes() const = 0; // Just for debugging
};

#endif

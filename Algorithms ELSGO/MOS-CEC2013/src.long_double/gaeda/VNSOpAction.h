#ifndef VNSOPACTION
#define VNSOPACTION

#include "genomes/GAGenome.h"
#include <string>

using namespace std;

/*
 * Used for storing (and reproducing) the actions of executing a VNSOpAction.
 * The initial idea for this class (based on the Command pattern) was to be able to
 * combine actions later and see if the combination of actions produced a better result.
 * Also this class can be used as a log of the best actions that improved the solution
 */
class VNSOpAction {
public:
  virtual void executeOver (GAGenome& gen) const = 0;
  virtual ~VNSOpAction() {}

  virtual bool         isCompatibleWith(const VNSOpAction& ot) const = 0;
  virtual bool         equal           (const VNSOpAction& ot) const = 0;
  virtual VNSOpAction* clone() const = 0;
  virtual string       name()  const = 0;
};

struct ActionDatum{
  virtual ~ActionDatum() {}

  virtual bool         isCompatibleWith(const ActionDatum& ot) const = 0;
  virtual bool         equal           (const ActionDatum& ot) const = 0;
  virtual ActionDatum* clone() const = 0;
};

class VNSShakerAction : public VNSOpAction {
protected:
  vector<ActionDatum*> data_;
public:

  VNSShakerAction(){}

  VNSShakerAction(const VNSShakerAction& other) ;

  virtual ~VNSShakerAction();

  virtual bool isCompatibleWith(const VNSOpAction& ot) const;

  virtual bool equal(const VNSOpAction& ot) const;
};



class NoVNSOpAction : public VNSOpAction {
  virtual ~NoVNSOpAction(){}
  void executeOver (GAGenome& gen) const {}

  virtual VNSOpAction* clone() const {
    return new NoVNSOpAction();
  }

  virtual bool isCompatibleWith(const VNSOpAction& ot) const { return true; }
  virtual bool equal           (const VNSOpAction& ot) const { return true; }


  virtual string name()  const { return "NoVNSOpAction"; }
};


#endif

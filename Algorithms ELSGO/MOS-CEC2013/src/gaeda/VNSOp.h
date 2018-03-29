#ifndef VNSOP
#define VNSOP

#include "VNSOpAction.h"
#include <string>
#include <iostream>
#include <sstream>

using namespace std;

class VNSOp {

protected:
  string name_;
  int    activations_;
  int    score_improvs_;
  int    cost_improvs_;
  int    nevals_;

  static int allops_nevals__;

public:
  VNSOp(string name) : name_(name), activations_(0), score_improvs_(0), cost_improvs_(0), nevals_(0) {}

  virtual ~VNSOp() {}

  void addResult(bool score_improved, bool cost_improved, int nevals);

  virtual string name() const { return name_; }

  virtual string getStatsResults();

  virtual VNSOp* clone() = 0;

  virtual VNSOpAction* operator() (GAGenome& gen) = 0;

};

class VNSShaker : public VNSOp{
protected:
  int    size_;

public:

  VNSShaker(string name, int sz) : VNSOp(name), size_(sz) {}

  virtual ~VNSShaker() {}

  int size() const { return size_; }

  virtual string name() const {
    stringstream msg; msg << name_ << " size: " << size_;
    return msg.str();
  }

  virtual VNSShaker& operator=(const VNSShaker& shaker) {
    size_ = shaker.size_;
    return *this;
  }

  virtual bool isStochastic() { return true; }

};

inline STD_OSTREAM& operator<< (STD_OSTREAM& os, const VNSShaker& v) {
  os << v.name() << " size: " << v.size();
  return (os);
}

#endif

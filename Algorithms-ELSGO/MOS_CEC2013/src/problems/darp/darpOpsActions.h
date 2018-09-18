#ifndef DAROPSACTIONS_H_
#define DAROPSACTIONS_H_

#include "aux.h"
#include <GAEDAlib.h>
#include <list>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>

using namespace std;


struct DARPActionDatum : public ActionDatum {
  int       orig_route;
  list<int> orig_vert_ids;

  DARPActionDatum(int o, list<int>& ids) : ActionDatum(), orig_route(o), orig_vert_ids(ids) {}
  virtual ~DARPActionDatum() {}

  virtual bool equal(const ActionDatum& ot) const {
    const DARPActionDatum& other = dynamic_cast<const DARPActionDatum&>(ot);
    return orig_route == other.orig_route and
           areReqsEqual(orig_vert_ids,other.orig_vert_ids);
  }

  virtual bool isCompatibleWith(const ActionDatum& ot) const {
    const DARPActionDatum& other = dynamic_cast<const DARPActionDatum&>(ot);

    return ! haveAnyReqInCommon(orig_vert_ids,other.orig_vert_ids);
  }
};

struct MoveActionDatum : public DARPActionDatum {
  int dest_route;

  MoveActionDatum(int o, list<int>& ids, int d) : DARPActionDatum(o,ids), dest_route(d) {}
  virtual ~MoveActionDatum() {}

  virtual ActionDatum* clone() const { return new MoveActionDatum(*this); }

  virtual bool equal(const ActionDatum& ot) const {
    const MoveActionDatum* other = dynamic_cast<const MoveActionDatum*>(&ot);
    return other and DARPActionDatum::equal(ot) and dest_route == other->dest_route;
  }
};

class MoveAction : public VNSShakerAction {
public:

  MoveAction(int o, list<int>& ids, int d) {
    addDatum(o,ids,d);
  }
  virtual ~MoveAction(){}

  virtual void addDatum (int orig, list<int>& ids, int dest) {
    data_.push_back(new MoveActionDatum(orig,ids,dest));
  }

  virtual void executeOver (GAGenome& gen)  const {
    for (int i=0; i<data_.size(); i++) { // Several move actions could be nested
      const MoveActionDatum& dat = dynamic_cast<const MoveActionDatum &>(* data_[i]);
      assert(dat.orig_vert_ids.size() > 0);

      dynamic_cast<DARPGenome&>(gen).moveVertices (dat.orig_route, dat.dest_route, dat.orig_vert_ids);
    }
  }

  virtual VNSOpAction* clone() const { return new MoveAction(*this); }

  string name() const {
    stringstream msg; msg << "MoveAction: List of actions: " << endl;
    for (int i=0; i<data_.size(); i++) {
      const MoveActionDatum& values = dynamic_cast<MoveActionDatum&>(* data_[i]);
      msg << "  Move from route: " << values.orig_route << " to route: " << values.dest_route << " these vertices: ";
      for (list<int>::const_iterator it=values.orig_vert_ids.begin(); it!=values.orig_vert_ids.end(); it++) {
        msg << *it << ",";
      }
      msg << endl;
    }

    return msg.str();
  }
};

struct SwapActionDatum : public MoveActionDatum {
  list<int> dest_vert_ids;

  SwapActionDatum(int o, list<int>& o_ids, int d, list<int>& d_ids) :
                 MoveActionDatum(o,o_ids,d), dest_vert_ids(d_ids) {}

  virtual ~SwapActionDatum() {}

  virtual bool isCompatibleWith(const ActionDatum& ot) const {
    if (!MoveActionDatum::isCompatibleWith(ot)) return false;

    const DARPActionDatum* other = dynamic_cast<const DARPActionDatum*>(&ot); assert(other!=0);
    if (haveAnyReqInCommon(dest_vert_ids,other->orig_vert_ids)) return false;


    const SwapActionDatum* swap_other = dynamic_cast<const SwapActionDatum*>(&ot);
    if (swap_other) {
      return ! haveAnyReqInCommon(orig_vert_ids,swap_other->orig_vert_ids) and
             ! haveAnyReqInCommon(orig_vert_ids,swap_other->dest_vert_ids) and
             ! haveAnyReqInCommon(dest_vert_ids,swap_other->orig_vert_ids) and
             ! haveAnyReqInCommon(dest_vert_ids,swap_other->dest_vert_ids);
    }

    return true;
  }

  virtual ActionDatum* clone() const { return new SwapActionDatum(*this); }

  virtual bool equal(const ActionDatum& ot) const {
    const SwapActionDatum* other = dynamic_cast<const SwapActionDatum*>(&ot);
    return other and
           MoveActionDatum::equal(ot) and
           areReqsEqual(dest_vert_ids,other->dest_vert_ids);
  }

};

class SwapAction : public VNSShakerAction {
public:
  SwapAction(int orig, list<int>& orig_ids, int dest, list<int>& dest_ids) {
    data_.push_back(new SwapActionDatum(orig,orig_ids,dest,dest_ids));
  }

  virtual ~SwapAction() {}

  virtual void executeOver (GAGenome& g) const {
    DARPGenome& gen = dynamic_cast<DARPGenome&>(g);

    assert(data_.size() == 1);
    const SwapActionDatum& datum = dynamic_cast<SwapActionDatum&>(* data_[0]);

    assert(datum.orig_route>=0 && datum.orig_route<gen.size()); assert(datum.dest_route>=0 && datum.dest_route<gen.size());
    assert(datum.orig_vert_ids.size() > 0);                      assert(datum.dest_vert_ids.size() > 0);

    gen.swapSeqs(datum.orig_route,datum.orig_vert_ids,datum.dest_route,datum.dest_vert_ids);
  }

  virtual VNSOpAction* clone() const { return new SwapAction(*this); }

  string name() const {
    stringstream msg; msg << "SwapAction: List of actions: " << endl;
    for (int i=0; i<data_.size(); i++) {
      const SwapActionDatum& values = dynamic_cast<const SwapActionDatum&>(* data_[i]);
      msg << "  Swap from route: " << values.orig_route << " these vertices: ";
      for (list<int>::const_iterator it=values.orig_vert_ids.begin(); it!=values.orig_vert_ids.end(); it++) {
        msg << *it << ",";
      }
      msg << endl;
      msg << "  with route: " << values.dest_route << " and its vertices:  ";
      for (list<int>::const_iterator it=values.dest_vert_ids.begin(); it!=values.dest_vert_ids.end(); it++) {
        msg << *it << ",";
      }
      msg << endl;
    }

    return msg.str();
  }
};

struct SwapWithInsertionPosDatum : public SwapActionDatum {
  int orig_ins_pos;
  int dest_ins_pos;

  SwapWithInsertionPosDatum(int or_route, int or_ins_pos, list<int>& or_ids,
                            int de_route, int de_ins_pos, list<int>& de_ids )
                           : SwapActionDatum(or_route,or_ids,de_route,de_ids),
                             orig_ins_pos(or_ins_pos), dest_ins_pos(de_ins_pos) {}

  virtual ActionDatum* clone() const { return new SwapWithInsertionPosDatum(*this); }
};


class SwapWithInsertionPosAction : public VNSShakerAction {
public:

  SwapWithInsertionPosAction() {}
  virtual ~SwapWithInsertionPosAction() {}

  SwapWithInsertionPosAction(SwapWithInsertionPosDatum& datum) {
    addDatum(datum);
  }

  void addDatum(SwapWithInsertionPosDatum& datum) {
    data_.push_back(datum.clone());
  }

  virtual void executeOver (GAGenome& g) const {
    DARPGenome& gen = dynamic_cast<DARPGenome&>(g);

    for (int i=0; i<data_.size(); i++) {
      const SwapWithInsertionPosDatum& dat = dynamic_cast<const SwapWithInsertionPosDatum&>(* data_[i]);

      int orig_ins_pos = dat.orig_ins_pos;
      int dest_ins_pos = dat.dest_ins_pos;
      int orig_max_ins_pos = gen.routeLength(dat.orig_route) - dat.orig_vert_ids.size();
      int dest_max_ins_pos = gen.routeLength(dat.dest_route) - dat.dest_vert_ids.size();

      // Minor correction for allowing chaining of complex actions. If the insertion pos is bigger than the
      // possible maximum insertion pos (taking into account that we are FIRST removing some vertices) we correct it
      // to the max pos

      if (orig_ins_pos >  orig_max_ins_pos) orig_ins_pos = orig_max_ins_pos;
      if (dest_ins_pos >  dest_max_ins_pos) dest_ins_pos = dest_max_ins_pos;

      gen.swapWithInsertPosSeqs(dat.orig_route, orig_ins_pos, dat.orig_vert_ids,
                                dat.dest_route, dest_ins_pos, dat.dest_vert_ids);

    }
  }

  virtual VNSOpAction* clone() const { return new SwapWithInsertionPosAction(*this); }

  string name() const {
    stringstream msg; msg << "SwapActionWithInsertionPos: List of actions: " << endl;
    for (int i=0; i<data_.size(); i++) {
      const SwapWithInsertionPosDatum& values = dynamic_cast<const SwapWithInsertionPosDatum&>(* data_[i]);
      msg << "  Swap from route: " << values.orig_route << " these vertices: ";
      for (list<int>::const_iterator it=values.orig_vert_ids.begin(); it!=values.orig_vert_ids.end(); it++) {
        msg << *it << ",";
      }
      msg << " and insertion pos: " << values.orig_ins_pos << endl;
      msg << "  with route: " << values.dest_route << " and its vertices:  ";
      for (list<int>::const_iterator it=values.dest_vert_ids.begin(); it!=values.dest_vert_ids.end(); it++) {
        msg << *it << ",";
      }
      msg << " and insertion pos: " << values.dest_ins_pos << endl;
      msg << endl;
    }

    return msg.str();
  }
};


#endif /* DAROPSACTIONS_H_ */

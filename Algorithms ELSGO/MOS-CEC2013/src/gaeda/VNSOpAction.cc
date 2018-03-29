#include "VNSOpAction.h"

VNSShakerAction::VNSShakerAction(const VNSShakerAction& other) {
  for (vector<ActionDatum*>::const_iterator it=other.data_.begin(); it!=other.data_.end(); it++) {
    data_.push_back( (*it)->clone() );
  }
}

VNSShakerAction::~VNSShakerAction() {
  for (vector<ActionDatum*>::iterator it=data_.begin(); it!=data_.end(); it++) delete *it;
}

bool VNSShakerAction::isCompatibleWith(const VNSOpAction& ot) const {
  const VNSShakerAction& other = dynamic_cast< const VNSShakerAction& >(ot);

  for (int i=0; i<data_.size(); i++) {
    for (int j=0;j<other.data_.size();j++) {
      if (! data_[i]->isCompatibleWith(* other.data_[j]) ) return false;
    }
  }
  return true;
}

bool VNSShakerAction::equal(const VNSOpAction& ot) const {
  const VNSShakerAction& other = dynamic_cast< const VNSShakerAction& >(ot);

  for (int i=0; i<data_.size(); i++) {
    for (int j=0;j<other.data_.size();j++) {
      if (! data_[i]->equal(* other.data_[j]) ) {
        return false;
      }
    }
  }

  return true;
}



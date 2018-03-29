#include "aux.h"
#include <GAEDAlib.h>
#include <algorithm>
#include <assert.h>

int findPosOfValue(int value, std::vector<int>& values, int startpos) {
  int pos = -1;

  if (startpos<0) startpos=0; // Used if we want to specify a starting position for efficiency reasons

  for (int i=startpos; i<values.size(); i++ ) {
    if (values[i] == value) {
      pos = i;
      break;
    }
  }

  return pos;
}

bool haveAnyReqInCommon(const list<int>& l1, const list<int>& l2) {
  // Since these are request values, there should be a positive and a negative value.
  // Therefore we should only check the positive value
  for (list<int>::const_iterator it=l1.begin(); it!=l1.end(); it++) {
    // Check that for each value, the -1 * equivalent is found
    assert( find(l1.begin(), l1.end(), -1 * (*it) ) != l1.end() );

    if (*it > 0) {
      if ( find( l2.begin(), l2.end(), *it) != l2.end() ) return true;
    }
  }

  return false;
}

bool areReqsEqual (const list<int>& l1, const list<int>& l2) {
  if (l1.size() != l2.size() ) return false;

  for (list<int>::const_iterator it=l1.begin(); it!=l1.end(); it++) {
    // Check that for each value, the -1 * equivalent is found
    assert( find(l1.begin(), l1.end(), -1 * (*it) ) != l1.end() );

    if (*it > 0) {
      if ( find( l2.begin(), l2.end(), *it) == l2.end() ) return false;
    }
  }
  return true;
}

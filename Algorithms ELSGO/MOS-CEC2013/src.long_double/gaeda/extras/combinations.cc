#include "combinations.h"
#include <algorithm>
#include <assert.h>

void getAllCombinations(const vector<int>& values, int size, list<int>& sol, int index, vector< list<int>* >& allsols) {
  assert(size >=0);

  if (size == 0) {
    allsols.push_back( new list<int>(sol) );
  }
  else {

    for (int pos=index; pos<values.size(); pos++) {
      sol.push_back(values[pos]);

      getAllCombinations(values,size-1,sol,pos+1, allsols);
      sol.pop_back();
    }

  }
}

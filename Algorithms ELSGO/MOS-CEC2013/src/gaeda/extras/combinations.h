#ifndef COMBI
#define COMBI

#include <vector>
#include <list>

using namespace std;

void getAllCombinations(const vector<int>& values, int size, list<int>& sol, int index, vector< list<int>* >& allsols);

#endif

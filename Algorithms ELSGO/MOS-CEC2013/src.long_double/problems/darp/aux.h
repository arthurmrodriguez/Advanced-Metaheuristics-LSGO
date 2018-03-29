#ifndef DARPAUX_
#define DARPAUX_

#include <vector>
#include <list>

using namespace std;

int findPosOfValue(int value, vector<int>& values, int startpos=-1);
bool haveAnyReqInCommon(const list<int>& l1, const list<int>& l2);
bool areReqsEqual      (const list<int>& l1, const list<int>& l2);

#endif

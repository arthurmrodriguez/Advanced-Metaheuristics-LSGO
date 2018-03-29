#include "RandomReqDistributor.h"
#include <garandom.h>
#include <algorithm>


// random generator function:
ptrdiff_t myrandom (ptrdiff_t i) { return rand()%i;}

// pointer object to it:
ptrdiff_t (*p_myrandom)(ptrdiff_t) = myrandom;


vector< vector<Request> > RandomReqDistributor::distributeReqs(vector<Request>& reqs, ReqDistMatrix& reqdistmatrix,  int nodes) {
  assert(reqs.size() >= nodes);

  srand( GAGetRandomNoRankSeed() );
  random_shuffle(reqs.begin(),reqs.end(), p_myrandom);

  int reqs_per_node = reqs.size()/nodes;

  vector< vector<Request> > dist_reqs(nodes);
  int reqspos = 0;
  for (int nodepos=0; nodepos<nodes; nodepos++) {
    for (int j=0; j<reqs_per_node; j++, reqspos++) {
      dist_reqs[nodepos].push_back( reqs[reqspos] );
    }
  }

  for (; reqspos<reqs.size(); reqspos++) {
    dist_reqs[ rand() % dist_reqs.size() ].push_back( reqs[reqspos] );
  }

#ifdef DEBUG  // Just for debugging
  int dist_reqs_size = 0;
  for (int i=0; i<nodes; i++) dist_reqs_size += dist_reqs[i].size();
  assert(dist_reqs_size == reqs.size());
#endif

  return dist_reqs;
}

#include "KMedoidsReqDistributor.h"
#include <algorithm>
#include <garandom.h>
#include <assert.h>

KMedoidsReqDistributor::KMedoidsReqDistributor(int npass, kmedoidsMode mode) : npass_(npass), mode_(mode) {
  assert(mode >= 0);
}

vector< vector<Request> > createDistReqs(vector<Request>& reqs, int nodes, int clusterassignment[]) {
  vector< vector<Request> > dist_reqs(nodes);

  for (int i=0; i<reqs.size(); i++) {
    int clusterid = clusterassignment[i];
    assert(clusterid >=0 && clusterid < nodes);
    dist_reqs[clusterid].push_back(reqs[i]);
  }

  return dist_reqs;
}

vector< vector<Request> > KMedoidsReqDistributor::distributeReqs(vector<Request>& reqs, ReqDistMatrix& reqdistmatrix, int nodes) {
  int    clusterids[reqs.size()]; // the ids of the cluster that each individual belong to
  long double error;                   // The error of the best solution
  int    ifound;                  // Not used (the number of times the optimal solution is found

  setClusterFuncsSeed(GAGetRandomNoRankSeed()); // all the islands must use the same seed

  kmedoids(nodes,reqs.size(),reqdistmatrix,npass_,clusterids,&error,&ifound,mode_);

  vector< vector<Request> > dist_reqs = createDistReqs(reqs,nodes,clusterids);

#ifdef DEBUG
  int assignedreqs = 0; for (int i=0; i<dist_reqs.size(); i++) assignedreqs += dist_reqs[i].size();
  assert(assignedreqs > 0);
#endif

  return dist_reqs;
}

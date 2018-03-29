#include "LazyReqDistMatrix.h"

LazyReqDistMatrix::LazyReqDistMatrix(vector<Request>& reqs, CostMatrix& costmatrix, DARPEvaluator& eval,
                                     double TWVWeight, double rideVWeight, double loadVWeight, bool useconstrpen,
                                     bool usewaitingpen, long waitingPenThreshold, double waitingPenConstant) :
                                      ReqDistMatrix(costmatrix,eval,TWVWeight,rideVWeight,loadVWeight,useconstrpen,usewaitingpen,waitingPenThreshold,waitingPenConstant) {
  // Due to the C++ limitation virtual methods of the child classes cannot be called from the constructor
  // Therefore, we repeat the steps that are already conducted in the DistMatrix class but, in this case,
  // For the specific parts of the this class

  setUpDistMatrix(reqs);
  initializeDistMatrix(reqs);
}

void LazyReqDistMatrix::clear() {
  ReqDistMatrix::clear();
  datastoragestatus_.clear();
}

void LazyReqDistMatrix::resize(int s) {
  DistMatrix::resize(s);
  datastoragestatus_.resize(s);
}

void LazyReqDistMatrix::resize(int pos, int s) {
  DistMatrix::resize(pos,s);
  datastoragestatus_[pos].resize(s,false);
}

void LazyReqDistMatrix::setUpDistMatrix(vector<Request>& reqs) {
  reqs_ = reqs;
  ReqDistMatrix::setUpDistMatrix(reqs);
}

void LazyReqDistMatrix::initializeDistMatrix(vector<Request>& reqs) {
  // Initialized as specified in cluster.h
  for (int i=0; i<reqs.size(); i++) {
    if (i>0) resize(i,i);
    reqid2matrixpos_[reqs_[i].pickup_vert->id_] = i;
  }
}

double LazyReqDistMatrix::getDist(int ori, int dest) const {
  assert(ori >=0 && ori <reqs_.size() );
  assert(dest >=0 && dest <reqs_.size() );

  // Hack to avoid the constant constraint (cannot remove const since it affects dist matrix which affects
  // more classes)
  LazyReqDistMatrix& mythis = const_cast<LazyReqDistMatrix&>(*this);
  mythis.computeDataIfNotStored(ori,dest);

  return DistMatrix::getDist(ori,dest);
}

void LazyReqDistMatrix::computeDataIfNotStored(int ori, int dest) {
  assert(ori  >=0 && ori  < datastoragestatus_.size() );
  assert(dest >=0 && dest < datastoragestatus_[ori].size() );

  if (! datastoragestatus_[ori][dest] ) {
    setDist(ori,dest,computeValue(ori,dest));
    datastoragestatus_[ori][dest] = true;
  }
}

double LazyReqDistMatrix::computeValue(int ori, int dest) {
  assert(dest<ori);
  assert(reqs_.size() > 0);
  assert(ori < reqs_.size() );
  assert(dest < reqs_.size() );
  assert( reqid2matrixpos_[reqs_[ori].pickup_vert->id_]  == ori);
  assert( reqid2matrixpos_[reqs_[dest].pickup_vert->id_] == dest);

  return bestScoreFromAllCombs(reqs_[ori],reqs_[dest]);
}

double LazyReqDistMatrix::computedValuesRatio() {
  int computedvalues = 0;
  int allvalues      = 0;
  for (int i=0; i<datastoragestatus_.size(); i++) {
    for (int j=0; j<i; j++) {
      if (datastoragestatus_[i][j]) computedvalues++;
      allvalues++;
    }
  }
  return (double) computedvalues/(double)allvalues;
}

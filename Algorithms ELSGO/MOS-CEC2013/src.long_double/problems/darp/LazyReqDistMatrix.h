#ifndef LAZYREQDISTMATRIX_H_
#define LAZYREQDISTMATRIX_H_

#include "VerticesList.h"
#include "ReqDistMatrix.h"
#include "CostMatrix.h"
#include "DARPEvaluator.h"
#include <vector>

using namespace std;

class LazyReqDistMatrix : public ReqDistMatrix  {
private:
  vector<Request>        reqs_;
  vector< vector<bool> > datastoragestatus_;

protected:

  virtual void clear();
  virtual void resize(int s);
  virtual void resize(int pos, int s);

  LazyReqDistMatrix(CostMatrix& costmatrix, DARPEvaluator& eval, bool useconstrpen,
                bool usewaitingpen, long waitingPenThreshold, long double waitingPenConstant) :
                  ReqDistMatrix(costmatrix,eval,useconstrpen,usewaitingpen,waitingPenThreshold,waitingPenConstant){}

public:

  LazyReqDistMatrix(vector<Request>& reqs, CostMatrix& costmatrix, DARPEvaluator& eval,
                    long double TWVWeight, long double rideVWeight, long double loadVWeight, bool useconstrpen,
                    bool usewaitingpen, long waitingPenThreshold, long double waitingPenConstant);

  virtual ~LazyReqDistMatrix(){}

  virtual void setUpDistMatrix(vector<Request>& reqs);
  virtual void initializeDistMatrix(vector<Request>& reqs);

  virtual long double getDist(int ori, int dest) const;

  virtual void computeDataIfNotStored(int ori, int dest);
  virtual long double computeValue(int ori, int dest);
  long double  computedValuesRatio();

};

#endif /* LAZYREQDISTMATRIX_H_ */

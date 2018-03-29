#ifndef REQDISTMATRIX_H_
#define REQDISTMATRIX_H_

#include "VerticesList.h"
#include "CostMatrix.h"
#include "DARPEvaluator.h"
#include <vector>

using namespace std;

class ReqDistMatrix : virtual public DistMatrix {
protected:

  map<int,int> reqid2matrixpos_; // stores the correspondence between a request id (in this case the pickup vertex id) and its position in distmatrix_

  // Needed for the efficient computation of the matrix
  CostMatrix&    costmatrix_;
  DARPEvaluator& eval_;

  // These are needed to compute the score with constant weights. Otherwise, if we use the evaluation stored in the
  // evaluation object and a lazy approach, these could change over time
  double TWVWeight_, rideVWeight_, loadVWeight_;

  // These three are used for computing the penalization of the waiting times
  bool   usewaitingpen_;
  bool   useconstrpen_;
  long   waitingPenThreshold_;
  double waitingPenConstant_;

  virtual void setUpDistMatrix     (vector<Request>& reqs);
  virtual void initializeDistMatrix(vector<Request>& reqs);
  virtual void clear();


  virtual double bestScoreFromAllCombs(Request& req1, Request& req2);
  double         reqSequenceScore     (Request& req1, Request& req2, int req1d, int req2p, int req2d);


  // For the tests, uses the weights of the evaluator
  ReqDistMatrix(CostMatrix& costmatrix, DARPEvaluator& eval, bool useconstrpen,
                bool usewaitingpen, long waitingPenThreshold, double waitingPenConstant);

public:

  /*
   * HACK due to the lack of being able to call virtual methods from the constructors in C++. This way,
   * the classes that inherit this class can call a constructor and then call the setup and initialize methods
   * so that the virtual methods from their classes are called and no code is duplicated (except for this constructor)
   */
  ReqDistMatrix(CostMatrix& costmatrix, DARPEvaluator& eval,
                double TWVWeight, double rideVWeight, double loadVWeight, bool useconstrpen,
                bool usewaitingpen, long waitingPenThreshold, double waitingPenConstant);


  ReqDistMatrix(vector<Request>& reqs, CostMatrix& costmatrix, DARPEvaluator& eval,
                double TWVWeight, double rideVWeight, double loadVWeight, bool useconstrpen,
                bool usewaitingpen, long waitingPenThreshold, double waitingPenConstant);

  virtual ~ReqDistMatrix() {}

  double getDistOfReqsPickupIds(const int reqid1, const int reqid2) {
    assert(Vertex::isVertIdPickup(reqid1) && Vertex::isVertIdPickup(reqid2) );
    assert(reqid2matrixpos_.find(reqid1) != reqid2matrixpos_.end());
    assert(reqid2matrixpos_.find(reqid2) != reqid2matrixpos_.end());

    int pos1 = reqid2matrixpos_[reqid1];
    int pos2 = reqid2matrixpos_[reqid2];
    assert(reqid1 == reqid2 || pos1 != pos2);

    if (pos1 < pos2) {  // The matrix stores the results of the combinations 1-0 2-0 2-1 3-0 3-1 3-2 ...
      int tmp = pos1;
      pos1    = pos2;
      pos2    = tmp;
    }

    return getDist(pos1,pos2);
  }

};

#endif /* REQDISTMATRIX_H_ */

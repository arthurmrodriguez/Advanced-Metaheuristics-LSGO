#ifndef PARREQDISTMATRIX_H_
#define PARREQDISTMATRIX_H_

#include "ReqDistMatrix.h"
#include <GAEDAlib.h>
#include <vector>

using namespace std;

class ParReqDistMatrix : virtual public ReqDistMatrix {

protected:

  void           initializeDistMatrix(vector<Request>& reqs, CommManager& comm);
  vector<long double> computeCorrespondingVectorData(vector<Request>& reqs, CommManager& comm);
  void           computeFirstRowColumnAndNElems(int matrix_size, CommManager& comm,
                                               /*out*/ int& matrix_firstpos_row, /*out*/ int& matrix_firstpos_col,/*out*/ int& nelems);
  void           parallelfillOfData(CommManager& comm, vector<long double>& data);
  void           fillMatrixWithSerData(vector<long double>& data);

public:

  ParReqDistMatrix(vector<Request>& reqs, CostMatrix& costmatrix, DARPEvaluator& eval,
                long double TWVWeight, long double rideVWeight, long double loadVWeight, bool useconstrpen,
                bool usewaitingpen, long waitingPenThreshold, long double waitingPenConstant,
                CommManager& comm);

  virtual ~ParReqDistMatrix() {}
};

#endif /* PARREQDISTMATRIX_H_ */

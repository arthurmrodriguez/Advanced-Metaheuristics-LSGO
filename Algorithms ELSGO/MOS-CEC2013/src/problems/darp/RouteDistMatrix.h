#ifndef ROUTEDISTMATRIX_H_
#define ROUTEDISTMATRIX_H_

#include "ReqDistMatrix.h"
#include "DistMatrix.h"
#include "DARPEvaluator.h"
#include <genomes/GAGenome.h>
#include <vector>

using namespace std;

class RouteDistMatrix : public DistMatrix {
public:

  RouteDistMatrix(vector<GAGenome*>& routes, ReqDistMatrix& reqdistmatrix);
  virtual ~RouteDistMatrix();
};

#endif /* ROUTEDISTMATRIX_H_ */

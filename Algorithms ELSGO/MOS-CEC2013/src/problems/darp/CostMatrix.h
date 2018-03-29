#ifndef COSTMATRIX_H
#define COSTMATRIX_H

#include "VerticesList.h"
#include "DistMatrix.h"
#include "DARPEvaluator.h"
#include <vector>
#include <assert.h>

class CostMatrix {
protected:
  DistMatrix& distM_;

public:
  static long M;

  CostMatrix(DistMatrix& distmatrix)  : distM_(distmatrix) {}

  virtual ~CostMatrix() {}

  virtual DistMatrix& distMatrix() const { return distM_; }

  virtual double getCost(Vertex& ori, Vertex& dest) const = 0;

  virtual double getCost(int idOri, int idDest) const = 0;

  virtual double getTravelTime(int idOri, int idDest) const = 0;

  virtual void pruneArcs(DARPEvaluator& patheval, VerticesList& reqlist) = 0;
};

#endif

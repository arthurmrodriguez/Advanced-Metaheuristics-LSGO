#ifndef ONLYDISTCOSTMATRIX_H
#define ONLYDISTCOSTMATRIX_H

#include "VerticesList.h"
#include "DistMatrix.h"
#include "DARPEvaluator.h"
#include "CostMatrix.h"
#include <vector>
#include <assert.h>

class OnlyDistCostMatrix : public CostMatrix {
protected:

  VerticesList& reqlist_;

public:

  OnlyDistCostMatrix(DistMatrix& distmatrix, VerticesList& reqlist) : CostMatrix(distmatrix), reqlist_(reqlist) {}

  virtual ~OnlyDistCostMatrix() {}

  virtual long double getCost(Vertex& ori, Vertex& dest) const {
    return distM_.getDist(ori.pos_,dest.pos_);
  }

  virtual long double getCost(int idOri, int idDest) const {
    Vertex& vori  = reqlist_.getVertex(idOri);
    Vertex& vdest = reqlist_.getVertex(idDest);

    assert(&vori != &vdest);

    return getCost(vori,vdest);
  }

  virtual long double getTravelTime(int idOri, int idDest) const {

    return getCost(idOri,idDest);
  }

  void pruneArcs(DARPEvaluator& patheval, VerticesList& reqlist) {}

};

#endif

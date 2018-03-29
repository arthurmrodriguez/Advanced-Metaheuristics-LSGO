#ifndef DARPCOSTMATRIX_H
#define DARPCOSTMATRIX_H

#include "VerticesList.h"
#include "DistMatrix.h"
#include "DARPEvaluator.h"
#include "CostMatrix.h"
#include <vector>
#include <assert.h>

class DARPCostMatrix : public CostMatrix {
protected: //For the tests we need to include the private here
  std::map<int, std::map<int, double> > costs_;
  std::map<int, std::map<int, double> > times_;

protected:

  void printMatrix(bool onlystats) const;

public:
  static long M;

  DARPCostMatrix(DistMatrix& distmatrix, VerticesList& reqlist);

  virtual ~DARPCostMatrix() {}

  virtual double getCost(Vertex& ori, Vertex& dest) const {
    return getCost(ori.id_,dest.id_);
  }

  virtual double getCost(int idOri, int idDest) const {
    assert(costs_.find(idOri) != costs_.end());
    assert(const_cast<DARPCostMatrix&>(*this).costs_[idOri].find(idDest) != const_cast<DARPCostMatrix&>(*this).costs_[idOri].end());

    return const_cast<DARPCostMatrix&>(*this).costs_[idOri][idDest];
  }

  virtual double getTravelTime(int idOri, int idDest) const {
    assert(times_.find(idOri) != times_.end());
    assert(const_cast<DARPCostMatrix&>(*this).times_[idOri].find(idDest) != const_cast<DARPCostMatrix&>(*this).times_[idOri].end());


    return const_cast<DARPCostMatrix&>(*this).times_[idOri][idDest];
  }

  void pruneArcs(DARPEvaluator& patheval, VerticesList& reqlist);
};

#endif

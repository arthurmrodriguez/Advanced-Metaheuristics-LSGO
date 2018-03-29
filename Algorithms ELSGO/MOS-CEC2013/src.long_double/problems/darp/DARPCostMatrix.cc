#include "DARPCostMatrix.h"

#include "DARPEvaluator.h"

long DARPCostMatrix::M = 999999;

DARPCostMatrix::DARPCostMatrix(DistMatrix& distmatrix, VerticesList& reqlist) : CostMatrix(distmatrix) {

  std::vector<int> ids = reqlist.getVerticesds();

  for (unsigned i = 0; i < ids.size(); i++) {
    for (unsigned j = 0; j < ids.size(); j++) {
      int pos_i = reqlist.getVertex(ids[i]).pos_;
      int pos_j = reqlist.getVertex(ids[j]).pos_;

      costs_[ids[i]][ids[j]] = times_[ids[i]][ids[j]] = distmatrix.getDist(pos_i, pos_j);
    }
  }

}


void DARPCostMatrix::pruneArcs(DARPEvaluator& patheval,VerticesList& reqlist) {

  std::vector<int> pickups    = reqlist.getPickUpIds  ();
  std::vector<int> deliveries = reqlist.getDeliveryIds();
  std::vector<int> all        = reqlist.getVerticesds ();

  // Discard arcs (n+i, i)
  for (unsigned i = 0; i < deliveries.size(); i++) {
    int delivery = deliveries[i];
    int pickup   = Vertex::getSiblingVertId(deliveries[i]);

    costs_[delivery][pickup] = DARPCostMatrix::M;
  }

  // Discard arcs (i, j) with: ei + di + tij > lj
  for (unsigned i = 0; i < all.size(); i++) {
    int vi = all[i];

    //std::cout << " =============================" << std::endl;

    for (unsigned j = 0; j < all.size(); j++) {
      int vj = all[j];

      long ei = reqlist.getVertex(vi).fbegin_;
      long lj = reqlist.getVertex(vj).fend_;

      if (ei + costs_[vi][vj] > lj) { // We assume that di = 0
        //std::cout << "vi: " << vi << "   vj: " << vj << "   pos_i: " << reqlist_.getRequest(vi).pos_ << "   pos_j: " << reqlist_.getRequest(vj).pos_ << "   ei: " << ei << "   tij: " << costs_[vi][vj] << "   lj: " << lj << std::endl;
        costs_[vi][vj] = DARPCostMatrix::M;
      }

    }
  }

  // Discard arcs (i, j) and (j, n+i) if: tij + dj + tjn+i > L
  for (unsigned i = 0; i < pickups.size(); i++) {
    int vi = pickups[i];
    int vni = Vertex::getSiblingVertId(vi);

//    std::cout << " =============================" << std::endl;

    for (unsigned j = 0; j < all.size(); j++) {
      int vj = all[j];


      if (costs_[vi][vj] + costs_[vj][vni] > Vertex::maxRideTime(vi,vni,*this) ) { // We assume that di = 0
//        if (costs_[vi][vj] != DARPCostMatrix::M and costs_[vj][vni] != DARPCostMatrix::M)
//          std::cout << "vi: " << vi << "   vni: " << vni << "   vj: " << vj << "   pos_i: " << reqlist_.getRequest(vi).pos_ << "   pos_j: " << reqlist_.getRequest(vj).pos_ << "   tij: " << costs_[vi][vj] << "   tjn+i: " << costs_[vj][vni] << std::endl;
        costs_[vi][vj] = DARPCostMatrix::M;
        costs_[vj][vni] = DARPCostMatrix::M;
      }

    }

  }

  // Discard arcs (i, n + j) if path P = {j, i, n + j, n + i} is infeasible
  for (unsigned i = 0; i < pickups.size(); i++) {
    int vi = pickups[i];
    int vni = Vertex::getSiblingVertId(vi);

    //std::cout << " =============================" << std::endl;

    for (unsigned j = 0; j < pickups.size(); j++) {
      int vj = pickups[j];
      int vnj = Vertex::getSiblingVertId(vj);

      std::vector<int> P;
      P.push_back(vj);
      P.push_back(vi);
      P.push_back(vnj);
      P.push_back(vni);

      long cost, loadViolation, TWViolation, rideViolation;
      long double pickupDelay, deliveryDelay;

      bool feasible = patheval.evalRoute(P, cost, loadViolation, TWViolation, rideViolation, pickupDelay, deliveryDelay);

      if (!feasible) {
        //std::cout << "  Path: {" << vj << ", " << vi << ", " << vnj << ", " << vni << "} is infeasible" << std::endl;
        costs_[vi][vnj] = DARPCostMatrix::M;
      }
      else {
        //std::cout << "  Path: {" << vj << ", " << vi << ", " << vnj << ", " << vni << "} is feasible" << std::endl;
      }

    }

  }

  // Discard arcs (n + i, j) if path P = {i, n + i, j, n + i} is infeasible
  for (unsigned i = 0; i < pickups.size(); i++) {
    int vi = pickups[i];
    int vni = Vertex::getSiblingVertId(vi);

    //std::cout << " =============================" << std::endl;

    for (unsigned j = 0; j < pickups.size(); j++) {
      int vj = pickups[j];
      int vnj = Vertex::getSiblingVertId(vj);

      std::vector<int> P;
      P.push_back(vi);
      P.push_back(vni);
      P.push_back(vj);
      P.push_back(vnj);

      long cost, loadViolation, TWViolation, rideViolation;
      long double pickupDelay, deliveryDelay;

      bool feasible = patheval.evalRoute(P, cost, loadViolation, TWViolation, rideViolation, pickupDelay, deliveryDelay);

      if (!feasible) {
        //std::cout << "  Path: {" << vi << ", " << vni << ", " << vj << ", " << vnj << "} is infeasible" << std::endl;
        costs_[vni][vj] = DARPCostMatrix::M;
      }
      else {
        //std::cout << "  Path: {" << vi << ", " << vni << ", " << vj << ", " << vnj << "} is feasible" << std::endl;
      }

    }

  }

  // Discard arcs (i, j) if paths P1 = {i, j, n + i, n + j} and P2 = {i, j, n + j, n + i} are both infeasible
  for (unsigned i = 0; i < pickups.size(); i++) {
    int vi = pickups[i];
    int vni = Vertex::getSiblingVertId(vi);

    //std::cout << " =============================" << std::endl;

    for (unsigned j = 0; j < pickups.size(); j++) {
      int vj = pickups[j];
      int vnj = Vertex::getSiblingVertId(vj);

      std::vector<int> P1;
      P1.push_back(vi);
      P1.push_back(vj);
      P1.push_back(vni);
      P1.push_back(vnj);

      long cost, loadViolation, TWViolation, rideViolation;
      long double pickupDelay, deliveryDelay;

      bool feasible1 = patheval.evalRoute(P1, cost, loadViolation, TWViolation, rideViolation, pickupDelay, deliveryDelay);

      std::vector<int> P2;
      P2.push_back(vi);
      P2.push_back(vj);
      P2.push_back(vnj);
      P2.push_back(vni);

      bool feasible2 = patheval.evalRoute(P2, cost, loadViolation, TWViolation, rideViolation, pickupDelay, deliveryDelay);

      //std::cout << "  Path: {" << vi << ", " << vj << ", " << vni << ", " << vnj << "} is " << (feasible1 ? "feasible" : "infeasible") << std::endl;
      //std::cout << "  Path: {" << vi << ", " << vj << ", " << vnj << ", " << vni << "} is " << (feasible2 ? "feasible" : "infeasible") << std::endl;

      if (!feasible1 and !feasible2) {
        costs_[vi][vj] = DARPCostMatrix::M;
        //std::cout << "  *** Discarding arc (" << vi << ", " << vj << ")" << std::endl;
      }

      //std::cout << "  +++++++++++++++++++++++++++++++++" << std::endl;

    }

  }

  // Discard arcs (n + i, n + j) if paths P1 = {i, j, n + i, n + j} and P2 = {j, i, n + i, n + j} are both infeasible
  for (unsigned i = 0; i < pickups.size(); i++) {
    int vi = pickups[i];
    int vni = Vertex::getSiblingVertId(vi);

    //std::cout << " =============================" << std::endl;

    for (unsigned j = 0; j < pickups.size(); j++) {
      int vj = pickups[j];
      int vnj = Vertex::getSiblingVertId(vj);

      std::vector<int> P1;
      P1.push_back(vi);
      P1.push_back(vj);
      P1.push_back(vni);
      P1.push_back(vnj);

      long cost, loadViolation, TWViolation, rideViolation;
      long double pickupDelay, deliveryDelay;

      bool feasible1 = patheval.evalRoute(P1, cost, loadViolation, TWViolation, rideViolation, pickupDelay, deliveryDelay);

      std::vector<int> P2;
      P2.push_back(vj);
      P2.push_back(vi);
      P2.push_back(vni);
      P2.push_back(vnj);

      bool feasible2 = patheval.evalRoute(P2, cost, loadViolation, TWViolation, rideViolation, pickupDelay, deliveryDelay);

      //std::cout << "  Path: {" << vi << ", " << vj << ", " << vni << ", " << vnj << "} is " << (feasible1 ? "feasible" : "infeasible") << std::endl;
      //std::cout << "  Path: {" << vj << ", " << vi << ", " << vni << ", " << vnj << "} is " << (feasible2 ? "feasible" : "infeasible") << std::endl;

      if (!feasible1 and !feasible2) {
        costs_[vni][vnj] = DARPCostMatrix::M;
        //std::cout << "  *** Discarding arc (" << vni << ", " << vnj << ")" << std::endl;
      }

      //std::cout << "  +++++++++++++++++++++++++++++++++" << std::endl;

    }

  }

  // Discard arcs {i, n + i} x {j, n + j} if no path among paths P1..P6 is feasible
  for (unsigned i = 0; i < pickups.size(); i++) {
    int vi = pickups[i];
    int vni = Vertex::getSiblingVertId(vi);

    //std::cout << " =============================" << std::endl;

    for (unsigned j = 0; j < pickups.size(); j++) {
      int vj = pickups[j];
      int vnj = Vertex::getSiblingVertId(vj);

      long cost, loadViolation, TWViolation, rideViolation;
      long double pickupDelay, deliveryDelay;

      std::vector<int> P1;
      P1.push_back(vi);
      P1.push_back(vj);
      P1.push_back(vni);
      P1.push_back(vnj);
      bool feasible1 = patheval.evalRoute(P1, cost, loadViolation, TWViolation, rideViolation, pickupDelay, deliveryDelay);

      std::vector<int> P2;
      P2.push_back(vi);
      P2.push_back(vj);
      P2.push_back(vnj);
      P2.push_back(vni);
      bool feasible2 = patheval.evalRoute(P2, cost, loadViolation, TWViolation, rideViolation, pickupDelay, deliveryDelay);

      std::vector<int> P3;
      P3.push_back(vj);
      P3.push_back(vi);
      P3.push_back(vni);
      P3.push_back(vnj);
      bool feasible3 = patheval.evalRoute(P3, cost, loadViolation, TWViolation, rideViolation, pickupDelay, deliveryDelay);

      std::vector<int> P4;
      P4.push_back(vj);
      P4.push_back(vi);
      P4.push_back(vnj);
      P4.push_back(vni);
      bool feasible4 = patheval.evalRoute(P4, cost, loadViolation, TWViolation, rideViolation, pickupDelay, deliveryDelay);

      std::vector<int> P5;
      P5.push_back(vi);
      P5.push_back(vni);
      P5.push_back(vj);
      P5.push_back(vnj);
      bool feasible5 = patheval.evalRoute(P5, cost, loadViolation, TWViolation, rideViolation, pickupDelay, deliveryDelay);

      std::vector<int> P6;
      P6.push_back(vj);
      P6.push_back(vnj);
      P6.push_back(vi);
      P6.push_back(vni);
      bool feasible6 = patheval.evalRoute(P6, cost, loadViolation, TWViolation, rideViolation, pickupDelay, deliveryDelay);

//       std::cout << "  Path: {" << vi << ", " << vj << ", " << vni << ", " << vnj << "} is " << (feasible1 ? "feasible" : "infeasible") << std::endl;
//       std::cout << "  Path: {" << vi << ", " << vj << ", " << vnj << ", " << vni << "} is " << (feasible2 ? "feasible" : "infeasible") << std::endl;
//       std::cout << "  Path: {" << vj << ", " << vi << ", " << vni << ", " << vnj << "} is " << (feasible3 ? "feasible" : "infeasible") << std::endl;
//       std::cout << "  Path: {" << vj << ", " << vi << ", " << vnj << ", " << vni << "} is " << (feasible4 ? "feasible" : "infeasible") << std::endl;
//       std::cout << "  Path: {" << vi << ", " << vni << ", " << vj << ", " << vnj << "} is " << (feasible5 ? "feasible" : "infeasible") << std::endl;
//       std::cout << "  Path: {" << vj << ", " << vnj << ", " << vi << ", " << vni << "} is " << (feasible6 ? "feasible" : "infeasible") << std::endl;

      if (!feasible1 and !feasible2 and !feasible3 and !feasible4 and !feasible5 and !feasible6) {
        costs_[vi][vj] = DARPCostMatrix::M;
        costs_[vi][vnj] = DARPCostMatrix::M;
        costs_[vj][vi] = DARPCostMatrix::M;
        costs_[vj][vni] = DARPCostMatrix::M;
        costs_[vni][vj] = DARPCostMatrix::M;
        costs_[vni][vnj] = DARPCostMatrix::M;
        costs_[vnj][vi] = DARPCostMatrix::M;
        costs_[vnj][vni] = DARPCostMatrix::M;
        //std::cout << "  *** Discarding arcs (" << vi << ", " << vj << "), (" << vi << ", " << vnj << "), (" << vj << ", " << vi << "), (" << vj << ", " << vni << "), (";
        //std::cout << vni << ", " << vj << "), (" << vni << ", " << vnj << "), (" << vnj << ", " << vi << "), (" << vnj << ", " << vni << ")" << std::endl;
      }

      //std::cout << "  +++++++++++++++++++++++++++++++++" << std::endl;

    }

  }

  printMatrix(true);

}


void DARPCostMatrix::printMatrix(bool onlystats) const {

  std::map<int, std::map<int, long double> >::const_iterator it1;
  std::map<int, long double>::const_iterator it2;

  int overallVertices = 0;
  int penalizedVertices = 0;

  for (it1 = costs_.begin(); it1 != costs_.end(); it1++) {
    for (it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
      overallVertices++;
      if (it2->second == DARPCostMatrix::M)
        penalizedVertices++;
      if (!onlystats)
        std::cout << it2->second << " ";
    }
    if (!onlystats)
      std::cout << std::endl;
  }

  std::cout << "  + Overall vertices: " << overallVertices << std::endl;
  std::cout << "  + Penalized vertices: " << penalizedVertices << std::endl;
  std::cout << "  + Ratio: " << (long double) penalizedVertices / (long double) overallVertices << std::endl;

}

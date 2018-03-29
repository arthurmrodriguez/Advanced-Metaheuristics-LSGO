#include "gtest/gtest.h"

// Hack for testing private members
#define protected public
#define private public
#include "../darpOps.cc"
#undef protected
#undef private

#include "../DARPGenome.cc"
#include "../DistMatrix.cc"
#include "../DARPCostMatrix.cc"
#include "../VerticesList.cc"
#include "../DARPEvaluator.cc"
#include <garandom.h>
#include <vector>
#include <list>
#include <math.h>

using namespace std;

// Set up for the tests
//
int MAXPOINTS       = 100;
int MAXUSERRIDETIME = 60;
int MAXUSERWAITTIME = 20;
int VEHICLECAPACITY = 18;
int SUMTIMEEXP      = 2;
int TESTLOADV       = 10;
int TESTTWV         = 10;
int TESTRIDEV       = 10;

struct DummyCostPoint {
  int   orig;
  int   dest;
  long double cost;

  DummyCostPoint(int o, int d, long double c) : orig(o), dest(d), cost(c) {}

};

struct DARPTestsVars {
  static VerticesList* verticeslist;
  static CostMatrix*   costmatrix;
  static DARPEvaluator* darpeval;

  static void initData(VerticesList* vert, CostMatrix* cost) ;
  static void initDataOriginalEval(VerticesList* vert, CostMatrix* cost, int vcapac, int sumtimeexp, int maxuserridetime, int maxuserwaittime, int twv, int loadv, int ridev, DARPEvaluator::darpOptCriterionType optcrit);
};

VerticesList* DARPTestsVars::verticeslist = 0;
CostMatrix* DARPTestsVars::costmatrix     = 0;
DARPEvaluator* DARPTestsVars::darpeval    = 0;

long double dummyRouteScore(const DARPGenome& gen, int route) {
  CostMatrix&   costmatrix   = * DARPTestsVars::costmatrix;
  VerticesList& requestslist = * DARPTestsVars::verticeslist;

  long double score = 0;
  for (int i=0; i<gen.routeLength(route)-1; i++) {
    Vertex& req1 = requestslist.getVertex( gen.gene(route,i) );
    Vertex& req2 = requestslist.getVertex( gen.gene(route,i+1) );
    score += costmatrix.getCost(req1,req2);
  }

  return score;
}

extern "C" long double DummyObjFunc(GAGenome& g) {
  DARPGenome& gen = dynamic_cast<DARPGenome&>(g);
  long double score = 0.0;

  for (int route=0; route<gen.length(); route++){
    score += dummyRouteScore(gen,route);
  }

  return score;
}

extern "C" long double DummyObjFuncWithCriticalVertexConstraint(GAGenome& g) {
  long double score = DummyObjFunc(g);
  DARPGenome& gen = dynamic_cast<DARPGenome&>(g);

  VerticesList& requestslist = *DARPTestsVars::verticeslist;

  const int NONFEASIBLE_penalization = 10;

  // For each infeasible arc (only that a non critical vertex comes before than a critical one for this case) a penalization is added
  for (int route=0; route<gen.length(); route++){
    for (int i=0; i<gen.routeLength(route); i++) {
      Vertex& req = requestslist.getVertex( gen.gene(route,i) );
      if (req.critic_) {
        int noncrit_vert_id  = req.getSiblingVertexId();
        int noncrit_vert_pos = gen.findPosOfVertex(route,noncrit_vert_id);
        if (noncrit_vert_pos < i ) score += NONFEASIBLE_penalization;
      }
    }
  }

  return score;
}

class DummyDARPEvaluator : public DARPEvaluator {
public:
  DummyDARPEvaluator(const long vehicleCapacity,
                     const VerticesList& vertlist, const CostMatrix& costMatrix,
                     int timeSumExp,
                     long double loadvW, long double TWVW, long double rideVW,
                     darpOptCriterionType optcrit) : DARPEvaluator(vehicleCapacity, vertlist, costMatrix,timeSumExp,loadvW,TWVW,rideVW,optcrit) {}

  long double nonPenalizedScore(const DARPGenome& gen, int route) const {
    return dummyRouteScore(gen,route);
  }

  long double score(DARPGenome& gen) {
    return DummyObjFunc(gen);
  }

};

void DARPTestsVars::initData(VerticesList* vert, CostMatrix* cost) {
  verticeslist = vert; costmatrix = cost;

  if (darpeval != 0) delete darpeval;
  darpeval = new DummyDARPEvaluator(VEHICLECAPACITY,*vert, *cost,SUMTIMEEXP,TESTLOADV, TESTTWV, TESTRIDEV, DARPEvaluator::TRAVELCOST);
  DARPGenome::setDARPEvaluator(darpeval);
  DARPVNSOp::setDARPEvaluator(darpeval);

  Vertex::setTimeConstants(false,0,MAXUSERRIDETIME,MAXUSERWAITTIME);
}

void DARPTestsVars::initDataOriginalEval(VerticesList* vert, CostMatrix* cost, int vcapac, int sumtimeexp, int maxuserridetime, int maxuserwaittime, int twv, int loadv, int ridev, DARPEvaluator::darpOptCriterionType optcrit) {
  verticeslist = vert; costmatrix = cost;

  if (darpeval != 0) delete darpeval;

  darpeval = new DARPEvaluator(vcapac,*vert, *cost,sumtimeexp,loadv, twv, ridev, optcrit);
  DARPGenome::setDARPEvaluator(darpeval);
  DARPVNSOp::setDARPEvaluator(darpeval);

  Vertex::setTimeConstants(false,0,maxuserridetime,maxuserwaittime);
}

extern "C" long double originalobjective(GAGenome& g) {
  DARPGenome& gen = dynamic_cast<DARPGenome&>(g);

  return DARPTestsVars::darpeval->score(gen);
}

struct DummyDistMatrix : public DistMatrix {

  DummyDistMatrix(int numpoints)   {
    data_.resize(numpoints);
    for (int i=0; i<numpoints; i++) data_[i].resize(numpoints);
  }

  DummyDistMatrix(int numpoints, vector<DummyCostPoint>& costinf) {
    data_.resize(numpoints);
    for (int i=0; i<numpoints; i++) data_[i].resize(numpoints);

    for (int i=0; i<costinf.size(); i++) {
      DummyCostPoint& inf = costinf[i];
      data_[inf.orig][inf.dest] = data_[inf.dest][inf.orig] = inf.cost;
    }
  }

  virtual ~DummyDistMatrix() {}

  int size() { return data_.size(); }

  void addPos(int pos) {
    if (pos >= data_.size()) {
      data_.resize(pos+1);
      for (int i=0; i<=pos; i++) data_[i].resize(pos+1);
    }
  }

  void setCost(int orig, int dest, long double cost) { 
    assert(orig >= 0 && orig < data_.size() );
    assert(dest >= 0 && dest < data_[orig].size() );
    data_[orig][dest] = cost;
    data_[dest][orig] = cost;
  }

  void setCostOneSide(int orig, int dest, long double cost) {
    assert(orig >= 0 && orig < data_.size() );
    assert(dest >= 0 && dest < data_[orig].size() );
    data_[orig][dest] = cost;
  }


  long double getCost(int orig, int dest) {
    return data_[orig][dest];
  }

  long double getDist(int orig, int dest) const {
    assert(orig < data_.size() && dest < data_[0].size() );
    return data_[orig][dest];
  }

};

struct DummyCostMatrix : public DARPCostMatrix {

  DummyCostMatrix(DummyDistMatrix* distmatrix, VerticesList& reqlist) : DARPCostMatrix(*distmatrix,reqlist) {}

  ~DummyCostMatrix() {}

  void setCost(int orig, int dest, long double cost) {
    costs_[orig][dest] = cost;
    costs_[dest][orig] = cost;
    /*LOG*/ assert(costs_.find(orig) != costs_.end());
    /*LOG*/ assert(const_cast<DummyCostMatrix&>(*this).costs_[orig].find(dest) != const_cast<DummyCostMatrix&>(*this).costs_[orig].end());

  }

  void setCostOneSide(int orig, int dest, long double cost) {
    costs_[orig][dest] = cost;
  }


};

void addRequestRoute(VerticesList& list, int id,int origin, int dest, int wbegin, int wend, vector<DummyCostPoint>& costinf, long double cost) {
  list.addVertex( new Vertex(id   ,origin,true,  1, wbegin,               wend,                 Vertex::PICKUP) );
  list.addVertex( new Vertex(id*-1,dest,  false, 1, wend+MAXUSERRIDETIME, wend+MAXUSERRIDETIME, Vertex::DELIVERY) );

  costinf.push_back( DummyCostPoint(origin,dest,cost) );
}

void addRequest(VerticesList& list, int id, int origin, int dest, int wbegin, DummyDistMatrix& distmatrix) {
  int wend = wbegin+MAXUSERWAITTIME;
  list.addVertex( new Vertex(id   ,origin,true,  1, wbegin,                                      wend                , Vertex::PICKUP) );
  list.addVertex( new Vertex(id*-1,dest,  false, 1, wbegin+(int)distmatrix.getCost(origin,dest), wend+MAXUSERRIDETIME, Vertex::DELIVERY) );
}

void printRoute(DARPGenome& gen, int route, VerticesList& list, CostMatrix& costmatrix) {
  cout << "Gen:  ";
  for (int i=0; i<gen.routeLength(route); i++) {
    int genv = gen.gene(route,i);
    int req_id = list.getVertex( genv ).pos_;
    cout << genv << "(" << req_id << ") ";
  }
  cout << endl;
  
  for (int i=0; i<gen.routeLength(route)-1; i++) {
    Vertex& req      = list.getVertex( gen.gene(route,i) );
    Vertex& next_req = list.getVertex( gen.gene(route,i+1) );
    long double cost     = costmatrix.getCost(req,next_req);

    cout << req.pos_ << " -" <<  cost << "- " << next_req.pos_ << "   ";
  }
  cout << endl << endl;
}

void printNatSeq(list<int>* natseq) {
  for (list<int>::iterator it=natseq->begin(); it!=natseq->end(); it++) {
    cout << *it << ",";
  }
  cout << " ";
}

void printNatSeq(vector< list<int>* >& route_natseqs, int natseqpos) {
  list<int>* natseq = route_natseqs[natseqpos];
  printNatSeq(natseq);
}

void printRouteNatSeqs( vector< list<int>* >& route_natseqs) {
  for (int j=0; j<route_natseqs.size(); j++) { printNatSeq(route_natseqs,j); cout << "   "; }
  cout << endl;
}

void printNatSeqs( vector< vector< list<int>* > >& natseqs) {
  for (int i=0; i<natseqs.size(); i++) {
    vector< list<int>* >& route_natseqs = natseqs[i];
    cout << "route " << i << " : "; 
    printRouteNatSeqs(route_natseqs);
  }
}

void printNatSeqsTimes( vector< vector< list<int>* > >& natseqs) {
  for (int i=0; i<natseqs.size(); i++) {
    vector< list<int>* >& route_natseqs = natseqs[i];
    cout << "route " << i << " : "; 
    for (int j=0; j<route_natseqs.size(); j++) { 
      list<int>* natseqs = route_natseqs[j];
      int id_front = natseqs->front();
      int id_back  = natseqs->back();
     
      VerticesList& requestslist = * DARPTestsVars::verticeslist;
      Vertex& vfront = requestslist.getVertex(id_front);
      Vertex& vback  = requestslist.getVertex(id_back);

      cout << "NatSeq " << j << ": " << vfront.fbegin_ << ":" << vfront.fend_  << " " << vback.fbegin_ << ":" << vback.fend_; cout << "   "; 
    }
    cout << endl;
  }
}


void recomputeWindowsNonCriticalVertex(VerticesList& list, DARPGenome& gen, int route, CostMatrix& costmatrix) {
  for (int i=0; i<gen.routeLength(route)-1; i++) {
    Vertex& req = list.getVertex( gen.gene(route,i) );

    if (!req.critic_) {

      Vertex* next_critic = 0; 
      for (int j=i+1; j<gen.routeLength(route); j++) {
        Vertex& tmp = list.getVertex( gen.gene(route,j) );
        if (tmp.critic_) {
          next_critic = &tmp;
          break;
        }
      }

      int bias_value1 = 8;
      int bias_value2 = 2;

      if ( next_critic ) { // found a next critical
        req.fend_   = next_critic->fend_ - (int) costmatrix.getCost(req,*next_critic) - bias_value1;
        req.fbegin_ = req.fend_ - bias_value2;
      }
    }

  }
}

void setUpTestingEnvironment() {


}




DARPGenome* createGen(int num_routes) {
  DARPGenome* gen = new DARPGenome(num_routes,0,DummyObjFunc);
  return gen;
}

DARPGenome* createGen(int num_routes,int route, vector<int>& seqv) {
  DARPGenome* gen = new DARPGenome(num_routes,0,DummyObjFunc);

  for (int i=0; i<seqv.size(); i++) {
    gen->pushBackVertex(route  ,seqv[i]);

    for (int j=1; j<num_routes; j++) {
      int next_route = (route+j) % num_routes;
      gen->pushBackVertex( next_route , (int) pow((long double)10,j) *seqv[i]);
    }
  }

  return gen;
}

DARPGenome* createGen(vector< vector<int> >& values) {
  int num_routes = values.size();

  DARPGenome* gen = new DARPGenome(num_routes,0,DummyObjFunc);
  for (int route=0; route<values.size(); route++) {
    vector<int>& route_values = values[route];
    for (int pos=0; pos<route_values.size(); pos++) {
      gen->pushBackVertex(route,route_values[pos]);
    }
  }

  return gen;
}

void setRoute(DARPGenome& gen, int route, vector<int>& seqv) {
  assert(route >=0 && route <gen.size());
  for (int i=0; i<seqv.size(); i++) {
    gen.pushBackVertex(route  ,seqv[i]);
  }
}

bool areRoutesEqual(DARPGenome& gen, int route1, int route2) {
  assert(route1 >= 0 and route1 < gen.size());
  assert(route2 >= 0 and route2 < gen.size());

  bool result = gen.routeLength(route1) == gen.routeLength(route2);
  if (result) {
    for (int i=0; i<gen.routeLength(route1); i++) {
      if (gen.gene(route1,i) != gen.gene(route2,i)) {
        result = false;
        break;
      }
    }
  }
  return result;
}

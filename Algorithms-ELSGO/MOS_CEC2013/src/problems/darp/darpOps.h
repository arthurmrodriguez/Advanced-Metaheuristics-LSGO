#ifndef DARPOPS_H
#define DARPOPS_H

#include "DARPEvaluator.h"
#include "darpOpsActions.h"
#include <VNS.h>
#include <limits>


class GAGenome;

// Little hack so both the DARPVNSShaker and DARPVNSLocalSearch inherit the darpeval class attribute
class DARPVNSOp  {
protected:
  static DARPEvaluator* darpeval__;
public:
  static void setDARPEvaluator(DARPEvaluator* darpeval) { darpeval__ = darpeval;}
};

class DARPVNSShaker : public VNSShaker, public DARPVNSOp {
protected:

  void getRandomNonEmptyRoutes(const DARPGenome& gen, int& route1, int& route2);
  void getRandomBadGoodRoutes (const DARPGenome& gen, /*out*/ int& badroute, /*out*/ int& goodroute);

  static void       getAllNaturalSeqs (const DARPGenome& gen, int route, vector< list<int>* >& natseqs) ;
  static list<int>* getVerticesIdsList(const DARPGenome& gen, int route, int start_pos, int end_pos) ;

  static list<int> extractRandomRouteSeqVertices(DARPGenome& gen, int route, int neighborhood_size) ;
  static void      extractMissingSiblings       (DARPGenome& gen, int route, list<int>& vertices_ids) ;

  static int  bestSeqPosToBeRemoved      (const DARPGenome& gen, int route, const vector< list<int>* >& allseqs);
  static int  bestSeqPosToBeInserted     (const DARPGenome& gen, int route, const vector< list<int>* >& allseqs);
  static int  worstSeqPosToBeRemoved     (const DARPGenome& gen, int route, const vector< list<int>* >& allseqs);
  static int  worstSeqPosToBeInserted    (const DARPGenome& gen, int route, const vector< list<int>* >& allseqs);
  static void bestAndWorstSeqToBeRemoved (const DARPGenome& gen, int route, const vector< list<int>* >& allseqs, /*out*/ int &bestpos, /*out*/ int &worstpos);
  static void bestAndWorstSeqToBeInserted(const DARPGenome& gen, int route, const vector< list<int>* >& allseqs, /*out*/ int &bestpos, /*out*/ int &worstpos);

  static vector< list<int>* >* getAllSequencesOfSize(const DARPGenome& gen, int route, int max_size, bool userandomsize=true);

  static int selectPos(vector<double>& probs);

  static string stringOfSeq(list<int>& seq);

  struct RouteAndScore { int route; double score;};

  static vector<RouteAndScore> getRouteAndScoresSorted          (const DARPGenome& gen, GAGenome::OptCriterion optcrit=GAGenome::MAXIMIZATION);
  static bool                  descendingRouteAndScoreOrderFunc (const RouteAndScore& a, const RouteAndScore& b);
  static vector<double>        getSelectionProbs                (const vector<RouteAndScore>& routescores, GAGenome::OptCriterion optcrit=GAGenome::MAXIMIZATION);

  // Instance methods
  list<int> bestSeqToBeRemoved (const DARPGenome& gen, int route);

public:
  static DARPVNSShaker* createNewShaker(string name, int size=0);

  DARPVNSShaker(string name, int size) : VNSShaker(name, size) {}

  virtual ~DARPVNSShaker() {}
};

class SwapNeighborhood : public DARPVNSShaker { 
public:

  SwapNeighborhood(int size) : DARPVNSShaker("SwapNeighborhood", size) {}

  virtual VNSOp* clone() {
    return new SwapNeighborhood(*this);
  }

  virtual VNSOpAction* operator() (GAGenome& gen);

};

class ChainNeighborhood : public DARPVNSShaker { 
public:

  ChainNeighborhood(int size) : DARPVNSShaker("ChainNeighborhood", size) {}

  virtual VNSOp* clone() {
    return new ChainNeighborhood(*this);
  }

  virtual VNSOpAction* operator() (GAGenome& gen);

};

class GreedyMoveNeighborhood : public DARPVNSShaker { 
protected:

  static int routeForExtraction       (const DARPGenome& gen);
  static int bestRoutePosForInsertion (const DARPGenome& gen, int orig_route, double orig_score, const list<int>& vertices);

public:

  GreedyMoveNeighborhood(int size) : DARPVNSShaker("GreedyMoveNeighborhood", size) {
    if (size % 2 == 1) cout << "GreedyMoveNeighborhood odd size=" << size << " will be round up to the previous even size" << endl;
  }

  virtual VNSOp* clone() {
    return new GreedyMoveNeighborhood(*this);
  }

  virtual VNSOpAction* operator() (GAGenome& gen);
};

class GreedyMoveNeighborhoodDestCentered : public DARPVNSShaker {
protected:

  list<int> bestSeqToBeRemovedAndInserted(const DARPGenome& gen, int fromroute, int toroute);

public:

  GreedyMoveNeighborhoodDestCentered(int size) : DARPVNSShaker("GreedyMoveNeighborhoodDestCentered", size) {
    if (size % 2 == 1) cout << "GreedyMoveNeighborhoodDestCentered odd size=" << size << " will be round up to the previous even size" << endl;
  }

  virtual VNSOp* clone() {
    return new GreedyMoveNeighborhoodDestCentered(*this);
  }

  virtual VNSOpAction* operator() (GAGenome& gen);

};

class GreedySwapNeighborhood : public DARPVNSShaker {
protected:

  enum swapopts {bestbest=0, bestworst, worstbest, worstworst, lastopt};

  vector<int> bestswaps_num_;

  virtual void getBestAndWorstSeqsToRemove(const DARPGenome& gen, int fromroute, int toroute, /*out*/ list<int>& best_seq, /*out*/ list<int>& worst_seq);

  virtual double scoreOfSwappingSeqs(const DARPGenome& gen, int route1, list<int>& seq1, int route2, list<int>& seq2);

public:

  GreedySwapNeighborhood(int size, string name="GreedySwapNeighborhood") : DARPVNSShaker(name, size),  bestswaps_num_(lastopt) {
    if (size % 2 == 1) cout << "GreedySwapNeighborhood odd size=" << size << " will be round up to the previous even size" << endl;
  }

  virtual VNSOp* clone() {
    return new GreedySwapNeighborhood(*this);
  }

  virtual string getStatsResults();

  virtual VNSOpAction* operator() (GAGenome& gen);

  void getBestSwap(const DARPGenome& gen,
                   int               route1, list<int>& best_seq1, list<int>& worst_seq1,
                   int               route2, list<int>& best_seq2, list<int>& worst_seq2,
                   /*out*/ list<int>& sel_seq1, /*out*/ list<int>& sel_seq2);

};

class GreedySwapNeighborhoodDestCentered : public GreedySwapNeighborhood {

protected:

  void getBestAndWorstSeqsToRemove(const DARPGenome& gen, int fromroute, int toroute, /*out*/ list<int>& best_seq, /*out*/ list<int>& worst_seq);

public:

  GreedySwapNeighborhoodDestCentered(int size) : GreedySwapNeighborhood(size,"GreedySwapNeighborhoodDestCentered") {
    if (size % 2 == 1) cout << "GreedySwapNeighborhoodDestCentered odd size=" << size << " will be round up to the previous even size" << endl;
  }

  virtual VNSOp* clone() {
    return new GreedySwapNeighborhoodDestCentered(*this);
  }
};


class ZeroSplitNeighborhood : public DARPVNSShaker { 
public:

  ZeroSplitNeighborhood(int size=0) : DARPVNSShaker("ZeroSplitNeighborhood", size) {}

  virtual VNSOp* clone() {
    return new ZeroSplitNeighborhood(*this);
  }

  virtual VNSOpAction* operator() (GAGenome& gen);
};

class CheckAllNaturalSeqsCombsNeighborhood : public DARPVNSShaker { 
protected: 

  double prev_score_;  

  static void bestSwapOfNatSeqs(const DARPGenome&                     gen,
                                const vector< vector< list<int>* > >& allnatseqs,
                                int                                   fromroute,
                                int                                   toroute,
                                int&    /*out*/                       bestfrom_pos,
                                int&    /*out*/                       bestto_pos,
                                double& /*out*/                       swap_score);
  
  static double scoreFromSwapNaturalSequences(const DARPGenome&               gen,
                                              vector< vector< list<int>* > >& allnatseqs,
                                              int                             fromroute,
                                              int                             from_natseq_pos,
                                              int                             toroute,
                                              int                             to_natseq_pos);

  static bool isSwapFeasible(vector< vector< list<int>* > >& allnatseqs,
                             int                             fromroute,
                             int                             from_natseq_pos,
                             int                             toroute,
                             int                             to_natseq_pos) ;

  static bool canNatSeqReplaceNatSeq(vector< vector< list<int>* > >& allnatseqs,
                                     int                             fromroute,
                                     int                             from_natseq_pos,
                                     int                             toroute,
                                     int                             to_natseq_pos);

  static void getPrevAndNextNatSeqs(vector< vector< list<int>* > >& allnatseqs, 
                                    int                             route, 
                                    int                             pos,
                                    list<int>*&                     prev_natseq,
                                    list<int>*&                     next_natseq) ;

  static SwapWithInsertionPosDatum swapNaturalSequences(DARPGenome&                     gen,
                                                        vector< vector< list<int>* > >& allnatseqs,
                                                        int                             fromroute,
                                                        int                             from_natseq_pos,
                                                        int                             toroute,
                                                        int                             to_natseq_pos);

  static int computeInsPosBasedOnTheNatSeqPos(const vector< list<int>* >& route_natseqs, int natseq_pos) ;

public:

  CheckAllNaturalSeqsCombsNeighborhood(int size=0) : DARPVNSShaker("CheckAllNaturalSeqsCombsNeighborhood", size) {
    prev_score_ = GAGenome::worstPossibleScore();
  }

  virtual VNSOp* clone() {
    return new CheckAllNaturalSeqsCombsNeighborhood(*this);
  }

  virtual VNSOpAction* operator() (GAGenome& gen);

  virtual bool isStochastic() { return false; }
};

class DARPLocalSearch : public VNSOp, public DARPVNSOp  {
protected:
  long maxevalspercall_;
  bool randomroutes_;

  int localSearch(DARPGenome& gen, int route, int maxevals=numeric_limits<int>::max());

  static void localSearchToCritVertexPos      (DARPGenome& gen, int route, int critvert_origpos, int critvert_searchpos,
                                               long maxevals,
                                               /*out*/ int& newpos, /*out*/ bool& unfinished);
  static int  firstCriticalVertexPos          (DARPGenome& gen, int route, int start_pos);
  static void insertNonCriticalVertexAtBestPos(DARPGenome& gen, int route, int critvert_pos);

public:

  DARPLocalSearch(long maxevals=numeric_limits<int>::max(), bool randomroutes=false) : VNSOp("DARPLocalSearch"), maxevalspercall_(maxevals),
                                                                                       randomroutes_(randomroutes) {
    assert(maxevalspercall_ > 0);
  }

  virtual VNSOp* clone() {
    return new DARPLocalSearch(*this);
  }

  virtual VNSOpAction* operator() (GAGenome& gen);
};

class DARPNOLocalSearch : public VNSOp, public DARPVNSOp  {

public:

  DARPNOLocalSearch() : VNSOp("DARPNOLocalSearch") {}

  virtual VNSOp* clone() {
    return new DARPNOLocalSearch(*this);
  }

  virtual VNSOpAction* operator() (GAGenome& gen) { return 0; }
};









#endif

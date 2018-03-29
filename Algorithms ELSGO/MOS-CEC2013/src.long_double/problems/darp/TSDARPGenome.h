#ifndef TSDARPGENOME_H
#define TSDARPGENOME_H

#include "DARPGenome.h"
#include <genomes/GAGenome.h>
#include <vector>
#include <utility>
#include <map>
#include <string>

using namespace std;

class TSDARPGenome : public DARPGenome {
private:

  long double evalWeightsUpdateValue_; // delta  in Cordeau's paper
  long double lackOfDivPenalization_;  // lambda in Cordeau's paper
  long double maxTabuLife_;            // teta   in Cordeau's paper

  map< pair<int,int> , int > tabuMem_;
  map< pair<int,int> , int > opFreq_; // rho in Cordeau's paper

  static pair<int,int> routeVertIdKey(int newroute, int vertid) {return pair<int,int> (newroute,vertid);}

  void   moveToNeighbor (int oldroute, int newroute, int vertid);
  bool   allowedMove    (int oldroute, int newroute, int vertid);
  long double divPenalization(int newroute, int vertid);

  // Adds maxTabuLife_ plus one since at the end of the step it is decremented by one and this cannot be avoided for the
  // first time the tabu life is set
  void addTabuMovement(pair<int,int>& key) { tabuMem_[key] = (int) round(maxTabuLife_); }
  void increaseFreq   (pair<int,int>& key) { opFreq_ [key] += 1; }

  void ageMemory();
  void updateTSMemsWithMov(int newroute,int vertid);

public:

  TSDARPGenome (unsigned x, GAGenome::Evaluator f = (GAGenome::Evaluator) 0, void* u = (void*) 0);

  TSDARPGenome (const TSDARPGenome& orig) : DARPGenome (orig) { copy(orig); }

  TSDARPGenome& operator= (const GAGenome& arr) {copy(arr); return *this;}

  GAGenome* clone (GAGenome::CloneMethod flag = GAGenome::CONTENTS) const {return new TSDARPGenome (*this);}

  virtual void copy (const GAGenome& orig);

  virtual int  write      (STD_OSTREAM& os) const;
  virtual void writeObject(ostream&     os) const;
  virtual void readObject (istream& is);


  // Tabu Algorithm Specific Methods

  virtual pair<RoutingGenome*,int> bestFromNeighborhood();

  virtual void updateDynParameters();
  virtual void updatePenalizations();

  virtual string memoryDataInformation();
};


#endif

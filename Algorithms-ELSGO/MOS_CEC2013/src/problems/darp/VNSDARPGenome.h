#ifndef VNSDARPGENOME_H
#define VNSDARPGENOME_H

#include "DARPGenome.h"
#include <vector>
#include <utility>
#include <map>
#include <string>

using namespace std;

class VNSDARPGenome : public DARPGenome {
private:
  double minEvalWeightsUpdateValue_; // Min possible value for evalWeightsUpdateValue_
  double maxEvalWeightsUpdateValue_; // Max possible value for evalWeightsUpdateValue_

public:

  VNSDARPGenome (unsigned x,
                 double minevalweightsvalue, double maxevalsweightsvalue,
                 GAGenome::Evaluator f = (GAGenome::Evaluator) 0, void* u = (void*) 0);

  VNSDARPGenome (const VNSDARPGenome& orig) : DARPGenome (orig) { copy(orig); }

  VNSDARPGenome& operator= (const GAGenome& arr) {copy(arr); return *this;}

  GAGenome* clone (GAGenome::CloneMethod flag = GAGenome::CONTENTS) const {return new VNSDARPGenome (*this);}

  virtual void copy (const GAGenome& orig);

  virtual int  write      (STD_OSTREAM& os) const;
  virtual void writeObject(ostream&     os) const;
  virtual void readObject (istream& is);

  virtual void updatePenalizations();
};


#endif

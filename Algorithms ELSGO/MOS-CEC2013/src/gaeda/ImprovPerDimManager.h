#ifndef IMPROVPERDIMMANAGER_H_
#define IMPROVPERDIMMANAGER_H_

#include <vector>

using namespace std;

class ImprovPerDimManager {
protected:

  static ImprovPerDimManager* instance__;

  int num_dims_;

  int numDims() const { return num_dims_; }

  vector<double>*         getProbs        (                     ) const;
  vector<double>*         getAccumProbs   (                     ) const;
  virtual vector<double>  dimPerfValues   (                     ) const;
  int                     posBasedOnProbs (vector<double>& probs) const;

public:

  ImprovPerDimManager(int num_dims);
  virtual ~ImprovPerDimManager();

  void getDims                      (vector<int>& dims                                                  ) const;
  void getDimsWithRep               (vector<int>& dims                                                  ) const;
  void getDimsWhichProbSumAThreshold(vector<int>& selected_dims, double prob_percent, int num_other_dims) const;

  //void         updateValues (vector<double>& new_values);
  virtual void updateValue  (int pos, double new_value );

  static ImprovPerDimManager* instance() { return instance__; }

};

#endif /* IMPROVPERDIMMANAGER_H_ */

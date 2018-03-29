#ifndef NSC_H
#define NSC_H

class GAPopulation;

class NSC {

 public:

  NSC (unsigned size);
  virtual ~NSC ();

  void   addPopulation  (const GAPopulation* pop);
  void   prepareNSCData (const GAPopulation* pop);
  double computeNSC     (const GAPopulation* pop, int tech);

 protected:

  //NSC structure node
  typedef struct {
    unsigned long int id;
    int island;
    double fit;
  } nscdata;

  unsigned _sz;
  unsigned _nscn;
  nscdata* _nscd;

};

#endif

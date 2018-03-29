#ifndef FSSGENOME_H
#define FSSGENOME_H

#include "GA1DArrayGenome.h"

template <class T> class FSSGenome : public GA1DArrayAlleleGenome<T> {

public:

  FSSGenome (unsigned x, const GAAlleleSet<T>& a,
             GAGenome::Evaluator f = (GAGenome::Evaluator) 0, void* u = (void*) 0)
            : _AUC (0), _LL (0), _LLS (0), _BS (0), _ACC (0),  GA1DArrayAlleleGenome<T> (x, a, f, u) {}

  FSSGenome (const GAAlleleSetArray<T>& a,
             GAGenome::Evaluator f = (GAGenome::Evaluator) 0, void* u = (void*) 0)
             : _AUC (0), _LL (0), _LLS (0), _BS (0), _ACC (0),  GA1DArrayAlleleGenome<T> (a, f, u) {}

  FSSGenome (const FSSGenome<T>& orig) : _AUC (0), _LL (0), _LLS (0), _BS (0), _ACC (0),  GA1DArrayAlleleGenome<T> (orig) {}

  FSSGenome<T>& operator= (const GAGenome& arr) {copy(arr); return *this;}

  GAGenome* clone (GAGenome::CloneMethod flag = GAGenome::CONTENTS) const {return new FSSGenome<T> (*this);}

  long double getAUC () const {return _AUC;}
  long double setAUC (long double auc) {return _AUC = auc;}

  long double getLL  () const {return _LL;}
  long double setLL  (long double ll) {return _LL = ll;}

  long double getLLS () const {return _LLS;}
  long double setLLS (long double lls) {return _LLS = lls;}

  long double getBS  () const {return _BS;}
  long double setBS  (long double bs) {return _BS = bs;}

  long double getACC () const {return _ACC;}
  long double setACC (long double acc) {return _ACC = acc;}

  virtual void copy (const GAGenome& orig) {

    GA1DArrayAlleleGenome<T>::copy (orig);

    const FSSGenome<T>& ori = DYN_CAST (const FSSGenome<T>&, orig);

    _AUC = ori._AUC;
    _LL  = ori._LL;
    _LLS = ori._LLS;
    _BS  = ori._BS;
    _ACC = ori._ACC;

  }

protected:

  long double _AUC;
  long double _LL;
  long double _LLS;
  long double _BS;
  long double _ACC;

};

#endif

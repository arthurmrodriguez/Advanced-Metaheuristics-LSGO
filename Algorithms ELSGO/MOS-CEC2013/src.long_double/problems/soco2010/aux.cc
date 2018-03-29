#include "aux.h"

#include <logger/GALogger.h>

#include <sstream>

static long double        m__;
static int           shared__;
static int           rest__;
static GA1DArrayAlleleGenome<long double>* ind_part1__;
static GA1DArrayAlleleGenome<long double>* ind_part2__;
static GA1DArrayAlleleGenome<long double>* ind_partrest__;

void setIndPart1AndPart2Values(GAGenome& ind);

void setInitValues(GAGenome& ind, long double param_m) {
  GA1DArrayAlleleGenome<long double>& ind_r = dynamic_cast<GA1DArrayAlleleGenome<long double>&>(ind);
  m__ = param_m;

  ind_part1__ = dynamic_cast<GA1DArrayAlleleGenome<long double>*>(ind_r.clone());
  ind_part2__ = dynamic_cast<GA1DArrayAlleleGenome<long double>*>(ind_r.clone());

  if (m__ <= 0.5) {
    ind_partrest__ = ind_part2__;
  } else {
    ind_partrest__ = ind_part1__;
    m__ = 1 - m__;
  }

  shared__ = (int) floor(ind_r.length() * m__);
  rest__   = 2 * shared__;

  int size1 = shared__;
  int size2 = ind_r.length() - size1;

  if (ind_partrest__ == ind_part1__) {
    int tmp = size1;
    size1 = size2;
    size2 = tmp;
  }

  ind_part1__->resize(size1);
  ind_part2__->resize(size2);
}

long double composeFuncs(GAGenome& ind, long double (*f1)(GAGenome&), long double bias1, long double (*f2)(GAGenome&), long double bias2 ) {
  setIndPart1AndPart2Values(ind);
  long double res_f1 = f1(*ind_part1__) - bias1;
  if (res_f1 < 0) {
      std::cout << "res_f1: " << res_f1 << std::endl;
      exit(0);
   }
   //assert(res_f1 >= 0.0);
  long double res_f2 = f2(*ind_part2__) - bias2; assert(res_f2 >= 0.0);
  ///*LOG*/ stringstream msg; msg << "Original ind: " << ind << endl;
  ///*LOG*/ msg << " ind1: " << *ind_part1__ << " score: " << res_f1 <<endl;
  ///*LOG*/ msg << " ind2: " << *ind_part2__ << " score: " << res_f2 << endl;
  ///*LOG*/ GALogger::instance()->appendLogMessage("Composing hybrid score functions",msg.str());
  return res_f1 + res_f2;
}

void setIndPart1AndPart2Values(GAGenome& ind) {
  GA1DArrayAlleleGenome<long double>& orig_ind = dynamic_cast<GA1DArrayAlleleGenome<long double>&>(ind);
  for (int i=0; i<shared__; i++) {
    ind_part1__->gene(i, orig_ind.gene(i*2)    );
    ind_part2__->gene(i, orig_ind.gene(i*2 + 1));
  }

  for (int  i=0; i<orig_ind.length()-rest__; i++) {
    ind_partrest__->gene(i + shared__, orig_ind.gene(i + rest__) );
  }
}

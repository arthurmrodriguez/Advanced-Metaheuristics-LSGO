#ifndef DISTAUX_H_
#define DISTAUX_H_

#include "../genomes/GAGenome.h"
#include <vector>

vector<long double> computeDistanceMatrix(vector<GAGenome*>& inds);

/*
 * Computes the D2 Method from Glover, F., C.C., Kuo and K.S. Dhir, 1998 Heuristic
 * algorithms fort he maximum diversity problem. Journal of Information and Optimization
 * Sciences
 */
void maxminD2Method(vector<GAGenome*>& inds, vector<GAGenome*>& slctd_inds);

long double computeDiversityAvg(vector<GAGenome*>& inds);

#endif /* DISTAUX_H_ */

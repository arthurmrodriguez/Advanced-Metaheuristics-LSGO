#ifndef MOSCONFIG_H
#define MOSCONFIG_H

#include <map>
#include <vector>

class MOSEA;

/* Definition of the ID type for techniques */
typedef int techIdType;

/* Definition of the associative participation vector */
typedef std::map<techIdType, long double> MOSProbVector;

/* Definition of the ID type for the encoding of solutions */
typedef int encodingType;

// mixProbVector: Computes the average value of two probability vectors (from parents)
void mixProbVector(const std::vector< const MOSProbVector* >& probs, const std::vector<long double>& scores, MOSProbVector& pvChild);

#endif

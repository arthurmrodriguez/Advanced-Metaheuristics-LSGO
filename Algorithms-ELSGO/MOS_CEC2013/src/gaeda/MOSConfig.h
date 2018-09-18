#ifndef MOSCONFIG_H
#define MOSCONFIG_H

#include <map>
#include <vector>

class MOSEA;

/* Definition of the ID type for techniques */
typedef int techIdType;

/* Definition of the associative participation vector */
typedef std::map<techIdType, double> MOSProbVector;

/* Definition of the ID type for the encoding of solutions */
typedef int encodingType;

/* Definition of the evolution type */
typedef enum {CentralEvolution, AutonomicEvolution} evolutionType;

// Definition of types
typedef void (*ParticipationFunction) (MOSEA& algorithm);

// mixProbVector: Computes the average value of two probability vectors (from parents)
void mixProbVector(const std::vector< const MOSProbVector* >& probs, const std::vector<double>& scores, MOSProbVector& pvChild);

#endif

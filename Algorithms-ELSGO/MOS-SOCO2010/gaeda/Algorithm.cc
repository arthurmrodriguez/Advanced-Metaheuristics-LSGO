#include "Algorithm.h"

Algorithm::Algorithm () {}

Algorithm::~Algorithm () {}

void Algorithm::initialize (unsigned seed) {}

void Algorithm::run (unsigned seed) {
  initialize(seed);
  evolve(seed);
}


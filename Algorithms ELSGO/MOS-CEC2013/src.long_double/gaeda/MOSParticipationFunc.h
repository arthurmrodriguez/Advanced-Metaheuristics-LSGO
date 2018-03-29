#ifndef MOSPARTICIPATIONFUNC_H
#define MOSPARTICIPATIONFUNC_H

class MOSEA;

/**
 * Constant participation function
 */
extern "C" void constantPF (MOSEA& algorithm);

/**
 * Alternating participation function which distributes a given number of
 * iterations among the available techniques.
 */
extern "C" void alternatingQualityPF (MOSEA& algorithm);

/**
 * Dynamic participation function. To calculate averages, it uses the population
 * quality.
 */
extern "C" void dynQualityMOSPF (MOSEA& alg);

/**
 * Dynamic participation function with maximum participation and decreasing
 * population size. To calculate averages, it uses the population quality.
 */
extern "C" void dynQualityMaxPartPF (MOSEA& alg);

extern "C" void dynQualityMOSPF2 (MOSEA& alg);

#endif

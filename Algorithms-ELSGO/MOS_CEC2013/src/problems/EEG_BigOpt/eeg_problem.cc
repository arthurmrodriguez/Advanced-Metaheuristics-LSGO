//
// Created by Arthur Rodriguez on 19/7/18.
// Adaptation from MOS_SOCO2010
//

#include "eeg_problem.h"
#include <GARealOps.h>
#include <MOSGenomeFactory.h>
#include <genomes/GA1DArrayGenome.h>
#include <GAEDAConfig.h>
#include <GAGeneticAlgorithm.h>


const float MIN_ALLELE_VALUE = -8;
const float MAX_ALLELE_VALUE =  8;

GA1DArrayAlleleGenome<double>* optimum_genome;

extern "C" double objective (GAGenome& g) {
    return eeg_Function(g);
}

extern "C" void individualInit (GAGenome& g) {
    return RealUniformInitializer (g);
}

extern "C" GAGenome* defineProblem (int size, char *data) {

    GAEDAConfig* cfg = GAEDAConfig::handle();

    int problemSize = ((cfg->getProblemSize()%2) == 0) ? cfg->getProblemSize() : cfg->getProblemSize()-1;
    GAAlleleSet<double> alleles (MIN_ALLELE_VALUE, MAX_ALLELE_VALUE);
    GA1DArrayAlleleGenome<double>* genome = new GA1DArrayAlleleGenome<double> (problemSize, alleles, objective);

    // Before everything else, load A, S and X according to problem size
    // so that objective function could be evaluated properly
    switch (cfg->getProblemSize()){
        case 1024:
            problemName = "D4";
            break;
        case 3072:
            problemName = "D12";
            break;
        case 4864:
            problemName = "D19";
            break;
        case 1025:
            problemName = "D4N";
            break;
        case 3073:
            problemName = "D12N";
            break;
        case 4865:
            problemName = "D19N";
            break;
        default:
            problemName = "D4";
            break;
    }

    loaddata(problemName);

    // Right after loading files, we establish the optimum
    optimum = getOptimum();

    // Common operators
    genome->initializer (RealUniformInitializer);
    genome->comparator  (RealEuclideanComparator);

    // Specific stuff for GAs
    genome->crossover   (RealBlendCrossover);
    genome->mutator     (RealGaussianMutator);

    // Specific stuff for DE
    genome->crossover   (RealExponentialCrossover);

    // Specific stuff for MOS
    MOSGenomeFactory::handle()->registerGenome (GAID::RealEncoding, genome);

    optimum_genome = new GA1DArrayAlleleGenome<double> (problemSize, alleles, objective);
    optimum_genome->comparator (RealEuclideanComparator);
    optimum_genome->crossover  (RealBlendCrossover);

    // Check ICAComponent (matrix S) size
    // cout<<"Tamaño ICAComponent " <<ICAcomponent.size();

    for (int i = 0; i < optimum_genome->length (); i++)
        optimum_genome->gene (i, optimum[i]);

    return genome;

}

// Get the scale values for SO scoreboard
extern "C" bool postprocess (GAPopulation* pop, int rank) {
    cout<< "-> SCALE: " <<f1min <<" " <<f1max <<" " <<f2min <<" " <<f2max <<endl;
    delete optimum_genome;
    return true;
}

extern "C" const char *describeProblem () {
    return ("Electroencephalography Big Opt - " + problemName).c_str();
}

extern "C" GAGenome* getOptimumGenome () {
    return optimum_genome;
}

extern "C" GAGenome::OptCriterion optCriterion(){
    return GAGenome::MINIMIZATION;
}

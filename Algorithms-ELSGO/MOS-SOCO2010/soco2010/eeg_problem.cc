//
// Created by Arthur Rodriguez on 5/7/18.
//

#include "eeg_problem.h"
#include <GARealOps.h>
#include <MOSGenomeFactory.h>
#include <genomes/GA1DArrayGenome.h>
#include <GAEDAConfig.h>
#include <GAGeneticAlgorithm.h>


const float MIN_ALLELE_VALUE = -8;
const float MAX_ALLELE_VALUE =  8;

GA1DArrayAlleleGenome<long double>* optimum_genome;

extern "C" long double objective (GAGenome& g) {
    return eeg_Function(g);
}

extern "C" void individualInit (GAGenome& g) {
    return RealUniformInitializer (g);
}

extern "C" GAGenome* defineProblem (int size, char *data) {

    GAEDAConfig* cfg = GAEDAConfig::handle();

    GAAlleleSet<long double> alleles (MIN_ALLELE_VALUE, MAX_ALLELE_VALUE);
    GA1DArrayAlleleGenome<long double>* genome = new GA1DArrayAlleleGenome<long double> (cfg->getProblemSize (), alleles, objective);

    // Before everything else, load A, S and X according to problem size
    // so that objective function could be evaluated properly
    switch (cfg->getProblemSize()){
        case 1024:
            problemName = "D4";
            loaddata("D4");
            break;
        case 3072:
            problemName = "D12";
            loaddata("D12");
            break;
        case 4864:
            problemName = "D19";
            loaddata("D19");
            break;
        default:
            problemName = "D4";
            loaddata("D4");
            break;
    }

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

    optimum_genome = new GA1DArrayAlleleGenome<long double> (cfg->getProblemSize (), alleles, objective);
    optimum_genome->comparator (RealEuclideanComparator);
    optimum_genome->crossover  (RealBlendCrossover);

    // Check ICAComponent (matrix S) size
    // cout<<"Tamaño ICAComponent " <<ICAcomponent.size();

    for (int i = 0; i < optimum_genome->length (); i++)
        optimum_genome->gene (i, optimum[i]);

    return genome;

}

// Get the scale values for SO optimization scoreboard
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


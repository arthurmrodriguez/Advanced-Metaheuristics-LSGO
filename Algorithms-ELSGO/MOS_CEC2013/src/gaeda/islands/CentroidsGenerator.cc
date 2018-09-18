#include "CentroidsGenerator.h"

CentroidsGenerator::CentroidsGenerator(unsigned int pcent_number, GAGenome& psample_gen,GAGenome::Initializer pinit_func) 
                                      : cent_number(pcent_number), sample_gen(psample_gen), init_func(pinit_func) {}

CentroidsGenerator::~CentroidsGenerator(){}

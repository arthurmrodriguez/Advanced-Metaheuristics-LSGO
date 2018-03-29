#include <iostream>

#include <GARealOps.h>

#include "ImageManager.h"
#include "VesicleGenome.h"

ImageManager* img = NULL;
std::vector<VesicleGenome<long double>*> synapses;

const int minWidth = 2, maxWidth = 2;   // Vesicle width range of iteration
const int minRadius = 4, maxRadius = 8; // Vesicle radius range of iteration

const long double alpha = 0.9;

long double penalty (VesicleGenome<long double>& genome) {

  long double val = RealEuclideanComparator(genome, *synapses[0]);

  for (unsigned i = 1; i < synapses.size(); i++) {
    long double aux = RealEuclideanComparator(genome, *synapses[i]);
    if (aux < val)
      val = aux;
  }
  std::cout << "penalty distance: " << val << std::endl;

  long double val2 = 0.0;

/*
  if (found_vesicles.size() > 0) {
    val2 = RealEuclideanComparator(genome, *found_vesicles[0]);

    for (unsigned i = 1; i < found_vesicles.size(); i++) {
      long double aux = RealEuclideanComparator(genome, *found_vesicles[i]);
      if (aux < val2)
        val2 = aux;
    }
  }
*/
  return val + val2;

}

long double fitness (VesicleGenome<long double>& genome) {

  long double res, aux;

  VesicleInfo info;

  for (int radius=minRadius; radius<=maxRadius; radius++) {
    for (int width=minWidth; width<=maxWidth; width++) {
      std::cout << "radius: " << radius << std::endl;
      std::cout << "width: " << width << std::endl;
      aux = img->objectiveFunction((int)genome.gene(0), (int)genome.gene(1), radius, width, &info);
      std::cout << "objective: " << aux << std::endl << std::endl;
      if ((radius==minRadius && width==minWidth) || (aux<res)) {
        res = aux;
        genome.radius  (radius);
        genome.width   (width );
      }
    }
  }

  return res;

}

extern "C" long double objective (GAGenome& g) {

  VesicleGenome<long double>& genome = DYN_CAST (VesicleGenome<long double>&, g);

  std::cout << "==================================" << std::endl;
  long double fit = fitness(genome);
  long double pen = penalty(genome);
  long double res = (alpha * fit) + ((1 - alpha) * pen);
  std::cout << "objective: " << fit << std::endl;
  std::cout << "res: " << res << std::endl;
  std::cout << "==================================" << std::endl;

  return res;

}

int main(int argc, char** argv) {

  img = new ImageManager(argv[1]);

  // Definition of the genome and the allele set
  GAAlleleSetArray<long double> alleles;
  alleles.add(0.0, img->getDimX());
  alleles.add(0.0, img->getDimY());
  alleles.add(0.0, 0.0           );

  // Discretized synapse (to compute penalties)
  VesicleGenome<long double>* genome2 = new VesicleGenome<long double> (alleles, objective);
  genome2->gene(0, 540);
  genome2->gene(1, 170);
  genome2->gene(2,   0);
  synapses.push_back(genome2);

  genome2 = new VesicleGenome<long double> (alleles, objective);
  genome2->gene(0, 520);
  genome2->gene(1, 200);
  genome2->gene(2,   0);
  synapses.push_back(genome2);

  genome2 = new VesicleGenome<long double> (alleles, objective);
  genome2->gene(0, 470);
  genome2->gene(1, 230);
  genome2->gene(2,   0);
  synapses.push_back(genome2);

  genome2 = new VesicleGenome<long double> (alleles, objective);
  genome2->gene(0, 360);
  genome2->gene(1, 200);
  genome2->gene(2,   0);
  synapses.push_back(genome2);


  genome2 = new VesicleGenome<long double> (alleles, objective);
  genome2->gene(0, 363);
  genome2->gene(1, 189);
  genome2->gene(2,   0);

  objective(*genome2);

  genome2 = new VesicleGenome<long double> (alleles, objective);
  genome2->gene(0, 519);
  genome2->gene(1, 201);
  genome2->gene(2,   0);

  objective(*genome2);

  genome2 = new VesicleGenome<long double> (alleles, objective);
  genome2->gene(0, 511);
  genome2->gene(1, 268);
  genome2->gene(2,   0);

  objective(*genome2);

  genome2 = new VesicleGenome<long double> (alleles, objective);
  genome2->gene(0, 382);
  genome2->gene(1, 130);
  genome2->gene(2,   0);

  objective(*genome2);


  return 0;
}

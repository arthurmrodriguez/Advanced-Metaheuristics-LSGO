#include <iostream>

#include "Vesicle.h"
#include "ImageManager.h"
#include "BruteForceAlgorithm.h"

int main(int argc, char** argv) {
  if (argc!=4) {
    std::cout << "Usage: ";
    std::cout << argv[0] << " input_image output_image number_of_vesicles" << std::endl;
    return -1;
  }
  else {
    BruteForceAlgorithm bfa(argv[1]);
    int n = atoi(argv[3]);
    for (int i=1; i<=n; i++) {
      bfa.run();
      bfa.getBestVesicle();
      std::cout << "Fin de iteracion " << i << std::endl;
    }
    bfa.printSelectedVesicles(argv[2]);
    return 0;
  }
}

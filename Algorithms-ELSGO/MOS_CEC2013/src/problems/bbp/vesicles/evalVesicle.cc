#include <iostream>

#include "ImageManager.h"
#include "VesicleInfo.h"

ImageManager* img = NULL;

int main (int argc, char** argv) {

  VesicleInfo info;
  img = new ImageManager(argv[1]);

  double fit = img->objectiveFunction(329, 167, 4, 2, &info);

  std::cout << "fit: " << fit << std::endl;
  return 0;

}

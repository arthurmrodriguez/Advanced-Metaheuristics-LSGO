#include "BruteForceAlgorithm.h"

#include "ImageManager.h"

// Class constructor
BruteForceAlgorithm::BruteForceAlgorithm(char* _path) {
  img = new ImageManager(_path);
  maxX = img->getDimX();
  maxY = img->getDimY();
}

BruteForceAlgorithm::~BruteForceAlgorithm() {
  if (img)
    delete img;
}

// Runs the vesicle brute force algorithm for vesicle detection
void BruteForceAlgorithm::run() {
  int xRange = maxX-maxRadius;
  int yRange = maxY-maxRadius;
  long double value, valueAux;
  Vesicle* vesicle = NULL;

  for (int x=maxRadius; x<xRange; x++) { // Iterates over the x axis
    for (int y=maxRadius; y<yRange; y++) { // Iterates over the y axis
      for (int radius=minRadius; radius<=maxRadius; radius++) { // Iterates over the given radius range
        for (int width=minWidth; width<=maxWidth; width++) { // Iterates over the given width range
          valueAux = img->objectiveFunction(x, y, radius, width);
          if ((radius==minRadius) || (valueAux<value)) {
            value = valueAux;
            if (vesicle)
              delete vesicle;
            vesicle = new Vesicle(x, y, radius, width, value);
          }
        } // width
      } // radius
      vesicles.push_back(*vesicle);
      delete vesicle;
      vesicle = NULL;
    } // y
  } // x
}

void BruteForceAlgorithm::getBestVesicle() {
  int len = vesicles.size();

  Vesicle vesicle = vesicles[0];
  long double valueAux=0, value=vesicle.getValue();

  for (int i=1; i<len; i++) {
    value = vesicles[i].getValue();
    if (value<valueAux) {
      valueAux = value;
      vesicle = vesicles[i];
    }
  }

  img->addVesicle(vesicle);
  vesicles.clear();
}

void BruteForceAlgorithm::printSelectedVesicles(char* _path) {
  img->printSelectedVesicles(_path);
}

#ifndef BRUTEFORCE
#define BRUTEFORCE

#include <vector>

class ImageManager;
class Vesicle;

class BruteForceAlgorithm {
  public:
    // Class constructor
    BruteForceAlgorithm(char* _path);
    ~BruteForceAlgorithm();
    void run();
    void getBestVesicle();
    void printSelectedVesicles(char* _path);

  private:
    ImageManager* img; // Data structure for image access
    int maxX, maxY;  // Limit coordinates of the image
    static const int minWidth = 2, maxWidth = 2;   // Vesicle width range of iteration
    static const int minRadius = 4, maxRadius = 6; // Vesicle radius range of iteration
    std::vector<Vesicle> vesicles; // Vesicles vector
};

#endif

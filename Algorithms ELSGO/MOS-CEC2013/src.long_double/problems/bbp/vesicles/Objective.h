#ifndef OBJECTIVE_H
#define OBJECTIVE_H

#include "ITKConfig.h"
#include "VesicleInfo.h"

class Objective {
  public:
    // Class constructor
    Objective (int _x, int _y, int _radius, int _width, NeighborhoodIteratorType* _nit, std::vector<bool>& _occupation, int _dimX, int _dimY);
    ~Objective();

    long double objectiveFunction(VesicleInfo* info=NULL);
    void markVesiclesMap();

  protected:
    void handlePixel(int pixelPosition, bool isEdge);
    void iterateCircle();
    void handlePixelDeviation(int pixelPosition, bool isEdge);
    void iterateCircleDeviation();
    void markPixel(int pixelPosition, bool isEdge);
    void iterateCircleMark();

    // Circunference variables
    int x, y;       // Circunference center
    int radius;     // Circunference radius
    int diameter;   // Circunference diameter
    int width;      // Circunference edge width

    // Data structure for image approaching
    NeighborhoodIteratorType* nit;

    // Pixel occupation matrix
    std::vector<bool>& occupation;
    int tamX;
    int tamY;

    // Color intensity values
    static int const black = 0;
    static int const white = 255;

    // Auxiliar variables for the objective function value computation
    int edgePx ;                // Number of pixels of the circunference edge
    int centerPx ;              // Number of pixels of the circunference surface
    int edgeIntensity;          // Edge intensities sumatory
    int centerIntensity;        // Center intensities sumatory
    long double edgeAvgIntensity;    // Average edge intensity
    long double centerAvgIntensity;  // Average center intensity
    long double edgeDeviation;       // Edge intensity deviation
    long double centerDeviation;     // Center intensity deviation

    // Weight constants
    long double kEdge, kCenter;
};

#endif

#ifndef DRAW
#define DRAW

#include "Magick++.h"

#include "ITKConfig.h"

class Vesicle;

class DrawVesicle {
  static Magick::Color currentColor;

  public:
    // Class constructor
    DrawVesicle (int _x, int _y, int _radius, int _width, NeighborhoodIteratorType* _nit, Magick::Image& _outImage, int _valid = 0);
    DrawVesicle (Vesicle _vesicle, NeighborhoodIteratorType* _nit, Magick::Image& _outImage);
    ~DrawVesicle();

    void draw() const;
    void draw2() const;

  private:
    void drawPixel(int pixelPosition, bool isEdge) const;

    // Circunference variables
    int x, y;       // Circunference center
    int radius;     // Circunference radius
    int diameter;   // Circunference diameter
    int width;      // Circunference edge width

    // Vesicle value
    long double value;

    // Is a valid vesicle?
    int valid;

    // Data structure for image approaching
    NeighborhoodIteratorType* nit;
    Magick::Image& outImage;

    // Circunference iterator variables
    int row;        // Row index
    int rowCount;   // Total number of rows (2*radius + central row)
    int initRow;    // Starting position of the row iterator range
    int endRow;     // Ending position of the row iterator range
    int gapOffset;  // Margin for excluding the pixels outside the circunference

};

#endif

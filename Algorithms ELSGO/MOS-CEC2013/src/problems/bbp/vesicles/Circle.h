#ifndef CIRCLE
#define CIRCLE

#include <vector>

class Circle {
  public:
    // Class constructor
    Circle(int _radius, int _width);
    ~Circle();

    int isMarked(int x, int y) const {return map[x][y];}

  protected:
    void circlePoints(int cx, int cy, int x, int y);
    void circleMidpoint(int xCenter, int yCenter, int radius);

    // Attributes
    int radius;  // Radius of the vesicle approaching circunference (in pixels)
    int width;   // Width of the vesicle approaching circunference outline (in pixels)
    std::vector< std::vector<int> > map;
};

#endif

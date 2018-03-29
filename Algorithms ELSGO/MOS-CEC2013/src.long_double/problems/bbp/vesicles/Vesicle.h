#ifndef VESICLE_H
#define VESICLE_H

class Vesicle {
  public:
    // Class constructor
    Vesicle(int _x, int _y, int _radius, int _width, long double _value, int _valid = 0);
    ~Vesicle();

    int    getX     () const {return x;     }
    int    getY     () const {return y;     }
    int    getRadius() const {return radius;}
    int    getWidth () const {return width; }
    long double getValue () const {return value; }
    int    isValid  () const {return valid; }

  protected:
    // Attributes
    int x, y;     // X and Y axis coordinates
    int radius;   // Radius of the vesicle approaching circunference (in pixels)
    int width;    // Width of the vesicle approaching circunference outline (in pixels)
    long double value; // Value of the vesicle fitness function
    int valid;    // Is this a valid (according to its parameters) vesicle?
};

#endif

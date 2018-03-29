#include "DoubleDip.h"

#include <math.h>

DoubleDip::DoubleDip()
  : UnitFunction1D()
{
}

DoubleDip::DoubleDip(long double centre, long double scale)
  : UnitFunction1D(centre, scale)
{
}

DoubleDip::~DoubleDip()
{
}

long double DoubleDip::getValue(long double point) {

  long double depth = 0.0;
  long double xs;
  long double x = (point - centre_) / scale_;

  if (x > -0.5 and x < 0.5) {
    xs = 4*x*x;               // scale to -1 to 1 -> (2x)^2
    depth = (-96*xs*xs*xs + 193*xs*xs - 98*xs +1) * scale_;
  }

  return depth;

}

long double DoubleDip::twist(long double x, long double y) {

  long double dx = 0;
  y = fmod (y, 1);
  //y = (y % 1);
  long double ys = y*y;

  if (y > 0)
    dx = 4 * (ys*ys - 2*ys*y + ys); // twisting quartic
  else
    dx = 4 * (ys*ys + 2*ys*y + ys); // twisting quartic

  return dx;

}

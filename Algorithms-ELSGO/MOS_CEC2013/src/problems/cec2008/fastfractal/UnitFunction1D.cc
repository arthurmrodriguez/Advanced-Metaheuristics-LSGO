#include "UnitFunction1D.h"

UnitFunction1D::UnitFunction1D()
  : centre_ (0),
    scale_  (1)
{
}

UnitFunction1D::UnitFunction1D(long double centre, long double scale)
  : centre_ (centre),
    scale_  (scale)
{
}

UnitFunction1D::~UnitFunction1D()
{
}

void UnitFunction1D::setParams(long double centre, long double scale) {
  centre_ = centre;
  scale_ = scale;
}
  
void UnitFunction1D::setCentre(long double centre) {
  centre_ = centre;
}
  
long double UnitFunction1D::getCentre() {
  return centre_;
}
  
void UnitFunction1D::setScale(long double scale) {
  scale_ = scale;
}

long double UnitFunction1D::getScale() {
  return scale_;
}

long double UnitFunction1D::getValue(long double point) {
  return 0.0;
}

long double UnitFunction1D::twist(long double x, long double y) {
  return 0.0;
}
 
std::string UnitFunction1D::getName() {
  return "UnitFunction1D";
}

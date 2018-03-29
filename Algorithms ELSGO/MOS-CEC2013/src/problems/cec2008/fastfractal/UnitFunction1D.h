#ifndef UNITFUNCTION1D_H_
#define UNITFUNCTION1D_H_

#include<string>

class UnitFunction1D {

 protected:

  long double centre_;
  long double scale_;

 public:

  UnitFunction1D();
  UnitFunction1D(long double centre, long double scale);

  virtual ~UnitFunction1D();
	
  void setParams(long double centre, long double scale);
  
  void        setCentre(long double centre);
  long double getCentre();
  
  void        setScale(long double scale);
  long double getScale();
  
  virtual long double getValue(long double point);
 
  virtual long double twist(long double x, long double y);
 
  std::string getName();
 
};

#endif /*UNITFUNCTION1D_H_*/

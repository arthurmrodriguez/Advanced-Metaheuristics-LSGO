#ifndef SERIALIZABLE_H_
#define SERIALIZABLE_H_

#include <iostream>

using namespace std;

/**
 * Used for serializing the genomes. The methods have not been declared pure virtual so 
 * not all the genomes are obliged to implement them (altough they should be), but for 
 * a quick compiling, this hack has been done. In a future all the genomes should be
 * changed in order to implement these methods
 */
class Serializable {
protected:
  virtual ~Serializable() {};
  virtual void readObject(istream& is)        = 0;
  virtual void writeObject(ostream& os) const = 0;
};

#endif /*SERIALIZABLE_H_*/

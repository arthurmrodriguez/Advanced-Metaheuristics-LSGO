#ifndef RANQD1_H_
#define RANQD1_H_

class RanQD1 {

 private:

  static long long   MASK_;    // lower order 32 bits of long
  static long double MAX_INT_; // 2^32-1 (MASK as a double)
  static long long   A_;       // suggested by Knuth
  static long long   C_;       // suggested by Lewis

  long long idum_;

 public:

  RanQD1();
  RanQD1(long seed);

  virtual ~RanQD1();

  void setSeed(long seed);

  long long   nextLong  ();
  long double nextDouble();
  int         nextInt   (int min, int max);
	
};

#endif /*RANQD1_H_*/

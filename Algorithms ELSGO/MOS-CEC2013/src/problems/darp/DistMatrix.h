#ifndef DISTMATRIX_H
#define DISTMATRIX_H

#include <string>
#include <vector>
#include <assert.h>
#include <string>

class DistMatrix {
private:
  std::vector< std::vector<double> > data_;

protected:
  virtual void clear()                { data_.clear();  }
  virtual void resize(int s)          { data_.resize(s);}
  virtual void resize(int pos, int s) { assert(pos < data_.size()); data_[pos].resize(s); }
  virtual int  numRows()              { return data_.size(); }
  virtual int  numColumns(int row)    { return data_[row].size(); }

public:

  DistMatrix() {}// Just for the tests
  virtual ~DistMatrix() {}

  DistMatrix(int size) : data_(size+1, std::vector<double>(size+1, 0.0)) {}

  DistMatrix(std::string distFile, int size);

  virtual void setDist(int ori, int dest, double value) {
    assert(data_.size() > 0);
    assert(ori >= 0);
    assert(dest >= 0);
    assert(ori < data_.size());
    assert(dest < data_[ori].size());
    data_[ori][dest] =value;
  }

  // Note: cannote be constant since the LazyEvalDistMatrix needs to redefine this method
  virtual double getDist(int ori, int dest) const {
    assert(data_.size() > 0);
    assert(ori >= 0);
    assert(dest >= 0);
    assert(ori < data_.size());
    assert(dest < data_[ori].size());
    return data_[ori][dest];
  }


};

#endif

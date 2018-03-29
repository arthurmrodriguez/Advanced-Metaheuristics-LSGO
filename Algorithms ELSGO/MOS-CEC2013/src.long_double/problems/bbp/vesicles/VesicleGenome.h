#ifndef VESICLE_GENOME
#define VESICLE_GENOME

#include <vector>

#include <genomes/GA1DArrayGenome.h>

#include "VesicleInfo.h"

template <class T> class VesicleGenome : public GA1DArrayAlleleGenome<T> {

   public:

      VesicleGenome (unsigned int x, const GAAlleleSet<T>& a, GAGenome::Evaluator f = (GAGenome::Evaluator) 0, void* u = (void*) 0)
        : GA1DArrayAlleleGenome<T> (x, a, f, u), radius_(0), width_(0), valid_(0) {}

      VesicleGenome (const GAAlleleSetArray<T>& a, GAGenome::Evaluator f = (GAGenome::Evaluator) 0, void* u = (void*) 0)
        : GA1DArrayAlleleGenome<T> (a, f, u), radius_(0), width_(0), valid_(0) {}

      VesicleGenome (const VesicleGenome<T>& v)
        : GA1DArrayAlleleGenome<T> (v), radius_(v.radius()), width_(v.width()), valid_(v.valid()) {fits_ = v.fits(); infos_ = v.infos();}

      VesicleGenome<T>& operator= (const GAGenome& arr) {copy(arr); return *this;}
      VesicleGenome<T>& operator= (const T array []) { // no err checks!

         GA1DArrayAlleleGenome<T>::operator= (array);
         return *this;

      }

      virtual ~VesicleGenome () {}

      virtual GAGenome* clone (GAGenome::CloneMethod flag = GAGenome::CONTENTS) const {return new VesicleGenome<T> (*this);}

      virtual void copy (const GAGenome& orig) {
        if (&orig == this)
          return;

        const VesicleGenome<T>* c = DYN_CAST (const VesicleGenome<T>*, &orig);

        if(c) {
          GA1DArrayAlleleGenome<T>::copy (*c);

          radius_ = c->radius_;
          width_  = c->width_;
          valid_  = c->valid_;
          fits_   = c->fits_;
          infos_  = c->infos_;
        }
      }

      int radius(           ) const {return radius_;}
      int radius(int _radius)       {return radius_ = _radius;}

      int width(          ) const {return width_;}
      int width(int _width)       {return width_ = _width;}

      int valid(          ) const {return valid_;}
      int valid(int _valid) {return valid_ = _valid;}

      std::vector<long double> fits() const {return fits_;}
      void addFit(long double fit) {fits_.push_back(fit);}
      void resetFits() {fits_.clear();}

      std::vector<VesicleInfo> infos() const {return infos_;}
      void addInfo(VesicleInfo info) {infos_.push_back(info);}
      void resetInfos() {infos_.clear();}

  protected:
    int radius_;
    int width_;
    int valid_;
    std::vector<long double> fits_;
    std::vector<VesicleInfo> infos_;
};

#endif

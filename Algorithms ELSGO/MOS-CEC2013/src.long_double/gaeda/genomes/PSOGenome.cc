#include "PSOGenome.h"

PSOGenome::PSOGenome(const PSOGenome& orig) : GARealGenome(orig){
  best_ = dynamic_cast<GARealGenome*> (orig.best()->clone());
  vel_  = dynamic_cast<GARealGenome*> (orig.velocity()->clone());
  copy(orig);
}

PSOGenome::PSOGenome(const GARealGenome& c) : GARealGenome(c){
  GARealGenome::copy(c);
  best_ = dynamic_cast <GARealGenome*> (c.clone());
  vel_  = dynamic_cast <GARealGenome*> (c.clone());
}

PSOGenome::~PSOGenome() {
 delete best_;
 delete vel_;
}

PSOGenome* PSOGenome::clone (GAGenome::CloneMethod flag) const {

   PSOGenome* cpy = new PSOGenome (*this);
   return cpy;
}

void PSOGenome::copy(const PSOGenome& orig){
  assert (best_ != NULL); assert(vel_ != NULL);
  
  GARealGenome::copy (orig);
   
  best_->copy(*(orig.best()));
  vel_->copy(*(orig.velocity()));  
}

GARealGenome* PSOGenome::best() const{
  return best_;
}

GARealGenome* PSOGenome::velocity() const{
  return vel_;
}

void PSOGenome::updateBest(){
  if(this->score() > best_->score()){
   best_->copy(* this);
  }
}

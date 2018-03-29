#include "ImprovPerDimManager.h"
#include "garandom.h"
#include "logger/GALogger.h"
#include <algorithm>
#include <stdexcept>
#include <sstream>

ImprovPerDimManager* ImprovPerDimManager::instance__ = NULL;

ImprovPerDimManager::ImprovPerDimManager(int num_dims) : num_dims_(num_dims) {
  instance__ = this;
}

ImprovPerDimManager::~ImprovPerDimManager() {}

#include "quicksort.h"

static void swap (vector<int>& dims, vector<long double>& values, int p1, int p2) {

   long double tmp_val = values [p1];
   int    tmp_dim = dims   [p1];

   values [p1] = values [p2];
   dims   [p1] = dims   [p2];

   values [p2] = tmp_val;
   dims   [p2] = tmp_dim;

}

//This function orders the vector from the 0 to the IND_SIZE-1 elements
//Way to call it: quicksort(keys, values, 0, IND_SIZE-1);
static void quicksort (vector<int>& keys, vector<long double>& values, int l, int r) {
  assert( keys.size() == values.size() );
   int i, j;
   long double v;

   if (r > l) {

      v = values [r];
      i = l - 1;
      j = r;

      for (;;) {

         while (values [++i] < v);
         while (values [--j] > v);

         if (i >= j) break;

         swap (keys, values, i, j);

      } // for

      swap (keys, values, i, r);

      quicksort (keys, values, l, i - 1);
      quicksort (keys, values, i + 1, r);

   } // if

}

void ImprovPerDimManager::getDimsWhichProbSumAThreshold(vector<int>& selected_dims, long double prob_percent, int num_other_dims) const {
  assert(num_other_dims >= 0 && num_other_dims < numDims() );

  selected_dims.clear();
  if (num_other_dims == 0) num_other_dims = 1;

  vector<long double>* probs = getProbs();
  vector<int>     dims(numDims()); for (int i=0; i<dims.size(); i++) dims[i] = i;
  assert(probs->size() == numDims());
  assert(probs->size() == dims.size());

  // /*LOG*/ stringstream msg; msg << "Before sorting:" << endl;
  // /*LOG*/ msg << "  Probs: "; for (int i=0; i<probs->size(); i++) msg << (*probs)[i] << ", "; msg << endl;
  // /*LOG*/ msg << "  Dims "  ; for (int i=0; i<dims.size(); i++)   msg << dims[i]     << ", "; msg << endl;
  quicksort(dims, *probs, 0, probs->size()-1);
  // /*LOG*/ msg << "After sorting:" << endl;
  // /*LOG*/ msg << "  Probs: "; for (int i=0; i<probs->size(); i++) msg << (*probs)[i] << ", "; msg << endl;
  // /*LOG*/ msg << "  Dims "  ; for (int i=0; i<dims.size(); i++)   msg << dims[i]     << ", "; msg << endl;

  long double sum_probs = 0.0;
  int pos;
  for (pos=probs->size()-1; pos>=0; pos--) {
    selected_dims.push_back(dims[pos]);
    sum_probs += (*probs)[pos];
    // /*LOG*/ std::cout << "sum para pos=" << pos << " = " << sum_probs << endl;
    if (sum_probs > prob_percent) break;
  }
  assert(pos >= 0 && pos < numDims() );

  // /*LOG*/ msg << "Pos alcanzada = " << pos << " sum probs=" << sum_probs << endl;
 
  // Now for the rest of the dim
  if (num_other_dims > 0 && pos > 0) {
    if (pos<num_other_dims) num_other_dims = pos;
    // /*LOG*/ msg << "num other dims fixed: " << num_other_dims << endl;
    // /*LOG*/ msg << "pos= " << pos << " dims before shuffling: "; for (int i=0; i<dims.size(); i++) msg << dims[i] << ", "; msg << endl;
    random_shuffle(dims.begin(), dims.begin()+pos);
    // /*LOG*/ msg << "dims after shuffling: "; for (int i=0; i<dims.size(); i++) msg << dims[i] << ", "; msg << endl;
    for (int i=0; i<num_other_dims; i++) selected_dims.push_back(dims[i]);
  }

  // /*LOG*/ msg << "Final selected dims: "; for (int i=0; i<selected_dims.size(); i++) msg << selected_dims[i] << ", "; msg << endl;

    
  // /*LOG*/ GALogger::instance()->appendLogMessage("ImprovePerDimManager:getDimsWhichProbSumAThreshold", msg.str());

  delete probs;
}

vector<long double>* ImprovPerDimManager::getProbs() const {
  vector<long double> dim_values = dimPerfValues();
  long double sum = 0.0; for (int i=0; i<dim_values.size(); i++) sum += dim_values[i];
  assert( sum >= 0.0 );

  vector<long double>* probs = new vector<long double>(dim_values.size(),0.0);

  if (sum > 0.0) { // If all performance values are = 0.0 nothing is done

    for (int i=0; i<probs->size(); i++) {
      assert(dim_values[i]>=0.0);

      (*probs)[i] = dim_values[i]/sum;

      assert((*probs)[i] >=-0.01 && (*probs)[i] <=1.01);
    }
  }
  else {

    for (int i=0; i<probs->size(); i++) {
      (*probs)[i] = (i+1) / probs->size();
      assert((*probs)[i] >=-0.01 && (*probs)[i] <=1.01);
    }
  }

  // /*LOG*/ stringstream msg; msg << "Probs computed: " << endl;
  // /*LOG*/ for (int i=0; i<probs->size(); i++) msg << (*probs)[i] << ", "; msg << endl;
  // /*LOG*/ GALogger::instance()->appendLogMessage("ImprovePerDimManager", msg.str());

  return probs;
}

vector<long double>* ImprovPerDimManager::getAccumProbs() const {
  vector<long double> dim_values = dimPerfValues();

  // /*LOG*/ stringstream msg; msg << "Dim Per values: ";
  // /*LOG*/ for (int i=0; i<dim_values.size(); i++) msg << dim_values[i] << ", "; msg << endl;

  long double sum = 0.0; for (int i=0; i<dim_values.size(); i++) sum += dim_values[i];
  assert( sum >= 0.0 );

  vector<long double>* probs = new vector<long double>(dim_values.size(),0.0);

  if (sum > 0.0) { // If all performance values are = 0.0 nothing is done
    //int    num_zeros = 0;
    //long double min_non_zero_prob_value = 1.1; // initialized with a bigger value than the ones reached

    for (int i=0; i<probs->size(); i++) {
      assert(dim_values[i]>=0.0);

      (*probs)[i] = dim_values[i]/sum;

      //if      ( (*probs)[i] == 0.0 )                    num_zeros++;
      //else if ( (*probs)[i] < min_non_zero_prob_value ) min_non_zero_prob_value = (*probs)[i];

      if (i > 0) (*probs)[i] += (*probs)[i-1];

      assert((*probs)[i] >=-0.01 && (*probs)[i] <=1.01);
    }
  }

  // COmmented since it was not a good idea, there are cases where the probability is so low that this correction does not change the scenario
  //// We always assign some probability to the less probable dimensions so they have a chance to recover
  //// The distribution is made so that the probability of the less probable dimension is divided among himself and the
  //// dimensions with zero probability
  //if (num_zeros > 0) {
  //  for (int i=0; i<probs->size(); i++) {
  //    if ( (*probs)[i] == 0.0 || (*probs)[i] == min_non_zero_prob_value ) {
  //      (*probs)[i] = min_non_zero_prob_value / (num_zeros + 1);
  //    }
  //  }
  //}

  // /*LOG*/ msg << "Probs computed: " << endl;
  // /*LOG*/ for (int i=0; i<probs->size(); i++) msg << (*probs)[i] << ", "; msg << endl;
  // /*LOG*/ GALogger::instance()->appendLogMessage("ImprovePerDimManager", msg.str());

  return probs;
}

vector<long double> ImprovPerDimManager::dimPerfValues() const {
  // Needed to be implemented in order to be called from the posBasedOnProbs method
  throw runtime_error("Method dimPerfValues from ImprovPerDimManager should not be called");
  vector<long double> foo;
  return foo;
}

const long double MIN_PROB = 0.01;

vector<int>* getDimsWithVSmallProb(vector<long double>& probs) {
  vector<int>* dims_with_vsmall_prob = new vector<int>();

  if (probs[0] < MIN_PROB) dims_with_vsmall_prob->push_back(0);
  for (int i=1; i<probs.size(); i++) {
    if ( probs[i] - probs[i-1] < MIN_PROB) dims_with_vsmall_prob->push_back(i);
  }
  random_shuffle(dims_with_vsmall_prob->begin(), dims_with_vsmall_prob->end());
  return dims_with_vsmall_prob;
}

void ImprovPerDimManager::getDimsWithRep (vector<int>& dims) const {
  assert(dims.size() <= numDims());
  // /*LOG*/ stringstream msg; msg << "Requested number of dims: " << dims.size() << " Selected Dims: ";

  vector<long double>* probs = getAccumProbs();

  // Hack needed to ensure that all the dimensions with small probabilities get selected
  // Done here and not in 
  // Every dimension with a probability less thatn MIN_PROB has the same chance (MIN_PROB)
  // of being selected
  vector<int>* dims_with_vsmall_prob = getDimsWithVSmallProb(*probs);
  // /*LOG*/ msg << "dims with small value set: "; for (int i=0;i<dims_with_vsmall_prob->size();i++) msg << (*dims_with_vsmall_prob)[i] << ", "; msg << endl;

  int dim;
  for (int i=0; i<dims.size(); i++) {
    if (GARandomDouble(0.0,1.0) < MIN_PROB && dims_with_vsmall_prob->size() > 0) {
      dim = dims_with_vsmall_prob->back();
      dims_with_vsmall_prob->pop_back();
        // /*LOG*/ msg << "dim=" << dim << " selected from the set of dims with small value set" << endl;
    }
    else {
      dim = posBasedOnProbs(*probs);
    }

    dims[i] = dim;
  }

  // /*LOG*/ msg << "Selected dims: "; for (int i=0; i<dims.size(); i++) msg << dims[i] << ", "; msg << endl;
  // /*LOG*/ GALogger::instance()->appendLogMessage("ImprovPerDimManager", msg.str());

  delete dims_with_vsmall_prob;
  delete probs;
}

void ImprovPerDimManager::getDims (vector<int>& dims) const {
  assert(dims.size() <= numDims());
  // /*LOG*/ stringstream msg; msg << "Requested number of dims: " << dims.size() << " Selected Dims: ";

  vector<bool>    used_dims (numDims(),false);
  vector<long double>* probs = getAccumProbs();

  // Hack needed to ensure that all the dimensions with small probabilities get selected
  // Done here and not in 
  // Every dimension with a probability less thatn MIN_PROB has the same chance (MIN_PROB)
  // of being selected
  vector<int>* dims_with_vsmall_prob = getDimsWithVSmallProb(*probs);
  // /*LOG*/ msg << "dims with small value set: "; for (int i=0;i<dims_with_vsmall_prob->size();i++) msg << (*dims_with_vsmall_prob)[i] << ", "; msg << endl;

  int dim;
  for (int i=0; i<dims.size(); i++) {
    do {
      if (GARandomDouble(0.0,1.0) < MIN_PROB && dims_with_vsmall_prob->size() > 0) {
        dim = dims_with_vsmall_prob->back();
        dims_with_vsmall_prob->pop_back();
        // /*LOG*/ msg << "dim=" << dim << " selected from the set of dims with small value set" << endl;
      }
      else {
        dim = posBasedOnProbs(*probs);
      }
    } while (used_dims[dim]);

    used_dims[dim]      = true;
    dims[i] = dim;
  }

  // /*LOG*/ msg << "Selected dims: "; for (int i=0; i<dims.size(); i++) msg << dims[i] << ", "; msg << endl;
  // /*LOG*/ GALogger::instance()->appendLogMessage("ImprovPerDimManager", msg.str());

  delete dims_with_vsmall_prob;
  delete probs;
}

int ImprovPerDimManager::posBasedOnProbs(vector<long double>& probs) const {
  int selected_pos = -1;
  
  if (probs.back() == 0.0) {
    // /*LOG*/ stringstream msg; msg << "probs in posBasedOnProbs: "; for (int i=0; i< probs.size(); i++) msg << probs[i] << ", "; msg << endl;
    // /*LOG*/ msg <<  "All probs are 0.0 " << endl;
    // /*LOG*/ GALogger::instance()->appendLogMessage("posBasedOnProbs", msg.str());
    selected_pos = GARandomInt(0,probs.size()-1);
  }
  else {
    long double rnd_value  = GARandomDouble(0.0,1.0);
    for (int pos=0; pos<probs.size(); pos++) {
      if (rnd_value <= probs[pos]) {
        selected_pos = pos;
        break;
      }
    }
  }

  assert(selected_pos != -1);

  return selected_pos;
}

//void ImprovPerDimManager::updateValues (vector<long double>& new_values) {
//  // /*LOG*/ vector<long double>  old_values = dimPerfValues();
//  // /*LOG*/ stringstream msg; msg << "old values:" << endl;
//  // /*LOG*/ for (int i=0; i<numDims(); i++) msg << old_values[i] << ", "; msg << endl;
//  // /*LOG*/ msg << "new passed values:" << endl;
//  // /*LOG*/ for (int i=0; i<numDims(); i++) msg << new_values[i] << ", "; msg << endl;
//
//  for (int i=0; i<numDims(); i++) updateValue(i,new_values[i]);
//
//  // /*LOG*/ msg << "new computedvalues:" << endl;
//  // /*LOG*/ vector<long double>  new_comp_values = dimPerfValues();
//  // /*LOG*/ for (int i=0; i<numDims(); i++) msg << new_comp_values[i] << ", "; msg << endl;
//
//  // /*LOG*/ GALogger::instance()->appendLogMessage("ImprovPerDimManager", msg.str());
//}

void ImprovPerDimManager::updateValue(int pos, long double new_value) {
  // Needed to be implemented in order to be called from the posBasedOnProbs method
  throw runtime_error("Method dimPerfValues from ImprovPerDimManager should not be called");
}

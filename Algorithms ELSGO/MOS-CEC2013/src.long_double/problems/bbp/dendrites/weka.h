/*
 * weka.h
 *
 *  Created on: 29/06/2010
 *      Author: jfernandez
 */

#ifndef WEKA_H_
#define WEKA_H_

#include <vector>
#include <map>
#include <string>

#include <boost/unordered_map.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/random.hpp>

#include <garandom.h>

using namespace std;

class LinearRegression {

public:
	LinearRegression(){};
	LinearRegression( char* fileName, vector<string>& all_attributes );
  LinearRegression(const LinearRegression& orig) {_model = orig._model;}
	long double Predict( boost::unordered_map<string, long double>& data );
	// Debug
	void Print( void );

  unsigned ModelSize() {return _model.size();}

  boost::unordered_map<string, long double>& getModel() {return _model;}

  void addModel(boost::unordered_map<string, long double>& avgModel) {
    for (boost::unordered_map<std::string, long double>::iterator it = avgModel.begin(); it != avgModel.end(); it++)
      it->second += _model[it->first];
  }

  void perturbModel(boost::unordered_map<std::string, long double>& perturbation) {
    for (boost::unordered_map<string, long double>::iterator it=_model.begin() ; it != _model.end(); it++ ){
//      std:cout << "Perturbing " << it->first << " with value: " << it->second << "  with: " << perturbation[it->first];
      it->second += perturbation[it->first];
//      std::cout << " and we obtain: " << it->second << std::endl;
    }
  }


private:
	boost::unordered_map<string, long double> _model;
};


class NormalDistribution {

public:
	NormalDistribution(){};
	NormalDistribution( long double mean, long double stdv, string attribute );
	NormalDistribution( char* fileName );
  NormalDistribution(const NormalDistribution& orig) {_terminationD = new boost::math::normal(*orig._terminationD);
                                                      _factorIndex = orig._factorIndex;
                                                      mean = orig.mean;
                                                      stdv = orig.stdv;}
	~NormalDistribution();

	long double Predict( boost::unordered_map<string, long double>& data );

  long double getMean() {return mean;}
  long double getStdDev() {return stdv;}

  void perturbModel(long double mean_pert, long double stddev_pert) {
    if (_terminationD)
      delete _terminationD;

    //std::cout << "Perturbing mean: " << mean << " with: " << mean_pert;

    mean += mean_pert;

    //std::cout << " and obtain: " << mean << std::endl;

    //std::cout << "Perturbing stdv: " << stdv << " with: " << stddev_pert;

    stdv += stddev_pert;

    //std::cout << " and obtain: " << stdv << std::endl;

  	_terminationD = new boost::math::normal(mean, stdv);
  }


private:
	string _factorIndex;
  long double mean;
  long double stdv;
	boost::math::normal* _terminationD;

};

#endif /* WEKA_H_ */

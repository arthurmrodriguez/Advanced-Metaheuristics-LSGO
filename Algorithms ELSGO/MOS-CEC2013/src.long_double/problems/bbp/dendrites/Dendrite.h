#ifndef _CONTROLDENDRITE_H_
#define _CONTROLDENDRITE_H_

//#include "attribute_names.h"

#include <list>
#include <string>
#include <iostream>
#include <fstream>
//#include <map>
#include <boost/unordered_map.hpp>

using namespace std;

#include "instance.h"
#include "weka.h"


class Dendrite{

public:
	// Se encarga de crear una dendrita completa
	Dendrite( NormalDistribution& terminationModel, LinearRegression& pathlengthModel, LinearRegression& eucDistModel);
	~Dendrite();

	void grow();
	// Almacena la dendrita en formato de instancias
	void storeInstances( char* path_file );
	// Añade en el fichero las 3 métricas
	void storeMetrics( char* path_file );
	// devuelve las métricas en un map
	boost::unordered_map<string, long double> Metrics( void );


	// Comunicación con el sistema de salida
	void printBranch( std::string mesg );
	// Comunicación con los modelos
	long double predict(boost::unordered_map<string, long double>& data, string param);
	// Cálculo de métricas
	void maxEucDist(long double dist);
	int incrementNumberOfBranches();
	void incrementLength(long double l);


private:

	Instance* _root; // La primera instancia
	std::ofstream _outputFile;

	// Modelos de predicción
	NormalDistribution& _terminationModel;
	LinearRegression& _pathlengthModel;
	LinearRegression& _eucDistModel;
	// Atributos de métricas
	long double _totalLength;
	long double _maxEucDist;
	int _numBranches;
  int _numRepetitions;

};

#endif //_CONTROLDENDRITE_H_

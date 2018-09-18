/*
 * weka.cc
 *
 *  Created on: 29/06/2010
 *      Author: jfernandez
 */


/*
=== Run information ===

Scheme:       weka.classifiers.functions.LinearRegression -S 0 -R 1.0E-8
Relation:     fold0.csv
Instances:    68569
Attributes:   12
              Tamano
              orden
              recorrido_Soma
              Distancia_euclidea
              tortuosidad
              ultima_ramificacion_dist_euc
              ultima_ramificacion_recorrido
              indice_Sholl
              termina
              numero_hijos
              tamano_del_hijo
              diferencia_entre_distancias_euclideas_padre-hijo
Test mode:    10-fold cross-validation

=== Classifier model (full training set) ===


Linear Regression Model

diferencia_entre_distancias_euclideas_padre-hijo =

      -0.0163 * Tamano +
       0.0173 * orden +
       0.0053 * recorrido_Soma +
      -0.0065 * Distancia_euclidea +
       1.4908 * tortuosidad +
       0.0127 * ultima_ramificacion_dist_euc +
      -0.0098 * ultima_ramificacion_recorrido +
       0.8487 * tamano_del_hijo +
      -1.3784


Time taken to build model: 3.35 seconds

=== Cross-validation ===
=== Summary ===

Mean absolute error                      0.4948
Root mean squared error                  0.8495
*/

#include <iostream>
#include <fstream>
#include <algorithm>
#include <boost/regex.hpp>

#include "weka.h"

boost::mt19937 randgen (static_cast<unsigned int>(std::time(0)));
boost::uniform_01<> u01;
boost::variate_generator<boost::mt19937&, boost::uniform_01<> > rand_unif01(randgen, u01);


LinearRegression::LinearRegression( char* fileName, vector<string>& all_attributes ){
	/*
	 * fileName: Fichero con la descipción del modelo que obtiene Weka
	 * all_attributes: vector de strings con todos los atributos que utiliza la regresión
	 *
	 * El atributo _model tendrá elementos con las claves que aparecen en all_attributes + error + intercept_factor
	 */

	fstream file;
	bool intoDesc = false;

	boost::regex number_re ("[-]{0,1}(\\d+.\\d+|\\d+)");
	boost::cmatch what, what2;
	double factor;

	vector<string>::iterator result;

	file.open( fileName, ios::in );
	if (! file)
		throw ("No se pudo abrir el fichero "+ string(fileName) );

	// Inicializacion del modelo
	for( vector<string>::const_iterator it = all_attributes.begin(); it != all_attributes.end(); ++it ) {
		_model[*it] = 0.0;
	}
	_model["intercept_factor"] = 0.0;

	// Cargar el modelo
	char input[4096]; // Auxiliar para cada linea del fichero
	while( file.getline( input, 4096 ) ){

		if ( boost::regex_search( input, boost::regex("Linear Regression Model")) ){
			intoDesc = true; // Estamos dentro de la especificación de los factores del modelo
		}
		else if( boost::regex_search( input, boost::regex("Time taken to build model")) ){
			intoDesc = false;
		}
		else if( intoDesc ){
			// Hay un factor asociado a un atributo
			if ( boost::regex_search( input , what, number_re ) ) {
				factor = atof(string(what[0].first, what[0].second).c_str());
				// Buscar nombre del atributo si lo tiene
				if (boost::regex_search( what.suffix().first , what2, boost::regex("\\w+"))){
					result = find( all_attributes.begin(), all_attributes.end(), what2[0] );
					// Si no existe error, no coinciden los atributos especificados con los encontrados
					if( result == all_attributes.end() )
						throw runtime_error("Error: '"+ what2[0] +"' attribute not recognized.");
					else if ( _model[what2[0]] != 0.0 )
						throw runtime_error("Error: '"+ what2[0] +"' factor already defined.");
					else
						_model[what2[0]] = factor;
				}
				else // Sino es el termino independiente
					_model["intercept_factor"] = factor;
			}
		}
		// El Error
		else if( boost::regex_search( input, what, boost::regex("Root mean squared error")) ){
			boost::regex_search( input, what, number_re);
			_model["error"] = atof(string(what[0].first, what[0].second).c_str());
		}
	}
	if ( _model.find("error") == _model.end() )
		throw runtime_error( "Error: 'Root mean squared error' value not found" );

}

// Print _model: DEBUG
void LinearRegression::Print( void ){
	boost::unordered_map<string, double>::iterator it;
	for ( it=_model.begin() ; it != _model.end(); it++ ){
		cout << (*it).first << " => " << (*it).second << endl;
	}
	cout << "Len:" << _model.size() << endl;
}


double LinearRegression::Predict( boost::unordered_map<string, double>& data ){
	/*
	 * data: map de string -> double, que representa la instancia a partir de la cual se quiere predecir la clase del modelo
	 *
	 * Devuelve la predicción realizada por el modelo a partir de la instancia "data".	 *
	 */

	double acumul = 0.0, error = 0.0, res;

	for(boost::unordered_map<string, double>::iterator it=data.begin() ; it != data.end(); it++ ){
		acumul += data[(*it).first] * _model[(*it).first];
	}
	acumul += _model["intercept_factor"];
	error = _model["error"];
	res = acumul + (-1*error - error)*rand_unif01() + (-1*error);

	// No puede ser negativo ni el inc. de la dist.euc. ni el del tamaño
	int loops = 10;
	boost::normal_distribution<double> paramND(0, error);
	boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > nD(randgen, paramND);
	while( res < 0 and loops > 0){
		res = acumul + nD();
		loops -= 1;
	}
	if( res < 0)
		res = abs(nD());

	return res;
}


NormalDistribution::NormalDistribution( double mean, double stdv, string attribute )
	: _factorIndex(attribute){

	_terminationD = new boost::math::normal(mean, stdv);
}


NormalDistribution::NormalDistribution( char* fileName )
	:_factorIndex("recorrido_soma")
{
	fstream file;

	boost::regex number_re ("[-]{0,1}(\\d+.\\d+|\\d+)");
	boost::cmatch what;

	file.open( fileName, ios::in );
	if (! file) {
		std::cerr << "No se pudo abrir el fichero " << fileName << std::endl;
		exit(-1);
	}
	char input[4096];
	while( file.getline( input, 4096 ) ){

		if ( boost::regex_search( input, boost::regex("mean")) ){
			if( boost::regex_search( input, what, number_re) )
				mean = atof(string(what[0].first, what[0].second).c_str());
			else
				throw ("mean value does not found");
		}
		else if( boost::regex_search( input, boost::regex("stdv")) ){
			if( boost::regex_search( input, what, number_re) )
				stdv = atof(string(what[0].first, what[0].second).c_str());
			else
				throw ("stdv value does not found");
		}
		else if( boost::regex_search( input, boost::regex("\\w")) ){
			_factorIndex = string(input);
		}
	}
	_terminationD = new boost::math::normal(mean, stdv);
}

NormalDistribution::~NormalDistribution(){
  if (_terminationD)
    delete _terminationD;
}


double NormalDistribution::Predict( boost::unordered_map<string, double>& data  ){

	double p = boost::math::cdf( *_terminationD, data[_factorIndex] );
	double rnd = rand_unif01();

	return ( rnd <= p ) ? 1.0 : 0.0;

}

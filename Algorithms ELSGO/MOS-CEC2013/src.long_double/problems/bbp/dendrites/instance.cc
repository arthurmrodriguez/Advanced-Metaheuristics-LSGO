/*
 * instance.cpp
 *
 *  Created on: 25/06/2010
 *      Author: jfernandez
 */
#include <iostream>
#include <sstream>

#include "attribute_names.h"
#include "Dendrite.h"


Instance::Instance( boost::unordered_map<string, long double>& data/*, ??? models*/, Dendrite* d )
	:_data(data),
	 _dendrite(d)
{
	_dendrite->incrementLength( data[str_recorrido_segmento] );
	_dendrite->maxEucDist( data[str_distancia_euclidea] );
}

Instance::~Instance(){
	for(unsigned i=0; i < _nexts.size(); i++)
		delete _nexts[i];
}

// crea las siguientes instancias a la actual
void Instance::createNexts( void ){

	/* Version SIN crecimiento en anchura (sin indice de Sholl:
	 * Crea los siguientes segmentos, haciéndolos crecer (llamando al
	 * método createNexts) hasta la finalización de la rama (la predicción
	 * de termina)
	 */
	boost::unordered_map<string, long double> dt;
	Instance* inst;
	long double incRec, incEuc;
	long double termina = _dendrite->predict( _data, str_termina );
	unsigned count = 0;

	// crear dos nuevas ramas si no es una terminacion
	if( termina == 0.0){

		// Si llega a mas de 20 hay que abortar la dendrita
		if( _dendrite->incrementNumberOfBranches() > 20 )
			return;

		for( unsigned i=0; i < NUM_HIJOS; i++){
			do {
				incRec = _dendrite->predict( _data, str_inc_recorrido );
				_data[str_inc_recorrido] = incRec;
				incEuc = _dendrite->predict( _data, str_inc_dist_euc );
			}while( (incRec < incEuc) and (count++ < 3));
			count = 0;
			if (incRec < incEuc){
				incEuc = incRec;
			}

			dt[str_recorrido_segmento] = incRec;
			dt[str_orden] = _data[str_orden] + 1;
			dt[str_recorrido_soma] = _data[str_recorrido_soma] + incRec;
			dt[str_distancia_euclidea] = _data[str_distancia_euclidea] + incEuc;
			dt[str_tortuosidad] = dt[str_distancia_euclidea] / dt[str_recorrido_soma];
			dt[str_distancia_euclidea_segmento] = incEuc;

			inst = new Instance(dt, _dendrite);
			_nexts.push_back(inst);
			inst->createNexts();
		}


	}

}

void Instance::storeInstance(){
	std::stringstream msg;

	// Se ha incluido el atributo str_inc_recorrido en el método de creación de las ramas hijas y no hay que mostrarlo
	for(unsigned i=0; i < LAST_ATT_I; i++)
		msg << _data[attributes[i]] << ", ";
	//for( boost::unordered_map<string, long double>::const_iterator it = _data.begin(); it != _data.end(); ++it )
	//	msg << (*it).second << ", ";

	if( _nexts.size() > 0){
		// Imprimir la info de los hijos
		for( unsigned i=0; i < _nexts.size(); i++){
			_dendrite->printBranch( msg.str() );
			_nexts[i]->printClassValuesForFather();
		}
	}
	else{ // Segmento de terminación
		msg << "1, 0.0, 0.0" << endl;;
		_dendrite->printBranch( msg.str() );
	}

	// Pintar las siguientes
	for( unsigned i=0; i < _nexts.size(); i++)
		_nexts[i]->storeInstance();

}

void Instance::printClassValuesForFather(){
	std::stringstream msg;
	msg << "0" << ", " << _data[str_recorrido_segmento] << ", " << _data[str_distancia_euclidea_segmento] << std::endl;
	_dendrite->printBranch( msg.str() );
}

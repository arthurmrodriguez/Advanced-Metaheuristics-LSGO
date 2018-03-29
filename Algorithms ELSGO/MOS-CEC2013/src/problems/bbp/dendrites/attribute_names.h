/*
 * attribute_names.h
 *
 *  Created on: 30/06/2010
 *      Author: jfernandez
 */

#ifndef __ATTRIBUTE_NAMES_H_
#define __ATTRIBUTE_NAMES_H_

#include <string>
#include <vector>

using namespace std;

// metricas
static const string str_numero_branches ("numero_ramificaciones");
static const string str_longitud_total ("longitud_total");
static const string str_max_distancia_euclidea ("max_distancia_euclidea");

// atributos y clases
static const string str_orden ("orden");
static const string str_recorrido_soma("recorrido_soma");
static const string str_distancia_euclidea("distancia_euclidea_soma");
static const string str_tortuosidad("tortuosidad");
static const string str_distancia_euclidea_segmento("distancia_euclidea_segmento");
static const string str_recorrido_segmento("recorrido_segmento");
static const string str_termina("termina");
static const string str_inc_recorrido("inc_recorrido");
static const string str_inc_dist_euc("inc_distancia_euclidea");

static const string NAMES[] = {//str_tamano,
						 str_orden,
						 str_recorrido_soma,
						 str_distancia_euclidea,
						 str_tortuosidad,
						 str_recorrido_segmento,
						 str_distancia_euclidea_segmento,
						 str_termina,
						 str_inc_recorrido,
						 str_inc_dist_euc
};

// vector con los atributos en orden
static vector<string> attributes( NAMES, NAMES + sizeof(NAMES)/sizeof(string));
// indice donde se encuentra el último atributo
static const unsigned LAST_ATT_I = 6;


#endif /* __ATTRIBUTE_NAMES_H_ */

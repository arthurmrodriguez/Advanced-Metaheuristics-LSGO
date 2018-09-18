#ifndef HEURISTICOINSERCION_H
#define HEURISTICOINSERCION_H

#include "vehiculo.h"
#include "peticionInsercion.h"
#include "../VerticesList.h"
#include "../DistMatrix.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <cstdlib>
#include <time.h>
#include <math.h>

using namespace std;

	/*void heuristicoAdicion(Datos &d,std::vector<PeticionInsercion> &candidatas, std::vector<Vehiculo> &vehiculos);
	int long long evaluaRuta(std::vector <OrdenacionPeticiones> &solucion, Datos &d);
	bool insertaPeticionPositiva(Datos &d,PeticionInsercion &candidata,std::vector<Vehiculo> &vehiculos);
	bool insertaPeticionNegativa(Datos &d,PeticionInsercion &candidata,std::vector<Vehiculo> &vehiculos);*/
	void heuristicoInicializacion(std::vector<Vehiculo> &vehiculos,long vehicleCapacity,long maxDelay, double DARP_ALPHA, double DARP_BETA, VerticesList& verticeslist, long nVehicles, DistMatrix& distmatrix);
	/*int buscarPeticionMasTempranaAleatoria(std::vector<PeticionInsercion> &c);
	void insertaEnSolucionPrimera(Vehiculo &v,std::vector<PeticionInsercion> &candidatas,int horaSubida,int horaBajada,int horaBajadaM,int indice,Datos &d);
	std::vector<OrdenacionPeticiones> insertaEnSolucionSubida(PeticionInsercion &candidata,std::vector<OrdenacionPeticiones> &ordenacion,
									int horaSubida,int indice,Datos &d,bool &viola);
	std::vector<OrdenacionPeticiones> insertaEnSolucionBajada(PeticionInsercion &candidata,std::vector<OrdenacionPeticiones> &ordenacion,
									int horaBajada,int indice,Datos &d, bool &viola);
	std::vector<PeticionInsercion> borrarCandidato(std::vector<PeticionInsercion> v,int borrado);*/



#endif

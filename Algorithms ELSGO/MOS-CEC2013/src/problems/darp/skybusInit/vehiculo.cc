#include "vehiculo.h"

Vehiculo::Vehiculo() {
	// TODO Auto-generated constructor stub

}

int Vehiculo::getPlazas() const{
    return nPlazas;
}

std::vector<OrdenacionPeticiones> Vehiculo::getSolucion() const{
    return solucion;
}

void Vehiculo::setPlazas(int nPlazas){
    this->nPlazas = nPlazas;
}

int Vehiculo::getIdVehiculo() const{
    return idVehiculo;
}

void Vehiculo::setAsignadas(std::vector<PeticionInsercion> asignadas){
    this->asignadas = asignadas;
}

std::vector<PeticionInsercion> Vehiculo::getAsignadas() const{
    return asignadas;
}

void Vehiculo::setIdVehiculo(int idVehiculo){
    this->idVehiculo = idVehiculo;
}

void Vehiculo::setSolucion(std::vector<OrdenacionPeticiones> solucion){
    this->solucion = solucion;
}

Vehiculo::~Vehiculo() {
	// TODO Auto-generated destructor stub
}

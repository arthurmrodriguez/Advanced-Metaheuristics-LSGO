#include "ordenacionPeticiones.h"

OrdenacionPeticiones::OrdenacionPeticiones() {
	// TODO Auto-generated constructor stub

}

long OrdenacionPeticiones::getHoraDeseada() const{
    return horaDeseada;
}

long OrdenacionPeticiones::getHoraRecogida() const{
    return horaRecogida;
}

int OrdenacionPeticiones::getIdPeticion() const{
    return idPeticion;
}

void OrdenacionPeticiones::setHoraDeseada(long horaDeseada){
    this->horaDeseada = horaDeseada;
}

void OrdenacionPeticiones::setHoraRecogida(long horaRecogida){
    this->horaRecogida = horaRecogida;
}

int OrdenacionPeticiones::getIdParada() const{
    return idParada;
}

int OrdenacionPeticiones::getIdParadaPareja() const{
    return idParadaPareja;
}

int OrdenacionPeticiones::getPlazasLibres() const{
    return nPlazasLibres;
}

int OrdenacionPeticiones::getNumPax() const{
    return numPax;
}

void OrdenacionPeticiones::setNumPax(int numPax){
    this->numPax = numPax;
}

void OrdenacionPeticiones::setPlazasLibres(int nPlazasLibres){
    this->nPlazasLibres = nPlazasLibres;
}

void OrdenacionPeticiones::setIdParadaPareja(int idParadaPareja){
    this->idParadaPareja = idParadaPareja;
}

void OrdenacionPeticiones::setIdParada(int idParada){
    this->idParada = idParada;
}

void OrdenacionPeticiones::setIdPeticion(int idPeticion){
    this->idPeticion = idPeticion;
}

OrdenacionPeticiones::~OrdenacionPeticiones() {
	// TODO Auto-generated destructor stub
}

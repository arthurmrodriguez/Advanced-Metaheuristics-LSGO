#include "peticionInsercion.h"

PeticionInsercion::PeticionInsercion() {
	// TODO Auto-generated constructor stub

}
long PeticionInsercion::getHoraRecogidaBajada() const{
    return horaRecogidaBajada;
}

long PeticionInsercion::getHoraRecogidaSubida() const{
    return horaRecogidaSubida;
}

int PeticionInsercion::getIdParadaBajada() const{
    return idParadaBajada;
}

int PeticionInsercion::getIdParadaSubida() const{
    return idParadaSubida;
}

int PeticionInsercion::getIdPeticion() const{
    return idPeticion;
}

int PeticionInsercion::getNumPax() const{
    return numPax;
}

void PeticionInsercion::setHoraRecogidaBajada(long horaRecogidaBajada){
    this->horaRecogidaBajada = horaRecogidaBajada;
}

void PeticionInsercion::setHoraRecogidaSubida(long horaRecogidaSubida){
    this->horaRecogidaSubida = horaRecogidaSubida;
}

void PeticionInsercion::setIdParadaBajada(int idParadaBajada){
    this->idParadaBajada = idParadaBajada;
}

void PeticionInsercion::setIdParadaSubida(int idParadaSubida){
    this->idParadaSubida = idParadaSubida;
}

void PeticionInsercion::setIdPeticion(int idPeticion){
    this->idPeticion = idPeticion;
}

void PeticionInsercion::setNumPax(int numPax){
    this->numPax = numPax;
}

PeticionInsercion::~PeticionInsercion() {
	// TODO Auto-generated destructor stub
}

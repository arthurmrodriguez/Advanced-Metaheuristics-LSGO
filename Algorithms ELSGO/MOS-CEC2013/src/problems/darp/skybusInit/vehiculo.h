#ifndef VEHICULO_H_
#define VEHICULO_H_


#include <vector>
#include "peticionInsercion.h"
#include "ordenacionPeticiones.h"


class Vehiculo {

private:
	int idVehiculo;
	std::vector<OrdenacionPeticiones> solucion;
	std::vector<PeticionInsercion> asignadas;
	int nPlazas;


public:
	Vehiculo();
	virtual ~Vehiculo();
    int getPlazas() const;
    int getPlazasLibres() const;
    std::vector<OrdenacionPeticiones> getSolucion() const;
    void setPlazas(int nPlazas);
    void setPlazasLibres(int nPlazasLibres);
    void setSolucion(std::vector<OrdenacionPeticiones> solucion);
    int getIdVehiculo() const;
    void setIdVehiculo(int idVehiculo);
    std::vector<PeticionInsercion> getAsignadas() const;
    void setAsignadas(std::vector<PeticionInsercion> asignadas);
};

#endif

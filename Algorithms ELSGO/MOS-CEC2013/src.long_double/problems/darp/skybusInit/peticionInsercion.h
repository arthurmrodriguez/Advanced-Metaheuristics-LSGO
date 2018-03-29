#ifndef PETICIONINSERCION_H_
#define PETICIONINSERCION_H_


class PeticionInsercion {

private:
	int idPeticion; //id que tiene la peticion en el array de peticiones cargado en la clase datos (desde 0 hasta nPeticiones-1)
	long horaRecogidaSubida; //hora en segundos a la que quiere ser recogida el cliente
	long horaRecogidaBajada; //hora en segundos a la que quiere ser recogida el cliente
	int idParadaSubida;
	int idParadaBajada;
	int numPax;

public:
	PeticionInsercion();
	virtual ~PeticionInsercion();
    long getHoraRecogidaBajada() const;
    long getHoraRecogidaSubida() const;
    int getIdParadaBajada() const;
    int getIdParadaSubida() const;
    int getIdPeticion() const;
    int getNumPax() const;
    void setHoraRecogidaBajada(long horaRecogidaBajada);
    void setHoraRecogidaSubida(long horaRecogidaSubida);
    void setIdParadaBajada(int idParadaBajada);
    void setIdParadaSubida(int idParadaSubida);
    void setIdPeticion(int idPeticion);
    void setNumPax(int numPax);

};

#endif

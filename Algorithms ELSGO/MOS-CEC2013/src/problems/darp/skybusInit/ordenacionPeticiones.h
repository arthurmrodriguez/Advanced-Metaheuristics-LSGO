#ifndef ORDENACIONPETICIONES_H_
#define ORDENACIONPETICIONES_H_

class OrdenacionPeticiones {

private:
	int idPeticion;//id que tiene la peticion en el array de peticiones cargado en la clase datos (desde 0 hasta nPeticiones-1)
	long horaDeseada;  //hora en segundos a la que quiere ser recogida el cliente
	long horaRecogida; //hora en segundos a la que es realmente recogido el peaton
	int idParada;
	int idParadaPareja;
	int nPlazasLibres;
	int numPax;
public:
	OrdenacionPeticiones();
	virtual ~OrdenacionPeticiones();
    long getHoraDeseada() const;
    long getHoraRecogida() const;
    int getIdPeticion() const;
    void setHoraDeseada(long horaDeseada);
    void setHoraRecogida(long horaRecogida);
    void setIdPeticion(int idPeticion);
    int getIdParada() const;
    void setIdParada(int idParada);
    int getIdParadaPareja() const;
    void setIdParadaPareja(int idParadaPareja);
    int getPlazasLibres() const;
    void setPlazasLibres(int nPlazasLibres);
    int getNumPax() const;
    void setNumPax(int numPax);

};

#endif

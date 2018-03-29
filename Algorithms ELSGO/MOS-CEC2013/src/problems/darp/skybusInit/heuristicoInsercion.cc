#include "heuristicoInsercion.h"
#include "../../../gaeda/garandom.h"

/*********************************************************************************************************************
 *  Funcion encargada de buscar la peticion con hora de recogida mas temprana
 */
int buscarPeticionMasTempranaAleatoria(std::vector<PeticionInsercion> &c){
	long minimo=10000000000;
	int indice=-1;
	long hora=0;
	std::vector<int> indices;
	//buscamos el minimo
	for(int i=0;i<c.size();i++){
		hora=c[i].getHoraRecogidaSubida();
		//cout<<"hora "<< hora<< endl;
		if(hora<minimo){
			minimo=hora;
		}
	}
	//metemos en un vector todas las horas minimas
	for(int i=0;i<c.size();i++){
		hora=c[i].getHoraRecogidaSubida();
		if(hora==minimo){
			indices.push_back(i);
		}
	}

	if(indices.size()>0){
	  srand(GAGetRandomNoRankSeed());
		int aleatorio=rand() % (indices.size());
		//cout<<"tamanio size"<<indices.size()<<endl;
		//cout<<"aleatorio"<<aleatorio<<endl;
		indice=indices[aleatorio];
	}

	return indice;
}
/***********************************************************************************************************
 * Funcion encargarda de eliminar una peticion de las peticiones asignadas debido a q ya ha sido enrutada
 */
std::vector<PeticionInsercion> borrarCandidato(std::vector<PeticionInsercion> v,int borrado){
	std::vector<PeticionInsercion> copia;
	for (int i = 0; i < v.size(); i++){
		if (i != borrado){
			copia.push_back(v[i]);
		}
	}
	//printf("size peticiones %d\n",v.size());
	//printf("borrado %d\n",borrado);
	//printf("size copia %d\n",copia.size());
	return copia;
}


/*****************************************************************************************************************
 *  Funcion encargada de insertar en la solucion del vehiculo una peticion elegida como la siguiente
 */
void insertaEnSolucionPrimera(Vehiculo &v,std::vector<PeticionInsercion> &candidatas,long horaSubida,long horaBajada,int indice){

			std::vector<OrdenacionPeticiones> ordenacion;
			ordenacion=v.getSolucion();
			//cout<<"size solucion antes "<<v.getSolucion().size()<<endl;

			///////////////////////////////////////insertamos la subida del elemento en la solucion
			OrdenacionPeticiones subida;
			subida.setHoraRecogida(horaSubida);
			subida.setHoraDeseada(candidatas[indice].getHoraRecogidaSubida());
			subida.setIdPeticion(candidatas[indice].getIdPeticion());
			subida.setIdParada(candidatas[indice].getIdParadaSubida());
			subida.setIdParadaPareja(candidatas[indice].getIdParadaBajada());
			//cout<<"plazas iniciales "<<v.getPlazasLibres()<<endl;
			subida.setPlazasLibres(v.getPlazas()-candidatas[indice].getNumPax());
			subida.setNumPax(candidatas[indice].getNumPax());
			///cout<<"plazas despues de la subida "<<subida.getPlazasLibres()<<endl;
			ordenacion.push_back(subida);


			//////////////////////////////////insertamos la bajada del elemento en la solucion
			OrdenacionPeticiones bajada;
			bajada.setHoraRecogida(horaBajada);
			bajada.setHoraDeseada(candidatas[indice].getHoraRecogidaBajada());
			bajada.setIdPeticion(candidatas[indice].getIdPeticion()*-1);
			bajada.setIdParada(candidatas[indice].getIdParadaBajada());
			bajada.setIdParadaPareja(candidatas[indice].getIdParadaSubida());
			//cout<<"plazas antes de la bajada "<<subida.getPlazasLibres()<<endl;
			bajada.setPlazasLibres(subida.getPlazasLibres()+candidatas[indice].getNumPax());
			bajada.setNumPax(candidatas[indice].getNumPax());
			//cout<<"plazas despues de la bajada "<<bajada.getPlazasLibres()<<endl;
			ordenacion.push_back(bajada);

			v.setSolucion(ordenacion);
			//cout<<"size solucion despues "<<v.getSolucion().size()<<endl;
			/////////////////////// anadimos la primera peticion
			std::vector<PeticionInsercion> a;
			a.push_back(candidatas[indice]);
			v.setAsignadas(a);

			candidatas=borrarCandidato(candidatas,indice);
			//printf("tamanio candidatas despues de borrar %d\n",candidatas.size());

}

/*****************************************************************************************************************
 *  Funcion encargada de insertar en la solucion del vehiculo una peticion elegida como la siguiente
 */
std::vector<OrdenacionPeticiones> insertaEnSolucionSubida(PeticionInsercion &candidata,std::vector<OrdenacionPeticiones> &ordenacion,
															int horaSubida,int indice,bool &viola,DistMatrix& distmatrix,long vehicleCapacity,long maxDelay){

	//cout<<"tamanio ordenacion "<<ordenacion.size()<<endl;

	//actualizamos el tiempo deseado de bajada
	int tiempo=distmatrix.getDist(candidata.getIdParadaSubida(),candidata.getIdParadaBajada());
	candidata.setHoraRecogidaBajada(horaSubida+tiempo);

	viola=false;
	std::vector<OrdenacionPeticiones> ord;

	//cout<<"--------------------- INSERCION SUBIDA --------------------"<<endl;

	OrdenacionPeticiones subida;
	subida.setHoraRecogida(horaSubida);
	subida.setHoraDeseada(candidata.getHoraRecogidaSubida());
	subida.setIdPeticion(candidata.getIdPeticion());
	subida.setIdParada(candidata.getIdParadaSubida());
	subida.setIdParadaPareja(candidata.getIdParadaBajada());
	subida.setPlazasLibres(ordenacion[indice].getPlazasLibres()-candidata.getNumPax());
	subida.setNumPax(candidata.getNumPax());

	//cout<<"indice "<<indice<<endl;

	//hacemos copia desde 0 hasta indice de lo que hay en la ruta hasta ahora
	for(int i=0;i<=indice;i++){
		ord.push_back(ordenacion[i]);
	}

	//introducimos la subida
	ord.push_back(subida);


	//copiamos el resto de la ruta cambiando tiempos y comprobando si no se jode
	int coste;
	int tiempoLlegada;
	int tiempoViolacion;

	for(int j=indice+1;j<ordenacion.size() && !viola;j++){
		//cout<<"j "<<j<<endl;
		OrdenacionPeticiones aux=ord.at(ord.size()-1);
		OrdenacionPeticiones aux2=ordenacion[j];

		int plazasLibres=aux2.getPlazasLibres()-candidata.getNumPax();
		if(plazasLibres>=0){

			//cout<<"vamos de "<<aux.getIdParada()<<" a "<<aux2.getIdParada()<<endl;

			coste=distmatrix.getDist(aux.getIdParada(),aux2.getIdParada());
			//cout<<"coste "<<coste<<endl;
			tiempoLlegada=aux.getHoraRecogida()+coste;
			//cout<<"tiempo llegada nuevo "<<tiempoLlegada<<endl;
			//cout<<"tiempo llegada antiguo"<<aux2.getHoraRecogida()<<endl;

			tiempoViolacion=tiempoLlegada-aux2.getHoraDeseada();
			//cout<<"tiempo violacion  "<<tiempoViolacion<<endl;

			//comprobamos la violacion
			if(aux2.getIdPeticion()>0){//subida
				//cout<<"es una subida "<<endl;
				if(tiempoViolacion<0){
					//cout<<"llegamos antes esperamos"<<endl;
					tiempoLlegada=aux2.getHoraDeseada();
				}

				if(tiempoViolacion>maxDelay){
					viola=true;
				}

			}else{//bajada
				//cout<<" es una bajada "<<endl;
				int tiempoDirecto=distmatrix.getDist(aux2.getIdParadaPareja(), aux2.getIdParada());
				int tiempoDuracionTrayecto=tiempoViolacion+tiempoDirecto;
				//cout<<"tiempo duracion trayecto "<<tiempoDuracionTrayecto<<endl;

				/////////////////////////////////////////////////////////////////
				//filtramos por tipo de violacion bajada

				//if(d.getTipoViolacionBajada()==1){	//violacion fija (60 min)
					if(tiempoDuracionTrayecto>Vertex::maxRideTime(tiempoDirecto)) {
						viola=true;
					}//no viola candidata bajada
				//}//violacion fija (60 min)


				/*if(d.getTipoViolacionBajada()==2){ //violacion bajada variable
					if(tiempoDuracionTrayecto>d.getViolacionVariableBajada()*tiempoDirecto){
						viola=true;
					}//no viola candidata bajada
				}//violacion bajada variable*/

				/////////////////////////////////////////////////////////////////
			}
		}else{
			viola=true;
		}

		if(!viola){// si no viola la metemos
			//cout<<"no viola"<<endl;
			//actualizamos el numero de pasajeros
			aux2.setPlazasLibres(plazasLibres);
			//acutalizamos la hora de recogida
			aux2.setHoraRecogida(tiempoLlegada);
			ord.push_back(aux2);
		}else{
			//cout<<"viola el resto de la ruta"<<endl;
		}
	}

	//cout<<"ord size final "<<ord.size()<<endl;
	return ord;



}

/*****************************************************************************************************************
 *  Funcion encargada de insertar en la solucion del vehiculo una peticion elegida como la siguiente
 */
std::vector<OrdenacionPeticiones> insertaEnSolucionBajada(PeticionInsercion &candidata,std::vector<OrdenacionPeticiones> &ordenacion,
 											int horaBajada,int indice,bool &viola, DistMatrix& distmatrix, long vehicleCapacity,long maxDelay){

 	viola=false;
 	std::vector<OrdenacionPeticiones> ord;
 	//cout<<"--------------------- INSERCION BAJADA --------------------"<<endl;

 	OrdenacionPeticiones bajada;
 	bajada.setHoraRecogida(horaBajada);
 	bajada.setHoraDeseada(candidata.getHoraRecogidaBajada());
 	bajada.setIdPeticion(candidata.getIdPeticion()*-1);
 	bajada.setIdParada(candidata.getIdParadaBajada());
 	bajada.setIdParadaPareja(candidata.getIdParadaSubida());
 	bajada.setPlazasLibres(ordenacion[indice].getPlazasLibres()+candidata.getNumPax());
 	bajada.setNumPax(candidata.getNumPax());


 	//hacemos copia desde 0 hasta indice de lo que hay en la ruta hasta ahora
 	for(int i=0;i<=indice;i++){
 			//cout<<"plazas libres hasta indice "<<ordenacion[i].getPlazasLibres()<<endl;
 			ord.push_back(ordenacion[i]);
 	}
 	//introducimos la subida
 	ord.push_back(bajada);
 	//cout<<"plazas libres despues de la bajada "<<bajada.getPlazasLibres()<<endl;

 	//copiamos el resto de la ruta cambiando tiempos y comprobando si no se jode
 	int coste;
 	int tiempoLlegada;
 	int tiempoViolacion;

 	for(int j=indice+1;j<ordenacion.size() && !viola;j++){
 			OrdenacionPeticiones aux=ord.at(ord.size()-1);
 			OrdenacionPeticiones aux2=ordenacion[j];

 			//cout<<"...........bajada "<<j<<endl;
 			int plazasLibres=aux2.getPlazasLibres()+candidata.getNumPax();
 			//cout<<"plazas libres resto ruta bajada bajada "<<plazasLibres<<endl;
 			//cout<<"plazas devueltas "<<candidata.getNumPax()<<endl;

 			//cout<<"vamos de "<<aux.getIdParada()<<" a "<<aux2.getIdParada()<<endl;

 			coste=distmatrix.getDist(aux.getIdParada(),aux2.getIdParada());
 			//cout<<"coste "<<coste<<endl;
 			tiempoLlegada=aux.getHoraRecogida()+coste;
 			tiempoViolacion=tiempoLlegada-aux2.getHoraDeseada();
 			//cout<<"tiempo violacion al insertar la subida candidata "<<tiempoViolacion<<endl;


 			//comprobamos la violacion
 			if(aux2.getIdPeticion()>0){//subida
 				//cout<<"estamos en una subida"<<endl;
 				if(tiempoViolacion<0){
 					//cout<<"llegamos antes esperamos"<<endl;
 					tiempoLlegada=aux2.getHoraDeseada();
 				}
 				//cout<<"violacion permitida subida "<<d.getMargenTiempoRecogida()<<endl;
 				if(tiempoViolacion>maxDelay){
 					viola=true;
 				}
 			}else{//bajada
 				int tiempoDirecto=distmatrix.getDist(aux2.getIdParadaPareja(), aux2.getIdParada());
 				int tiempoDuracionTrayecto=tiempoViolacion+tiempoDirecto;

 				/////////////////////////////////////////////////////////////////
 				//filtramos por tipo de violacion bajada

 				//if(d.getTipoViolacionBajada()==1){	//violacion fija (60 min)
 					if(tiempoDuracionTrayecto>Vertex::maxRideTime(tiempoDirecto)){
 						viola=true;
 					}//no viola candidata bajada
 				//}//violacion fija (60 min)

 				/*if(d.getTipoViolacionBajada()==2){ //violacion bajada variable
 					if(tiempoDuracionTrayecto>d.getViolacionVariableBajada()*tiempoDirecto){
 						viola=true;
 					}//no viola candidata bajada
 				}//violacion bajada variable*/

 				/////////////////////////////////////////////////////////////////

 			}
 			if(!viola){// si no viola la metemos
 				//cout<<"no viola"<<endl;
 				aux2.setHoraRecogida(tiempoLlegada);
 				//acutalizamos el numero de pasajeros
 				aux2.setPlazasLibres(aux2.getPlazasLibres()+candidata.getNumPax());
 				ord.push_back(aux2);
 			}else{
 				//cout<<"viola el resto de la ruta"<<endl;
 			}
 		}

 	return ord;

 }

/*****************************************************************************************************************
  *  Funcion encargada de insertar en la solucion del vehiculo una peticion elegida como la siguiente
  */
void transformaPeticiones(std::vector<PeticionInsercion> &candidatas, VerticesList& verticeslist){
	 std::vector<Request> requests=verticeslist.getReqsList();
	//realizamos la transformacion de las peticiones en subidas y bajadas
	for(int j=0;j<requests.size();j++){
		////////////////////////////////////////////////primero hago las subidas
		PeticionInsercion candidata;
		Request peticion =requests[j];
		candidata.setIdPeticion(peticion.pickup_vert->id_);
		candidata.setNumPax(peticion.pickup_vert->load_);
		candidata.setIdParadaSubida(peticion.pickup_vert->pos_);
		candidata.setHoraRecogidaSubida(peticion.pickup_vert->fbegin_);
		//////////////////////////////////////////////ahora hago las bajadas
		candidata.setIdParadaBajada(peticion.delivery_vert->pos_);
		candidata.setHoraRecogidaBajada(peticion.delivery_vert->fbegin_);

		//std:cout<<"id peticion "<<candidata.getIdPeticion() <<" numpax "<<candidata.getNumPax()<<" parada subida "<<candidata.getIdParadaSubida()<<
			//	" parada bajada "<<candidata.getIdParadaBajada()<<" hora subida "<<candidata.getHoraRecogidaSubida()<<" hora bajada "<<candidata.getHoraRecogidaBajada()<<std::endl;
		candidatas.push_back(candidata);



	}

}

/*******************************************************************************************
 * Funcion encargada de realizar el fitness de la ruta generada,
 * inicialmente
 ******************************************************************************************/
//long int 
 int long long evaluaRuta(std::vector <OrdenacionPeticiones> &solucion,double DARP_ALPHA,double DARP_BETA){

    int long long coste=0;
    int long long retrasoSubida=0;
    int long long retrasoBajada=0;
	for(int i=0;i<solucion.size();i++){
		int retraso=solucion[i].getHoraRecogida()-solucion[i].getHoraDeseada();
		//cout<<"i "<<i<<" retraso "<<retraso<<"---------";
		if(retraso<0){
			retraso=0;
		}
		//vemos si son subidas o bajadas y en funcion de eso realizamos la penalizacion
		if(solucion[i].getIdPeticion()>0){//subidas
			int long long retrasocuartaSubida=pow((double)retraso,4);
            // cout<<" retraso cuarta subida"<<retrasocuartaSubida<<"---------";
			retrasoSubida+=retrasocuartaSubida;

		}else{//bajadas
			int long long retrasocuartaBajada=pow((double)retraso,4);
            //cout<<" retraso cuarta bajada"<<retrasocuartaBajada<<"---------";
			retrasoBajada+=retrasocuartaBajada;
		}
	}

	//divido por el numero de peticiones
	retrasoSubida=retrasoSubida;//d.getNpeticiones();
	retrasoBajada=retrasoBajada;//d.getNpeticiones();

	//hacemos la raiz
	int long long raizSubida=sqrt(retrasoSubida);
	int long long raizBajada=sqrt(retrasoBajada);

	raizSubida=sqrt(raizSubida);
	raizBajada=sqrt(raizBajada);

	//cout<<" raiz subida "<<raizSubida<<"---------";
	//cout<<" raiz bajada "<<raizBajada<<endl;

	coste=(DARP_ALPHA*raizSubida)+(DARP_BETA*raizBajada);

	return coste;
}


 /**************************************************************************************************
  *
  */
  bool insertaPeticionNegativa(PeticionInsercion &candidata,std::vector<Vehiculo> &vehiculos,double DARP_ALPHA,double DARP_BETA,
		  	  	  	  	  	  	  long vehicleCapacity,long maxDelay, DistMatrix& distmatrix){
 	 std::vector <OrdenacionPeticiones> solucionVirtualMinima;

 	 PeticionInsercion insertada;
 	 int long long costeMinimo=1000000;
 	 bool encajada=false;
 	 int idVehiculo;


 	 //vemos en que vehiculo es mejor insertarla
 	 for(int j=0;j<vehiculos.size();j++){
 		 std::vector <OrdenacionPeticiones> solucion=vehiculos[j].getSolucion();
 		 //evaluamos la ruta
 		 int long long coste=evaluaRuta(solucion,DARP_ALPHA, DARP_BETA);
 		// cout<< "coste inicial"<<coste<<endl;

 		 long horaDeseada=candidata.getHoraRecogidaSubida();
 		 bool seguir=true;
 		 int costeSubida;
 		 long tiempoLlegadaSubida;
 		 long tiempoViolacionSubida;
 		 int costeBajada;
 		 long tiempoLlegadaBajada;
 		 long tiempoViolacionBajada;

 		// cout<<"size solucion "<<solucion.size()<<endl;
 		 for(int k=0;k<solucion.size();k++){
 			//vemos si esta dentro de la ventana de tiempo, si lo esta, vemos si se puede insertar sin incumplir las restricciones en el resto de la ruta
 			 if(solucion[k].getHoraRecogida()<=horaDeseada+maxDelay){
 				//coste de ir desde la seleccionada hasta la parada de subida de la candidata
 				costeSubida=distmatrix.getDist(solucion[k].getIdParada(),candidata.getIdParadaSubida());
 				//cout<<"coste subida "<<costeSubida<<endl;
 				//hora de llegada a la subida de la peticion a encajar
 				tiempoLlegadaSubida=solucion[k].getHoraRecogida()+costeSubida;
 				//tiempo de violacion de subida (con respecto al deseado)
 				tiempoViolacionSubida=tiempoLlegadaSubida-candidata.getHoraRecogidaSubida();
 				//cout<<" tiempo violacion subida "<<tiempoViolacionSubida<<endl;
 				//si llegamos antes de la hora que ha pedido el cliente esperamos hasta la hora pedida
 				if(tiempoViolacionSubida<0){
 					//cout<<"llegamos antes esperamos "<<endl;
 					tiempoLlegadaSubida=candidata.getHoraRecogidaSubida();
 				}


 				int nPlazasLibres=solucion[k].getPlazasLibres()-candidata.getNumPax();
 				//cout<<"n plazas libres "<<nPlazasLibres<<endl;
 				if(tiempoViolacionSubida<=maxDelay && nPlazasLibres>=0){ //no viola el tiempo en cuanto a subida y la capacidad seguimos
 				//	cout<<"no viola el tiempo en cuanto a subida y capacidad"<<endl;
 					//inserto la subida y veo si el insertar esa subida no me jode el resto de peticiones con respecto a las restricciones
 					bool violaSubida;
 					std::vector <OrdenacionPeticiones> solucionVirtual=insertaEnSolucionSubida(candidata,solucion,tiempoLlegadaSubida,k,violaSubida,distmatrix,vehicleCapacity,maxDelay);

 				//	cout<<" size solucion virtual "<<solucionVirtual.size()<<endl;
 					//si el insertar la subida no me ha jodido la ruta pruebo a insertar la bajada
 					if(!violaSubida){
 					//	cout<<" subida insertada --> insertamos la bajada de la candidata "<<endl;
 						//calculo lo que se tardaria en ir directamente de la subida a la bajada de la candidata para calcular violacion luego
 						int tiempoDirecto=distmatrix.getDist(candidata.getIdParadaSubida(),candidata.getIdParadaBajada());

 						for(int l=k+1;l<solucionVirtual.size(); l++){
 							//insertamos ntra bajada detras de la peticion k+1 (siendo k+1 la subida correspondiente a ntra candidata)
 							//cout<<"vamos de ultima "<<solucionVirtual[l].getIdParada()<<" a "<<candidata.getIdParadaBajada()<<endl;
 							costeBajada=distmatrix.getDist(solucionVirtual[l].getIdParada(),candidata.getIdParadaBajada());
 							tiempoLlegadaBajada=solucionVirtual[l].getHoraRecogida()+costeBajada;
 							tiempoViolacionBajada=tiempoLlegadaBajada-candidata.getHoraRecogidaBajada();
 							//calculamos la duracion total del trayecto
 							int tiempoDuracionTrayecto=tiempoViolacionBajada+tiempoDirecto;

 							/////////////////////////////////////////////////////////////////
 							//filtramos por tipo de violacion bajada

 							//if(d.getTipoViolacionBajada()==1){	//violacion fija (60 min)

								if(tiempoDuracionTrayecto<=Vertex::maxRideTime(tiempoDirecto)){
									//cout<<"no viola el tiempo en cuanto a bajada"<<endl;
									bool violaBajada;
									std::vector <OrdenacionPeticiones> solucionVirtual2=insertaEnSolucionBajada(candidata,solucionVirtual,tiempoLlegadaBajada,l,violaBajada,distmatrix,vehicleCapacity,maxDelay);
									//si no viola las restricciones al encajar la bajada
									if(!violaBajada){
										//cout<<" bajada insertada" <<endl;
										int long long costeVirtual=evaluaRuta(solucionVirtual2,DARP_ALPHA,DARP_BETA);

										//cout<<" coste virtual "<<costeVirtual<<endl;
										int long long resta=costeVirtual-coste;
										if(resta<costeMinimo){/////////////////ojoooooooooooooo
											costeMinimo=resta;
											solucionVirtualMinima=solucionVirtual2;
											encajada=true;
											insertada=candidata;
											idVehiculo=j;
										}
									}//no viola resto paradas la bajada
								}//no viola candidata bajada
 							//}//violacion fija (60 min)


 							/*if(d.getTipoViolacionBajada()==2){ //violacion bajada variable

 								if(tiempoDuracionTrayecto<=d.getViolacionVariableBajada()*tiempoDirecto){
 									bool violaBajada;
 									std::vector <OrdenacionPeticiones> solucionVirtual2=insertaEnSolucionBajada(candidata,solucionVirtual,tiempoLlegadaBajada,l,d,violaBajada);
 									if(!violaBajada){
 										int long long costeVirtual=evaluaRuta(solucionVirtual2,d);

 										int id=candidata.getIdPeticionBD();
 										int long long resta=costeVirtual-coste;
 										if(resta<costeMinimo){/////////////////ojoooooooooooooo
 												costeMinimo=resta;
 												solucionVirtualMinima=solucionVirtual2;
 												encajada=true;
 												insertada=candidata;
 												idVehiculo=j;
 											}

									}//no viola resto paradas la bajada
 								}//no viola candidata bajada
 							}//violacion bajada variable*/

 							/////////////////////////////////////////////////////////////////

 						}//colocamos bajada*/
 					}//no viola resto paradas la subida
 				}//no viola candidata subida
 			 }else{ //la peticion ya no se encuentra dentro de la ventana de tiempo
 				 seguir=false;
 			 }
 		 }//for solucion
 	 }//for vehiculos



 	 if(encajada){
 		vehiculos[idVehiculo].setSolucion(solucionVirtualMinima);
 		std::vector<PeticionInsercion> a= vehiculos[idVehiculo].getAsignadas();
 		a.push_back(insertada);
 		vehiculos[idVehiculo].setAsignadas(a);
 		std::vector<OrdenacionPeticiones> sol=vehiculos[idVehiculo].getSolucion();
 		int long long coste=evaluaRuta(sol,DARP_ALPHA,DARP_BETA);
 		//int id=insertada.getIdPeticionBD();
 		//cout<<"coste final"<< coste <<endl;
 	 }

 	 return encajada;

  }


  /**************************************************************************************************
    *
    */
   void heuristicoInicializacion(std::vector<Vehiculo> &vehiculos,long vehicleCapacity,long maxDelay, double DARP_ALPHA, double DARP_BETA, VerticesList& verticeslist, long nVehicles, DistMatrix& distmatrix) {


  	 //1.- generamos tantos vehiculos como nos dice el fichero parametros e inicializamos la informacion
  	//std::vector<Vehiculo> vehiculos;
  	for(int i=0;i<nVehicles;i++){
  		//printf("vehiculo %d\n", i);
  	 	Vehiculo v;
  	 	v.setPlazas(vehicleCapacity);
  	 	v.setIdVehiculo(i);
  	 	vehiculos.push_back(v);
  	}

  	 std::vector<PeticionInsercion> candidatas;
  	 std::vector<PeticionInsercion> resto;

   	 //2.- Transformamos las peticiones a nuestra representacion del problema
   	 transformaPeticiones(candidatas,verticeslist);

   	 //3.- Buscamos la peticion mas temprana y la insertamos en cada vehiculo
   	int idVehiculo=0;

   	while(idVehiculo<nVehicles){
   		//Buscamos la peticion mas temprana
   		int indice=buscarPeticionMasTempranaAleatoria(candidatas);
   		//std::cout<<"indice "<<indice<<std::endl;
   		//La metemos como peticion inicial en la solucion y la eliminamos de candidatas

   		if(indice!=-1){
   			int horaSubida=candidatas[indice].getHoraRecogidaSubida();
   			int horaBajada=horaSubida+distmatrix.getDist(candidatas[indice].getIdParadaSubida(),candidatas[indice].getIdParadaBajada());;
   			insertaEnSolucionPrimera(vehiculos[idVehiculo],candidatas,horaSubida,horaBajada,indice);
   			std::vector<OrdenacionPeticiones> sol=vehiculos[idVehiculo].getSolucion();
   			//int long long coste=evaluaRuta(sol,d);
   			//cout<<"coste "<< coste << " insertada "<<id<< " idvehiculo "<<idVehiculo<<endl;
   		}
   		idVehiculo++;
   	}
   	while(candidatas.size()>0){
   		//cogemos la candidata mas temprana de todas y la quitamos de candidatos
   		int indiceCandidata=buscarPeticionMasTempranaAleatoria(candidatas);

   		PeticionInsercion candidata=candidatas[indiceCandidata];
   		//la eliminamos
   		candidatas=borrarCandidato(candidatas,indiceCandidata);
   		bool insertado=insertaPeticionNegativa(candidata,vehiculos,DARP_ALPHA,DARP_BETA,vehicleCapacity,maxDelay,distmatrix);

   		if(!insertado){
   			resto.push_back(candidata);
   		}
   	}


   	/**********************************************/
   	//hacemos el print de las rutas

   	/*ofstream ficheroSalida;
   	ficheroSalida.open("-BI.txt");
   	ofstream ficheroSalida1;
   	ficheroSalida1.open("-solucion.txt");

   	for (int i = 0; i < vehiculos.size(); i++) {
   				ficheroSalida1 << "vehiculo ";
   				ficheroSalida1 << vehiculos[i].getIdVehiculo();
   				ficheroSalida1 << "\n";
   				for (int j = 0; j < vehiculos[i].getSolucion().size(); j++) {
   					ficheroSalida << vehiculos[i].getIdVehiculo();
   					ficheroSalida << ",";
   					ficheroSalida << vehiculos[i].getSolucion()[j].getHoraRecogida();
   					ficheroSalida << ",";
   					ficheroSalida << vehiculos[i].getSolucion()[j].getIdParada();
   					ficheroSalida << ",";
   					ficheroSalida << vehiculos[i].getSolucion()[j].getIdPeticion();
   					ficheroSalida << ",";

   					if (vehiculos[i].getSolucion()[j].getIdPeticion()>0) {
   						ficheroSalida << "subida";
   					} else {
   						ficheroSalida << "bajada";
   					}
   					ficheroSalida << "\n";

   					if (vehiculos[i].getSolucion()[j].getIdPeticion()>0) {
   						int retraso = vehiculos[i].getSolucion()[j].getHoraRecogida()- vehiculos[i].getSolucion()[j].getHoraDeseada();
   						ficheroSalida1 << "retraso subida ";
   						ficheroSalida1 << retraso / 60;
   						ficheroSalida1 << "\n";
   					} else {
   						int retraso = vehiculos[i].getSolucion()[j].getHoraRecogida()- vehiculos[i].getSolucion()[j].getHoraDeseada();
   						ficheroSalida1 << "retraso bajada ";
   						ficheroSalida1 << retraso / 60;
   						ficheroSalida1 << "\n";
   					}
   				}
   			}

   	ficheroSalida.close();
   	ficheroSalida1.close();*/
   /***************************************************************************/

   	//cout<<" size resto "<<resto.size()<<endl;
   	//return resto;

   }

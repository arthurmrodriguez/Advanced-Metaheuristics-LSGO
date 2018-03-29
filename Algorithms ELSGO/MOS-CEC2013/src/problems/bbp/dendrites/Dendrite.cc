#include "Dendrite.h"
#include "attribute_names.h"
#include <boost/random.hpp>
#include <stdexcept>

static boost::mt19937 randgen (static_cast<unsigned int>(std::time(0)));
static boost::normal_distribution<double> NDInitialLength (13.058, 15.341);
static boost::normal_distribution<double> NDInitialTort(0.91, 0.118);
static boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > nDInitialLength(randgen, NDInitialLength);
static boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > nDInitialTort(randgen, NDInitialTort);

// Determnina el número de reintentos que se genera una dendrita no-real (num.ramas <=1, etc)
static unsigned NUM_REPETITIONS = 3;

Dendrite::Dendrite( NormalDistribution& terminationModel, LinearRegression& pathlengthModel, LinearRegression& eucDistModel )
	:_root(NULL),
	 _totalLength (0.0),
	 _maxEucDist  (0.0),
	 _numBranches (0),
	 _terminationModel ( terminationModel ),
	 _pathlengthModel ( pathlengthModel ),
	 _eucDistModel ( eucDistModel ),
   _numRepetitions (NUM_REPETITIONS)

{

	// informaciÃ³n comÃºn para la inicializaciÃ³n de los segmentos iniciales
/*	boost::mt19937 randgen (static_cast<unsigned int>(std::time(0)));
	boost::normal_distribution<double> NDInitialLength (13.058, 15.341);
	//boost::normal_distribution<double> NDInitialEucDist(11.149, 11.917);
	boost::normal_distribution<double> NDInitialTort(0.91, 0.118);
	boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > nDInitialLength(randgen, NDInitialLength);
	boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > nDInitialTort(randgen, NDInitialTort);
	//_nDInitialLength = new boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > (randgen, NDInitialLength);
	//_nDInitialTort = new boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > (randgen, NDInitialTort);
*/

}


Dendrite::~Dendrite(){
  if (_root)
    delete _root;
}

// Método para la generación de una dendrita. Si puede llamar repetidas veces, pero esto destruirá
// los datos generados anteriormente.
void Dendrite::grow(){
	_totalLength = 0.0;
	_maxEucDist = 0.0;
	_numBranches = 0;

	// datos de la raÃ­z
	boost::unordered_map<string, double> data;
  // Inicializa el recorrido inicial mientas que no sea menor o igual que 0.0.
  // Lo reintenta 3 veces si no va bien lo inicializa con la media de los reales
	unsigned count = 0;
	do{
		data[str_recorrido_segmento] = nDInitialLength();
	}while( data[str_recorrido_segmento] <= 0.0 and (count++ < 3) );
	if( data[str_recorrido_segmento] <= 0.0 ){
		data[str_recorrido_segmento] = 13.058;
	}

	data[str_orden] = 0;
	data[str_recorrido_soma] = data[str_recorrido_segmento];

  // Inicializa la tortuosidad mientas que no sea menor o igual que 0.0.
  // Lo reintenta 3 veces si no va bien lo inicializa con la media de los reales
	count = 0;
	do{
		data[str_tortuosidad] = nDInitialTort();
	}while( (data[str_tortuosidad] > 1 or data[str_tortuosidad] < 0) and (count++ < 3) );
	if( data[str_tortuosidad] > 1 or data[str_tortuosidad] < 0){
		data[str_tortuosidad] = 0.91;
	}

	data[str_distancia_euclidea] = data[str_tortuosidad] * data[str_recorrido_soma];
	data[str_distancia_euclidea_segmento] = data[str_distancia_euclidea];


  if (_root)
    delete _root;

	_root = new Instance(data, this);

	// crear dendrita completa
	_root->createNexts();

  if( _numBranches <= 1 and _numRepetitions > 0){
        _numRepetitions--;
        this->grow();
  }
  else // Por si se reutiliza con el mismo objeto, hay que inicializarlo con el número de repeticiones inicial
        _numRepetitions = NUM_REPETITIONS;

}

void Dendrite::storeInstances( char* path_file ) {
	if( _root ){
		_outputFile.open( path_file, std::ofstream::out | std::ofstream::trunc );
		_root->storeInstance();
		_outputFile.close();
	}
	else{
		std::cerr << "Error: Trying to store Dendrite before grow"<< std::endl;
		exit(-1);
	}
}

void Dendrite::storeMetrics( char* path_file ){

	std::ofstream outputMetrics ( path_file, std::ofstream::out | std::ofstream::app );
	outputMetrics << setprecision(9) << _numBranches << ", " << _maxEucDist << ", " << _totalLength << std::endl;
	outputMetrics.close();
}


boost::unordered_map<string, double> Dendrite::Metrics( void ){
	// Devuelve un map con los valores de cada una de las mÃ©tricas
	boost::unordered_map<string,double> metrics;

	if( _numBranches > 20 ){
		_maxEucDist = 0.0;
		_totalLength = 0.0;
		_numBranches = 0;
	}

	metrics[str_numero_branches] = (double)_numBranches;
	metrics[str_longitud_total] = _totalLength;
	metrics[str_max_distancia_euclidea] = _maxEucDist;

	return metrics;
}

void Dendrite::printBranch( std::string mesg ){
	//std::cout << mesg;
	_outputFile << mesg;
}

double Dendrite::predict(boost::unordered_map<string, double>& data, string param){
	if( param.compare(str_termina) == 0){
		return _terminationModel.Predict( data );
	}
	else if( param.compare(str_inc_recorrido) == 0 ){
		return _pathlengthModel.Predict( data );
	}
	else if( param.compare(str_inc_dist_euc) == 0){
		return _eucDistModel.Predict( data );
	}
  else
    throw runtime_error("[Dendrite] Error: wrong submodel in 'predict' method.");
}


void Dendrite::maxEucDist(double dist){
	if( _maxEucDist < dist)
		_maxEucDist = dist;
}

int Dendrite::incrementNumberOfBranches(){
	return ++_numBranches;
}

void Dendrite::incrementLength(double l){
	_totalLength += l;
}

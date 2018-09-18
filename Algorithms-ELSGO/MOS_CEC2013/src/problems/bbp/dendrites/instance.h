#ifndef _INSTANCE_
#define _INSTANCE_

//#include <map>
#include <vector>
#include <string>
#include <boost/unordered_map.hpp>

using namespace std;

#define NUM_HIJOS 2

class Dendrite;

class Instance{
	
public:
	Instance(){};
	Instance( boost::unordered_map<string, double>& data, Dendrite* d);
	~Instance();
	void createNexts();
	void storeInstance();
	void printClassValuesForFather();
	
private:
	vector<Instance*> _nexts; // Las dos instancias siguientes
	boost::unordered_map<string, double> _data; // conjunto de datos que define una rama
	Dendrite* _dendrite; // Acceso a la dendrita para tener acceso a los modelos y a las salida
};

#endif //_INSTANCE_

#include <GAEDAConfig.h>
#include <genomes/GA1DArrayGenome.h>
#include <MOSGenomeFactory.h>
#include <GARealOps.h>
#include <libconfig.h++>

#include <algorithm>
#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <time.h>

/*
 * This problem implements an optimized version of the inverse problem for detecting the sources
 */

using namespace boost::numeric::ublas;
using namespace std;

static string          transfmat_data_file  = "";
static string          optsources_data_file = "";
static string          signal_data_file     = "";
static matrix<long double>* transfmat            = 0;
static matrix<long double>* optsourcesvector     = 0;
static matrix<long double>* signalvector         = 0;
static matrix<long double>* solvector            = 0; // Used for storing the values when computing the objective value
static double          minvalue;
static double          maxvalue;

bool parse_config (char* fich) {

   libconfig::Config cfg;
   cfg.setAutoConvert (true);

   try {
      cfg.readFile (fich);
   }
   catch (libconfig::ParseException& e) {
     std::cerr << "Error: parsing the configuration file in line (" << e.getLine () << "): " << e.getError() << ". File name: " << fich << std::endl;
      exit (-1);
   }
   catch (libconfig::FileIOException& e) {
      std::cerr << "Error: file not found or could not be open." << std::endl;
      exit (-1);
   }

   // First, we read mandatory parameters
   try {
      transfmat_data_file  = (const char*) cfg.lookup("transfmat_data_file");
      signal_data_file     = (const char*) cfg.lookup("signal_data_file");
      minvalue             = (double)      cfg.lookup("minvalue");
      maxvalue             = (double)      cfg.lookup("maxvalue");
   }
   catch (libconfig::SettingNotFoundException& e) {
      std::cerr << "[DARP] Error: one or more mandatory parameter(s) were not found." << std::endl;
      exit (-1);
   }

   // Optional Stats
   if (cfg.exists("optsources_data_file")) optsources_data_file = (const char*) cfg.lookup("optsources_data_file");

   return true;
}

//extern "C" void populationInit (GAPopulation& pop, double perc) {
//  throw runtime_error("[populationInit] Error: initialization method not defined.");
//}

void fillSolVectorWithValues(GA1DArrayAlleleGenome<double>& genome) {
  assert(genome.size() == solvector->size1() && solvector->size2() == 1);

  for (int i=0; i<genome.size(); i++) (*solvector)(i,0) = genome.gene(i);
}

#include <sys/time.h>

extern "C" double objective(GAGenome& g) {
  GA1DArrayAlleleGenome<double>& gen = dynamic_cast< GA1DArrayAlleleGenome<double>& > (g);


  fillSolVectorWithValues(gen);
  ///*LOG*/ timeval startTime;
  ///*LOG*/ gettimeofday(&startTime, NULL);

  matrix<long double> resmatrix = prod(*transfmat, *solvector);
  ///*LOG*/ timeval endTime;
  ///*LOG*/ long seconds, useconds;
  ///*LOG*/ double duration;
  ///*LOG*/ gettimeofday(&endTime, NULL);
  ///*LOG*/ seconds  = endTime.tv_sec  - startTime.tv_sec;
  ///*LOG*/ useconds = endTime.tv_usec - startTime.tv_usec;

  ///*LOG*/ duration = seconds + useconds/1000000.0;
  ///*LOG*/ cout << "Evaluation time=" << duration << endl;

  assert(resmatrix.size2() == 1);

  resmatrix = (*signalvector) - resmatrix;

  long double res = 0;
  for (int i=0; i<resmatrix.size1(); i++) {
    res += resmatrix(i,0) * resmatrix(i,0);
  }


  return res;
}

matrix<long double>* readMatrixData(string filename) {
  assert(!filename.empty());

  ifstream fin(filename.c_str());

  // First we read the data into a temporary variable
  string line;
  std::vector< std::vector<long double> > tmp_matrix;
  while ( true ) {
    getline(fin,line);
    if (!fin.good() ) break;

    std::vector<long double> tmp_values;
    stringstream l(line);
    string s;

    while(getline(l,s,',')) { // The fields need to be separated by commas
      stringstream p(s);
      long double t;
      p >> t;
      tmp_values.push_back(t);
    }

    tmp_matrix.push_back(tmp_values);
  }

  // Then we create the matrix with the appropriate size and we fill the matrix variable with the data
  matrix<long double>* resmat = new matrix<long double>( tmp_matrix.size(), tmp_matrix[0].size() );
  for (int i=0; i<resmat->size1(); i++) {
    for (int j=0; j<resmat->size2(); j++) {
      (*resmat)(i,j) = tmp_matrix[i][j];
    }
  }

  return resmat;
}

matrix<long double>* readVectorData(string filename) {
  assert(!filename.empty());

  ifstream fin(filename.c_str());

  // First we read the data into a temporary variable
  string line;
  std::vector< long double > tmp_vector;
  while ( true ) {
    getline(fin,line);
    if (!fin.good() ) break;

    stringstream l(line);
    long double t;
    l >> t;
    tmp_vector.push_back(t);
  }

  // Then we create the vector with the appropriate size and we fill the matrix variable with the data

  matrix<long double>* sourcesvector = new matrix<long double>( tmp_vector.size(), 1 );

  // Then we fill the vector variable with the data
  for (int i=0; i<sourcesvector->size1(); i++) {
      (*sourcesvector)(i,0) = tmp_vector[i];
  }

  return sourcesvector;
}

void readData() {

  transfmat = readMatrixData(transfmat_data_file.c_str());

  // Since we now know the dimensions of the vector used in the objective value, we create it
  solvector = new matrix<long double>( transfmat->size2() , 1 );

  if (!optsources_data_file.empty() ) {
    optsourcesvector = readVectorData(optsources_data_file.c_str());
    assert(transfmat->size2() == optsourcesvector->size1());
    assert(transfmat->size2() != optsourcesvector->size2());
  }

  signalvector = readVectorData(signal_data_file.c_str());
  assert(transfmat->size1() == signalvector->size1());  assert(transfmat->size2() != signalvector->size2());
}

extern "C" void individualInit (GAGenome& g) {
   return RealUniformInitializer (g);
}

extern "C" GAGenome* defineProblem () {
  // Parse the config file for the problem
  if (!parse_config (GAEDAConfig::handle()->getProblemData ())) {
    throw runtime_error ("Error: A valid configuration file must be provided for the TSP problem.");
  }

  readData();

  GAEDAConfig* cfg = GAEDAConfig::handle();
  if (cfg->getProblemSize() != transfmat->size2() ) {
    cerr << "Problem size (" << cfg->getProblemSize() << ") != transfmat size2=" <<  transfmat->size2() << endl;
    exit(-1);
  }

  GAAlleleSet<double> alleles (minvalue, maxvalue);
  GA1DArrayAlleleGenome<double>* genome = new GA1DArrayAlleleGenome<double> (cfg->getProblemSize(), alleles, objective);

  genome->initializer (RealUniformInitializer);
  genome->comparator  (RealEuclideanComparator);

  // Specific stuff for MOS
  MOSGenomeFactory::handle()->registerGenome (GAID::RealEncoding, genome);


  cout << endl;
  cout << "Neural Sources specific parameter values: " << endl;
  cout << "Transfmatrix data file    : " << transfmat_data_file << endl;
  cout << "Optimum sources data file : " << optsources_data_file << endl;
  cout << "Measured Signal data file : " << signal_data_file << endl;
  cout << "Min allowed value         : " << minvalue            << endl;
  cout << "Max allowed value         : " << maxvalue           << endl;

  cout << endl;

  /*LOG*/ if (optsourcesvector != 0) {
  /*LOG*/  assert(optsourcesvector->size1() == genome->size() && optsourcesvector->size2() == 1);
  /*LOG*/  for (int i=0; i<genome->size(); i++) {
  /*LOG*/    genome->gene(i, (*optsourcesvector)(i,0) );
  /*LOG*/  }
  /*LOG*/  cout << "El valor devuelto por la funcion objetivo para el optimo es " << objective(*genome) << endl;
  /*LOG*/ }

  return genome;
}


extern "C" const char *describeProblem (void) {
   return "Neural Sources Localization Problem";
}


extern "C" GAGenome::OptCriterion optCriterion(){
  return GAGenome::MINIMIZATION;
}

extern "C" bool postprocess (GAPopulation* pop, int rank, int repetition) {
  return true;
}

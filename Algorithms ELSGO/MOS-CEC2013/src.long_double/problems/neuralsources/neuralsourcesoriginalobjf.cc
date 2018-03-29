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
 * This problem implements the original version of the inverse problem for detecting the sources.
 * Therefore, it implements the procedure described in the document presentacion-detec-fuentes-Cesvima.pdf
 */



using namespace boost::numeric::ublas;
using namespace std;

// Constants
const static int NUMSOURCES  = 2988;
const static int NUMCOORDS   = 3;
const static int NUMCHANNELS = 306;


static string                               transfmat_data_file  = "";
static string                               optsources_data_file = "";
static string                               signal_data_file     = "";
static std::vector< matrix<long double>* >* leadfield_matrices   = 0;
static std::vector< matrix<long double>* >* sources_sol          = 0; // Used for storing the values when computing the objective value
matrix<long double>*                        opt_sol              = 0;


static matrix<long double>* signalvector         = 0;
static long double               minvalue;
static long double               maxvalue;

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
      minvalue             = (long double)      cfg.lookup("minvalue");
      maxvalue             = (long double)      cfg.lookup("maxvalue");
   }
   catch (libconfig::SettingNotFoundException& e) {
      std::cerr << "[DARP] Error: one or more mandatory parameter(s) were not found." << std::endl;
      exit (-1);
   }

   // Optional Stats
   if (cfg.exists("optsources_data_file")) optsources_data_file = (const char*) cfg.lookup("optsources_data_file");

   return true;
}

//extern "C" void populationInit (GAPopulation& pop, long double perc) {
//  throw runtime_error("[populationInit] Error: initialization method not defined.");
//}

void fillSolVectorWithValues(GA1DArrayAlleleGenome<long double>& genome) {
  assert(genome.size() == sources_sol->size() * ( (*sources_sol)[0]->size1() ) );
  assert(sources_sol->size() == NUMSOURCES);

  for (int i=0; i<sources_sol->size(); i++) {          // it should be == NUMSOURCES

    matrix<long double>* sources_sol_i = (*sources_sol)[i];
    assert(sources_sol_i->size1() == NUMCOORDS);
    assert(sources_sol_i->size2() == 1);

    for (int j=0; j<sources_sol_i->size1(); j++) {     // it should be == NUMCOORDS
      int pos = i*sources_sol_i->size1()+j;
      assert(pos >= 0 && pos < genome.size() );
      (*sources_sol_i)(j,0) = genome.gene(pos);
    }
  }

}

#include <sys/time.h>

extern "C" long double objective(GAGenome& g) {
  GA1DArrayAlleleGenome<long double>& gen = dynamic_cast< GA1DArrayAlleleGenome<long double>& > (g);

  fillSolVectorWithValues(gen);
  ///*LOG*/ timeval startTime;
  ///*LOG*/ gettimeofday(&startTime, NULL);

  matrix<long double> tmpmatrix(NUMCHANNELS,1,0);

  for (int i=0; i<NUMSOURCES; i++) {
    matrix<long double>* leadfield_matrix_i = (*leadfield_matrices)[i];
    matrix<long double>* sources_sol_i      = (*sources_sol)[i];

    tmpmatrix = tmpmatrix + prod(*leadfield_matrix_i,*sources_sol_i);
  }

  matrix<long double> resmatrix = (*signalvector) - tmpmatrix;
  ///*LOG*/ timeval endTime;
  ///*LOG*/ long seconds, useconds;
  ///*LOG*/ long double duration;
  ///*LOG*/ gettimeofday(&endTime, NULL);
  ///*LOG*/ seconds  = endTime.tv_sec  - startTime.tv_sec;
  ///*LOG*/ useconds = endTime.tv_usec - startTime.tv_usec;

  ///*LOG*/ duration = seconds + useconds/1000000.0;
  ///*LOG*/ cout << "Evaluation time=" << duration << endl;

  assert(resmatrix.size2() == 1);

  long double res = 0;
  for (int i=0; i<resmatrix.size1(); i++) {
    res += resmatrix(i,0) * resmatrix(i,0);
  }


  return res;
}

std::vector< std::vector<long double>* >* readMatrixFromFile(string filename) {
  assert(!filename.empty());

  ifstream fin(filename.c_str());

  // First we read the data into a temporary variable
  string line;
  std::vector< std::vector<long double>* >* tmp_matrix = new std::vector< std::vector<long double>* >();
  while ( true ) {
    getline(fin,line);
    if (!fin.good() ) break;

    std::vector<long double>* tmp_values = new std::vector<long double>();
    stringstream l(line);
    string s;

    while(getline(l,s,',')) { // The fields need to be separated by commas
      stringstream p(s);
      long double t;
      p >> t;
      tmp_values->push_back(t);
    }

    tmp_matrix->push_back(tmp_values);
  }

  return tmp_matrix;
}

void deleteTmpMatrix(std::vector< std::vector<long double>* >* matrix) {
  for (int i=0; i<matrix->size(); i++) {
      delete (*matrix)[i];
  }
  delete matrix;
}

/*
 * Reads the set of leadfield matrices. It is assumed that the file contains the data. The structure
 * of this data is the following NUMCHANNELS rows x ( NUMCOORDS * NUMSOURCES) columns. Each
 * NUMCHANNELS x (NUMCOORDS) matrix corresponds to the matrix associated to a source.
 */
std::vector< matrix<long double>* >* readLeadfieldMatrix(string filename) {
  std::vector< std::vector<long double>* >* tmp_matrix = readMatrixFromFile(filename);

  assert(tmp_matrix->size() == NUMCHANNELS);
  for (int i=0; i<tmp_matrix->size(); i++) assert((*tmp_matrix)[i]->size() == NUMCOORDS*NUMSOURCES);

  // Then we create the matrix with the appropriate size and we fill the matrix variable with the data
  std::vector< matrix<long double>* >* leadfield_matrices = new std::vector< matrix<long double>* > (NUMSOURCES);

  for (int i=0; i<NUMSOURCES; i++) {

    (*leadfield_matrices)[i] = new matrix<long double>(NUMCHANNELS,NUMCOORDS);
    matrix<long double>* matrix_i = (*leadfield_matrices)[i];

    for (int j=0; j<NUMCHANNELS; j++) {

      for (int k=0; k<NUMCOORDS; k++) {
        (*matrix_i) (j,k) = ( * (*tmp_matrix)[j] )[i*NUMCOORDS+k];
      }

    }
  }

  deleteTmpMatrix(tmp_matrix);

  return leadfield_matrices;
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

void readAndInitData() {

  leadfield_matrices = readLeadfieldMatrix(transfmat_data_file.c_str());

  if (!optsources_data_file.empty() ) {
    opt_sol = readVectorData(optsources_data_file.c_str());
    assert( opt_sol->size1() == NUMSOURCES*NUMCOORDS);
  }

  signalvector = readVectorData(signal_data_file.c_str());
  assert(signalvector->size1() == NUMCHANNELS);
  assert(signalvector->size2() == 1);

  // This variable is going to hold the values of the solution. We reserve the space that is going to use
  sources_sol = new std::vector< matrix<long double>* >(NUMSOURCES);
  for (int i=0; i<sources_sol->size(); i++) (*sources_sol)[i] = new matrix<long double>(NUMCOORDS,1);
}

extern "C" void individualInit (GAGenome& g) {
   return RealUniformInitializer (g);
}

extern "C" GAGenome* defineProblem () {
  // Parse the config file for the problem
  if (!parse_config (GAEDAConfig::handle()->getProblemData ())) {
    throw runtime_error ("Error: A valid configuration file must be provided for the TSP problem.");
  }

  readAndInitData();

  GAEDAConfig* cfg = GAEDAConfig::handle();
  if (cfg->getProblemSize() != leadfield_matrices->size() * (*leadfield_matrices)[0]->size2() ) {
    cerr << "Problem size (" << cfg->getProblemSize() << ") != transfmat num sources=" <<  leadfield_matrices->size() << " size 2 of matrix 0=" << (*leadfield_matrices)[0]->size2() << endl;
    exit(-1);
  }

  GAAlleleSet<long double> alleles (minvalue, maxvalue);
  GA1DArrayAlleleGenome<long double>* genome = new GA1DArrayAlleleGenome<long double> (cfg->getProblemSize(), alleles, objective);

  genome->initializer (RealUniformInitializer);
  genome->comparator  (RealEuclideanComparator);

  // Specific stuff for MOS
  MOSGenomeFactory::handle()->registerGenome (GAID::RealEncoding, genome);


  cout << endl;
  cout << "Neural Sources specific parameter values: " << endl;
  cout << "Transfmatrix data fil e   : " << transfmat_data_file << endl;
  cout << "Optimum sources data file : " << optsources_data_file << endl;
  cout << "Measured Signal data file : " << signal_data_file << endl;
  cout << "Min allowed value         : " << minvalue            << endl;
  cout << "Max allowed value         : " << maxvalue           << endl;

  cout << endl;

  /*LOG*/ if (opt_sol) {
  /*LOG*/  assert(opt_sol->size1() == genome->size() );
  /*LOG*/  for (int i=0; i<genome->size(); i++) {
  /*LOG*/    genome->gene(i, (*opt_sol)(i,0) );
  /*LOG*/  }
  /*LOG*/  cout << "El valor devuelto por la funcion objetivo para el optimo es " << objective(*genome) << endl;
  /*LOG*/ }

  return genome;
}


extern "C" const char *describeProblem (void) {
   return "Neural Sources Localization Problem (Exact version)";
}


extern "C" GAGenome::OptCriterion optCriterion(){
  return GAGenome::MINIMIZATION;
}

extern "C" bool postprocess (GAPopulation* pop, int rank, int repetition) {
  return true;
}

// old code
//std::vector< matrix<long double>* >* readSourcesData(string filename) {
//  assert(!filename.empty());
//
//  std::vector< std::vector<long double>* >* tmp_matrix = readMatrixFromFile(filename);
//
//  assert(tmp_matrix->size() == NUMSOURCES);
//  for (int i=0; i<tmp_matrix->size(); i++) assert( (*tmp_matrix)[i]->size() == NUMCOORDS );
//
//  // Then we create the vector with the appropriate size and we fill the matrix variable with the data
//  std::vector< matrix<long double>* >* sources = new std::vector< matrix<long double>* >(NUMSOURCES);
//
//  for (int i=0; i<NUMSOURCES; i++) {
//
//    (*sources)[i]                 = new matrix<long double>(NUMCOORDS,1);
//    matrix<long double>& source_i = (*sources)[i];
//
//    for (int j=0; j<NUMCOORDS; j++) {
//      source_i (j,0) = ( * (*tmp_matrix)[i] )[j];
//    }
//  }
//
//  deleteTmpMatrix(tmp_matrix);
//
//  return sources;
//}

#include <algorithm>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include <boost/unordered_map.hpp>
#include <studentttests.h>
#include <libconfig.h++>

#include <GARealOps.h>
#include <GAPopulation.h>
#include <MOSGenomeFactory.h>
#include <genomes/GA1DArrayGenome.h>
#include <GAEDAConfig.h>

#include "attribute_names.h"
#include "Dendrite.h"
#include "kolmogorov.h"
#include "weka.h"

typedef struct ModelsT {
  LinearRegression*   euclidean;
  LinearRegression*   pathway;
  NormalDistribution* avgDev;
} models;

typedef struct MetricsT {
  std::vector<long double> numBranches;
  std::vector<long double> maxEucDistance;
  std::vector<long double> totalLength;
} metrics;

typedef std::pair<models, metrics> bsDataPerIter;

unsigned bsIters      = 2;
unsigned valIters     = 0;
int      currentFold  = 0;
std::string directory = "";

std::vector<bsDataPerIter> bsData;
bsDataPerIter bsResustitutionData;

// Maps to compute the average values of each component of the models and the standard deviations
boost::unordered_map<std::string, long double> incRec_avg, incEuc_avg, incRec_dev, incEuc_dev;

long double avgDev_mean_avg, avgDev_stddev_avg, avgDev_mean_dev, avgDev_stddev_dev;

bool pprocess = false;

bool parse_config (char* fich) {

  libconfig::Config cfg;
  cfg.setAutoConvert (true);

  try {
    cfg.readFile (fich);
  }
  catch (libconfig::ParseException e) {
    std::cerr << "[Dendrites] Error: parsing the configuration file in line (" << e.getLine () << "): " << e.getError() << "." << std::endl;
    exit (-1);
  }
  catch (libconfig::FileIOException e) {
    std::cerr << "[Dendrites] Error: file not found or could not be open." << std::endl;
    exit (-1);
  }

  // First, we read mandatory parameters
  try {
    directory   = (const char*) cfg.lookup ("directory"  );
    //bsIters     =               cfg.lookup ("bsIters"    );
    valIters    =               cfg.lookup ("valIters"   );
    currentFold =               cfg.lookup ("currentFold");

    if (bsIters <= 0) {
      std::cerr << "[Dendrites] Error: bsIters (number of iterations for the bootstrap) must be greater than zero." << std::endl;
      exit (-1);
    }

    if (currentFold < 0 || currentFold > 9) {
      std::cerr << "[Dendrites] Error: currentFold (current iteration of the 10-fold cross validation) must be within the interval [0,9]." << std::endl;
      exit (-1);
    }
  }
  catch (libconfig::SettingNotFoundException e) {
    std::cerr << "[Dendrites] Error: one or more mandatory parameter(s) were not found." << std::endl;
    exit (-1);
  }

  // Then, we read optional vars (the default value is not modified if the
  // parameter is not provided in the config file)

  return true;

}

void csvline_populate(std::vector<std::string> &record, const std::string& line, char delimiter) {
  int linepos=0;
  int inquotes=false;
  char c;
  int i;
  int linemax=line.length();
  std::string curstring;

  record.clear();

  while(line[linepos] != 0 && linepos < linemax) {

    c = line[linepos];

    if (!inquotes && curstring.length() == 0 && c == '"') { //beginquotechar
      inquotes=true;
    }
    else if (inquotes && c == '"') { //quotechar
      if ((linepos+1 < linemax) && (line[linepos+1] == '"')) { //encountered 2 long double quotes in a row (resolves to 1 long double quote)
        curstring.push_back(c);
        linepos++;
      }
      else { //endquotechar
        inquotes=false;
      }
    }
    else if (!inquotes && c == delimiter) { //end of field
      record.push_back(curstring);
      curstring="";
    }
    else if (!inquotes && (c == '\r' || c == '\n')) {
      record.push_back( curstring );
      return;
    }
    else {
      curstring.push_back(c);
    }

    linepos++;

  }

  record.push_back( curstring );

  return;
}

bool parse_data (void) {

  /* *********************************************** */
  /* Initialization of variables and data structures */
  /* *********************************************** */

  // Allocate memory for the data of the bootstrap
  // We allocate soace for one more iteration for the resustitution information
  bsData.reserve(bsIters + 1);

  // Vectors with the attributes for both regression models
  std::vector<std::string> attrib_incRec(NAMES, NAMES + LAST_ATT_I);

  // This model has one more attribute
  std::vector<std::string> attrib_incEuc(attrib_incRec);
  attrib_incEuc.push_back(str_inc_recorrido);

  // Initialization of the global variables that will store the average
  // and standard deviation of both Linear Regression models
  for (unsigned i = 0; i < attrib_incRec.size(); i++) {
    incRec_avg[attrib_incRec[i]] = 0.0;
    incRec_dev[attrib_incRec[i]] = 0.0;
  }

  incRec_avg["intercept_factor"] = 0.0;
  incRec_dev["intercept_factor"] = 0.0;

  incRec_avg["error"] = 0.0;
  incRec_dev["error"] = 0.0;

  for (unsigned i = 0; i < attrib_incEuc.size(); i++) {
    incEuc_avg[attrib_incEuc[i]] = 0.0;
    incEuc_dev[attrib_incEuc[i]] = 0.0;
  }

  incEuc_avg["intercept_factor"] = 0.0;
  incEuc_dev["intercept_factor"] = 0.0;

  incEuc_avg["error"] = 0.0;
  incEuc_dev["error"] = 0.0;

  // Initialization of the global variables that will store the average
  // and standard deviation of the Normal Distribution model
  avgDev_mean_avg = avgDev_stddev_avg = avgDev_mean_dev = avgDev_stddev_dev = 0.0;

  // Aux variables to store regressional and normal models
  LinearRegression* regModel;
  NormalDistribution* normalDist;

  /* ********************************************* */
  /* Main loop to parse the data for the bootstrap */
  /* ********************************************* */

  for (unsigned i = 0; i < bsIters; i++) {

    bsDataPerIter bsDataIter;

    /* *************************** */
    /* Load models (training data) */
    /* *************************** */

    stringstream filenamebuf;
    filenamebuf << "/bootStrap/bootStrap" << currentFold << "/Train/models_b_train" << i << ".csv/";

    // Load the path model and store it in current's BS iteration structure
    regModel = new LinearRegression((char*) (directory + filenamebuf.str() + "iRecorrido.model").c_str(), attrib_incRec);
    bsDataIter.first.pathway = regModel;

    // Add the values of this model to the previous ones into the incRec_avg vector tu compute the mean afterwards
    regModel->addModel(incRec_avg);

    // Load the euclidean model and store it in current's BS iteration structure
    regModel = new LinearRegression((char*) (directory + filenamebuf.str() + "iEuclideo.model").c_str(), attrib_incEuc);
    bsDataIter.first.euclidean = regModel;

    // Add the values of this model to the previous ones into the incEuc_avg vector tu compute the mean afterwards
    regModel->addModel(incEuc_avg);

    // Load the termination model and store it in current's BS iteration structure
    normalDist = new NormalDistribution((char*) (directory + filenamebuf.str() + "iMediaDesv.model").c_str());
    bsDataIter.first.avgDev = normalDist;

    // Add the values of this model to the previous ones into the avgDev_mean_avg and avgDev_stddev_avg variables tu compute the mean afterwards
    avgDev_mean_avg   += normalDist->getMean();
    avgDev_stddev_avg += normalDist->getStdDev();


    /* ****************************** */
    /* Load metrics (validation data) */
    /* ****************************** */

    stringstream filenamebuf2;
    filenamebuf2 << "/bootStrap/bootStrap" << currentFold << "/Test/metrics_b_test" << i << ".csv";

    // Open the file containing the metrics
    std::ifstream metrics_file((directory + filenamebuf2.str()).c_str());

    // If file does no exist, echo an error message and exit
    if (metrics_file.fail()) {
      std::cerr << "[Dendrites] Metrics file '" << directory << filenamebuf2.str() << "' not found" << std::endl;
      exit(0);
    }

    // Read the file and parse the metrics
    std::string line;

    while (getline(metrics_file, line) && metrics_file.good()) {

      // Parse current line and retrieve metrics
      std::vector<std::string> res;
      csvline_populate(res, line, ',');

      // Check for the correct number of metrics
      assert(res.size() == 3);

      // Push metrics to each of the vectors associated with its respective model
      bsDataIter.second.numBranches.push_back   (strtod(res[0].c_str(), NULL));
      bsDataIter.second.maxEucDistance.push_back(strtod(res[1].c_str(), NULL));
      bsDataIter.second.totalLength.push_back   (strtod(res[2].c_str(), NULL));

    }

    // Sort the three vectors (vectors must be sorted for the KS Test)
    std::sort(bsDataIter.second.numBranches.begin(),    bsDataIter.second.numBranches.end()   );
    std::sort(bsDataIter.second.maxEucDistance.begin(), bsDataIter.second.maxEucDistance.end());
    std::sort(bsDataIter.second.totalLength.begin(),    bsDataIter.second.totalLength.end()   );

    // Ensure that the same number of metrics have been parsed for each model
    assert(bsDataIter.second.numBranches.size()    == bsDataIter.second.maxEucDistance.size() &&
           bsDataIter.second.numBranches.size()    == bsDataIter.second.totalLength.size()    &&
           bsDataIter.second.maxEucDistance.size() == bsDataIter.second.totalLength.size()      );

    // Push back the data of this iteration into the bsData structure
    bsData.push_back(bsDataIter);

  }

  /* ************************************************************************************** */
  /* Computation of average values and standard deviations for each attribute of each model */
  /* ************************************************************************************** */

  boost::unordered_map<std::string, long double>::iterator it;

  // Compute the avg and the stddev for each attribute of the path model
  for (it = incRec_avg.begin(); it != incRec_avg.end(); it++) {
    // Divide the overall value for this attribute (we only added them before) by the number of iterations
    it->second /= bsIters;

    // Compute the cummulated differences for the stddev
    long double diffs_acum = 0.0;
    for (unsigned i = 0; i < bsIters; i++)
      diffs_acum += pow(bsData[i].first.pathway->getModel()[it->first] - it->second, 2);

    // Compute the stddev for that attribute
    incRec_dev[it->first] = sqrt(diffs_acum / bsIters);
  }


  // Compute the avg and the stddev for each attribute of the euclidean model
  for (it = incEuc_avg.begin(); it != incEuc_avg.end(); it++) {
    // Divide the overall value for this attribute (we only added them before) by the number of iterations
    it->second /= bsIters;

    // Compute the cummulated differences for the stddev
    long double diffs_acum = 0.0;
    for (unsigned i = 0; i < bsIters; i++)
      diffs_acum += pow(bsData[i].first.euclidean->getModel()[it->first] - it->second, 2);

    // Compute the stddev for that attribute
    incEuc_dev[it->first] = sqrt(diffs_acum / bsIters);
  }


  // Compute the avg and the stddev for both the mean and the stddev of the termination model
  avgDev_mean_avg /= bsIters;
  avgDev_stddev_avg /= bsIters;

  // Compute the cummulated differences for both attributes
  for (unsigned i = 0; i < bsIters; i++) {
    avgDev_mean_dev   += pow(bsData[i].first.avgDev->getMean()   - avgDev_mean_avg,   2);
    avgDev_stddev_dev += pow(bsData[i].first.avgDev->getStdDev() - avgDev_stddev_avg, 2);
  }

  // Compute the stddev for both attributes
  avgDev_mean_dev   = sqrt(avgDev_mean_dev   / bsIters);
  avgDev_stddev_dev = sqrt(avgDev_stddev_dev / bsIters);

  /* ******************************************* */
  /* Parsing of the resustitution data of the BS */
  /* ******************************************* */

  stringstream resustmodels, resustmetrics;
  bsDataPerIter bsDataResust;

  /* *************************** */
  /* Load models (training data) */
  /* *************************** */

  resustmodels << "/resustitution/models_train_fold" << currentFold << ".csv/";

  // Load the path model and store it in last's BS iteration structure
  regModel = new LinearRegression((char*) (directory + resustmodels.str() + "iRecorrido.model").c_str(), attrib_incRec);
  bsDataResust.first.pathway = regModel;

  // Load the euclidean model and store it in last's BS iteration structure
  regModel = new LinearRegression((char*) (directory + resustmodels.str() + "iEuclideo.model").c_str(), attrib_incEuc);
  bsDataResust.first.euclidean = regModel;

  // Load the termination model and store it in last's BS iteration structure
  normalDist = new NormalDistribution((char*) (directory + resustmodels.str() + "iMediaDesv.model").c_str());
  bsDataResust.first.avgDev = normalDist;


  /* ****************************** */
  /* Load metrics (validation data) */
  /* ****************************** */

  resustmetrics << "/resustitution/metrics_train_fold" << currentFold << ".csv";

  // Open the file containing the metrics
  std::ifstream resust_metrics_file((directory + resustmetrics.str()).c_str());

  // If file does no exist, echo an error message and exit
  if (resust_metrics_file.fail()) {
    std::cerr << "[Dendrites] Metrics file for resustitution '" << directory << resustmetrics.str() << "' not found" << std::endl;
    exit(0);
  }

  // Read the file and parse the metrics
  std::string line;

  while (getline(resust_metrics_file, line) && resust_metrics_file.good()) {

    // Parse current line and retrieve metrics
    std::vector<std::string> res;
    csvline_populate(res, line, ',');

    // Check for the correct number of metrics
    assert(res.size() == 3);

    // Push metrics to each of the vectors associated with its respective model
    bsDataResust.second.numBranches.push_back   (strtod(res[0].c_str(), NULL));
    bsDataResust.second.maxEucDistance.push_back(strtod(res[1].c_str(), NULL));
    bsDataResust.second.totalLength.push_back   (strtod(res[2].c_str(), NULL));

  }

  // Sort the three vectors (vectors must be sorted for the KS Test)
  std::sort(bsDataResust.second.numBranches.begin(),    bsDataResust.second.numBranches.end()   );
  std::sort(bsDataResust.second.maxEucDistance.begin(), bsDataResust.second.maxEucDistance.end());
  std::sort(bsDataResust.second.totalLength.begin(),    bsDataResust.second.totalLength.end()   );

  // Ensure that the same number of metrics have been parsed for each model
  assert(bsDataResust.second.numBranches.size()    == bsDataResust.second.maxEucDistance.size() &&
         bsDataResust.second.numBranches.size()    == bsDataResust.second.totalLength.size()    &&
         bsDataResust.second.maxEucDistance.size() == bsDataResust.second.totalLength.size()      );

  bsData.push_back(bsDataResust);

  return true;

}


extern "C" long double objective (GAGenome& g) {

  static unsigned indiv = 0;

  /* *********************************************** */
  /* Initialization of variables and data structures */
  /* *********************************************** */

  GA1DArrayAlleleGenome<long double>& genome = DYN_CAST (GA1DArrayAlleleGenome<long double> &, g);

  // Aux vars to store the copies of the actual models that will be perturbated
  LinearRegression*   eucModel;
  LinearRegression*   pathModel;
  NormalDistribution* avgDevModel;

  // Aux vars to store partial p-values (per model and for resustitution/no-resustitution)
  long double pvalue_euc = 0, pvalue_path = 0, pvalue_avg = 0, pvalue_euc_res = 0, pvalue_path_res = 0, pvalue_avg_res = 0;

  // Aux var to store the resustitution p-value (this phase of the BS is executed
  // in the loop, with the info in the last position of the vector)
  long double pvalue_res = 0.0;

  // Aux var to store the no-resustitution p-value
  long double pvalue_nores = 0.0;


  /* *********************************************************************** */
  /* Main loop to conduct the bootstrap.                                     */
  /* Recall: the resustitution is done in the last iteration of this loop.   */
  /*         For this reason there is an extra iteration in the loop header. */
  /* *********************************************************************** */

  for (unsigned i = 0; i < bsIters + 1; i++) {

    /* ****************************************************** */
    /* Perturbation of the models of this iteration of the BS */
    /* ****************************************************** */

    // Create copies of the models of the current iteration for the perturbation take place
    eucModel    = new LinearRegression  (*bsData[i].first.euclidean);
    pathModel   = new LinearRegression  (*bsData[i].first.pathway  );
    avgDevModel = new NormalDistribution(*bsData[i].first.avgDev   );

    // Perturb the three copies of the models with the values obtained by the EA
    boost::unordered_map<std::string, long double> perturbation;
    boost::unordered_map<std::string, long double>::iterator it;
    unsigned pos = 0;


    // Construct the perturbation vector
    for (; pos < LAST_ATT_I; pos++) {
      perturbation[NAMES[pos]] = genome.gene(pos) * 3 * incRec_dev[NAMES[pos]];
      //std::cout << NAMES[pos] << " : " << perturbation[NAMES[pos]] << std::endl;
    }

    perturbation["intercept_factor"] = genome.gene(pos) * 3 * incRec_dev["intercept_factor"];
    //std::cout << "intercept_factor: " << perturbation["intercept_factor"] << std::endl;
    pos++;

    perturbation["error"] = genome.gene(pos) * 3 * incRec_dev["error"];
    //std::cout << "error: " << perturbation["error"] << std::endl;
    pos++;

    //std::cout << std::endl;

    // Perturb the path model
    pathModel->perturbModel(perturbation);

    //std::cout << std::endl;


    // Construct the perturbation vector
    for (unsigned pos2 = 0; pos2 < LAST_ATT_I; pos2++, pos++) {
      perturbation[NAMES[pos2]] = genome.gene(pos) * 3 * incEuc_dev[NAMES[pos2]];
      //std::cout << NAMES[pos2] << " : " << incEuc_dev[NAMES[pos2]] << std::endl;
    }

    perturbation[str_inc_recorrido] = genome.gene(pos) * 3 * incRec_dev[str_inc_recorrido];
    //std::cout << str_inc_recorrido << " : " << genome.gene(pos) << std::endl << std::endl;
    pos++;

    perturbation["intercept_factor"] = genome.gene(pos) * 3 * incRec_dev["intercept_factor"];
    //std::cout << "intercept_factor: " << perturbation["intercept_factor"] << std::endl;
    pos++;

    perturbation["error"] = genome.gene(pos) * 3 * incRec_dev["error"];
    //std::cout << "error: " << perturbation["error"] << std::endl;
    pos++;

    // Perturb the euclidean model
    eucModel->perturbModel(perturbation);

    //std::cout << std::endl;


    // Perturb the termination model
    avgDevModel->perturbModel(genome.gene(pos++) * 3 * avgDev_mean_dev,
                              genome.gene(pos++) * 3 * avgDev_stddev_dev);

    //std::cout << std::endl;


    /* ********************************************* */
    /* Sampling of new metrics from perturbed models */
    /* ********************************************* */

    long double pvalueBranches = 0.0;
    long double pvalueEuc      = 0.0;
    long double pvalueLength   = 0.0;

    // To obtain more stable p-values, we repeat the sampling/validating process 'valIters' times
    for (unsigned k = 0; k < valIters; k++) {

      // Then, these models are used to obtain as many metrics as needed
      // (as many as the available metrics for validation)
      unsigned nMetrics = bsData[i].second.numBranches.size();

      std::vector<long double> newNumBranches(nMetrics);
      std::vector<long double> newMaxEucDist (nMetrics);
      std::vector<long double> newTotalLength(nMetrics);

      for (unsigned j = 0; j < nMetrics; j++) {
        boost::unordered_map<std::string, long double> metrics;

        // Now we create a sampler with the three models and grow as many dendritic trees
        // as needed. From each tree, we obtain its associated metrics
        if (i < bsIters) {
          Dendrite sampler(*avgDevModel, *pathModel, *eucModel);
          sampler.grow();
          metrics = sampler.Metrics();
        }
        else {
          Dendrite sampler(*bsData[i].first.avgDev, *bsData[i].first.pathway, *bsData[i].first.euclidean);
          sampler.grow();
          metrics = sampler.Metrics();
        }

        /*
        if (metrics[str_numero_branches] == 0) {
          j--;
        }
        else {
        */
          newNumBranches[j] = metrics[str_numero_branches];
          newMaxEucDist [j] = metrics[str_max_distancia_euclidea];
          newTotalLength[j] = metrics[str_longitud_total];
        //}

      }

      // The metrics must be sorted in order to the KS Test to work
      std::sort(newNumBranches.begin(), newNumBranches.end());
      std::sort( newMaxEucDist.begin(),  newMaxEucDist.end());
      std::sort(newTotalLength.begin(), newTotalLength.end());

      // Now we obtain the p-value by means of K-S test
      pvalueBranches += KolmogorovTest(bsData[i].second.numBranches,    newNumBranches, "");
      pvalueEuc      += KolmogorovTest(bsData[i].second.maxEucDistance, newMaxEucDist,  "");
      pvalueLength   += KolmogorovTest(bsData[i].second.totalLength,    newTotalLength, "");

    }

    // Compute the mean of the valIters iterations
    pvalueBranches /= valIters;
    pvalueEuc      /= valIters;
    pvalueLength   /= valIters;

    // We compute the geometric mean of these p-values and add it to the cummulative vars pvalue_nores or pvalue_res
    // (depending on the iteration, the last one is for the resustitution phase). We also update the aux vars that store
    // the p-values of the individual models (for debugging purposes, so that we can print them individually as needed)
    if (i < bsIters) {
      pvalue_nores += pow(pvalueBranches * pvalueEuc * pvalueLength, 1. / 3.);
      pvalue_euc += pvalueEuc;
      pvalue_path += pvalueLength;
      pvalue_avg += pvalueBranches;
    }
    // In this case, we are in the last iteration (the resustitution phase). For this reason we do not add the current p-values.
    else {
      pvalue_res = pow(pvalueBranches * pvalueEuc * pvalueLength, 1. / 3.);
      pvalue_euc_res = pvalueEuc;
      pvalue_path_res = pvalueLength;
      pvalue_avg_res = pvalueBranches;
    }

    // Free memory
    delete eucModel;
    delete pathModel;
    delete avgDevModel;

  }

  // Compute the average p-value of the no-resustitution step
  pvalue_nores /= bsIters;

  // Compute the average individual p-values (for debugging purposes)
  pvalue_euc   /= bsIters;
  pvalue_path  /= bsIters;
  pvalue_avg   /= bsIters;

  // If we should print debugging information, do so
  if (pprocess) {
    std::cout << "==================================================================" << std::endl;
    std::cout << "pvalue Euc: "          << pvalue_euc      << std::endl;
    std::cout << "pvalue Length: "       << pvalue_path     << std::endl;
    std::cout << "pvalue Branches: "     << pvalue_avg      << std::endl;
    std::cout << std::endl;
    std::cout << "pvalue Euc Res: "      << pvalue_euc_res  << std::endl;
    std::cout << "pvalue Length Res: "   << pvalue_path_res << std::endl;
    std::cout << "pvalue Branches Res: " << pvalue_avg_res  << std::endl;
    std::cout << std::endl;
    std::cout << "res: "    << pvalue_res << std::endl;
    std::cout << "nores: "  << pvalue_nores << std::endl;
    std::cout << "pvalue: " << (0.368 * pvalue_res) + (0.632 * pvalue_nores) << std::endl;
    std::cout << "==================================================================" << std::endl;
  }

  return (0.368 * pvalue_res) + (0.632 * pvalue_nores);

}

extern "C" void individualInit (GAGenome& g) {
  return RealUniformInitializer (g);
}

extern "C" GAGenome* defineProblem () {

  GAEDAConfig* cfg = GAEDAConfig::handle();

  // Parse the config file for the problem
  if (!parse_config (cfg->getProblemData ()))
    throw runtime_error ("[Dendrites] Error: A valid configuration file must be provided.");

  if (!parse_data ())
    throw runtime_error ("[Dendrites] Error: Either the learning or the validation data could not be loaded.");

  // Perturbations take values in the interval [-1, 1]
  GAAlleleSet<long double> alleles (-1, 1);

  // Compute problem size from the parsed models
  unsigned sz = bsData[0].first.euclidean->ModelSize() + bsData[0].first.pathway->ModelSize() + 2;
  GA1DArrayAlleleGenome<long double>* genome = new GA1DArrayAlleleGenome<long double> (sz, alleles, objective);

  // Common operators
  genome->initializer (RealUniformInitializer);
  genome->comparator  (RealEuclideanComparator);

  // Specific stuff for GAs
  genome->crossover   (RealBlendCrossover);
  genome->mutator     (RealGaussianMutator);

  // Specific stuff for DE
  genome->crossover   (RealExponentialCrossover);

  // Specific stuff for MOS
  MOSGenomeFactory::handle()->registerGenome (GAID::RealEncoding, genome);

  return genome;

}

extern "C" const char *describeProblem (void) {
  return "Dendrites Grow";
}


extern "C" bool postprocess (GAPopulation* pop, int rank, int repetition) {
  pprocess = true;
  pop->best().evaluate(gaTrue);

  for (unsigned i = 0; i < bsIters + 1; i++) {
    delete bsData[i].first.euclidean;
    delete bsData[i].first.pathway;
    delete bsData[i].first.avgDev;
  }
  return true;
}




/*
    ap::real_1d_array br      = ap::real_1d_array();
    ap::real_1d_array br_new  = ap::real_1d_array();
    ap::real_1d_array euc     = ap::real_1d_array();
    ap::real_1d_array euc_new = ap::real_1d_array();
    ap::real_1d_array len     = ap::real_1d_array();
    ap::real_1d_array len_new = ap::real_1d_array();

    br.setcontent     (0, bsData[i].second.numBranches.size    () -1, &(bsData[i].second.numBranches   [0]));
    br_new.setcontent (0, newNumBranches.size                  () -1, &(newNumBranches                 [0]));
    euc.setcontent    (0, bsData[i].second.maxEucDistance.size () -1, &(bsData[i].second.maxEucDistance[0]));
    euc_new.setcontent(0, newMaxEucDist.size                   () -1, &(newMaxEucDist                  [0]));
    len.setcontent    (0, bsData[i].second.totalLength.size    () -1, &(bsData[i].second.totalLength   [0]));
    len_new.setcontent(0, newTotalLength.size                  () -1, &(newTotalLength                 [0]));

    long double bothtails, lefttail, righttail;

    studentttest2(br , bsData[i].second.numBranches.size    (), br_new , newNumBranches.size(), bothtails, lefttail, righttail);
    long double pvalueBranches = bothtails;

    studentttest2(euc, bsData[i].second.maxEucDistance.size (), euc_new, newMaxEucDist.size (), bothtails, lefttail, righttail);
    long double pvalueEuc = bothtails;

    studentttest2(len, bsData[i].second.totalLength.size    (), len_new, newTotalLength.size(), bothtails, lefttail, righttail);
    long double pvalueLength = bothtails;
*/

#include "attribute_names.h"
#include "Dendrite.h"
#include "kolmogorov.h"
#include "weka.h"

#include <libconfig.h++>

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

typedef std::pair<models, metrics> dataPerIter;

std::vector<dataPerIter> data;

// Maps to compute the average values of each component of the models and the standard deviations
boost::unordered_map<std::string, long double> incRec_avg, incEuc_avg, incRec_dev, incEuc_dev;

long double avgDev_mean_avg, avgDev_stddev_avg, avgDev_mean_dev, avgDev_stddev_dev;

unsigned valIters     = 0;
unsigned iters        = 0;
int      currentFold  = 0;
std::string directory = "";

std::vector<long double> best;

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
    iters       =               cfg.lookup ("iters"      );
    valIters    =               cfg.lookup ("valIters"   );
    currentFold =               cfg.lookup ("currentFold");

    libconfig::Setting& Best = cfg.lookup ("best");

    for (unsigned i = 0; i < (unsigned) Best.getLength (); i++)
       best.push_back ((long double) Best [i]);

    if (iters <= 0) {
      std::cerr << "[Dendrites] Error: iters (number of iterations for the k-fold) must be greater than zero." << std::endl;
      exit (-1);
    }

    if (currentFold < 0 || currentFold > iters) {
      std::cerr << "[Dendrites] Error: currentFold (current iteration of the k-fold cross validation) must be within the interval [0,iters]." << std::endl;
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

bool parse_data () {

  // Allocate memory for the data of the k-fold
  data.reserve(iters);

  // Vectors with the attributes for both regression models
  std::vector<std::string> attrib_incRec(NAMES, NAMES + LAST_ATT_I);

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

  LinearRegression* regModel;
  NormalDistribution* normalDist;

  // Main loop to parse the data for the k-fold
  for (unsigned i = 0; i < iters; i++) {

    dataPerIter dataIter;

    // for each file in 'directory' do (parse)
    stringstream filenamebuf;

    filenamebuf << "/models_train_fold" << i << ".csv/";

    regModel = new LinearRegression((char*) (directory + filenamebuf.str() + "iRecorrido.model").c_str(), attrib_incRec);
    dataIter.first.pathway = regModel;

    // Add the values of this model to the previous ones into the incRec_avg vector tu compute the mean afterwards
    regModel->addModel(incRec_avg);

    regModel = new LinearRegression((char*) (directory + filenamebuf.str() + "iEuclideo.model").c_str(), attrib_incEuc);
    dataIter.first.euclidean = regModel;

    // Add the values of this model to the previous ones into the incEuc_avg vector tu compute the mean afterwards
    regModel->addModel(incEuc_avg);

    normalDist = new NormalDistribution((char*) (directory + filenamebuf.str() + "iMediaDesv.model").c_str());
    dataIter.first.avgDev = normalDist;

    // Add the values of this model to the previous ones into the avgDev_mean_avg and avgDev_stddev_avg variables tu compute the mean afterwards
    avgDev_mean_avg   += normalDist->getMean();
    avgDev_stddev_avg += normalDist->getStdDev();

    stringstream filenamebuf2;
    filenamebuf2 << "/metrics_fold" << i << ".csv";

    std::ifstream metrics_file((directory + filenamebuf2.str()).c_str());

    if (metrics_file.fail()) {
      std::cerr << "[Dendrites] Metrics file '" << directory << filenamebuf2.str() << "' not found" << std::endl;
      exit(0);
    }

    std::string line;

    while (getline(metrics_file, line) && metrics_file.good()) {

      std::vector<std::string> res;
      csvline_populate(res, line, ',');

      // Check for the correct number of metrics
      assert(res.size() == 3);

      dataIter.second.numBranches.push_back   (strtod(res[0].c_str(), NULL));
      dataIter.second.maxEucDistance.push_back(strtod(res[1].c_str(), NULL));
      dataIter.second.totalLength.push_back   (strtod(res[2].c_str(), NULL));

    }

    std::sort(dataIter.second.numBranches.begin(),    dataIter.second.numBranches.end()   );
    std::sort(dataIter.second.maxEucDistance.begin(), dataIter.second.maxEucDistance.end());
    std::sort(dataIter.second.totalLength.begin(),    dataIter.second.totalLength.end()   );

    assert(dataIter.second.numBranches.size()    == dataIter.second.maxEucDistance.size() &&
           dataIter.second.numBranches.size()    == dataIter.second.totalLength.size()    &&
           dataIter.second.maxEucDistance.size() == dataIter.second.totalLength.size()      );

    data.push_back(dataIter);

  }

  /* ************************************************************************************** */
  /* Computation of average values and standard deviations for each attribute of each model */
  /* ************************************************************************************** */

  boost::unordered_map<std::string, long double>::iterator it;

  // Compute the avg and the stddev for each attribute of the path model
  for (it = incRec_avg.begin(); it != incRec_avg.end(); it++) {
    // Divide the overall value for this attribute (we only added them before) by the number of iterations
    it->second /= iters;

    // Compute the cummulated differences for the stddev
    long double diffs_acum = 0.0;
    for (unsigned i = 0; i < iters; i++)
      diffs_acum += pow(data[i].first.pathway->getModel()[it->first] - it->second, 2);

    // Compute the stddev for that attribute
    incRec_dev[it->first] = sqrt(diffs_acum / iters);
  }


  // Compute the avg and the stddev for each attribute of the euclidean model
  for (it = incEuc_avg.begin(); it != incEuc_avg.end(); it++) {
    // Divide the overall value for this attribute (we only added them before) by the number of iterations
    it->second /= iters;

    // Compute the cummulated differences for the stddev
    long double diffs_acum = 0.0;
    for (unsigned i = 0; i < iters; i++)
      diffs_acum += pow(data[i].first.euclidean->getModel()[it->first] - it->second, 2);

    // Compute the stddev for that attribute
    incEuc_dev[it->first] = sqrt(diffs_acum / iters);
  }


  // Compute the avg and the stddev for both the mean and the stddev of the termination model
  avgDev_mean_avg /= iters;
  avgDev_stddev_avg /= iters;

  // Compute the cummulated differences for both attributes
  for (unsigned i = 0; i < iters; i++) {
    avgDev_mean_dev   += pow(data[i].first.avgDev->getMean()   - avgDev_mean_avg,   2);
    avgDev_stddev_dev += pow(data[i].first.avgDev->getStdDev() - avgDev_stddev_avg, 2);
  }

  // Compute the stddev for both attributes
  avgDev_mean_dev   = sqrt(avgDev_mean_dev   / iters);
  avgDev_stddev_dev = sqrt(avgDev_stddev_dev / iters);

  return true;

}


int main (int argc, char** argv) {

  parse_config(argv[1]);
  parse_data();

  std::cout << "==================================================================" << std::endl;
  std::cout << " Fold: " << currentFold << "====================================================" << std::endl;
  std::cout << "==================================================================" << std::endl;

  long double branches_acum = 0, euc_acum = 0, len_acum = 0;

  // Create copies of these models for the perturbation take place
  LinearRegression*   eucModel    = new LinearRegression  (*data[currentFold].first.euclidean);
  LinearRegression*   pathModel   = new LinearRegression  (*data[currentFold].first.pathway  );
  NormalDistribution* avgDevModel = new NormalDistribution(*data[currentFold].first.avgDev   );

  boost::unordered_map<std::string, long double> perturbation;
  boost::unordered_map<std::string, long double>::iterator it;
  unsigned pos = 0;

  for (; pos < LAST_ATT_I; pos++)
    perturbation[NAMES[pos]] = best[pos] * 3 * incRec_dev[NAMES[pos]];

  perturbation["intercept_factor"] = best[pos] * 3 * incRec_dev["intercept_factor"];
  pos++;

  perturbation["error"] = best[pos] * 3 * incRec_dev["error"];
  pos++;

  pathModel->perturbModel(perturbation);

  for (unsigned pos2 = 0; pos2 < LAST_ATT_I; pos2++, pos++)
    perturbation[NAMES[pos2]] = best[pos] * 3 * incRec_dev[NAMES[pos2]];

  perturbation[str_inc_recorrido] = best[pos] * 3 * incRec_dev[str_inc_recorrido];
  pos++;

  perturbation["intercept_factor"] = best[pos] * 3 * incRec_dev["intercept_factor"];
  pos++;

  perturbation["error"] = best[pos] * 3 * incRec_dev["error"];
  pos++;

  eucModel->perturbModel(perturbation);

  avgDevModel->perturbModel(best[pos++] * 3 * avgDev_mean_dev,
                            best[pos++] * 3 * avgDev_stddev_dev);

  long double pvalueBranches = 0.0;
  long double pvalueEuc      = 0.0;
  long double pvalueLength   = 0.0;

  // To obtain more stable p-values, we repeat the sampling/validating process 'valIters' times
  for (unsigned k = 0; k < valIters; k++) {

    // Then, these models are used to obtain as many metrics as needed
    // (as many as the available metrics for validation)
    unsigned nMetrics = data[currentFold].second.numBranches.size();

    std::vector<long double> newNumBranches(nMetrics);
    std::vector<long double> newMaxEucDist (nMetrics);
    std::vector<long double> newTotalLength(nMetrics);

    for (unsigned j = 0; j < nMetrics; j++) {
      boost::unordered_map<std::string, long double> metrics;

      // Now we create a sampler with the three models and grow as many dendritic trees
      // as needed. From each tree, we obtain its associated metrics
      Dendrite sampler(*avgDevModel, *pathModel, *eucModel);

      sampler.grow();
      metrics = sampler.Metrics();

      newNumBranches[j] = metrics[str_numero_branches];
      newMaxEucDist [j] = metrics[str_max_distancia_euclidea];
      newTotalLength[j] = metrics[str_longitud_total];
    }

    // Now we obtain the p-value by means of K-S test
    std::sort(newNumBranches.begin(), newNumBranches.end());
    std::sort( newMaxEucDist.begin(),  newMaxEucDist.end());
    std::sort(newTotalLength.begin(), newTotalLength.end());

    pvalueBranches += KolmogorovTest(data[currentFold].second.numBranches,    newNumBranches, "");
    pvalueEuc      += KolmogorovTest(data[currentFold].second.maxEucDistance, newMaxEucDist,  "");
    pvalueLength   += KolmogorovTest(data[currentFold].second.totalLength,    newTotalLength, "");

  }

  // Compute the mean of the valIters iterations
  pvalueBranches /= valIters;
  pvalueEuc      /= valIters;
  pvalueLength   /= valIters;

  long double pvalue = pow(pvalueBranches * pvalueEuc * pvalueLength, 1. / 3.);

  std::cout << "  ============" << std::endl;
  std::cout << "  = p-values =" << std::endl;
  std::cout << "  ============" << std::endl << std::endl;

  std::cout << "  branches: " << pvalueBranches << std::endl;
  std::cout << "  euc: " << pvalueEuc << std::endl;
  std::cout << "  length: " << pvalueLength << std::endl << std::endl;
  std::cout << "  p-value: " << pvalue << std::endl;

  return 0;

}

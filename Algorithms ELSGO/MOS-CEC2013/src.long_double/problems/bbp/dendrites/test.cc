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

typedef std::pair<models, metrics> dataPerIter;

std::vector<dataPerIter> data;


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

bool parse_data (unsigned iters, std::string directory) {

  // Allocate memory for the data of the k-fold
  data.reserve(iters);

  // Vectors with the attributes for both regression models
  std::vector<std::string> attrib_incRec(NAMES, NAMES + LAST_ATT_I);

  std::vector<std::string> attrib_incEuc(attrib_incRec);
  attrib_incEuc.push_back(str_inc_recorrido);

  LinearRegression* regModel;
  NormalDistribution* normalDist;

  // Main loop to parse the data for the k-fold
  for (unsigned i = 0; i < iters; i++) {

    dataPerIter dataIter;

    // for each file in 'directory' do (parse)
    stringstream filenamebuf;

    filenamebuf << "/models_train_fold" << i << ".csv/";

//     std::cout << "iEuclideo: "  << directory << filenamebuf.str() << "iEuclideo.model"  << std::endl;
//     std::cout << "iMediaDesv: " << directory << filenamebuf.str() << "iMediaDesv.model" << std::endl;
//     std::cout << "iRecorrido: " << directory << filenamebuf.str() << "iRecorrido.model" << std::endl;

    regModel = new LinearRegression((char*) (directory + filenamebuf.str() + "iRecorrido.model").c_str(), attrib_incRec);
    dataIter.first.pathway = regModel;
    // std::cout << "regModel: " << regModel->ModelSize() << std::endl;

    regModel = new LinearRegression((char*) (directory + filenamebuf.str() + "iEuclideo.model").c_str(), attrib_incEuc);
    dataIter.first.euclidean = regModel;
    // std::cout << "regModel: " << regModel->ModelSize() << std::endl;

    normalDist = new NormalDistribution((char*) (directory + filenamebuf.str() + "iMediaDesv.model").c_str());
    dataIter.first.avgDev = normalDist;

    stringstream filenamebuf2;
    filenamebuf2 << "/metrics_fold" << i << ".csv";

//     std::cout << "metrics: " << directory  << filenamebuf2.str() << std::endl;

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

//    std::cout << "Metrics length: " << dataIter.second.maxEucDistance.size() << std::endl;

    data.push_back(dataIter);

  }

  return true;

}


int main (int argc, char** argv) {

  unsigned iters = 10;
  int valIters = 10;
  std::string dir = "/scratch/atorre/7picos-4/dendrites/resustitution/";

  parse_data(iters, dir);

  for (unsigned i = 0; i < iters; i++) {

    std::cout << "==================================================================" << std::endl;
    std::cout << " Fold: " << i << "====================================================" << std::endl;
    std::cout << "==================================================================" << std::endl;

    // Create copies of these models for the perturbation take place
    LinearRegression*   eucModel    = new LinearRegression  (*data[i].first.euclidean);
    LinearRegression*   pathModel   = new LinearRegression  (*data[i].first.pathway  );
    NormalDistribution* avgDevModel = new NormalDistribution(*data[i].first.avgDev   );

    long double pvalueBranches = 0.0;
    long double pvalueEuc      = 0.0;
    long double pvalueLength   = 0.0;

    // To obtain more stable p-values, we repeat the sampling/validating process 'valIters' times
    for (unsigned k = 0; k < valIters; k++) {

      // Then, these models are used to obtain as many metrics as needed
      // (as many as the available metrics for validation)
      unsigned nMetrics = data[i].second.numBranches.size();

      std::vector<long double> newNumBranches(nMetrics);
      std::vector<long double> newMaxEucDist (nMetrics);
      std::vector<long double> newTotalLength(nMetrics);

      for (unsigned j = 0; j < nMetrics; j++) {
        boost::unordered_map<std::string, long double> metrics;

        // Now we create a sampler with the three models and grow as many dendritic trees
        // as needed. From each tree, we obtain its associated metrics
        Dendrite sampler(*data[i].first.avgDev, *data[i].first.pathway, *data[i].first.euclidean);

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

      pvalueBranches += KolmogorovTest(data[i].second.numBranches,    newNumBranches, "");
      pvalueEuc      += KolmogorovTest(data[i].second.maxEucDistance, newMaxEucDist,  "");
      pvalueLength   += KolmogorovTest(data[i].second.totalLength,    newTotalLength, "");

    }

    // Compute the mean of the valIters iterations
    pvalueBranches /= valIters;
    pvalueEuc      /= valIters;
    pvalueLength   /= valIters;

    std::cout << "  ============" << std::endl;
    std::cout << "  = p-values =" << std::endl;
    std::cout << "  ============" << std::endl << std::endl;

    std::cout << "  branches: " << pvalueBranches << std::endl;
    std::cout << "  euc: "      << pvalueEuc      << std::endl;
    std::cout << "  length: "   << pvalueLength   << std::endl << std::endl;

  }

  return 0;

}


int main_test (int argc, char** argv) {

  unsigned iters = 10;
  int nIters = 1;
  std::string dir = "/scratch/atorre/7picos-4/dendrites/resustitution/";

  parse_data(iters, dir);

  for (unsigned i = 0; i < iters; i++) {

    std::cout << "==================================================================" << std::endl;
    std::cout << " Fold: " << i << "====================================================" << std::endl;
    std::cout << "==================================================================" << std::endl;

    for (int p = 1; p <= nIters; p++) {

      // Create copies of these models for the perturbation take place
      LinearRegression*   eucModel    = new LinearRegression  (*data[i].first.euclidean);
      LinearRegression*   pathModel   = new LinearRegression  (*data[i].first.pathway  );
      NormalDistribution* avgDevModel = new NormalDistribution(*data[i].first.avgDev   );

      // Then, these models are used to obtain as many metrics as needed
      // (as many as the available metrics for validation)
      unsigned nMetrics = pow(10.0, p) * data[i].second.numBranches.size();

      std::vector<long double> newNumBranches(nMetrics);
      std::vector<long double> newMaxEucDist (nMetrics);
      std::vector<long double> newTotalLength(nMetrics);

      for (unsigned j = 0; j < nMetrics; j++) {
        boost::unordered_map<std::string, long double> metrics;

        // Now we create a sampler with the three models and grow as many dendritic trees
        // as needed. From each tree, we obtain its associated metrics
        Dendrite sampler(*data[i].first.avgDev, *data[i].first.pathway, *data[i].first.euclidean);

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

      long double pvalueBranches = KolmogorovTest(data[i].second.numBranches,    newNumBranches, "");
      long double pvalueEuc      = KolmogorovTest(data[i].second.maxEucDistance, newMaxEucDist,  "");
      long double pvalueLength   = KolmogorovTest(data[i].second.totalLength,    newTotalLength, "");

      std::cout << "  ============" << std::endl;
      std::cout << "  = p-values =" << std::endl;
      std::cout << "  ============" << std::endl << std::endl;

      std::cout << "  branches: " << pvalueBranches << std::endl;
      std::cout << "  euc: "      << pvalueEuc      << std::endl;
      std::cout << "  length: "   << pvalueLength   << std::endl << std::endl;

      unsigned rep = 100;

      std::cout << "  ======================================================" << std::endl;
      std::cout << "  = Tolerance to sample replacement (" << rep << " iters) =" << std::endl;
      std::cout << "  ======================================================" << std::endl;

      for (unsigned z = 0; z < rep; z++) {
        boost::unordered_map<std::string, long double> metrics;

        // Now we create a sampler with the three models and grow as many dendritic trees
        // as needed. From each tree, we obtain its associated metrics
        Dendrite sampler(*data[i].first.avgDev, *data[i].first.pathway, *data[i].first.euclidean);

        sampler.grow();
        metrics = sampler.Metrics();

        unsigned randPos = GARandomInt(0, nMetrics);

        newNumBranches[randPos] = metrics[str_numero_branches];
        newMaxEucDist [randPos] = metrics[str_max_distancia_euclidea];
        newTotalLength[randPos] = metrics[str_longitud_total];

        // Now we obtain the p-value by means of K-S test
        std::sort(newNumBranches.begin(), newNumBranches.end());
        std::sort( newMaxEucDist.begin(),  newMaxEucDist.end());
        std::sort(newTotalLength.begin(), newTotalLength.end());

        pvalueBranches = KolmogorovTest(data[i].second.numBranches,    newNumBranches, "");
        pvalueEuc      = KolmogorovTest(data[i].second.maxEucDistance, newMaxEucDist,  "");
        pvalueLength   = KolmogorovTest(data[i].second.totalLength,    newTotalLength, "");

        std::cout << "  branches: " << pvalueBranches << std::endl;
        std::cout << "  euc: " << pvalueEuc << std::endl;
        std::cout << "  length: " << pvalueLength << std::endl << std::endl;
      }

    }

  }

  return 0;

}


int main_other (int argc, char** argv) {

  unsigned iters = 1;
  std::string dir = "/scratch/atorre/7picos-4/dendrites/resustitution/";

  parse_data(iters, dir);

  for (unsigned i = 0; i < iters; i++) {

    std::cout << "==================================================================" << std::endl;
    std::cout << " Fold: " << i << "====================================================" << std::endl;
    std::cout << "==================================================================" << std::endl;

    long double branches_acum = 0, euc_acum = 0, len_acum = 0;

    long double best[] = {0.05972187829191, 0.29779727538014, 0.59706682227415, -0.80947874943236, -0.42551609847020, -0.11565566580540, 0.73890077031166, 0.49946019961474, 0.78102556792136, -0.95666897667417, 0.43557121764662, 0.57921600974128, -0.60280865412336, -0.38899560709903, -0.12810965493699, -0.26231747645061, 0.41435924145130, 0.34445060572794, 0.26786990615906, -0.41291502742698};

    int nIters = 1;

    // Maps to compute the average values of each component of the models and the standard deviations
    boost::unordered_map<std::string, long double> incRec_avg, incEuc_avg, incRec_dev, incEuc_dev;

    // Vectors with the attributes for both regression models
    std::vector<std::string> attrib_incRec(NAMES, NAMES + LAST_ATT_I);

    std::vector<std::string> attrib_incEuc(attrib_incRec);
    attrib_incEuc.push_back(str_inc_recorrido);

    for (unsigned k = 0; k < attrib_incRec.size(); k++) {
      incRec_avg[attrib_incRec[k]] = 0.0;
      incRec_dev[attrib_incRec[k]] = 0.0;
      std::cout << attrib_incRec[k] << " : " << incRec_dev[attrib_incRec[k]] << std::endl;
    }

    std::cout << std::endl;

    for (unsigned k = 0; k < attrib_incEuc.size(); k++) {
      incEuc_avg[attrib_incEuc[k]] = 0.0;
      incEuc_dev[attrib_incEuc[k]] = 0.0;
      std::cout << attrib_incEuc[k] << " : " << incEuc_dev[attrib_incEuc[k]] << std::endl;
    }

    std::cout << std::endl;

    long double avgDev_mean_avg, avgDev_stddev_avg, avgDev_mean_dev, avgDev_stddev_dev;
    avgDev_mean_avg = avgDev_stddev_avg = avgDev_mean_dev = avgDev_stddev_dev = 0.0;

    for (int p = 1; p <= nIters; p++) {

      // Create copies of these models for the perturbation take place
      LinearRegression*   eucModel    = new LinearRegression  (*data[i].first.euclidean);
      LinearRegression*   pathModel   = new LinearRegression  (*data[i].first.pathway  );
      NormalDistribution* avgDevModel = new NormalDistribution(*data[i].first.avgDev   );

      boost::unordered_map<std::string, long double> perturbation;
      boost::unordered_map<std::string, long double>::iterator it;
      unsigned pos = 0;

      for (; pos < LAST_ATT_I; pos++) {
        perturbation[NAMES[pos]] = best[pos];
        std::cout << NAMES[pos] << " : " << best[pos] << std::endl;
      }

      std::cout << std::endl;

      pathModel->perturbModel(perturbation);

      std::cout << std::endl;

      for (unsigned pos2 = 0; pos2 < LAST_ATT_I; pos2++, pos++) {
        perturbation[NAMES[pos2]] = best[pos];
        std::cout << NAMES[pos2] << " : " << best[pos] << std::endl;
      }

      perturbation[str_inc_recorrido] = best[pos];
      std::cout << str_inc_recorrido << " : " << best[pos] << std::endl << std::endl;
      pos++;

      eucModel->perturbModel(perturbation);

      std::cout << std::endl;

      avgDevModel->perturbModel(best[pos++],
                                best[pos++]);

      std::cout << std::endl;

      // Then, these models are used to obtain as many metrics as needed
      // (as many as the available metrics for validation)
      unsigned nMetrics = data[i].second.numBranches.size();

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

      long double pvalueBranches = KolmogorovTest(data[i].second.numBranches,    newNumBranches, "");
      long double pvalueEuc      = KolmogorovTest(data[i].second.maxEucDistance, newMaxEucDist,  "");
      long double pvalueLength   = KolmogorovTest(data[i].second.totalLength,    newTotalLength, "");

      branches_acum += pvalueBranches;
      euc_acum += pvalueEuc;
      len_acum += pvalueLength;

      std::cout << "Old branches: ";
      for (unsigned kkk = 0; kkk < data[i].second.numBranches.size(); kkk++)
        std::cout << data[i].second.numBranches[kkk] << " ";
      std::cout << std::endl << std::endl;

      std::cout << "New branches: ";
      for (unsigned kkk = 0; kkk < newNumBranches.size(); kkk++)
        std::cout << newNumBranches[kkk] << " ";
      std::cout << std::endl << std::endl;

      std::cout << "Old euc: ";
      for (unsigned kkk = 0; kkk < data[i].second.maxEucDistance.size(); kkk++)
        std::cout << data[i].second.maxEucDistance[kkk] << " ";
      std::cout << std::endl << std::endl;

      std::cout << "New euc: ";
      for (unsigned kkk = 0; kkk < newMaxEucDist.size(); kkk++)
        std::cout << newMaxEucDist[kkk] << " ";
      std::cout << std::endl << std::endl;

      std::cout << "Old length: ";
      for (unsigned kkk = 0; kkk < data[i].second.totalLength.size(); kkk++)
        std::cout << data[i].second.totalLength[kkk] << " ";
      std::cout << std::endl << std::endl;

      std::cout << "New length: ";
      for (unsigned kkk = 0; kkk < newTotalLength.size(); kkk++)
        std::cout << newTotalLength[kkk] << " ";
      std::cout << std::endl << std::endl;


      std::cout << "  branches: " << pvalueBranches << std::endl;
      std::cout << "  euc: " << pvalueEuc << std::endl;
      std::cout << "  length: " << pvalueLength << std::endl << std::endl;

  }

  std::cout << "  ============" << std::endl;
  std::cout << "  = p-values =" << std::endl;
  std::cout << "  ============" << std::endl << std::endl;

  std::cout << "  branches: " << branches_acum / nIters << std::endl;
  std::cout << "  euc: " << euc_acum / nIters << std::endl;
  std::cout << "  length: " << len_acum / nIters << std::endl << std::endl;

  }

  return 0;

}

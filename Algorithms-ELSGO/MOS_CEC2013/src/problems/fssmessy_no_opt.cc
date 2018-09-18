#include <stdexcept>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <math.h>

#include <libconfig.h++>

#include <GAIntOps.h>
#include <MOSGenomeFactory.h>
#include <genomes/FSSGenome.h>
#include <GAEDAConfig.h>
#include <quicksort.h>

/*
/////////////////////////////////////////////////////////////
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
*/

enum FitnessMeasures {AUC, AUC_PIV, LL, LLS, BS, ACC};

// Number of neighbours for knn
int nNeighbors = 7;

// Number of iterations for Bootstrap
unsigned iters = 200;

// Dataset (name and actual data)
std::string dataset;
std::vector< std::vector<double> > data;

long pivot = 0;

unsigned fit_measure = AUC;

/*
/////////////////////////////////////////////////////////////
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
*/

bool parse_config (char* fich) {

   libconfig::Config cfg;
   cfg.setAutoConvert (true);

   try {
      cfg.readFile (fich);
   }
   catch (libconfig::ParseException e) {
      std::cerr << "Error: parsing the configuration file in line (" << e.getLine () << "): " << e.getError() << "." << std::endl;
      exit (-1);
   }
   catch (libconfig::FileIOException e) {
      std::cerr << "Error: file not found or could not be open." << std::endl;
      exit (-1);
   }

   // First, we read mandatory parameters
   try {
      dataset    = (const char*) cfg.lookup ("dataset");
      nNeighbors =               cfg.lookup ("neighbors");
      iters      =               cfg.lookup ("iters");
      pivot      =               cfg.lookup ("pivot");

      std::string measure;
      measure    = (const char*) cfg.lookup ("measure");

      if (measure == "auc")
        if (pivot == 0)
          fit_measure = AUC;
        else
          fit_measure = AUC_PIV;
      else if (measure == "ll")
        fit_measure = LL;
      else if (measure == "lls")
        fit_measure = LLS;
      else if (measure == "bs")
        fit_measure = BS;
      else if (measure == "acc")
        fit_measure = ACC;

   }
   catch (libconfig::SettingNotFoundException e) {
      std::cerr << "Error: one or more mandatory parameter(s) were not found." << std::endl;
      exit (-1);
   }

   return true;

}


void strSplitDouble (string str, string delim, vector<double>& results) {

   std::string::size_type cutAt;

   while ((cutAt = str.find_first_of (delim)) != str.npos) {

      if (cutAt > 0)
         results.push_back (strtod ((str.substr (0, cutAt)).c_str (), NULL));

      str = str.substr (cutAt + 1);

   }

   if (str.length () > 0)
      results.push_back (strtod (str.c_str (), NULL));

}


void strSplitInt (string str, string delim, vector<int>& results) {

   std::string::size_type cutAt;

   while ((cutAt = str.find_first_of (delim)) != str.npos) {

      if (cutAt > 0)
         results.push_back (atoi ((str.substr (0, cutAt)).c_str ()));

      str = str.substr (cutAt + 1);

   }

   if (str.length () > 0)
      results.push_back (atoi (str.c_str ()));

}


bool parse_data (void) {

   ifstream arch (dataset.c_str ());
   string line;

   while (arch.bad ())
      std::cerr << "Error!!" << endl;

   // Discard headers
   getline (arch, line);

   while (!arch.eof ()) {

      vector<double> columns;

      getline (arch, line);
      strSplitDouble (line, ",", columns);

      if (columns.size () > 0)
         data.push_back (columns);

   }

   arch.close ();

   return true;

}


/*
/////////////////////////////////////////////////////////////
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
*/


// Distance functions
double EuclideanDistance (const double* huesped, const double* vecino, unsigned size, const int* genes) {

  double num = 0.0;

  for (unsigned i = 0; i < size; i++)
    if (genes [i])
      num += (huesped [i] - vecino [i]) * (huesped [i] - vecino [i]);

  num = sqrt (num);

  return num;

}


double ManhattanDistance (const double* huesped, const double* vecino, unsigned size) {

   double num = 0.0;

   for (unsigned i = 0; i < size; i++)
      num += fabs (huesped [i] - vecino [i]);

  return num;

}


double ChevyshevDistance (const double* huesped, const double* vecino, unsigned size) {

   double max = 0.0;

   for (unsigned i = 0; i < size; i++)
      if (fabs (huesped [i] - vecino [i]) > max)
         max = fabs (huesped [i] - vecino [i]);

  return max;

}


double distance (const double* huesped, const double* vecino, unsigned size, const int* genes) {

  return EuclideanDistance (huesped, vecino, size, genes);
   //return ManhattanDistance (huesped, vecino, size);
   //return ChevyshevDistance (huesped, vecino, size);

}


// Random number between an interval
double rand2 (double min, double max) {

   double diff = max - min;
   double n = (rand() / (RAND_MAX + 1.0));

   return (min + n * diff);

}


// Function to select instances for learning
int select_instances (vector<int>& res, unsigned sz) {

   for (unsigned i = 0; i < sz; i++) {

      int pos = (int) rand2 (0, sz - 1);
      res [pos]++;

   }

   return 0;

}


// Function to calculate the AUC
double ROCArea (vector<double>& YNew, vector<double>& Y) {

   int seguir = 1;
   int numUnos = 0;
   int pE, fE, index;
   int len = Y.size ();

   double aux = 0.0;
   double area = 0.0;
   double indiceMedio = 0.0;

   double rango [len];


   double sumY = 0.0, sumYNew = 0.0;

   for (unsigned i = 0; i < Y.size (); i++) {

      sumY+= Y [i];
      sumYNew+= YNew [i];

   }

   if ((sumY    == (double) Y.size    () || sumY    == 0.0) ||
       (sumYNew == (double) YNew.size () || sumYNew == 0.0)    ) {

      if (sumY == sumYNew)
         return 1.0;
      else
         return 0.0;

   }

   // Ordenar las listas de menor a mayor en result
   // La lista class debe tener
   while (seguir) {

      seguir = 0;

      for (int i = 0; i < len - 1 ; i++) {

         if (YNew [i + 1] < YNew [i]) {

            aux = YNew [i];
            YNew [i] = YNew [i + 1];
            YNew [i + 1] = aux;

            aux = Y [i];
            Y [i] = Y [i + 1];
            Y [i + 1] = aux;

            seguir = 1;
         }

      }

   }


   pE = fE = 0;
   indiceMedio = 0.0;
   index = 1;

   while (fE < len - 1) {

      if (YNew [fE + 1] > YNew [fE]) {

         if (pE == fE) {

            rango [fE] = index;
            pE++;
            fE++;
            index++;

         }
         else { // Rellenar todas las posiciones intermedias
                  // Asegurar que pE y fE son los mismos de nuevo

            for (int i = pE; i <= fE ; i++)
               rango [i]= (double) (indiceMedio + index) / (double) (fE - pE + 1);

            fE++;
            pE = fE;
            index++;
            indiceMedio = 0;

         }

      }
      else { // hay empate

         // poner rangos al final del empate
         // pE=fE;
         fE++;
         indiceMedio += index;
         index++;

      }

   }

   if (pE != fE)
      for (int i = pE; i <= fE ; i++)
         rango [i]= (double) (indiceMedio + index) / (double) (fE - pE + 1);
   else
      rango [pE]= index;

   numUnos = 0;
   indiceMedio = 0.0;

   for (int i = 0 ; i < len; i++)
      if (Y [i] == 1) {

         numUnos++;
         indiceMedio += rango [i];

      }

//    // This avoids division by zero
//    if ((numUnos == len) || (numUnos == 0))
//       return 0;

   indiceMedio /= (double) numUnos;

   area = (indiceMedio - ((double) (numUnos + 1) / 2.0)) / (double) (len - numUnos);

   return area;

}


double linear_interpolation (double distance) {

   // F(X) = F(X_MAX) + (((F(X_MAX) - F(X_MIN)) / (X_MAX - X_MIN)) * (X - X_MIN)
   //return max_w + ((((max_w - min_w)/(max_distance-min_distance))
   //* (distance - min_distance));

   if (distance == 0.0)
      return 10000.0;

   return 1.0 / distance;

}


// Function to calculate the LogL
double logL (vector<double>& YNew, vector<double>& Y) {

   double logl = 0.0;
   double prob;

   for (unsigned i = 0; i < YNew.size (); i++) {

      if (YNew [i] == 0.0)
         prob = 0.00001;
      else if (YNew [i] == 1.0)
         prob = 0.99999;
      else
         prob = YNew [i];

      logl += (Y [i] * log (prob)) + ((1 - Y [i]) * log (1 - prob));

   }

   return logl;

}


// Function to calculate the Brier Score
double brierScore (vector<double>& YNew, vector<double>& Y) {

   double bs = 0.0;

   for (unsigned i = 0; i < YNew.size (); i++)
      bs += (Y [i] - YNew [i]) * (Y [i] - YNew [i]);

   return bs;

}


// Function to calculate the accuracy
double accuracy (vector<double>& YNew, vector<double>& Y) {

   double acum = 0.0;

   for (unsigned i = 0; i < YNew.size (); i++)
      if ((Y [i] <= 0.5 && YNew [i] <= 0.5) || (Y [i] >  0.5 && YNew [i] >  0.5))
         acum++;

   return acum / YNew.size ();

}


// Algoritmo KNN
bool knn (const vector<int>& learn, const vector<int>& test, unsigned learnSize, unsigned testSize, const int* genes,
          double& ll_val, double& auc_val, double& bs_val, double& acc_val             ) {

   int    keys [learnSize];
   int   keys2 [learnSize];
   double dist [learnSize];

   vector<double> real      (testSize, 0.0);
   vector<double> predicted (testSize, 0.0);

   nNeighbors = ((unsigned) nNeighbors <= learnSize) ? nNeighbors : learnSize;

   unsigned it = 0;

   for (unsigned i = 0; i < learnSize; i++) {

      if (test [i] > 0) {

         for (unsigned j = 0, k = 0; j < learnSize && k < learnSize; j++) {

            int n = learn [j];

            while (n > 0) {

               keys [k] = k;
               keys2 [k] = j;

               // We ignore first column from each instance because it is the class
               dist [k] = distance (&(data [i][1]), &(data [j][1]), data [i].size () - 1, genes);

               n--;
               k++;

            }

         }

         // Sort all the values
         quicksort (keys, dist, 0, learnSize - 1);

         double num_zeros = 0;
         double num_ones  = 0;

         for (int k = 0; k < nNeighbors; k++) {

            if (data [keys2[keys [k]]][0] == 0.0)
               num_zeros += linear_interpolation (dist [k]);
            else
               num_ones += linear_interpolation (dist [k]);

         }

         // We calculate the probability of being 1.0
         predicted[it] = (num_ones / (num_ones + num_zeros));

         // First column of each instance stores the class
         real[it] = (data [i][0]);

         it++;

      }

   }

   ll_val  = logL       (predicted, real) / testSize;
   auc_val = ROCArea    (predicted, real);
   bs_val  = brierScore (predicted, real);
   acc_val = accuracy   (predicted, real);

   return true;

}


extern "C" double objective (GAGenome& g) {

   FSSGenome<int>& genome = DYN_CAST (FSSGenome<int>&, g);

   // Score initialization
   double ll = 0.0, auc = 0.0, bs = 0.0, acc = 0.0;

   unsigned size = genome.size ();
   int genes[size];

   // Count selected genes
   unsigned nGenes = 0;
   for (unsigned i = 0; i < size; i++) {
      int gene = genome.gene (i);
      genes [i] = gene;
      if (gene)
         nGenes++;
   }

   // We copy all the instances to a tmp vector
   unsigned totalInstances = data.size ();
   std::vector<int> instances (totalInstances, 1);

   // Resustitution
   double ll_res, auc_res, bs_res, acc_res;
   knn (instances, instances, totalInstances, totalInstances, genes, ll_res, auc_res, bs_res, acc_res);

   // No resustitution
   std::vector<double> ll_nores (iters, 0), auc_nores (iters, 0), bs_nores (iters, 0), acc_nores (iters, 0);

   for (unsigned i = 0; i < iters; i++) {

     std::vector<int> learn (totalInstances, 0), test (totalInstances, 0);

     // Select instances to learn
     select_instances (learn, totalInstances);

     unsigned nTestInstances = 0;

     for (unsigned j = 0; j < totalInstances; j++)
       if (learn [j] == 0) { // Instance selected to validate
          test [j] = 1;
          nTestInstances++;
       }

      double r1, r2, r3, r4;

      knn (learn, test, totalInstances, nTestInstances, genes, r1, r2, r3, r4);

      ll_nores [i] = r1;
      auc_nores [i] = r2;
      bs_nores [i] = r3;
      acc_nores [i] = r4;

   }

   double ll_avg = 0.0, auc_avg = 0.0, bs_avg = 0.0, acc_avg = 0.0;

   for (unsigned i = 0; i < ll_nores.size (); i++) {

      ll_avg  += ll_nores  [i];
      auc_avg += auc_nores [i];
      bs_avg  += bs_nores  [i];
      acc_avg += acc_nores [i];

   }

   ll_avg  /= ll_nores.size  ();
   auc_avg /= auc_nores.size ();
   bs_avg  /= bs_nores.size  ();
   acc_avg /= acc_nores.size ();

   ll  = (0.368 * ll_res ) + (0.632 * ll_avg );
   auc = (0.368 * auc_res) + (0.632 * auc_avg);
   bs  = (0.368 * bs_res ) + (0.632 * bs_avg );
   acc = (0.368 * acc_res) + (0.632 * acc_avg);

   ll *= totalInstances;

   genome.setLL  (ll );
   genome.setLLS (ll - (2 * nGenes));
   genome.setAUC (auc);
   genome.setBS  (bs );
   genome.setACC (acc);

   double res = 0.0;

   switch (fit_measure) {
   case AUC:
     res = auc;
     break;
   case AUC_PIV:
     res =  auc * (1.0 /  (double) (abs ((long) nGenes - pivot) + 1));
     break;
   case LL:
     res = 1 / - ll;
     break;
   case LLS:
     res = 1 / - (ll - (2 * abs ((long) nGenes - pivot)));
     break;
   case BS:
     res = bs;
     break;
   case ACC:
     res = acc;
     break;
   default:
     throw runtime_error ("Error: Invalid fitness measure.");
     break;
   }

   //res = (auc + (1.0 - ((double) nGenes / (double) genome.length ()))) / 2.0; // AUC
   //res =  auc -  0.5 * ((double) nGenes / (double) genome.length ()); // AUC

   if (res < 0.0)
      return 0.0;
   else
      return res;

}


/*
/////////////////////////////////////////////////////////////
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
*/

// Probability function for destructive mutation
//vector<double> prob_f;

   // NOTE: for the Messy encoding
/*
extern "C" double destProb (int n) {

   if (n < 0)
      return prob_f.front ();

   if (n > (int) prob_f.size ())
      return prob_f.back ();

   return prob_f [n];

}
// Parse of probability function for destructive mutation

   double x, y;
   ifstream arch ("/home/alum/atorre/datasets/FSS/ecml07/mutationProb.dat");

   while (arch.bad ())
      std::cerr << "Error: config file not ready." << std::endl;

   while (!arch.eof ()) {

      arch >> x;
      arch >> y;

      prob_f.push_back (y);

   }

   GADualIntMessyGenome::destructiveProb = destProb;
*/

  // NOTE: ????
//   GADualIntMessyGenome::setMaxK (data.size ());

/*
/////////////////////////////////////////////////////////////
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
*/

extern "C" void individualInit (GAGenome& g) {
   return IntegerUniformInitializer (g);
}

extern "C" GAGenome* defineProblem () {

   GAEDAConfig* cfg = GAEDAConfig::handle ();

   // Parse the config file for the problem
   if (!parse_config (cfg->getProblemData ()))
     throw runtime_error ("Error: A valid configuration file must be provided for the FSS problem.");

   if (!parse_data ())
     throw runtime_error ("Error: A dataset must be provided for the FSS problem.");

   unsigned nGens = data [0].size () - 1;

   // Create the genome and the allele set
   GAAlleleSet<int> alleles;
   alleles.add (0);
   alleles.add (1);

   FSSGenome<int>* genome = new FSSGenome<int> (nGens, alleles, objective);

   // Common operators
   genome->initializer (IntegerUniformInitializer);
   genome->comparator  (IntegerElementComparator);

   // Specific stuff for GAs
   genome->crossover   (IntegerTwoPointsCrossover);
   genome->mutator     (IntegerSwapMutator);

   // Specific stuff for DE
   genome->crossover   (IntegerExponentialCrossover);

   // Specific stuff for MOS
   MOSGenomeFactory::handle()->registerGenome (GAID::IntegerEncoding, genome);

   return genome;

}

extern "C" bool postprocess (GAPopulation* pop, int rank, int repetition) {
  return true;
}

extern "C" const char* describeProblem (void) {
   return "Feature Subset Selection";
}


#include <stdexcept>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include <libconfig.h++>

#include <GAEDAConfig.h>
#include <GARealOps.h>
#include <GAIntOps.h>
#include <MOSGenomeFactory.h>
#include <MOSConversion.h>
#include <MOSConversionFunc.h>
#include <genomes/GA1DArrayGenome.h>
#include <quicksort.h>

/*
/////////////////////////////////////////////////////////////
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
*/

// Dataset file name
std::string dataset;

// Optimum value for this dataset
long double optimum;

// Number of local search iterations
unsigned lsiters = 0;

// Structure to store the problem distance matrix
typedef struct {
   long double** matrix;
   int n;
} tsp_problem_t;

// Problem dataset
tsp_problem_t problem;

/*
/////////////////////////////////////////////////////////////
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
*/

void print_matrix (void) {

   for (int i = 0; i < problem.n; i++) {

      for (int j = 0; j < problem.n; j++) {

         long double weight = problem.matrix [i][j];

         if (weight > 999)
            std::cout << weight << " ";
         else if (weight > 99)
            std::cout << " " << weight << " ";
         else if (weight > 9)
            std::cout << "  " << weight << " ";
         else
            std::cout << "   " << weight << " ";

      }

      std::cout << std::endl;

   }

}


bool parse_data (void) {

   enum EdgeFormat {FullMatrix, LowerDiagRow, UpperDiagRow, UpperRow};

   std::ifstream arch (dataset.c_str ());

   std::string name;
   std::string type;
   std::string comment;

   int dimension;
   int edgeF = FullMatrix;

   char buffer [1024];

   while (!arch.eof ()) {

      std::string token;

      arch >> token;

      if (token == "NAME:") {

         arch >> name;

      }
      else if (token == "TYPE:") {

         arch >> type;
      }
      else if (token == "COMMENT:") {

         arch.getline (buffer, 1024);

         comment = buffer;

      }
      else if (token == "DIMENSION:") {

         arch >> dimension;

         problem.n = dimension;

         problem.matrix = (long double**) calloc ((size_t) sizeof (long double) * problem.n, (size_t) sizeof (long double));

         for (int i = 0; i < problem.n; i++)
            problem.matrix [i] = (long double*) calloc ((size_t) sizeof (long double) * problem.n, (size_t) sizeof (long double));

      }
      else if (token == "EDGE_WEIGHT_TYPE:") {

         arch >> token;

         if (token != "EXPLICIT") {

            std::cerr << "Error. Edge weight type not recognized." << std::endl;
            return false;

         }

      }
      else if (token == "EDGE_WEIGHT_FORMAT:") {

         arch >> token;

         if (token == "FULL_MATRIX")
            edgeF = FullMatrix;
         else if (token == "LOWER_DIAG_ROW")
            edgeF = LowerDiagRow;
         else if (token == "UPPER_DIAG_ROW")
            edgeF = UpperDiagRow;
         else if (token == "UPPER_ROW")
            edgeF = UpperRow;
         else {

            std::cerr << "Error: Edge weight format not recognized." << std::endl;
            return false;

         }

      }
      else if (token == "DISPLAY_DATA_TYPE:") {

         arch.getline (buffer, 1024);

      }
      else if (token == "NODE_COORD_TYPE:") {

         arch.getline (buffer, 1024);

      }
      else if (token == "EDGE_WEIGHT_SECTION") {

         switch (edgeF) {

            case FullMatrix:

               for (int i = 0; i < problem.n; i++)
                  for (int j = 0; j < problem.n; j++) {

                     long double weight;

                     arch >> weight;

                     problem.matrix [i][j] = weight;

                  }

               break;

            case LowerDiagRow:

               for (int i = 0; i < problem.n; i++)
                  for (int j = 0; j < i + 1; j++) {

                     long double weight;

                     arch >> weight;

                     problem.matrix [i][j] = weight;
                     problem.matrix [j][i] = weight;

                  }

               break;

            case UpperDiagRow:

               for (int i = 0; i < problem.n; i++)
                  for (int j = i; j < problem.n; j++) {

                     long double weight;

                     arch >> weight;

                     problem.matrix [i][j] = weight;
                     problem.matrix [j][i] = weight;

                  }

               break;

            case UpperRow:

               for (int i = 0; i < problem.n; i++)
                  for (int j = i+1; j < problem.n; j++) {

                     long double weight;

                     arch >> weight;

                     problem.matrix [i][j] = weight;
                     problem.matrix [j][i] = weight;

                  }

               break;

            default:

               break;

         }

         arch >> token;

         if (token != "EOF" && token != "DISPLAY_DATA_SECTION")
            std::cerr << "=> Error: File was not succesfully parsed." << std::endl;

      }

   }

   arch.close ();

   return true;

}


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
      dataset = (const char*) cfg.lookup ("dataset");
      optimum =               cfg.lookup ("optimum");
   }
   catch (libconfig::SettingNotFoundException e) {
      std::cerr << "Error: one or more mandatory parameter(s) were not found." << std::endl;
      exit (-1);
   }

   // Then, we read optional vars (the default value is not modified if the
   // parameter is not provided in the config file)
   cfg.lookupValue ("localsearch", lsiters);

   return true;

}


extern "C" long double objectiveInt (GAGenome& g) {

   GA1DArrayAlleleGenome<int>& gen = DYN_CAST (GA1DArrayAlleleGenome<int>&, g);

   int city_orig = gen.gene (0);
   int city_dest = gen.gene (0);

   long double totaldist = problem.matrix [0][gen.gene (0)];

   for (register int i = 1; i < problem.n - 1; i++) {
      city_orig = gen.gene (i - 1);
      city_dest = gen.gene (i    );
      totaldist += problem.matrix [city_orig][city_dest];
   }

   totaldist += problem.matrix [city_dest][0];

   return optimum / totaldist;

}

extern "C" long double objectiveReal (GAGenome& g) {

   GA1DArrayAlleleGenome<long double>& gen = DYN_CAST (GA1DArrayAlleleGenome<long double>&, g);

   unsigned len = gen.length ();
   int      keys   [len];
   long double   values [len];

   // Prepare data
   for (register unsigned i = 0; i < len; i++) {
      keys   [i] = i;
      values [i] = gen.gene (i);
   }

   // Sort all the values
   quicksort (keys, values, 0, len - 1);

   int city_orig = keys [0];
   int city_dest = keys [0];

   long double totaldist = problem.matrix [0][keys[0]];

   for (register int i = 1; i < problem.n - 1; i++) {
     city_orig = keys[i - 1];
     city_dest = keys[i];
     totaldist += problem.matrix [city_orig][city_dest];
   }

   totaldist += problem.matrix [city_dest][0];

   return optimum / totaldist;

}

extern "C" void individualInit (GAGenome& g) {
   return IntegerOrderedInitializer (g);
}

extern "C" GAGenome* defineProblem () {

   GAEDAConfig* cfg = GAEDAConfig::handle ();

   // Parse the config file for the problem
   if (!parse_config (cfg->getProblemData ()))
     throw runtime_error ("Error: A valid configuration file must be provided for the TSP problem.");

   if (!parse_data ())
     throw runtime_error ("Error: A dataset must be provided for the TSP problem.");

   // Definition of the genomes and the allele sets
   GAAlleleSet<int>    int_alleles;
   GAAlleleSet<long double> real_alleles (0, 1);

   for (int i = 1; i < problem.n; i++)
      int_alleles.add (i);

   GA1DArrayAlleleGenome<long double>* genomeReal = new GA1DArrayAlleleGenome<long double> (problem.n - 1, real_alleles, objectiveReal);
   GA1DArrayAlleleGenome<int>*    genomeInt  = new GA1DArrayAlleleGenome<int>    (problem.n - 1, int_alleles,  objectiveInt );

   // Common operators
   genomeInt->initializer (IntegerOrderedInitializer);
   genomeInt->comparator  (IntegerElementComparator);

   // Specific stuff for GAs
   genomeInt->crossover   (IntegerAlternativeOrderCrossover);
   genomeInt->mutator     (IntegerRepeatedExchangeMutator);

   // Specific stuff for DE
   genomeInt->crossover   (IntegerExponentialCrossover);

   // Specific stuff for MOS
   MOSGenomeFactory::handle ()->registerGenome (GAID::RealEncoding,    genomeReal);
   MOSGenomeFactory::handle ()->registerGenome (GAID::IntegerEncoding, genomeInt );

   MOSConversion::handle ()->registerConvFunction (GAID::IntegerEncoding, GAID::RealEncoding,    convertIntToReal);
   MOSConversion::handle ()->registerConvFunction (GAID::RealEncoding,    GAID::IntegerEncoding, convertRealToInt);

   // No need for this genome anymore
   delete genomeReal;

   return genomeInt;

}

extern "C" bool postprocess (GAPopulation* pop, int rank, int repetition) {
  for (int i = 0; i < problem.n; i++)
    free (problem.matrix[i]);
  free (problem.matrix);
  return true;
}

extern "C" const char *describeProblem (void) {
   return "TSP";
}

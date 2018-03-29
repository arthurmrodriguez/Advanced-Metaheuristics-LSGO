#include <iostream>

#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/toplev.h> /* do_octave_atexit */

int main (const int argc, char** argv) {
   const char* argvv[] = {"" /* name of program, not relevant */, "--silent"};

   octave_main(2, (char**) argvv, true /* embedded */);

   // Initialize path of Matlab code
   octave_value_list pathArguments;
   pathArguments(0) = argv[1];
   const octave_value_list resultPath = feval("addpath", pathArguments, 1);

   // First, initialize Benchmark Global Functions
   octave_value_list initializeArguments;
   initializeArguments(0) = argv[1];
   const octave_value_list result = feval("initializeBenchmark", initializeArguments, 1);

   // Now, call the evaluation function
   octave_value_list evalArguments;
   Matrix inMatrix (1, 1000);

   for (unsigned i = 0; i < 1000; i++)
      inMatrix(0, i) = 1;

   evalArguments(0) = inMatrix;
   evalArguments(1) = atoi(argv[2]);

   const octave_value_list fitness = feval("benchmark_func", evalArguments, 1);

   std::cout << "fitness is " << fitness(0).scalar_value() << std::endl;

   do_octave_atexit ();
}

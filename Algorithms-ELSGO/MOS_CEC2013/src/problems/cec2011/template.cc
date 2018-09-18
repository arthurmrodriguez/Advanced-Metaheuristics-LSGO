#ifndef CEC2011_H_
#define CEC2011_H_

#include <stdexcept>
#include <limits>

#include <GARealOps.h>
#include <MOSGenomeFactory.h>
#include <GAEDAConfig.h>
#include <GAGeneticAlgorithm.h>
#include <genomes/GA1DArrayGenome.h>

#include <mCEC_Function.h>

double GD_Max[6][21] = { {0.2170, 0.0240, 0.0760, 0.8920, 0.1280, 0.2500, 0.0580, 0.1120, 0.0620, 0.0820, 0.0350, 0.0900, 0.0320, 0.0950, 0.0220, 0.1750, 0.0320, 0.0870, 0.0350, 0.0240, 0.1060},
                         {0.2170, 0.0240, 0.0260, 0.4910, 0.2280, 0.3000, 0.0580, 0.1120, 0.0620, 0.0820, 0.0350, 0.0900, 0.0320, 0.0950, 0.0220, 0.1750, 0.0320, 0.0870, 0.0350, 0.0240, 0.1060},
                         {0.2160, 0.0240, 0.0760, 0.2160, 0.2160, 0.2160, 0.0580, 0.1120, 0.0620, 0.0820, 0.0350, 0.0900, 0.0320, 0.0950, 0.0220, 0.1750, 0.0320, 0.0870, 0.0350, 0.0240, 0.0810},
                         {0.2170, 0.0240, 0.0760, 0.2280, 0.2280, 0.2280, 0.0580, 0.1120, 0.0620, 0.0820, 0.0350, 0.0900, 0.0320, 0.0950, 0.0220, 0.0250, 0.0320, 0.0870, 0.0350, 0.0240, 0.0810},
                         {0.1240, 0.0240, 0.0760, 0.1240, 0.1240, 0.1240, 0.0580, 0.1120, 0.0620, 0.0820, 0.0350, 0.0650, 0.0320, 0.0950, 0.0220, 0.1240, 0.0320, 0.0870, 0.0350, 0.0240, 0.1060},
                         {0.1160, 0.0240, 0.0760, 0.1160, 0.1160, 0.1160, 0.0580, 0.0870, 0.0620, 0.0820, 0.0350, 0.0900, 0.0320, 0.0950, 0.0220, 0.1160, 0.0320, 0.0870, 0.0350, 0.0240, 0.1060}
                       };

double ELD6  [  6][2] = { {100, 500}, {50, 200}, {80, 300}, {50, 150}, { 50, 200}, {50, 120} };
double ELD13 [ 13][2] = { {0, 680}, {0, 360}, {0, 360}, {60, 180}, {60, 180}, {60, 180}, {60, 180}, { 60, 180}, {60, 180}, {40, 120}, {40, 120}, {55, 120}, {55, 120} };
double ELD15 [ 15][2] = { {150, 455}, {150, 455}, {20, 130}, {20, 130}, {150, 470}, { 135, 460}, {135, 465}, {60, 300}, {25, 162}, {25, 160}, { 20, 80}, {20, 80}, {25, 85}, {15, 55}, {15, 55} };
double ELD40 [ 40][2] = { {36, 114}, {36, 114}, {60, 120}, {80, 190}, {47, 97}, {68, 140}, {110, 300}, {135, 300}, {135, 300}, {130, 300}, {94, 375}, {94, 375}, {125, 500}, {125, 500}, { 125, 500}, {125, 500}, {220, 500}, {220, 500}, {242, 550}, {242, 550}, {254, 550}, {254, 550}, {254, 550}, {254, 550}, {254, 550}, {254, 550}, {10, 150}, {10, 150 }, {10, 150}, {47, 97}, {60, 190}, {60, 190}, {60, 190}, {90, 200}, {90, 200}, {90, 200}, { 25, 110}, { 25, 110}, {25, 110}, {242, 550} };
double ELD140[140][2] = { {71, 119}, {120, 189}, {125, 190}, {125, 190}, {90, 190}, {90, 190}, {280, 490}, { 280, 490}, {260, 496}, {260, 496}, {260, 496}, {260, 496}, {260, 506}, {260, 509}, {260, 506}, {260, 505}, {260, 506}, {260, 506}, {260, 505}, {260, 505}, { 260, 505}, {60, 505}, {260, 505}, {260, 505}, {280, 537}, {280, 537}, {280, 549}, { 280, 549}, {260, 501}, {260, 501}, {260, 506}, {260, 506}, {260, 506}, { 260, 506}, {260, 500}, {260, 500}, {120, 241}, {120, 241}, {423, 774}, { 423, 769}, {3, 19}, {3, 28}, {160, 250}, {160, 250}, {160, 250}, {160, 250}, { 160, 250}, {160, 250}, {160, 250}, {160, 250}, {165, 504}, {165, 504}, { 165, 504}, {165, 504}, {180, 471}, {180, 561}, {103, 341}, {198, 617}, { 100, 312}, {153, 471}, {163, 500}, {95, 302}, {160, 511}, {160, 511}, { 196, 490}, {196, 490}, {196, 490}, {196, 490}, {130, 432}, {130, 432}, { 137, 455}, {137, 455}, {195, 541}, {175, 536}, {175, 540}, {175, 538}, {175, 540}, {330, 574}, {160, 531}, {160, 531}, {200, 542}, {56, 132}, { 115, 245}, {115, 245}, {115, 245}, {207, 307}, {207, 307}, {175, 345}, { 175, 345}, {175, 345}, {175, 345}, {360, 580}, {415, 645}, {795, 984}, { 795, 978}, {578, 682}, {615, 720}, {612, 718}, {612, 720}, {758, 964}, { 755, 958}, {750, 1007}, {750, 1006}, {713, 1013}, {718, 1020}, {791, 954}, { 786, 952}, {795, 1006}, {795, 1013}, {795, 1021}, {795, 1015}, {94, 203}, { 94, 203}, {94, 203}, {244, 379}, {244, 379}, {244, 379}, {95, 190}, {95, 189}, { 116, 194}, {175, 321}, {2, 19}, {4, 59}, {15, 83}, {9, 53}, {12, 37}, {10, 34}, { 112, 373}, {4, 20}, {5, 38}, {5, 19}, {50, 98}, {5, 10}, {42, 74}, {42, 74}, { 41, 105}, {17, 51}, {7, 19}, {7, 19}, {26, 40} };

GAAlleleSetArray<double>* alleles = 0;

std::string pname;

FIELD_TYPE sol[1000];
FIELD_TYPE obj_val[10];

void initializeAlleles (unsigned func);

extern "C" double objective (GAGenome& g) {
   GA1DArrayAlleleGenome<double>& gen = dynamic_cast<GA1DArrayAlleleGenome<double>&> (g);

   for (unsigned i = 0; i < gen.size(); i++)
     sol[i] = gen.gene(i);

   cost_functionXXX(sol, obj_val);

   if (isnan(obj_val[0]))
     return numeric_limits<double>::max();
   else
     return obj_val[0];
}

extern "C" void RealUniformInitializerF4 (GAGenome& g) {

   GA1DArrayAlleleGenome<double>& gen = DYN_CAST (GA1DArrayAlleleGenome<double>&, g);

   gen.resize (GAGenome::ANY_SIZE);

   for (int i = 0; i < gen.length (); i++)
      gen.gene (i, GARandomDouble (0, 5));

}

extern "C" void individualInit (GAGenome& g) {
   if (YYY == 4)
      return RealUniformInitializerF4 (g);
   else
      return RealUniformInitializer (g);
}

extern "C" GAGenome* defineProblem (int size, char *data) {

   if (Initial_CEC2011_Cost_Function() != _NO_ERROR)
      throw runtime_error("Error: could not initialize Matlab wrapper.");

   GAEDAConfig* cfg = GAEDAConfig::handle();

   initializeAlleles(YYY);
   GA1DArrayAlleleGenome<double>* genome = new GA1DArrayAlleleGenome<double> (*alleles, objective);

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

extern "C" bool postprocess (GAPopulation* pop, int rank, int repetition) {
   delete alleles;
   Terminate_CEC2011_Cost_Function();
   return true;
}

extern "C" const char *describeProblem () {
   return pname.c_str();
}

extern "C" GAGenome::OptCriterion optCriterion(){
  return GAGenome::MINIMIZATION;
}

void initializeAlleles (unsigned func) {
  alleles = new GAAlleleSetArray<double>();

  switch (func) {
    case 1:
      for (unsigned i = 0; i < 6; i++)
        alleles->add(-6.40, 6.35);
      break;

    case 2:
    case 5:
    case 6:
      alleles->add(0, 4 );
      alleles->add(0, 4 );
      alleles->add(0, M_PI);
      for (unsigned i = 4; i <= 30; i++) {
        double fl = floor(((double)(i - 4.0)) / 3.0);
        double lb = -4.0 - ((1.0 / 4.0) * fl);
        double ub =  4.0 + ((1.0 / 4.0) * fl);
        alleles->add(lb, ub);
      }
      break;

    case 3:
      alleles->add(0.6, 0.9);
      break;

    case 4:
      alleles->add(-1000, 1000);
      break;

    case 7:
      for (unsigned i = 0; i < 20; i++)
        alleles->add(0, 2*M_PI);
      break;

    case 8:
      for (unsigned i = 0; i < 7; i++)
        alleles->add(0, 15, GAAllele::EXCLUSIVE, GAAllele::INCLUSIVE);
      break;

    case 9:
      for (unsigned i = 0; i < 6; i++)
        for (unsigned j = 0; j < 21; j++)
          alleles->add(0, GD_Max[i][j]);
      break;

    case 10:
      for (unsigned i = 0; i < 6; i++)
        alleles->add(0.2, 1);
      for (unsigned i = 0; i < 6; i++)
        alleles->add(-180, 180);
      break;

    case 115:
      for (unsigned i = 0; i < 24; i++) {
        alleles->add(10,  75);
        alleles->add(20, 125);
        alleles->add(30, 175);
        alleles->add(40, 250);
        alleles->add(50, 300);
      }
      break;

    case 1110:
      for (unsigned i = 0; i < 24; i++)
        alleles->add( 55,  55);
      for (unsigned i = 0; i < 24; i++) {
        alleles->add(150, 470);
        alleles->add(135, 460);
        alleles->add( 73, 340);
        alleles->add( 60, 300);
        alleles->add( 73, 243);
        alleles->add( 57, 160);
        alleles->add( 20, 130);
        alleles->add( 47, 120);
        alleles->add( 20,  80);
      }
      break;

    case 126:
      for (unsigned i = 0; i < 6; i++)
        alleles->add(ELD6[i][0], ELD6[i][1]);
      break;

    case 1213:
      for (unsigned i = 0; i < 13; i++)
        alleles->add(ELD13[i][0], ELD13[i][1]);
      break;

    case 1215:
      for (unsigned i = 0; i < 15; i++)
        alleles->add(ELD15[i][0], ELD15[i][1]);
      break;

    case 1240:
      for (unsigned i = 0; i < 40; i++)
        alleles->add(ELD40[i][0], ELD40[i][1]);
      break;

    case 12140:
      for (unsigned i = 0; i < 140; i++)
        alleles->add(ELD140[i][0], ELD140[i][1]);
      break;

    case 131:
    case 132:
    case 133:
      for (unsigned i = 0; i < 24; i++) {
        alleles->add( 5, 15);
        alleles->add( 6, 15);
        alleles->add(10, 30);
        alleles->add(13, 25);
      }
      break;

    case 14:
      alleles->add( 1900.0 , 2300.0); // Range for dimension 1 MJD2000 units
      alleles->add( 2.5 , 4.05); // Range for dimension 2 Km/sec
      alleles->add( 0.0 , 1.0); // Range for dimension 3 n/a
      alleles->add( 0.0 , 1.0); // Range for dimension 4 n/a
      alleles->add( 100.0 , 500.0); // Range for dimension 5 days
      alleles->add( 100.0 , 500.0); // Range for dimension 6 days
      alleles->add( 100.0 , 500.0); // Range for dimension 7 days
      alleles->add( 100.0 , 500.0); // Range for dimension 8 days
      alleles->add( 100.0 , 500.0); // Range for dimension 9 days
      alleles->add( 100.0 , 600.0); // Range for dimension 10 days
      alleles->add( 0.01, 0.99); // Range for dimension 11 n/a
      alleles->add( 0.01, 0.99); // Range for dimension 12 n/a
      alleles->add( 0.01, 0.99); // Range for dimension 13 n/a
      alleles->add( 0.01, 0.99); // Range for dimension 14 n/a
      alleles->add( 0.01, 0.99); // Range for dimension 15 n/a
      alleles->add( 0.01, 0.99); // Range for dimension 16 n/a
      alleles->add( 1.1, 6.0); // Range for dimension 17 n/a
      alleles->add( 1.1, 6.0); // Range for dimension 18 n/a
      alleles->add( 1.05, 6.0); // Range for dimension 19 n/a
      alleles->add( 1.05, 6.0); // Range for dimension 20 n/a
      alleles->add( 1.05, 6.0); // Range for dimension 21 n/a
      alleles->add( -M_PI , M_PI); // Range for dimension 22 n/a
      alleles->add( -M_PI , M_PI); // Range for dimension 23 n/a
      alleles->add( -M_PI , M_PI); // Range for dimension 24 n/a
      alleles->add( -M_PI , M_PI); // Range for dimension 25 n/a
      alleles->add( -M_PI , M_PI); // Range for dimension 26 n/a
      break;

    case 15:
      alleles->add(-1000.0 , 0.0); // Range for dimension 1 MJD2000 units
      alleles->add( 3.0 , 5.0); // Range for dimension 2 Km/sec
      alleles->add( 0.0 , 1.0); // Range for dimension 3 n/a
      alleles->add( 0.0 , 1.0); // Range for dimension 4 n/a
      alleles->add( 100.0 , 400.0); // Range for dimension 5 days
      alleles->add( 100.0 , 500.0); // Range for dimension 6 days
      alleles->add( 30.0 , 300.0); // Range for dimension 7 days
      alleles->add( 400.0 , 1600.0); // Range for dimension 8 days
      alleles->add( 800.0 , 2200.0); // Range for dimension 9 days
      alleles->add( 0.01, 0.9); // Range for dimension 10 n/a
      alleles->add( 0.01, 0.9); // Range for dimension 11 n/a
      alleles->add( 0.01, 0.9); // Range for dimension 12 n/a
      alleles->add( 0.01, 0.9); // Range for dimension 13 n/a
      alleles->add( 0.01, 0.9); // Range for dimension 14 n/a
      alleles->add( 1.05, 6.0); // Range for dimension 15 n/a
      alleles->add( 1.05, 6.0); // Range for dimension 16 n/a
      alleles->add( 1.15, 6.5); // Range for dimension 17 n/a
      alleles->add( 1.7 , 291.0); // Range for dimension 18 n/a
      alleles->add( -M_PI , M_PI); // Range for dimension 19 rad
      alleles->add( -M_PI , M_PI); // Range for dimension 20 rad
      alleles->add( -M_PI , M_PI); // Range for dimension 21 rad
      alleles->add( -M_PI , M_PI); // Range for dimension 22 rad
      break;

    default:
      throw runtime_error("Error: alleles not defined for this function");
      break;
  }

  return;
}

#endif /* CEC2011_H_ */

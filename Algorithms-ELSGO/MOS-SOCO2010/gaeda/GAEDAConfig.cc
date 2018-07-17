/*****************************************************************************
 * GAEDAlib: A C++ GA library with EDA and multiprocessor (MPI) support      *
 *                                                                           *
 * (C) 2005 Pedro Diaz (pdiaz@laurel.datsi.fi.upm.es)                        *
 *                                                                           *
 * GAEDAlib is distributed under the terms of the BSD software license       *
 *                                                                           *
 * GAEDAlib is heavily based on GAlib, a C++ GA library by Mathew Wall:      *
 * Copyright (c) 1995-1996 Massachusetts Institute of Technology (MIT)       *
 * Copyright (c) 1996-2000 Matthew Wall (author of GAlib)                    *
 *                                                                           *
 * Some portions of GAEDAlib's source code come from the GNU C++ compiler    *
 * library and therefore are covered under the terms of a different license, *
 * the GNU Public License.                                                   *
 *                                                                           *
 * You should have received a file named LICENSE along with this software.   *
 * This file contains more information about the licensing conditions of     *
 * GAEDAlib as well as the full text of each license involved.               *
 *                                                                           *
 * The file AUTHORS lists the people who have contributed (directly or       *
 * indirectly) to GAEDAlib                                                   *
 *****************************************************************************/

#include <stdexcept>
#include <assert.h>

#include <boost/program_options.hpp>

#include "GAEDAConfig.h"

#include "GARealOps.h"
#include "MOSQualityFunction.h"
#include "islands/CommManager.h"

namespace po = boost::program_options;

/* CLASS IMPLEMENTATION */

po::options_description desc("Available options");

GAEDAConfig* GAEDAConfig::ptrSelf = NULL;


GAEDAConfig* GAEDAConfig::handle(CommManager& comm_manager) {
  if (ptrSelf) {
    std::cerr << "[GAEDAConfig] Warning: trying to construct again the GAEDAConfig Singleton. Nothing is done..." << std::endl;
    return ptrSelf;
  }
  else {
    ptrSelf = new GAEDAConfig (comm_manager);
    return ptrSelf;
  }
}


GAEDAConfig* GAEDAConfig::handle() {
  if (!GAEDAConfig::ptrSelf) {
    std::cerr << "[GAEDAConfig] Fatal error: the GAEDAConfig Singleton object has not been previously created. Aborting..." << std::endl;
   exit(-1);
  }
  return GAEDAConfig::ptrSelf;
}


GAEDAConfig::GAEDAConfig(CommManager& comm_manager) : comm_manager_(comm_manager) {
   // Set default values
   mIters = 100;
   mEvals = 5000;
   mConvCrit = EVALS;
   mPopSize = 100;

   mProblemStruct = NULL;
   mLogFile       = NULL;
   mProblemData   = NULL;
}


// Destructor
GAEDAConfig::~GAEDAConfig() {
}

void GAEDAConfig::destroy() {
   if (GAEDAConfig::ptrSelf)
      delete GAEDAConfig::ptrSelf;
   GAEDAConfig::ptrSelf = NULL;
}


// Parse command line arguments
bool GAEDAConfig::parse(int argc, char **argv, bool embedded) {

   int tmp, size;
   string cnvGens_msg, mutProb_msg, crossProb_msg;
   stringstream converter;

   try {

     po::options_description general_opts("General options:");
     general_opts.add_options()

       // Help message
       ("help,h", "Show this help message.")
       ("quiet,q", "Turns output off (to embed the algorithm in other software.")

       // Common options (for all evolutionary algorithms and models)
       ("cfg", po::value<string>()->default_value("./config.cfg"), "Use the config file. The command line options have higher priority than the config file. "
                                                                   "To use it you must put the large option name of the options.")

       ("log-file,T", po::value< string >(), "Turns logging to the provided file on. The given filename will be appended with the rank number (identifier for this island).")

       ("log-level,L", po::value< string >()->default_value("only_stats"), "Log verbosity (none, only_stats, normal, debug_lite, debug).")

       ("log-option,b", po::value< vector<string> >()->composing(), "Set the statistics that should be traced among: score, fitness, participation, quality, etc.")

       ("stats-frequency,S", po::value<int>()->default_value(1), "Set the frequency to display the statistics.")

       ("scaling", po::value<string>()->default_value("linear"), "Set the scaling scheme to compute fitness from score (none, standard, linear, sigma_truncation, power_law).")

       ("algorithm,a", po::value< string >()->default_value("mos"), "Set the algorithm to be used (mos).")

       ("problem-size,s", po::value<int>(), "Set the size of the problem (not honored by all problems).")

       ("additional-data,A", po::value< string >(), "Additional data to be sent to the problem module.")

       ;

     po::options_description convergence_opts("Convergence options:");
     convergence_opts.add_options()

       ("iterations,i", po::value<int>(), "Set maximum number of generations. In combination with 'convergence' sets a mixed convergence criterion of "
                                          "maximum number of iterations and fitness convergence.")

       ("convergence,C", po::value<long double>(), "Set the fitness convergence ratio. Use 'convergence-iters' to set the number of generations considered to compute the convergence"
                                              "ratio. In combination with 'iterations' sets a mixed convergence criterion of fitness convergence and maximum number of iterations.")

       ("convergence-iters", po::value<int>()->default_value(20), "Set the number of iterations used to calculate finess convergence.")

       ("pop-convergence,U", po::value<long double>(), "Set the population convergence ratio.")

       ("partial-pop-convergence,V", po::value<long double>(), "Set a partial population convergence criterion.")

       ("pop-or-gen-convergence,y",po::value< string >(),"Set a population or generation convergence criterion.")

       ("evaluations", po::value<unsigned long>(), "Set maximum number of fitness evaluations.")

       ;

     po::options_description parallel_opts("Parallel model options:");
     parallel_opts.add_options()

       ("pop-island,p", po::value<int>(), "Set the population size per island.")

       ("pop-total,P", po::value<int>(), "Set the population size shared among all islands.")

       ;

     po::options_description ga_opts("Genetic Algorithms options:");
     ga_opts.add_options()

       ("mut-prob,M", po::value<long double>()->default_value(0.01), "Set mutation probability.")

       ("cross-prob,X", po::value<long double>()->default_value(0.9), "Set crossover probability.")

       ("selection-scheme,R", po::value<string>()->default_value("tournament:2"), "Selection scheme (tournament:n,roulette,random,rank,stochastic,deterministic).")

       ;

     po::options_description mos_opts("MOS options:");
     mos_opts.add_options()

       ("techs-repository", po::value<string>(), "Set the repository with all the MOS techniques that can be used.")

       ("use-tech", po::value< vector<string> >()->composing(), "Set the techniques to use in MOS algorithm.")

       ("minimum-part", po::value<long double>()->default_value(0.05) , "Set the minimum participation of a MOS technique.")

       ("mos-genome-printing", po::value<string>()->default_value("full"), "Sets the way that MOS genomes should be printed ('default' for default encoding and 'full' for every encoding).")

       ("quality-measure", po::value<string>()->default_value("fAvg"), "Set type of quality measure used for dynamic adjustment of participation (fAvg, fitIncAvg, dAvg, Compass).")

       ("mos-bonus", po::value<long double>()->default_value(0.01), "Sets the reward bonus parameter.")

       ("mos-shared-evals-num", po::value<int>()->default_value(1500), "Sets the number of shared evaluations to distribute in a MOS HRH algorithm.")

       ("mos-shared-evals-percent", po::value<long double>()->default_value(0.01), "Sets the percentage of shared evaluations to distribute in a MOS HRH algorithm.")

       ;

     po::options_description debug_opts("Debug options:");
     debug_opts.add_options()

       ("random-seed,q", po::value<unsigned>()->default_value(0), "Seed for random number generation subsystem. 0 for random seed.")

       ("debug-flag,Z", po::value<int>()->default_value(0), "Set a sleep point in the main program to allow easy GDB attaching and debugging.")

       ;

     // Join all options groups together
     desc.add(general_opts).add(convergence_opts).add(parallel_opts).add(ga_opts).add(mos_opts).add(debug_opts);

     po::options_description hidden("Hidden Options");
     hidden.add_options()
       ("problem_module", po::value< string >(), "")
       ;

     // For options that don't have any values
     po::positional_options_description problem_module;
     problem_module.add("problem_module", -1);

     // To get cmdline and config file options
     po::options_description cmdline_options;
     cmdline_options.add(desc).add(hidden);

     po::options_description config_file_options;
     config_file_options.add(desc).add(hidden);

     // Variables map containing all the parsed options
     po::variables_map vm;
     store(po::command_line_parser(argc, argv).
           options(cmdline_options).positional(problem_module).run(), vm);

     // Beginning of the options parsing

     /*
      * General Options
      */

     // Read configuration file (if any)
     if (vm.count("cfg")) {
       ifstream ifs(vm["cfg"].as<string>().c_str());
       store(parse_config_file(ifs, config_file_options), vm);
     }

     notify(vm);


     // Show help message
     if (vm.count("help")) {
       showUsage (argv [0]);
       exit(0);
     }


     if (vm.count("quiet"))
       mQuiet = true;
     else
       mQuiet = false;


     if (vm.count("log-file")) {
       mLogFile = strdup((char*)vm["log-file"].as< string >().c_str());
     }


     if (vm.count("log-level")) {
       if (vm["log-level"].as< string >() == "none")
         mLogLevel = GALogger::none;
       else if (vm["log-level"].as< string >() == "only_stats")
         mLogLevel = GALogger::only_stats;
       else if (vm["log-level"].as< string >() == "normal")
         mLogLevel = GALogger::normal;
         else if (vm["log-level"].as< string >() == "debug_lite")
         mLogLevel = GALogger::debug_lite;
       else if (vm["log-level"].as< string >() == "debug")
         mLogLevel = GALogger::debug;
       else
         std::cerr << "[GAEDAConfig] Warning: wrong log level selected. Using default value (only_stats)." << std::endl;
     }


     if (vm.count("log-option")) {
       for (unsigned int i=0; i < vm["log-option"].as< vector<string> >().size(); i++){
         if (vm["log-option"].as< vector<string> >()[i] == "score")
           logstats_names.push_back(Score);
         if (vm["log-option"].as< vector<string> >()[i] == "fitness")
           logstats_names.push_back(Fitness);
         if (vm["log-option"].as< vector<string> >()[i] == "fit_inc")
           logstats_names.push_back(FitInc);
         if (vm["log-option"].as< vector<string> >()[i] == "participation")
           logstats_names.push_back(Participation);
         if (vm["log-option"].as< vector<string> >()[i] == "quality")
           logstats_names.push_back(Quality);
         if (vm["log-option"].as< vector<string> >()[i] == "quality_func_active")
           logstats_names.push_back(QualityFuncActive);
         if (vm["log-option"].as< vector<string> >()[i] == "evaluations")
           logstats_names.push_back(Evals);
         if (vm["log-option"].as< vector<string> >()[i] == "improvements")
           logstats_names.push_back(Improve);
       }
     }


     if (vm.count("stats-frequency")) {
       mStatsDispFreq = vm["stats-frequency"].as<int>();

       if (mStatsDispFreq <= 0){
         std::cerr << "[GAEDAConfig] Error: the stats display frequency specified: " << mStatsDispFreq << " is <= 0." << std::endl;
         return false;
       }
     }


     if (vm.count("scaling")) {
       if (vm["scaling"].as<string>() == "none")
         scaling_ = None;
       else if (vm["scaling"].as<string>() == "standard")
         scaling_ = Standard;
       else if (vm["scaling"].as<string>() == "linear")
         scaling_ = Linear;
       else if (vm["scaling"].as<string>() == "sigma_truncation")
         scaling_ = Sigma_Truncation;
       else if (vm["scaling"].as<string>() == "power_law")
         scaling_ = Power_Law;
       else {
         cerr << endl << "[GAEDAConfig] Error: unrecognized scaling method." << endl;
         return false;
       }
     }


     if (vm.count("problem-size")) {
       mProbSize = vm["problem-size"].as<int>();

       if (mProbSize <= 0) {
         std::cerr << "[GAEDAConfig] Error: Problem size must be greater than 0." << std::endl;
         return false;
       }
     }


     if (vm.count("additional-data")) {
       mProblemData = strdup((char*)vm["additional-data"].as< string >().c_str());
     }


     /*
      * Convergence Options
      */

     if (vm.count("iterations")) {
       mConvCrit = ITERS;
       mIters = vm["iterations"].as<int>();

       if (mIters <= 0) {
         std::cerr << "[GAEDAConfig] Error: Wrong number of iterations." << std::endl;
         return false;
       }
     }


     if (vm.count("convergence")) {
       // If both 'iterations' and 'convergence' are defined, we use the mixed criterion of convergence
       if (vm.count("iterations"))
         mConvCrit =  SCOREORGEN;
       else
         mConvCrit = CONV;

       mConv = vm["convergence"].as<long double>();

       if (mConv < 0 || mConv > 1){
         std::cerr << "[GAEDAConfig] Error: convergence must be within [0, 1] interval." << std::endl;
         return false;
       }
     }


     if (vm.count("convergence-iters")) {
       mConvIters = vm["convergence-iters"].as<int>();

       if (mConvIters <= 0) {
         std::cerr << "[GAEDAConfig] Error: convergence-iters must be greater than 0." << std::endl;
         return false;
       }
     }


     if (vm.count("pop-convergence")) {
       mConvCrit = POPCONV;
       mConv = vm["pop-convergence"].as<long double>();

       if (mConv < 0 || mConv > 1) {
         std::cerr << "[GAEDAConfig] Error: Population convergence ratio  must be within [0, 1] interval." << std::endl;
         return false;
       }
     }


     if (vm.count("partial-pop-convergence")) {
       mConvCrit = PARTIALPOPCONV;
       mConv = vm["partial-pop-convergence"].as<long double>();

       if (mConv < 0 || mConv > 1) {
         std::cerr << "[GAEDAConfig] Error: The partial pop convergence ratio must be within [0, 1] interval. " << std::endl;
         return false;
       }
     }


     if (vm.count("pop-or-gen-convergence")) {
       if (sscanf(vm["y"].as<string>().c_str(),"%lf:%d",&mConv,&mIters) ) {
         mConvCrit = POPCONVITERS;
       }
       else {
         cerr << "[GAEDAConfig] Error: Options for convergence upon population or ngens not recognised." << endl;
         return false;
       }
     }


     // The algorithm should be processed after convergence options to allow
     // the error checking in the case of the MOS algorithm
     if (vm.count("algorithm")) {
       if (vm["algorithm"].as< string >() == "mos") {
         mAlg = GAEDAConfig::MOS;
         if (mConvCrit != EVALS) {
           if ( comm_manager_.isIslandMaster() )
             std::cerr << "[GAEDAConfig] Error: You can not use the 'mos' algorithm with a convergence criterion different to 'evaluations'." << std::endl;
           return false;
         }
       }
       else{
         cerr << "[GAEDAConfig] Error: unrecognised algorithm." << endl;
         return false;
       }
     }


     /*
      * Parallel Model Options
      */

     if (vm.count("pop-island")) {
       mPopSize = vm["pop-island"].as<int>();

       if (mPopSize <= 0) {
         std::cerr << "[GAEDAConfig] Error: Population size must be greater than 0." << std::endl;
         return false;
       }
     }


     if (vm.count("pop-total")) {
       tmp = vm["pop-total"].as<int>();
       size = comm_manager_.getNumIslands();

       if (tmp <= 0 || tmp < size ){
         std::cerr << "[GAEDAConfig] Error: Population size must be greater than 0 and the number of islands being used." << std::endl;
         return false;
       }

       mPopSize = tmp / size;
     }


     if (vm.count("evaluations")) {
       mConvCrit = EVALS;
       mEvals = vm["evaluations"].as<unsigned long>();

       if (mEvals < mPopSize) {
         std::cerr << "[GAEDAConfig] Error: Wrong number of fitness evaluations: " << mEvals << std::endl;
         return false;
       }
     }


     /*
      * Genetic Algorithms Options
      */

     if (vm.count("mut-prob")) {
       mMutProb = vm["mut-prob"].as<long double>();

       if (mMutProb < 0 || mMutProb > 1) {
         std::cerr << "[GAEDAConfig] Error: mutation probability out of range [0, 1]. " << std::endl;
         return false;
       }
     }


     if (vm.count("cross-prob")) {
       mCrossProb = vm["cross-prob"].as< long double>();

       if (mCrossProb < 0 || mCrossProb > 1) {
         cerr << "[GAEDAConfig] Error: crossover probability out of range [0, 1]." << endl;
         return false;
       }
     }


     if (vm.count("selection-scheme")) {
       if (sscanf(vm["selection-scheme"].as<string>().c_str(), "tournament:%d",&selector_degree_) ) {
         assert(selector_degree_>1);
         selector_ = Tournament;
       }
       else if (sscanf(vm["selection-scheme"].as<string>().c_str(), "unif_tournament:%d",&selector_degree_) ) {
         assert(selector_degree_>1);
         selector_ = Unif_Tournament;
       }
       else if (vm["selection-scheme"].as<string>() == "roulette"){
         selector_ = Roulette;
       }
       else if (vm["selection-scheme"].as<string>() == "random"){
         selector_ = RandomSelec;
           }
       else if (vm["selection-scheme"].as<string>() == "rank"){
         selector_ = Rank;
       }
       else if (vm["selection-scheme"].as<string>() == "stochastic"){
         selector_ = Stochastic;
       }
       else if (vm["selection-scheme"].as<string>() == "deterministic"){
         selector_ = Deterministic;
       }
       else {
         std::cerr << "[GAEDAConfig] Error: unrecognized selection method." << std::endl;
         return false;
       }
     }


     /*
      * MOS Options
      */

     if (vm.count("techs-repository")){
       techsRepository = strdup( (char*) vm["techs-repository"].as< string >().c_str() );
     }
     else if (mAlg == GAEDAConfig::MOS) {
       std::cerr << "[GAEDAConfig] Error: you must specify a repository of techniques." << std::endl;
       return false;
     }


     if (vm.count("use-tech")) {
       techsToUse = vm["use-tech"].as< vector<std::string> >();
     }


     if (vm.count("minimum-part")) {
       if (vm["minimum-part"].as<long double>() < 0 || vm["minimum-part"].as<long double>() > 1 ){
         std::cerr << "[GAEDAConfig] Error: Wrong minimum participation for MOS techniques: " << vm["minimum-part"].as<long double>() << ". Please select a value within [0, 1]. " << std::endl;
         return false;
       }
       minPart = vm["minimum-part"].as<long double>();
     }


     if (vm.count("mos-genome-printing")) {
       if (vm["mos-genome-printing"].as<string>() == "full")
         mosGenomePrinting = FullPrinting;
       else if (vm["mos-genome-printing"].as<string>() == "default")
         mosGenomePrinting = DefaultPrinting;
       else  {
         std::cerr << "[GAEDAConfig] Error: MOS Genome printing type not recognised." << std::endl;
         return false;
       }
     }


     if (vm.count("quality-measure")) {
       string quality_measure = vm["quality-measure"].as<string>();
       if ( quality_measure == "fAvg"      ||
            quality_measure == "fitIncAvg" ||
            quality_measure == "dAvg"      ||
            quality_measure == "Compass"      ) {
         mQualityMeasure = quality_measure;
       }
       else{
         std::cerr << "[GAEDAConfig] Error: Wrong quality measure selected." << std::endl;
         return false;
       }
     }


     if (vm.count("mos-bonus")) {
       if (vm["mos-bonus"].as<long double>() < 0 || vm["mos-bonus"].as<long double>() > 1 ){
         std::cerr << "[GAEDAConfig] Error: Wrong reward bonus for MOS: " << vm["mos-bonus"].as<long double>() << ". Please select a value within [0, 1]. " << std::endl;
         return false;
       }
       mBonus = vm["mos-bonus"].as<long double>();
     }


     if (vm.count("mos-shared-evals-num")) {
       if (vm["mos-shared-evals-num"].as<int>() < 1){
         std::cerr << "[GAEDAConfig] Error: Wrong number of shared evaluations for MOS HRH: " << vm["mos-shared-evals-num"].as<int>() << ". Please select a value within [1, +inf). " << std::endl;
         return false;
       }
       mSharedEvals = vm["mos-shared-evals-num"].as<int>();
     }


     if (vm.count("mos-shared-evals-percent")) {
       if (vm["mos-shared-evals-percent"].as<long double>() < 0 || vm["mos-shared-evals-percent"].as<long double>() > 1){
         std::cerr << "[GAEDAConfig] Error: Wrong number of shared evaluations for MOS HRH: " << vm["mos-shared-evals-percent"].as<long double>() << ". Please select a value within (0, 1]. " << std::endl;
         return false;
       }
       mSharedEvalsPercent = vm["mos-shared-evals-percent"].as<long double>();
     }


     /*
      * Debug Options
      */

     if (vm.count("random-seed")) {
       seed = vm["random-seed"].as<unsigned>();
     }


     if (vm.count("debug-flag")) {
       loop = vm["debug-flag"].as<int>();
     }


     /* Get the problem module */
     if (!embedded) {
       if (!vm.count("problem_module")) {
         std::cerr << "[GAEDAConfig] Error: problem module not provided." << std::endl;
         return false; // No problem module provided
       }

       mProblemFile   = strdup((char*)vm["problem_module"].as<string>().c_str());
       mProblemStruct = GAloadProblem (mProblemFile);

       if (mProblemStruct == NULL) {
         if( comm_manager_.isIslandMaster() )
           std::cerr << "[GAEDAConfig] ERROR: Couldn't load problem " << vm["problem_module"].as<string>() << std::endl;
         return false;
       }
     }

   }
   catch(exception& e){
     std::cerr << "[GAEDAConfig] Exception: " << e.what() << std::endl;
     return false;
   }

   return true;

}


// Show command line arguments
void GAEDAConfig::showUsage( char *prgname ) {

  std::cerr << "Usage: mpirun [MPI OPTIONS] " << prgname << " [OPTIONS] problem_module" << std::endl
            << desc
            <<  "Notes:" << std::endl
            <<  "  * If you use -c, you have to also use -C" << std::endl
            <<  "  * problem_module is a path to a shared library defining the problem" << std::endl
            <<  std::endl;
}


GAGenome::OptCriterion GAEDAConfig::getOptCriterion() const{
  int struct_opt_crit = mProblemStruct->optCriterion;
  return (struct_opt_crit >= 0) ? (GAGenome::OptCriterion) struct_opt_crit : GAGenome::MAXIMIZATION;
}


GASelectionScheme* GAEDAConfig::getSelector () const{
  GASelectionScheme *selector = NULL, *sel_tmp = NULL;
  switch(selector_) {
  case Tournament:
    sel_tmp  = new GARouletteWheelSelector();
    selector = new GATournamentSelector(selector_degree_, *sel_tmp);
    delete sel_tmp;
    break;
  case Unif_Tournament:
    sel_tmp  = new GAUniformSelector();
    selector = new GATournamentSelector(selector_degree_, *sel_tmp);
    delete sel_tmp;
    break;
  case Roulette:
    selector = new GARouletteWheelSelector();
    break;
  case RandomSelec:
    selector = new GAUniformSelector();
    break;
  case Rank:
    selector= new GARankSelector();
    break;
  case Stochastic:
    selector = new GASRSSelector();
    break;
  case Deterministic:
    selector = new GADSSelector();
    break;
  default:
    throw runtime_error("unrecognized selector");
  }
  return selector;
}


GAScalingScheme* GAEDAConfig::getScaling () const{
  GAScalingScheme* scaling = NULL;
  switch(scaling_) {
  case None:
    scaling = new GANoScaling();
    break;
  case Standard:
    scaling = new GAStandardScaling();
    break;
  case Linear:
    scaling = new GALinearScaling();
    break;
  case Sigma_Truncation:
    scaling = new GASigmaTruncationScaling();
    break;
  case Power_Law:
    scaling = new GAPowerLawScaling();
    break;
  default:
    throw runtime_error("unrecognized scaling method");
  }
  return scaling;
}


GAGenome* GAEDAConfig::getOptimum() const{
  return (mProblemStruct->optimum == NULL) ? NULL : mProblemStruct->optimum();
}


MOSQuality* GAEDAConfig::getQualityFunction (void) const {
  if      (mQualityMeasure == "fAvg")      return new AverageFitnessQuality();
  else if (mQualityMeasure == "fitIncAvg") return new AverageFitnessIncrementQuality();
  else if (mQualityMeasure == "dAvg")      return new AverageDiversityQuality();
  else if (mQualityMeasure == "Compass")   return new CompassQuality();
}

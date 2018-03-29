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

#include "Recombinator.h"
#include "ClassicElitism.h"
#include "DEElitism.h"
#include "PopElitism.h"
#include "PercentElitism.h"
#include "GARealOps.h"
#include "GAGenealogy.h"
#include "GAGenealogyMemory.h"
#include "GAGenealogyTracer.h"
#include "MOSParticipationFunc.h"
#include "MOSParticipationFunction.h"
#include "MOSQualityFunction.h"
#include "VNS.h"
#include "islands/CommManager.h"
#include "islands/CentroidsGenerator.h"
#include "islands/RandomHighSeparatedCentGenerator.h"
#include "islands/EAIslandsTopology.h"
#include "islands/GAIslandsTopologyA2A.h"
#include "islands/GAIslandsTopologyRandom.h"
#include "islands/GAIslandsTopologyRing.h"
#include "islands/GAIslandsTopologyRing2.h"
#include "islands/GAIslandsTopologyHyperCube.h"
#include "islands/GAIslandsTopologyNearestNeighbor.h"
#include "islands/GAIslandsTopologyFurthestNeighbor.h"
#include "islands/EAEmigrantsSelector.h"
#include "islands/EABestEmigrantsSelector.h"
#include "islands/EARandomEmigrantsSelector.h"

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
   mEvalRouteCalls = 10;
   mConvCrit = EVALS;
   mPopSize = 100;
   CrossOp = 0;
   mPrintBestSol = true;

   mUsePrecission = false;

   mUseCentroids  = false;
   mProblemStruct = NULL;
   mLogFile       = NULL;
   mProblemData   = NULL;

   // Default values
   vnstinitratio_              = 0;
   vnstinitprob_               = 0;
   vnslsratio1_                = 0;
   vnslsprob_                  = 0;
   vnslsratio2_                = 0;
   vnsmaxevalsshaker_          = 0;
   vnsuseallshakersineachstep_ = false;
   vnsmaxshakersres_           = 1;
   vnsupdatepenalizations_     = false;
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

       ("log-option,b", po::value< vector<string> >()->composing(), "Set the statistics that should be traced among: score, fitness, participation, quality, age, etc.")

       ("stats-frequency,S", po::value<int>()->default_value(1), "Set the frequency to display the statistics.")

       ("elitism-percentage", po::value<long double>()->default_value(1.0), "Set the percentage of the old population that is brought into the new one.")

       ("recombination,B", po::value<string>()->default_value("full_elitism"), "Set the recombination scheme (single_elitism, full_elitism, percent_elitism, de_elitism).")

       ("scaling", po::value<string>()->default_value("linear"), "Set the scaling scheme to compute fitness from score (none, standard, linear, sigma_truncation, power_law).")

       ("algorithm,a", po::value< string >()->default_value("ga"), "Set the algorithm to be used (ga, ssga, eda, de, mos, mos2, mosmultideme, mosrl).")

       ("problem-size,s", po::value<int>()->default_value(1), "Set the size of the problem (not honored by all problems).")

       ("additional-data,A", po::value< string >(), "Additional data to be sent to the problem module.")

       ("repetitions,r", po::value<int>()->default_value(1), "Number of repetitions of the execution of the algorithm.")

       ("function", po::value<int>()->default_value(1), "Number of function for those benchmarks in which all the functions are encoded in the same problem module.")

       ("print-best-sol",po::value<bool>()->default_value(true), "Set if the best solution is going to be print at the end of the execution")

       ("cmaes-max-pop-size", po::value<int>()->default_value(100), "Maximum population size for the CMAES algorithm (when using population size increase).")

       ("cmaes-initials", po::value<string>()->default_value("./initials.par"), "Initial config for the CMAES algorithm.")

       ("bbob-output-path", po::value<string>()->default_value("./"), "Output directory for the BBOB problems.")

       ("bbob-config-name", po::value<string>()->default_value("BBOB_EXEC"), "Config name for this BBOB execution.")

       ("reset-percentage", po::value<long double>()->default_value(1.0), "Percentage of the population that will be re-initialized after a population restart.")

       ("rerun-percentage", po::value<long double>()->default_value(1.0), "Percentage of the population that will be re-initialized in multiple runs.")

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

       ("pop-or-gen-convergence",po::value< string >(),"Set a population or generation convergence criterion.")

       ("pop-or-evals-convergence",po::value< string >(),"Set a population or evaluations convergence criterion.")

       ("evaluations", po::value<unsigned long>(), "Set maximum number of fitness evaluations.")

       ("maxEvalRouteCalls", po::value<int>(), "Set maximum number of calls to the routeEval method (Routing Algs only).")

       ("timeConv", po::value<long double>(), "Set maximum number of seconds for the algorithm to execute.")

       ("precission", po::value<long double>()->default_value(1e-20), "Set required precission for fitness to consider convergence (in combination with evaluations).")

       ("use-precission", "If present, it activates the use of the given precission to stop the algorithm before the convergence cirterion is met.")

       ;

     po::options_description parallel_opts("Parallel model options:");
     parallel_opts.add_options()

       ("pop-island,p", po::value<int>(), "Set the population size per island.")

       ("pop-total,P", po::value<int>(), "Set the population size shared among all islands.")

       ("topology,t", po::value< string >()->default_value("ring"), "Set the topology (mesh, a2a, ring, ring2, pairs, nearestneighbor, furthestneighbor, random, hypercube).")

       ("sync-model,m", po::value< string >()->default_value("seq"), "Set the synchronization model (sync, async, asyncfast, seq).")

       ("frequency,f", po::value<int>()->default_value(10), "Set the migration frequency.")

       ("emigrant-size,e", po::value<int>()->default_value(1), "Set the size of emigrant populations.")

       ("emigrant-percentage,E", po::value<int>(), "Set the percentage of the island population to migrate.")

       ("emigrant-selector,w", po::value<string>()->default_value("best"), "Set emigrant selector (random, best).")

       ("degree,d", po::value<int>()->default_value(0), "Degree of Topology")

       ("n,n", po::value<string>()/*->default_value("gabriel")*/, "Set the neighbor condition for the nearest neighbor topology (gabriel,relative_neighbors,square).")
       ("Y,Y", "Set dynamic migration type used with the nearest neighbor topology (sender_oriented, receivers_oriented). Default: senders_oriented")
       ("N,N", po::value<int>()/*->default_value(1)*/, "Set the minimum number of neighbors that each medoid must have.")
       ("G,G", po::value<int>()/*->default_value(number of islands)*/,"Set the maximum number of neighbors that each medoid must have. Default: number of islands")
       ("u,u", po::value<int>()/*->default_value(mig_period)*/,"Set islands clustering period . Default: mig_period")
       ("c2", po::value<string>()->default_value("random"), "Use centroid initializers. Actual types: random, hs-random (high separated random). with random the user can optionally specify the maximum number ot tries to find a good distribution. Default: Centroids not used")
       ("z,z", po::value<string>()->default_value("random"), "Set islands initializer (random|voronoi).")

       ("c,c", po::value< string >(), "Use centroid initializers. Prob is the mutation probability."
        "used to generate individuals from the centroid. min is the minimum "
        "distance between centroids. iters is the iteration limit set on the "
        "centroid generation procedure.")

       ;

     po::options_description dist_routealg_options(" distributed routing algorithm options:");
     dist_routealg_options.add_options()

       ("distrouting-freq", po::value<int>(), "Set the frequency of the sync model of the distributed routing algorithm.")

       ("routedistributor", po::value<string>()->default_value("random"), "Set the name of the route distributor used.")

       ("distrouting-emigrant", po::value< string >()->default_value("best"), "Set the type of the emigrant route selected for the exchanging of information")
     ;

     po::options_description ga_opts("Genetic Algorithms options:");
     ga_opts.add_options()

       ("mut-prob,M", po::value<long double>()->default_value(0.01), "Set mutation probability.")

       ("cross-prob,X", po::value<long double>()->default_value(0.9), "Set crossover probability.")

       ("selection-scheme,R", po::value<string>()->default_value("tournament:2"), "Selection scheme (tournament:n,roulette,random,rank,stochastic,deterministic).")

       ;

     po::options_description de_opts("Differential Evolution options:");
     de_opts.add_options()

       ("F,F", po::value<long double>(), "Set the DE mutation factor Phi.")

       ("CR", po::value<long double>(), "Set the DE crossover probability.")

       ("de-crossover,J", po::value<string>(), "Set the DE corssover operator.")

       ;

     po::options_description es_opts("Evolution Strategies options:");
     es_opts.add_options()

       ("es-mu", po::value<unsigned>()->default_value(10), "Set the ES population size.")

       ("es-ro", po::value<unsigned>()->default_value(2), "Set the ES mixing degree (number of parents used in recombination).")

       ("es-lambda", po::value<unsigned>()->default_value(100), "Set the ES offspring size.")

       ;


     po::options_description vns_opts("VNS options:");
     vns_opts.add_options()
       ("vnstinitratio", po::value<long double>()->default_value(0.005),"Set the VNS t init ratio")

       ("vnstinitprob",  po::value<long double>()->default_value(0.2),  "Set the VNS t init probability")

       ("vnslsratio1",   po::value<long double>()->default_value(1.02), "Set the VNS 1st worsening ratio for executing the LS")

       ("vnslsprob",     po::value<long double>()->default_value(0.021), "Set the VNS probability of executing the LS")

       ("vnslsratio2",   po::value<long double>()->default_value(1.05), "Set the VNS 2nd worsening ratio for executing the LS")

       ("vnsmaxevalsshaker",po::value<int>()->default_value(0), "Set the VNS maximum #evals per shaker")

       ("vnsuseallshakersineachstep",po::value<bool>()->default_value(false), "Set if the VNS is using all the shakers at each step")

       ("vnsmaxshakersres",po::value<int>()->default_value(1), "Set the maximum number of results used for combining the results")

       ("vnsupdatepenalizations",po::value<bool>()->default_value(false), "Set if the VNS is going to update the penalizations from the objective function at each step")

       ;

     po::options_description mts_opts("MTS options:");
     mts_opts.add_options()

       ("mts-adjustFailed", po::value<long double>()->default_value(0.5), "Set the adjustment factor (*) to be used after one iteration of the LS without any improvement")

       ("mts-adjustMin",    po::value<long double>()->default_value(0.4), "Set the adjustment factor (*) to be used if SR reaches its minimum vlue (1e-14)")

       ("mts-moveLeft",     po::value<long double>()->default_value(1.0), "Set the movement factor (*) to be used when moving a solution left")

       ("mts-moveRight",    po::value<long double>()->default_value(0.5), "Set the movement factor (*) to be used when moving a solution right")

     ;

     po::options_description mtsred_opts("MTS Reduced options:");
     mtsred_opts.add_options()

       ("mtsred-searchProb", po::value<long double>()->default_value(0.9), "Set the cummulated percentage of dimensions that will be searched (instead of exploring all of them , we just explore those cummulating this percentage")

       ("mtsred-minProb", po::value<long double>()->default_value(0.05), "Set the minimum percentage of dimensions to explore");

     po::options_description sw_opts("Solis-Wets options:");
     sw_opts.add_options()

       ("sw-maxSuccess",    po::value<int>()->default_value(5), "Set the maximum number of improvements with the current delta value, before increasing delta")

       ("sw-maxFailed",     po::value<int>()->default_value(3), "Set the maximum number of failures with the current delta value, before decreasing delta")

       ("sw-adjustSuccess", po::value<long double>()->default_value(2.0), "Set the adjustment factor (*) to be used after 'sw-maxSuccess' iterations of the LS with improvements")

       ("sw-adjustFailed",  po::value<long double>()->default_value(0.5), "Set the adjustment factor (*) to be used after 'sw-maxFailed' iterations of the LS without improvements")

       ("sw-delta",         po::value<long double>()->default_value(1.2), "Set the delta factor of Solis-Wets algorithm")

     ;

     po::options_description adapde_opts("Adaptive DE options:");
     adapde_opts.add_options()

       ("adapde-Fl",    po::value<long double>()->default_value(0.1), "Set the lower bound for the adjustment interval of the F parameter of a DE algorithm")

       ("adapde-Fu",    po::value<long double>()->default_value(0.9), "Set the upper bound for the adjustment interval of the F parameter of a DE algorithm")

       ("adapde-tauF",     po::value<long double>()->default_value(0.1), "Set the threshold to choose if an adjustment on the F parameter should be done")

       ("adapde-tauCR",     po::value<long double>()->default_value(0.1), "Set the threshold to choose if an adjustment on the CR parameter should be done")

     ;

     po::options_description stsde_opts("STSDE options:");
     stsde_opts.add_options()

       ("stsde-prob",    po::value<long double>()->default_value(0.4), "Set the probability to apply STSDE")

     ;

     po::options_description mos_opts("MOS options:");
     mos_opts.add_options()

       ("techs-repository", po::value<string>(), "Set the repository with all the MOS techniques that can be used.")

       ("use-tech", po::value< vector<string> >()->composing(), "Set the techniques to use in MOS algorithm.")

       ("part-function", po::value<string>()->default_value("dynamic"), "Set the participation function of a MOS algorithm (constant, dynamic, dynamicMaxPart, dynamic2, alternating_qual).")

       ("minimum-part", po::value<long double>()->default_value(0.05) , "Set the minimum participation of a MOS technique.")

       ("evolutive-approach", po::value<string>()->default_value("central"), "Set the evolutive approach of a MOS algorithm (central, autonomic).")

       ("mos-genome-printing", po::value<string>()->default_value("full"), "Sets the way that MOS genomes should be printed ('default' for default encoding and 'full' for every encoding).")

       ("mosql-pmax", po::value<long double>()->default_value(0.9), "Sets the selection pressure for actions in MOSEA-QL.")

       ("mosrl-alfa", po::value<long double>()->default_value(0.7), "Sets the learning rate for MOSEA-RL.")

       ("mosrl-gamma", po::value<long double>()->default_value(0.8), "Sets the discount rate for MOSEA-RL.")

       ("mosrl-beta", po::value<long double>()->default_value(0.3), "Sets the reward rate for MOSEA-RL.")

       ("mosrl-delta", po::value<long double>()->default_value(0.1), "Sets the Delta rate for MOSEA-RL (HPC only).")

       ("mosrl-deltal", po::value<long double>()->default_value(0.2), "Sets the Delta Lose for MOSEA-RL (WoLF only).")

       ("mosrl-deltaw", po::value<long double>()->default_value(0.05), "Sets the Delta Win for MOSEA-RL (WoLF only).")

       ("mosrl-policy", po::value<string>()->default_value("pmax"), "Sets the Learnign Policy for MOSEA-RL (pmax, phc, wolf).")

       ("genealogy,g", po::value<string>()->default_value("none"), "Set type of genealogy (none, tracer, memory, all, full, statistics).")

       ("print-genealogy", po::value<string>()->default_value("no"), "Activate/deactivate genealogy printing (yes, no).")

       ("quality-measure", po::value<string>()->default_value("fAvg"), "Set type of quality measure used for dynamic adjustment of participation (fAvg, fitIncAvg, dAvg, Compass, NSC).")

       ("weighted-autonomic", po::value<string>()->default_value("no"), "Weighted (yes) or Arithmetic (no) average in MOS autonomic.")

       ("mos-bonus", po::value<long double>()->default_value(0.01), "Sets the reward bonus parameter.")

       ("het-mos-cfg", po::value<string>(), "Set the configuration file for a heterogeneous MOS algorithm. If non file is given, the MOS algorithm will be homogeneous (the same techniques will be used in all the islands).")

       ("mos-shared-evals-num", po::value<int>()->default_value(1500), "Sets the number of shared evaluations to distribute in a MOS HRH algorithm.")

       ("mos-shared-evals-percent", po::value<long double>()->default_value(0.01), "Sets the percentage of shared evaluations to distribute in a MOS HRH algorithm.")

       ("mos-shared-evals-factor", po::value<long double>(), "Sets the factor to compute the number of shared evaluations to distribute in a MOS HRH algorithm from the population size.")

       ;

     po::options_description debug_opts("Debug options:");
     debug_opts.add_options()

       ("random-seed,q", po::value<unsigned>()->default_value(0), "Seed for random number generation subsystem. 0 for random seed.")

       ("debug-flag,Z", po::value<int>()->default_value(0), "Set a sleep point in the main program to allow easy GDB attaching and debugging.")

       ;

     // Join all options groups together
     desc.add(general_opts).add(convergence_opts).add(parallel_opts).add(dist_routealg_options).add(ga_opts).add(de_opts).add(es_opts).add(vns_opts).add(mts_opts).add(mtsred_opts).add(sw_opts).add(adapde_opts).add(stsde_opts).add(mos_opts).add(debug_opts);

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
         if (vm["log-option"].as< vector<string> >()[i] == "score")               logstats_names.push_back(Score);
         if (vm["log-option"].as< vector<string> >()[i] == "fitness")             logstats_names.push_back(Fitness);
         if (vm["log-option"].as< vector<string> >()[i] == "age")                 logstats_names.push_back(Age);
         if (vm["log-option"].as< vector<string> >()[i] == "opt_dist")            logstats_names.push_back(OptDist);
         if (vm["log-option"].as< vector<string> >()[i] == "opt_ncomps")          logstats_names.push_back(OptNcomps);
         if (vm["log-option"].as< vector<string> >()[i] == "gen_div")             logstats_names.push_back(GenDiv);
         if (vm["log-option"].as< vector<string> >()[i] == "fit_inc")             logstats_names.push_back(FitInc);
         if (vm["log-option"].as< vector<string> >()[i] == "native_prcnt")        logstats_names.push_back(NativePrcnt);
         if (vm["log-option"].as< vector<string> >()[i] == "entropy")             logstats_names.push_back(Entropy);
         if (vm["log-option"].as< vector<string> >()[i] == "gref_bias")           logstats_names.push_back(GrefBias);
         if (vm["log-option"].as< vector<string> >()[i] == "participation")       logstats_names.push_back(Participation);
         if (vm["log-option"].as< vector<string> >()[i] == "quality")             logstats_names.push_back(Quality);
         if (vm["log-option"].as< vector<string> >()[i] == "evaluations")         logstats_names.push_back(Evals);
         if (vm["log-option"].as< vector<string> >()[i] == "fss-avg")             logstats_names.push_back(FSSAvg);
         if (vm["log-option"].as< vector<string> >()[i] == "fss-best")            logstats_names.push_back(FSSBest);
         if (vm["log-option"].as< vector<string> >()[i] == "fss-worst")           logstats_names.push_back(FSSWorst);
         if (vm["log-option"].as< vector<string> >()[i] == "improvements")        logstats_names.push_back(Improve);
         if (vm["log-option"].as< vector<string> >()[i] == "routingalgscore")     logstats_names.push_back(RoutingAlgScore);
         if (vm["log-option"].as< vector<string> >()[i] == "vnssuccessfulshaker") logstats_names.push_back(VNSSuccessfulShaker);
         if (vm["log-option"].as< vector<string> >()[i] == "vnsshakerpos")        logstats_names.push_back(VNSShakerpos);
         if (vm["log-option"].as< vector<string> >()[i] == "bestsolpickuptime")   logstats_names.push_back(RoutingAlgBestSolPickupTime);
         if (vm["log-option"].as< vector<string> >()[i] == "bestsoldeliverytime") logstats_names.push_back(RoutingAlgBestSolDeliveryTime);
         if (vm["log-option"].as< vector<string> >()[i] == "bestsolnonpenscore")  logstats_names.push_back(RoutingAlgBestSolNonPenScore);
         if (vm["log-option"].as< vector<string> >()[i] == "bestsoltwv")          logstats_names.push_back(RoutingAlgBestSolTWV);
         if (vm["log-option"].as< vector<string> >()[i] == "bestsolloadv")        logstats_names.push_back(RoutingAlgBestSolLoadV);
         if (vm["log-option"].as< vector<string> >()[i] == "bestsolridev")        logstats_names.push_back(RoutingAlgBestSolRideV);

         if (vm["log-option"].as< vector<string> >()[i] == "evalroutecalls")      logstats_names.push_back(EvalRouteCalls);
         if (vm["log-option"].as< vector<string> >()[i] == "runtime")             logstats_names.push_back(Runtime);
         if (vm["log-option"].as< vector<string> >()[i] == "freemem")             logstats_names.push_back(Freemem);

       }
     }


     if (vm.count("stats-frequency")) {
       mStatsDispFreq = vm["stats-frequency"].as<int>();

       if (mStatsDispFreq <= 0){
         std::cerr << "[GAEDAConfig] Error: the stats display frequency specified: " << mStatsDispFreq << " is <= 0." << std::endl;
         return false;
       }
     }


     if( vm.count("elitism-percentage") ){
       elitism_percentage = vm["elitism-percentage"].as<long double>();
       if (elitism_percentage > 1.0 || elitism_percentage < 0.0){
         std::cerr << "[GAEDAConfig] Error: the elisitm percentage value must be within [0, 1] interval." << std::endl;
         return false;
       }
     }


     if (vm.count("recombination")) {
       if (vm["recombination"].as<string>() == "single_elitism")
         recombinator_ = SingleElitism;
       else if (vm["recombination"].as<string>() == "full_elitism")
         recombinator_ = FullElitism;
       else if (vm["recombination"].as<string>() == "de_elitism")
         recombinator_ = DE_Elitism;
       else if (vm["recombination"].as<string>() == "percent_elitism")
         recombinator_ = Percent_Elitism;
       else {
         cerr << endl << "[GAEDAConfig] Error: unrecognized recombination method." << endl;
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


     if (vm.count("repetitions")) {
       mRepetitions = vm["repetitions"].as<int>();

       if (mRepetitions < 1) {
         std::cerr << "[GAEDAConfig] Error: The number of repetitions must be greater than 0." << std::endl;
         return false;
       }
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
       if (sscanf(vm["pop-or-gen-convergence"].as<string>().c_str(),"%lf:%d",&mConv,&mIters) ) {
         mConvCrit = POPCONVITERS;
       }
       else {
         cerr << "[GAEDAConfig] Error: Options for convergence upon population or ngens not recognised." << endl;
         return false;
       }
     }


     if (vm.count("pop-or-evals-convergence")) {
       if (sscanf(vm["pop-or-evals-convergence"].as<string>().c_str(),"%lf:%d",&mConv,&mEvals) ) {
         mConvCrit = POPCONVEVALS;
       }
       else {
         cerr << "[GAEDAConfig] Error: Options for convergence upon population or evaluations not recognised." << endl;
         return false;
       }
     }


     if (vm.count("evaluations")) {
       mConvCrit = EVALS;
       mEvals = vm["evaluations"].as<unsigned long>();

       if (mEvals < mPopSize) {
         std::cerr << "[GAEDAConfig] Error: Wrong number of fitness evaluations: " << mEvals << std::endl;
         return false;
       }
     }


     if (vm.count("maxEvalRouteCalls")) {
       mConvCrit   = EVALROUTECALLS;
       mEvalRouteCalls = vm["maxEvalRouteCalls"].as<int>();

       if (mEvalRouteCalls < 0) {
         std::cerr << "[GAEDAConfig] Error: Wrong number of maximum route calls: " << mEvals << std::endl;
         return false;
       }
     }


     if (vm.count("timeConv")) {
       mConvCrit      = TIMECONV;
       maxExecSeconds = vm["timeConv"].as<long double>();

       if (maxExecSeconds < 0) {
         std::cerr << "[GAEDAConfig] Error: Wrong number of seconds: " << maxExecSeconds << std::endl;
         return false;
       }
     }



     if (vm.count("precission")) {
       mPrecission = vm["precission"].as<long double>();
     }


     if (vm.count("use-precission")) {
       mUsePrecission = true;
     }


     // The algorithm should be processed after convergence options to allow
     // the error checking in the case of the MOS2 algorithm
     if (vm.count("algorithm")) {
       if (vm["algorithm"].as< string >() == "ga")
         mAlg = GAEDAConfig::GA;
       else if (vm["algorithm"].as< string >() == "ssga")
         mAlg = GAEDAConfig::SSGA;
       else if (vm["algorithm"].as< string >() == "eda")
         mAlg = GAEDAConfig::EDA;
       else if (vm["algorithm"].as< string >() == "de")
         mAlg = GAEDAConfig::DE;
       else if (vm["algorithm"].as< string >() == "es")
         mAlg = GAEDAConfig::ES;
       else if (vm["algorithm"].as< string >() == "vns")
         mAlg = GAEDAConfig::VNS;
       else if (vm["algorithm"].as< string >() == "tscordeau")
         mAlg = GAEDAConfig::TSCordeau;
       else if (vm["algorithm"].as< string >() == "mos")
         mAlg = GAEDAConfig::MOS;
       else if (vm["algorithm"].as< string >() == "mos2") {
         mAlg = GAEDAConfig::MOS2;
         if (mConvCrit != EVALS && mConvCrit != POPCONVEVALS) {
           if ( comm_manager_.isIslandMaster() )
             std::cerr << "[GAEDAConfig] Error: You can not use the 'mos2' algorithm with a convergence criterion different to 'evaluations' or 'pop-or-evals-convergence'." << std::endl;
             return false;
         }
       }
       else if (vm["algorithm"].as< string >() == "mosmultideme")
         mAlg = GAEDAConfig::MOSMultiDeme;
       else if (vm["algorithm"].as< string >() == "mosrl")
         mAlg = GAEDAConfig::MOSRL;
       else{
         cerr << "[GAEDAConfig] Error: unrecognised algorithm." << endl;
         return false;
       }
     }


     if (vm.count("function")) {
       mFunction = vm["function"].as<int>();
     }

     if (vm.count("print-best-sol")) {
       mPrintBestSol = vm["print-best-sol"].as<bool>();
     }

     if (vm.count("cmaes-max-pop-size")) {
       mCMAESMaxPopSize = vm["cmaes-max-pop-size"].as<int>();
     }


     if (vm.count("cmaes-initials")) {
       mCMAESInitials = vm["cmaes-initials"].as<string>();
     }


     if (vm.count("bbob-output-path")) {
       mBBOBOutputPath = vm["bbob-output-path"].as<string>();
     }


     if (vm.count("bbob-config-name")) {
       mBBOBConfigName = vm["bbob-config-name"].as<string>();
     }


     if (vm.count("reset-percentage")) {
       mResetPerc = vm["reset-percentage"].as<long double>();

       if (mResetPerc < 0.0 || mResetPerc > 1.0) {
         std::cerr << "[GAEDAConfig] Error: Wrong percentage for the reset-percentage parameter (values must be within the [0, 1] interval): " << mResetPerc << std::endl;
         return false;
       }
     }


     if (vm.count("rerun-percentage")) {
       mRerunPerc = vm["rerun-percentage"].as<long double>();

       if (mRerunPerc < 0.0 || mRerunPerc > 1.0) {
         std::cerr << "[GAEDAConfig] Error: Wrong percentage for the rerun-percentage parameter (values must be within the [0, 1] interval): " << mRerunPerc << std::endl;
         return false;
       }
     }


     /*
      * Parallel Model Options
      */

     if (vm.count("topology")) {
       if (vm["topology"].as< string >() == "mesh")
         mTop = GAEDAConfig::Mesh;
       else if (vm["topology"].as< string >() == "a2a")
         mTop = GAEDAConfig::a2a;
       else if (vm["topology"].as< string >() == "ring")
         mTop = GAEDAConfig::Ring;
       else if (vm["topology"].as< string >() == "ring2")
         mTop = GAEDAConfig::Ring2;
       else if (vm["topology"].as< string >() == "pairs")
         mTop = GAEDAConfig::Pairs;
       else if (vm["topology"].as< string >() == "nearestneighbor")
         mTop = GAEDAConfig::NearestNeighbor;
       else if (vm["topology"].as< string >() == "furthestneighbor")
         mTop = GAEDAConfig::FurthestNeighbor;
       else if (vm["topology"].as< string >() == "random")
         mTop = GAEDAConfig::RandomTopology;
       else if (vm["topology"].as< string >() == "hypercube")
         mTop = GAEDAConfig::HybridTopology;
       else{
         std::cerr << "[GAEDAConfig] Error: Topology not recognised." << std::endl;
         return false;
       }
     }


     if (vm.count("sync-model")) {
       if (vm["sync-model"].as< string >() == "sync")
         mModel = GAEDAConfig::SYNC;
       else if (vm["sync-model"].as< string >() == "async")
         mModel = GAEDAConfig::ASYNC;
       else if (vm["sync-model"].as< string >() == "asyncfast")
         mModel = GAEDAConfig::ASYNCFAST;
       else if (vm["sync-model"].as< string >() == "seq")
         mModel = GAEDAConfig::SEQ;
       else if (vm["sync-model"].as< string >() == "clust")
         mModel = GAEDAConfig::CLUST;
       else if (vm["sync-model"].as< string >() == "routingdistsync")
         mModel = GAEDAConfig::ROUTINGDISTSYNC;
       else{
         std::cerr << "[GAEDAConfig] Error: Synchronization model not recognized." << std::endl;
         return false;
       }

       comm_manager_.initialize(mModel);

       if (comm_manager_.getNumIslands() > 1 && vm["sync-model"].as< string >() == "seq") {
         if ( comm_manager_.isIslandMaster() ) {
           std::cerr << "[GAEDAConfig] Error: You can not use 'seq'(sequential) mode with more than 1 node. Aborting..." << std::endl;
         }
         return false;
       }
       else if (comm_manager_.getNumIslands() == 1 && vm["sync-model"].as< string >() != "seq") {
         if ( comm_manager_.isIslandMaster() ) {
           std::cerr << "[GAEDAConfig] Warning: You can not use '" << vm["sync-model"].as< string >() << "' mode with just 1 node. Aborting..." << std::endl;
         }
         return false;
       }
     }


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


     if (vm.count("frequency")) {
       mMF = vm["frequency"].as<int>();

       if (mMF <= 0) {
         std::cerr << "[GAEDAConfig] Error: Migration frequency must be greater than 0." << std::endl;
         return false;
       }
     }


     if (vm.count("emigrant-size")) {
       mMigPopSize = vm["emigrant-size"].as<int>();

       if (mMigPopSize <= 0) {
         std::cerr << "[GAEDAConfig] Error: emigrant size must be greater than 0." << std::endl;
         return false;
       }
     }


     if (vm.count("emigrant-percentage")) {
       tmp = vm["emigrant-percentage"].as<int>();

       if (tmp < 0 || tmp > 100){
         std::cerr << "[GAEDAConfig] Error: emigrant percentage out of range [0, 100]." << std::endl;
         return false;
       }

       mMigPopSize = tmp * mPopSize / 100;

       if (mMigPopSize == 0)
         mMigPopSize = 1;
     }


     if (vm.count("emigrant-selector")) {
       if (vm["emigrant-selector"].as<string>() == "random")
         mEmigrantsSelector = RandomEmigrantsSelector;
       else if (vm["emigrant-selector"].as<string>() == "best")
         mEmigrantsSelector = BestEmigrantsSelector;
       else{
         std::cerr << "[GAEDAConfig] Error: Unrecognized emigrants selector." << std::endl;
         return false;
       }
     }


     if (vm.count("degree")) { // Degree of Topology
       mTopDegree = vm["degree"].as<int>();
       mMinNeighborsNum = mTopDegree;
       mMaxNeighborsNum = mTopDegree;
       if (mTopDegree <0){
         std::cerr << "[GAEDAConfig] Error: the degree of topology  must be greater than or equal to zero." << std::endl;
         return false;
       }
     }


     if (vm.count("n")) {
       if (vm["n"].as<string>() == "gabriel")
         mNeighbor_cond = GAEDAConfig::Gabriel;
       else if (vm["n"].as<string>() == "relative_neighbors")
         mNeighbor_cond = GAEDAConfig::Rel_neighbors;
       else if (vm["n"].as<string>() == "square")
         mNeighbor_cond = GAEDAConfig::Square;
       else{
         std::cerr << "[GAEDAConfig] Error: neighborhood condition unrecognised." << std::endl;
         return false;
       }
     }


     if (vm.count("Y")) { // TODO...
     }


     if (vm.count("N")) {
       mMinNeighborsNum = vm["N"].as<int>();

       if (mTopDegree != 0){
         std::cerr << "[GAEDAConfig] Error: topology degree has been stablished with option -d to value: " << mTopDegree << std::endl;
         return false;
       }

       if (mMinNeighborsNum < 0) {
         std::cerr << "[GAEDAConfig] Error: minimum number of neighbors must be greater than or equal to zero." << std::endl;
         return false;
       }
     }


     if (vm.count("G")) {
       mMaxNeighborsNum = vm["G"].as<int>();

       if (mTopDegree != 0){
         std::cerr << "[GAEDAConfig] Error: topology degree has been stablished with option -d to value: " << mTopDegree << endl;
         return false;
       }

       if (mMaxNeighborsNum < 0) {
         std::cerr << "[GAEDAConfig] Error: maximum number of neighbors must be greater than or equal to zero." << std::endl;
         return false;
       }
     }


     if (vm.count("u")) {
       clustering_period_ = vm["u"].as<int>();

       if (clustering_period_ <= 0){
         cerr << "[GAEDAConfig] Error: the clustering period must be greater than 0." << endl;
         return false;
       }
     }


     if (vm.count("c2")) {
       if (vm["c2"].as<string>() == "random")
         mCentGen = Random;
       else if (vm["c2"].as<string>() == "hs-random")
         mCentGen = HighSeparatedRandom;
       else {
         std::cerr << "[GAEDAConfig] Error: the centroid initializer specified: " << vm["c2"].as<string>() << " is not valid." << std::endl;
         return false;
       }
     }


     if (vm.count("z")) {
       if (vm["z"].as<string>() == "random")
         mIslandsInit =  RandomIslandInit;
       else if (vm["z"].as<string>() == "voronoi")
         mIslandsInit = VoronoiIslandInit;
       else{
         std::cerr << "[GAEDAConfig] Error: islands initializer not recognaised." << std::endl;
         return false;
       }
     }


     if (vm.count("c")) {
       if (sscanf (vm["c"].as< string >().c_str(), "%lf:%lf:%d", &mCentProb, &mMinCentDist, &mMaxCentIters) != 3) {
         std::cerr << "[GAEDAConfig] Error: Wrong centroids initializer." << std::endl;
         return false;
       }
       if (mCentProb < 0 || mCentProb > 1 || mMinCentDist< 0 || mMaxCentIters <= 0){
         std::cerr << "[GAEDAConfig] Error: Wrong centroids initializer." << std::endl;
         return false;
       }
       mUseCentroids = true;
     }

     /*
      * Distributed Routing Algorithm options
      */
     if (vm.count("distrouting-freq")) {        // Using the same attribute as the migration frequency. Need to see if
       mMF = vm["distrouting-freq"].as<int>();  // it is harmful to be so redundant

       if (mMF <= 0) {
         std::cerr << "[GAEDAConfig] Error: Migration frequency must be greater than 0." << std::endl;
         return false;
       }
     }

     if (vm.count("distrouting-emigrant")) {
       string name = vm["distrouting-emigrant"].as<string>();
       if      (name == "best")   emigrant_route_t_ = DistRoutingAlg::BEST;
       else if (name == "actual") emigrant_route_t_ = DistRoutingAlg::ACTUAL;
       else {
         std::cerr << "[GAEDAConfig] Error: invalid distributed routing emigrant option" << std::endl;
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
      * Differential Evolution Options
      */

     if (vm.count("F")) {
       if (vm["F"].as<long double>() <= 0 || vm["F"].as<long double>() >= 2 ){
         std::cerr << "[GAEDAConfig] Error: Wrong Phi constant value: " << vm["F"].as<long double>() << ". Please select a value within (0, 2)." << std::endl;
         return false;
       }
       DE_F_factor = vm["F"].as<long double>();
     }


     if (vm.count("CR")) {
       if (vm["CR"].as<long double>() <= 0 || vm["CR"].as<long double>() >= 1 ){
         std::cerr << "[GAEDAConfig] Error: Wrong crossover probability: " << vm["CR"].as<long double>() << ". Please select a value within (0, 1) " << std::endl;
         return false;
       }
       DE_CR_factor = vm["CR"].as<long double>();
     }


     if (vm.count("de-crossover")) {
       if (vm["de-crossover"].as<string>() == "binomial")
         CrossOp = RealBinomialCrossover;
       else if (vm["de-crossover"].as<string>() == "exponential")
         CrossOp = RealExponentialCrossover;
       else  {
         std::cerr << "[GAEDAConfig] Error: DE crossover operator not recognised." << std::endl;
         return false;
       }
     }


     /*
      * Evolution Strategies Options
      */

     if (vm.count("es-mu")) {
       if (vm["es-mu"].as<unsigned>() <= 0){
         std::cerr << "[GAEDAConfig] Error: Wrong Mu constant value: " << vm["es-mu"].as<unsigned>() << ". Please select a value greater than 0." << std::endl;
         return false;
       }
       ES_mu = vm["es-mu"].as<unsigned>();
     }


     if (vm.count("es-ro")) {
       if (vm["es-ro"].as<unsigned>() <= 0){
         std::cerr << "[GAEDAConfig] Error: Wrong Ro constant value: " << vm["es-ro"].as<unsigned>() << ". Please select a value greater than 0." << std::endl;
         return false;
       }
       ES_ro = vm["es-ro"].as<unsigned>();
     }


     if (vm.count("es-lambda")) {
       if (vm["es-lambda"].as<unsigned>() <= 0){
         std::cerr << "[GAEDAConfig] Error: Wrong Lambda constant value: " << vm["es-lambda"].as<unsigned>() << ". Please select a value greater than 0." << std::endl;
         return false;
       }
       ES_lambda = vm["es-lambda"].as<unsigned>();
     }

     /*
      * VNS Options
      */

#define VNSOPT(optname,type) \
     if (vm.count(#optname)) { \
       type value = vm[#optname].as<type>(); \
       if (value < 0.0) { \
         std::cerr << "[GAEDAConfig] Error: Wrong " << #optname << " value: " << value <<". Please select a value greater or equals than 0" << std::endl;\
         return false; \
       } \
       optname##_ = value;\
     }

     VNSOPT(vnstinitratio,long double);
     VNSOPT(vnstinitprob,long double);
     VNSOPT(vnslsratio1,long double);
     VNSOPT(vnslsprob,long double);
     VNSOPT(vnslsratio2,long double);
     VNSOPT(vnsmaxevalsshaker,int);
     VNSOPT(vnsuseallshakersineachstep,bool);
     VNSOPT(vnsmaxshakersres,int);
     VNSOPT(vnsupdatepenalizations,bool);

     /*
      * MTS Options
      */

     if (vm.count("mts-adjustFailed")) {
       mts_adjustFailed = vm["mts-adjustFailed"].as<long double>();
     }

     if (vm.count("mts-adjustMin")) {
       mts_adjustMin = vm["mts-adjustMin"].as<long double>();
     }

     if (vm.count("mts-moveLeft")) {
       mts_moveLeft = vm["mts-moveLeft"].as<long double>();
     }

     if (vm.count("mts-moveRight")) {
       mts_moveRight = vm["mts-moveRight"].as<long double>();
     }

     /*
      * MTS Reduced Options
      */

     if (vm.count("mtsred-searchProb")) {
       mtsred_searchProb = vm["mtsred-searchProb"].as<long double>();
     }

     if (vm.count("mtsred-minProb")) {
       mtsred_minProb = vm["mtsred-minProb"].as<long double>();
     }

     /*
      * Solis-Wets Options
      */

     if (vm.count("sw-maxSuccess")) {
       sw_maxSuccess = vm["sw-maxSuccess"].as<int>();
     }

     if (vm.count("sw-maxFailed")) {
       sw_maxFailed = vm["sw-maxFailed"].as<int>();
     }

     if (vm.count("sw-adjustSuccess")) {
       sw_adjustSuccess = vm["sw-adjustSuccess"].as<long double>();
     }

     if (vm.count("sw-adjustFailed")) {
       sw_adjustFailed = vm["sw-adjustFailed"].as<long double>();
     }

     if (vm.count("sw-delta")) {
       sw_delta = vm["sw-delta"].as<long double>();
     }

     /*
      * Adaptive DE Options
      */

     if (vm.count("adapde-Fl")) {
       adapde_Fl = vm["adapde-Fl"].as<long double>();
     }

     if (vm.count("adapde-Fu")) {
       adapde_Fu = vm["adapde-Fu"].as<long double>();
     }

     if (vm.count("adapde-tauF")) {
       adapde_tauF = vm["adapde-tauF"].as<long double>();
     }

     if (vm.count("adapde-tauCR")) {
       adapde_tauCR = vm["adapde-tauCR"].as<long double>();
     }

     /*
      * STS DE Options
      */

     if (vm.count("stsde-prob")) {
       stsde_prob = vm["stsde-prob"].as<long double>();
     }

     /*
      * MOS Options
      */

     if (vm.count("techs-repository")){
       techsRepository = strdup( (char*) vm["techs-repository"].as< string >().c_str() );
     }
     else if (mAlg == GAEDAConfig::MOS || mAlg == GAEDAConfig::MOSRL) {
       std::cerr << "[GAEDAConfig] Error: you must specify a repository of techniques." << std::endl;
       return false;
     }


     if (vm.count("use-tech")) {
       techsToUse = vm["use-tech"].as< vector<std::string> >();
     }


     if (vm.count("part-function")) {
       if (vm["part-function"].as<string>() == "constant")
         partFunction = constantPF;
       else if (vm["part-function"].as<string>() == "dynamic")
         partFunction = dynQualityMOSPF;
       else if (vm["part-function"].as<string>() == "dynamicMaxPart")
         partFunction = dynQualityMaxPartPF;
       else if (vm["part-function"].as<string>() == "dynamic2")
         partFunction = dynQualityMOSPF2;
       else if (vm["part-function"].as<string>() == "alternating_qual")
         partFunction = alternatingQualityPF;
       else  {
         std::cerr << "[GAEDAConfig] Error: MOS participation function not recognised." << std::endl;
         return false;
       }
     }


     if (vm.count("minimum-part")) {
       if (vm["minimum-part"].as<long double>() < 0 || vm["minimum-part"].as<long double>() > 1 ){
         std::cerr << "[GAEDAConfig] Error: Wrong minimum participation for MOS techniques: " << vm["minimum-part"].as<long double>() << ". Please select a value within [0, 1]. " << std::endl;
         return false;
       }
       minPart = vm["minimum-part"].as<long double>();
     }


     if (vm.count("evolutive-approach")) {
       if (vm["evolutive-approach"].as<string>() == "central")
         evolutiveApproach = CentralEvolution;
       else if (vm["evolutive-approach"].as<string>() == "autonomic")
         evolutiveApproach = AutonomicEvolution;
       else  {
         std::cerr << "[GAEDAConfig] Error: MOS evolutive approach not recognised." << std::endl;
         return false;
       }
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


     if (vm.count("mosql-pmax")) {
       if (vm["mosql-pmax"].as<long double>() < 0 || vm["mosql-pmax"].as<long double>() > 1 ){
         std::cerr << "[GAEDAConfig] Error: Wrong PMax value for MOSQL: " << vm["mosql-pmax"].as<long double>() << ". Please select a value within [0, 1]. " << std::endl;
         return false;
       }
       mMOSQLPmax = vm["mosql-pmax"].as<long double>();
     }


     if (vm.count("mosrl-alfa")) {
       if (vm["mosrl-alfa"].as<long double>() < 0 || vm["mosrl-alfa"].as<long double>() > 1 ){
         std::cerr << "[GAEDAConfig] Error: Wrong learning rate for MOSRL: " << vm["mosrl-alfa"].as<long double>() << ". Please select a value within [0, 1]. " << std::endl;
         return false;
       }
       mMOSRLAlfa = vm["mosrl-alfa"].as<long double>();
     }


     if (vm.count("mosrl-gamma")) {
       if (vm["mosrl-gamma"].as<long double>() < 0 || vm["mosrl-gamma"].as<long double>() > 1 ){
         std::cerr << "[GAEDAConfig] Error: Wrong discount rate for MOSRL: " << vm["mosrl-gamma"].as<long double>() << ". Please select a value within [0, 1]. " << std::endl;
         return false;
       }
       mMOSRLGamma = vm["mosrl-gamma"].as<long double>();
     }


     if (vm.count("mosrl-beta")) {
       if (vm["mosrl-beta"].as<long double>() < 0 || vm["mosrl-beta"].as<long double>() > 1 ){
         std::cerr << "[GAEDAConfig] Error: Wrong reward rate for MOSRL: " << vm["mosrl-beta"].as<long double>() << ". Please select a value within [0, 1]. " << std::endl;
         return false;
       }
       mMOSRLBeta = vm["mosrl-beta"].as<long double>();
     }


     if (vm.count("mosrl-delta")) {
       if (vm["mosrl-delta"].as<long double>() < 0 || vm["mosrl-delta"].as<long double>() > 1 ){
         std::cerr << "[GAEDAConfig] Error: Wrong Delta rate for MOSRL (HPC): " << vm["mosrl-delta"].as<long double>() << ". Please select a value within [0, 1]. " << std::endl;
         return false;
       }
       mMOSRLDelta = vm["mosrl-delta"].as<long double>();
     }


     if (vm.count("mosrl-deltal")) {
       if (vm["mosrl-deltal"].as<long double>() < 0 || vm["mosrl-deltal"].as<long double>() > 1 ){
         std::cerr << "[GAEDAConfig] Error: Wrong Delta Lose rate for MOSRL (WoLF): " << vm["mosrl-deltal"].as<long double>() << ". Please select a value within [0, 1]. " << std::endl;
         return false;
       }
       mMOSRLDeltaL = vm["mosrl-deltal"].as<long double>();
     }


     if (vm.count("mosrl-deltaw")) {
       if (vm["mosrl-deltaw"].as<long double>() < 0 || vm["mosrl-deltaw"].as<long double>() > 1 ){
         std::cerr << "[GAEDAConfig] Error: Wrong Delta Win rate for MOSRL (WoLF): " << vm["mosrl-deltaw"].as<long double>() << ". Please select a value within [0, 1]. " << std::endl;
         return false;
       }
       mMOSRLDeltaW = vm["mosrl-deltaw"].as<long double>();
     }


     if (vm.count("mosrl-policy")) {
       if (vm["mosrl-policy"].as<string>() == "pmax")
         mMOSRLPolicy = GAEDAConfig::PMAX;
       else if (vm["mosrl-policy"].as<string>() == "phc")
         mMOSRLPolicy = GAEDAConfig::PHC;
       else if (vm["mosrl-policy"].as<string>() == "wolf")
         mMOSRLPolicy = GAEDAConfig::WOLF;
       else  {
         std::cerr << "[GAEDAConfig] Error: Wrong MOSRL policy: " << vm["mosrl-policy"].as<string>() << "." << std::endl;
         return false;
       }
     }


     if (vm.count("genealogy")) {
       if (vm["genealogy"].as<string>() == "memory")
         geneal = GAEDAConfig::Memory;
       else if (vm["genealogy"].as<string>() ==  "statistic")
         { // Do genealogy and get age statistics
           geneal = GAEDAConfig::Memory;
           logstats_names.push_back(Genealogy);
         }
       else if (vm["genealogy"].as<string>() == "all")
         { // Do genealogy and get age statistics
           geneal = GAEDAConfig::All;
           logstats_names.push_back(Genealogy);
         }
       else if (vm["genealogy"].as<string>() == "full")
         { // Do genealogy and get age statistics
           geneal = GAEDAConfig::Full;
           logstats_names.push_back(Genealogy);
         }
       else if (vm["genealogy"].as<string>() == "tracer")
         geneal = GAEDAConfig::Tracer;
       else if (vm["genealogy"].as<string>() == "none")
         geneal = GAEDAConfig::NoGenealogy;
       else{
         std::cerr << "[GAEDAConfig] Error: Wrong type of genealogy selected." << std::endl;
         return false;
       }
     }


     if (vm.count("print-genealogy")) {
       if (vm["print-genealogy"].as<string>() == "yes")
         mPrintGenealogy = true;
       else if (vm["print-genealogy"].as<string>() ==  "no")
         mPrintGenealogy = false;
       else {
         std::cerr << "[GAEDAConfig] Error: Select if the genealogy must be printed (yes) or not (no)." << std::endl;
         return false;
       }
     }


     if (vm.count("quality-measure")) {
       string quality_measure = vm["quality-measure"].as<string>();
       if ( quality_measure == "fAvg"      ||
            quality_measure == "fitIncAvg" ||
            quality_measure == "dAvg"      ||
            quality_measure == "Compass"   ||
            quality_measure == "NSC"         ) {
         mQualityMeasure = quality_measure;

         if (quality_measure == "NSC" && (geneal == NoGenealogy || geneal == Tracer) ) {
           std::cerr << "[GAEDAConfig] Warning: changing genealogy to 'memory' to be able to compute NSC." << std::endl;
           geneal = Memory;
         }
       }
       else{
         std::cerr << "[GAEDAConfig] Error: Wrong quality measure selected." << std::endl;
         return false;
       }
     }


     if (vm.count("weighted-autonomic")) {
       if (vm["weighted-autonomic"].as<string>() == "yes")
         mWeightedAutonomic = true;
       else if (vm["weighted-autonomic"].as<string>() ==  "no")
         mWeightedAutonomic = false;
       else {
         std::cerr << "[GAEDAConfig] Error: Select if MOS autonomic should compute weighted (yes) or arithmetic (no) averages." << std::endl;
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


     if (vm.count("het-mos-cfg")) {
       if (vm["het-mos-cfg"].as<string>() != "")
          mHetMOSCfgFile = vm["het-mos-cfg"].as<string>();
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


     if (vm.count("mos-shared-evals-factor")) {
       if (vm["mos-shared-evals-factor"].as<long double>() < 1.0){
         std::cerr << "[GAEDAConfig] Warning: The value for the factor to compute the shared evaluations for MOS HRH is smaller than the population size: " << vm["mos-shared-evals-factor"].as<long double>() << ". Be sure that this is what you actually want. " << std::endl;
       }
       mSharedEvalsFactor = vm["mos-shared-evals-factor"].as<long double>();
       mSharedEvals       = mPopSize * mSharedEvalsFactor;
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
   catch ( const boost::program_options::multiple_occurrences& e ) {
      std::cerr << e.what() << " from option: " << e.get_option_name() << std::endl;
      return false;
   }

   catch(exception& e){
     std::cerr << "[GAEDAConfig] Exception: " << e.what() << std::endl;
     return false;
   }

   return true;

}

VNSOp* GAEDAConfig::getRoutingLocalSearch() const {
  assert(mProblemStruct->VNSLocalSearch);
  return mProblemStruct->VNSLocalSearch();
}

std::vector< VNSShaker* > GAEDAConfig::getVNSShakersList() const {
  assert(mProblemStruct->VNSShakers);
  return mProblemStruct->VNSShakers();
}

long double GAEDAConfig::getVNSTinitRatio()              const { return vnstinitratio_;              }
long double GAEDAConfig::getVNSTinitProb()               const { return vnstinitprob_;               }
long double GAEDAConfig::getVNSLSRatio1()                const { return vnslsratio1_;                }
long double GAEDAConfig::getVNSLSProb()                  const { return vnslsprob_;                  }
long double GAEDAConfig::getVNSLSRatio2()                const { return vnslsratio2_;                }
int    GAEDAConfig::getVNSMaxEvalsShaker()          const { return vnsmaxevalsshaker_;          }
bool   GAEDAConfig::getVNSUseAllShakersInEachStep() const { return vnsuseallshakersineachstep_; }
int    GAEDAConfig::getVNSMaxShakersRes()           const { return vnsmaxshakersres_;           }
bool   GAEDAConfig::getVNSUpdatePenalizations()     const { return vnsupdatepenalizations_;     }


// Show command line arguments
void GAEDAConfig::showUsage( char *prgname ) {

  std::cerr << "Usage: mpirun [MPI OPTIONS] " << prgname << " [OPTIONS] problem_module" << std::endl
            << desc
            <<  "Notes:" << std::endl
            <<  "  * If you use -c, you have to also use -C" << std::endl
            <<  "  * problem_module is a path to a shared library defining the problem" << std::endl
            <<  std::endl;
}


Recombinator* GAEDAConfig::getRecombinator() const {
  Recombinator* recomb = NULL;
  switch(recombinator_){
  case SingleElitism:
    recomb = new ClassicElitism();
    break;
  case DE_Elitism:
    recomb = new DEElitism();
    break;
  case FullElitism:
    recomb = new PopElitism();
    break;
  case Percent_Elitism:
    recomb = new PercentElitism(elitism_percentage);
    break;
  default:
    throw runtime_error("Error unexpected recombination method");
  }
  return recomb;
}


GAGenome::OptCriterion GAEDAConfig::getOptCriterion() const{
  int struct_opt_crit = mProblemStruct->optCriterion;
  return (struct_opt_crit >= 0) ? (GAGenome::OptCriterion) struct_opt_crit : GAGenome::MAXIMIZATION;
}


EAIslandsTopology* GAEDAConfig::getTopology () const {

  int count = comm_manager_.getNumIslands();
  int rank  = comm_manager_.getMyRank();

  EAIslandsTopology* top = NULL;

   switch (mTop) {

      case GAEDAConfig::Mesh:

        //return new GAIslandsTopologyMesh (count);
         throw runtime_error("Error mesh topology is currently not implemented");
         return NULL;
         break;

      case GAEDAConfig::a2a:

        top = new GAIslandsTopologyA2A (count, rank);
        break;

      case GAEDAConfig::Ring2:

        top = new GAIslandsTopologyRing2 (count, rank);
        break;

      case GAEDAConfig::Pairs:

        //top = new GAIslandsTopologyPairs (count, rank);
        throw runtime_error("Error pairs topology is currently not implemented");
        return NULL;
        break;

      case GAEDAConfig::NearestNeighbor:
        top = new GAIslandsTopologyNearestNeighbor ( count,getNeighborConditionFunc(), getMinNeighborsNum(),getMaxNeighborsNum(), rank, comm_manager_);
        break;

      case GAEDAConfig::FurthestNeighbor:
        top = new GAIslandsTopologyFurthestNeighbor ( count, getMinNeighborsNum(),getMaxNeighborsNum(), rank, comm_manager_ );
        break;

      case GAEDAConfig::RandomTopology:
        top = new GAIslandsTopologyRandom ( count, rank, getTopDegree());
        break;

      case GAEDAConfig::HyperCubeTopology:
        top = new GAIslandsTopologyHyperCube ( count, rank, getTopDegree() );
        break;

        //case GAEDAConfig::HybridTopology:
        //top = new GAIslandsTopologyHybrid( count, getNeighborConditionFunc(), comm_manager_.getMyRank(), getTopDegree() , 0.5, comm_manager_ );
        //break;

      case GAEDAConfig::Ring:
      default:
        top = new GAIslandsTopologyRing (count, rank);
        break;
   }
   return top;

}


EAEmigrantsSelector* GAEDAConfig::getEmigrantsSelector(GAGeneticAlgorithm& ea) const{
  EAEmigrantsSelector* emigrants_selector = NULL;
  switch(mEmigrantsSelector){
  case BestEmigrantsSelector:
    emigrants_selector = new EABestEmigrantsSelector( ea, getMigrationPopSize() );
    break;
  case RandomEmigrantsSelector:
    emigrants_selector = new EARandomEmigrantsSelector( ea, getMigrationPopSize() );
    break;
  }
  return emigrants_selector;
}


unsigned int GAEDAConfig::getTopDegree () const{
   if (mTopDegree < 2) throw runtime_error("Error topologyDegree is < 2");
   if (mTopDegree > comm_manager_.getNumIslands()) throw runtime_error("Error mTopDegree is > nislands");
   return mTopDegree;
}


neighbor_cond_func GAEDAConfig::getNeighborConditionFunc ()const{
  switch (mNeighbor_cond){
    case GAEDAConfig::Square:
      return squareNeighborCond;
      break;
    case GAEDAConfig::Rel_neighbors:
      return relativeNeighborCond;
      break;
    case GAEDAConfig::Gabriel:
    default:
      return gabrielNeighborCond;
      break;
  }
}


unsigned int GAEDAConfig::getMinNeighborsNum () const{
   if (mMinNeighborsNum < 1) throw runtime_error("Error min neighbors is < 1");
   if ((int) mMinNeighborsNum > comm_manager_.getNumIslands()) throw runtime_error("Error min neighbors is > nislands");
   return mMinNeighborsNum;
}


unsigned int GAEDAConfig::getMaxNeighborsNum () const{
   if (mMaxNeighborsNum < 1) throw runtime_error("Error min neighbors is < 1");
   if ((int) mMaxNeighborsNum > comm_manager_.getNumIslands()) throw runtime_error("Error max neighbors is > nislands");
   return mMaxNeighborsNum;
}


CentroidsGenerator* GAEDAConfig::getCentroidsGenerator(GAGenome& sample_gen) const{
  CentroidsGenerator* cent_gen;

  switch(mCentGen){
  case HighSeparatedRandom:
    cent_gen = new RandomHighSeparatedCentGenerator(comm_manager_.getNumIslands(), sample_gen, mProblemStruct->individualInit);
    break;
  case Random:
    cent_gen = new RandomCentGenerator(comm_manager_.getNumIslands(), sample_gen, mProblemStruct->individualInit);
    break;
  default:
    throw runtime_error("Error: centroids generator not recognized");
  }
  return cent_gen;
}


int GAEDAConfig::getIslandInit() const{
  if (isSequential()) return RandomIslandInit;
  else                return mIslandsInit;
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


GAGenealogy* GAEDAConfig::getGenealogy() const{
   int rank = comm_manager_.getMyRank();

   switch (geneal) {

      case GAEDAConfig::Memory:

	return GAGenealogyMemory::create(rank, GAGenealogyMemory::NORMAL);
	break;

      case GAEDAConfig::All:

	 return GAGenealogyMemory::create(rank, GAGenealogyMemory::ALL);
	 break;

      case GAEDAConfig::Full:

	 return GAGenealogyMemory::create(rank, GAGenealogyMemory::FULL);
	 break;

      case GAEDAConfig::Tracer:

         return GAGenealogyTracer::create(rank);
         break;

      default:

         return NULL;
         break;

   }
}


GAGenome* GAEDAConfig::getOptimum() const{
  return (mProblemStruct->optimum == NULL) ? NULL : mProblemStruct->optimum();
}


MOSParticipation* GAEDAConfig::getParticipationFunction (void) const {
  if (partFunction == constantPF)
    return new ConstantParticipation(minPart);
  else if (partFunction == dynQualityMOSPF)
    return new DynamicParticipation(0.05, minPart);
  else
    throw runtime_error("ERROR: Wrong value in partFunction");
}


MOSQuality* GAEDAConfig::getQualityFunction (void) const {
  if      (mQualityMeasure == "fAvg")      return new AverageFitnessQuality();
  else if (mQualityMeasure == "fitIncAvg") return new AverageFitnessIncrementQuality();
  else if (mQualityMeasure == "dAvg")      return new AverageDiversityQuality();
  else if (mQualityMeasure == "Compass")   return new CompassQuality();
  else throw runtime_error("ERROR: Wrong value in mQualityMeasure");
}

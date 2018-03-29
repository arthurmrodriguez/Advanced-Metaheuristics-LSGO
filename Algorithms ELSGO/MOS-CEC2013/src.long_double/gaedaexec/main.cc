#include <assert.h>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <memory>
#include <execinfo.h>
#include <signal.h>

#include <GAEDAlib.h>
#include <GAEDAConfig.h>
#include <ConfigFileParser.h>
#include <MOSDefineProblem.h>
#include <MOSConfig.h>
#include <MOSConversion.h>
#include <RLImprovPerDimManager.h>

GAGenome::Initializer getInitializer (const GAEDAConfig& cfg, CommManager& comm_manager, GAGenome& gen) {

   GAGenome::Initializer initializer = NULL;
   GAGenome::Initializer ind_init    = cfg.getIndInitFunction ();

   if (ind_init == NULL)
      throw runtime_error("Error: Problem module does not have an initializer for the genome (function individualInit)");

   auto_ptr<CentroidsGenerator> cent_gen (cfg.getCentroidsGenerator (gen)); // only neccesary for Voronoi but must be placed
                                                                            // here (crosses initialization error)

   switch (cfg.getIslandInit ()) {

      case GAEDAConfig::VoronoiIslandInit:
         if (cfg.isSequential ())
            throw runtime_error("Error: voronoi initialization called with a sequential model");
         VoronoiIndInit::setVoronoiIndInitData (comm_manager, *cent_gen, ind_init);
         initializer = VoronoiIndInit::initialize;
         break;

      case GAEDAConfig::RandomIslandInit:
      default:
         initializer = ind_init;
         break;

   }

   return initializer;

}

static GALogger* generateLogger (const GAEDAConfig&   cfg,
                                 const int            repetition,
                                 const CommManager&   comm_manager) {

   char* baseName = cfg.getLogFile ();
   char* logFile;

   if (baseName == NULL) {
     return new GANullLogger (logFile, cfg.getStatsDisplayFreq (), cfg.getLogLevel ());
   }
   else {
     if (GAEDAConfig::handle()->getRepetitions() > 1) {
        logFile = (char*) malloc (strlen (baseName) + 32);
        sprintf (logFile, "%s-r%d-%d", baseName, repetition, comm_manager.getMyRank ());
     }
     else {
        logFile = (char*) malloc (strlen (baseName) + 32);
        sprintf (logFile, "%s-%d", baseName, comm_manager.getMyRank ());
     }
   }

   return new GAFileLogger (logFile, cfg.getStatsDisplayFreq (), cfg.getLogLevel ());

}

static void setConvergenceCriterion (const GAEDAConfig& cfg, GAGeneticAlgorithm& ga, CommManager& comm) {

   const GAEDAConfig::convergencecriterion& criterion = cfg.getConvCrit ();
   int                                      conviters = 0;

   RoutingAlg* ralg;

   switch (criterion) {

      case GAEDAConfig::CONV:
         ga.setTerminator (GAGeneticAlgorithm::TerminateUponConvergence);
         ga.pConvergence  (cfg.getConv ());

         conviters = cfg.getConvIters ();

//          if (conviters > 0)
//             ga.nConvergence (conviters);
         break;

      case GAEDAConfig::POPCONV:
         ga.setTerminator (GAGeneticAlgorithm::TerminateUponPopConvergence);
         ga.pConvergence  (cfg.getConv ());

         conviters = cfg.getConvIters ();

//          if (conviters > 0)
//             ga.nConvergence (conviters);
         break;

      case GAEDAConfig::PARTIALPOPCONV:
         ga.setTerminator (GAGeneticAlgorithm::TerminateUponPartialPopConvergence);
         ga.pConvergence  (cfg.getConv ());

         break;

      case GAEDAConfig::POPCONVITERS:
         ga.setTerminator (GAGeneticAlgorithm::TerminateUponPopConvergenceOrGen);
         ga.pConvergence  (cfg.getConv  ());
         ga.nGenerations  (cfg.getIters ());
         break;

      case GAEDAConfig::POPCONVEVALS:
         ga.setTerminator (GAGeneticAlgorithm::TerminateUponPopConvergenceOrEvals);
         ga.pConvergence  (cfg.getConv  ());
         ga.nEvaluations  (cfg.getEvals ());
         break;

      case GAEDAConfig::SCOREORGEN:
         ga.setTerminator (GAGeneticAlgorithm::TerminateUponScoreOrGen);
         ga.pConvergence (cfg.getConv  ());
         ga.nGenerations (cfg.getIters ());

      case GAEDAConfig::EVALS:
         ga.setTerminator (GAGeneticAlgorithm::TerminateUponFitnessEvaluations);
         ga.nEvaluations  (cfg.getEvals ());
         break;

      case GAEDAConfig::EVALROUTECALLS:
        ralg = dynamic_cast<RoutingAlg*>(&ga); assert(ralg);
        ralg->maxEvalRouteCalls(cfg.getMaxEvalRouteCalls());
        ralg->setTerminator (RoutingAlg::TerminateUponRouteCalls);
        break;

      case GAEDAConfig::TIMECONV:
        TimeConvergence::setMaxTime(cfg.getMaxExecSeconds());
        TimeConvergence::setCommManager(&comm);
        ga.setTerminator( TimeConvergence::TerminateUponTime);
        break;

      case GAEDAConfig::ITERS:
      default:
         ga.setTerminator (GAGeneticAlgorithm::TerminateUponGeneration);
         ga.nGenerations  (cfg.getIters ());
         break;

   }

}

static GAGeneticAlgorithm* createAndConfigureGA (const GAEDAConfig& cfg,
                                                       GAGenome*    gen,
                                                       CommManager& comm_manager,
                                                       GAPopulation* oldPop,
                                                       int repetition) {

   GAGeneticAlgorithm* ga  = NULL;
//   int                 val = 0 ;

   switch (cfg.getAlgorithm ()) {

      case GAEDAConfig::GA:
         ga = new GASimpleGA (*gen);
         break;

      case GAEDAConfig::SSGA:
         ga = new GASteadyStateGA (*gen);
         break;

      case GAEDAConfig::EDA:
         ga  = new GASimpleEDA (*gen);

// TODO: Parametrizar la configuracion de los EDA (y el resto de algoritmos)
//       para poder seleccionar el metodo de simulacion adecuado.
//         val = GABayesianNetwork::PLS_CORRECT;
//         ga->setptr (gaNdiscreteSimulationMethod, &val);

         break;

      case GAEDAConfig::DE:
        //std::cerr << "Warning: Remind that you should have defined suitable F and CR factors." << std::endl;
         if (cfg.getDECrossover() != 0) gen->crossover( cfg.getDECrossover() );
         ga = new DE (*gen, cfg.getFFactor (), cfg.getCRFactor ());
         break;

      case GAEDAConfig::ES:
        ga = new EvolutionStrategy (*gen);
        break;

      case GAEDAConfig::VNS:
        ga = new VNS (*gen,
                      cfg.getVNSTinitRatio(), cfg.getVNSTinitProb(),
                      cfg.getVNSLSRatio1(), cfg.getVNSLSProb(), cfg.getVNSLSRatio2(),
                      cfg.getVNSUpdatePenalizations(),
                      cfg.getVNSMaxEvalsShaker(),
                      cfg.getVNSUseAllShakersInEachStep(),
                      cfg.getVNSMaxShakersRes(),
                      cfg.getVNSShakersList(),
                      cfg.getRoutingLocalSearch() );
        break;
      case GAEDAConfig::TSCordeau:
        ga = new TSCordeau (*gen, cfg.getRoutingLocalSearch() );
        break;

      case GAEDAConfig::MOS:

        if (cfg.getEvolutiveApproach() == CentralEvolution)
          ga = genericDefineMOSProblemCentral   (cfg.getPopSize ());
        else
          ga = genericDefineMOSProblemAutonomic (cfg.getPopSize ());

         break;

      case GAEDAConfig::MOS2:

         if (repetition == 0) {
           ga = genericDefineMOS2ProblemCentral (cfg.getPopSize ());
           std::cout << "It is first run, we create the algorithm as usually (full initialization)..." << std::endl;
         }
         else {
           ga = genericDefineMOS2ProblemCentralReusing (cfg.getPopSize(), oldPop);
           std::cout << "It is run: " << repetition << ". We create the algorithm with partial initialization..." << std::endl;
         }
         break;

      case GAEDAConfig::MOSMultiDeme:

        if (cfg.getEvolutiveApproach() == CentralEvolution)
          ga = genericDefineMOSMultiDemeProblemCentral   (cfg.getPopSize ());
        else
          throw runtime_error("Error: the MOSMultiDeme algorithm can not be executed with the Autonomic Evolution approach.");

         break;

      case GAEDAConfig::MOSRL:
         ga = genericDefineMOSProblemRL (cfg.getPopSize ());
         break;

      default:
         throw runtime_error("Error: Unexpected algorithm.");

   }

   GASelectionScheme* sel = cfg.getSelector();
   GAScalingScheme*  scal = cfg.getScaling ();

   ga->populationSize           ( cfg.getPopSize         ());
   ga->pMutation                ( cfg.getMutProb         ());
   ga->pCrossover               ( cfg.getCrossProb       ());
   ga->recombinator             ( cfg.getRecombinator    ());
   ga->setRank                  ( comm_manager.getMyRank ());
   ga->selector                 (*sel );
   ga->scaling                  (*scal);

   setConvergenceCriterion (cfg,*ga, comm_manager);

   delete sel;
   delete scal;

   return ga;

}

vector<LogStat*>* createLogStats (const GAEDAConfig& cfg, const Algorithm& alg, CommManager& comm, long double startruntime) {

   vector<GAEDAConfig::logstat_t>logsnames = cfg.getLogStatsNames ();
   vector<LogStat*>*             logstats  = new vector<LogStat*> ();
   LogStat*                      stat;

   for (unsigned i=0; i<logsnames.size(); i++){

      if      (logsnames[i] == GAEDAConfig::Age)                           stat = new AgeLogStat                           (alg);
      else if (logsnames[i] == GAEDAConfig::DistAvgPoint)                  stat = new DistAvgPoint                         (alg);
      else if (logsnames[i] == GAEDAConfig::Entropy)                       stat = new EntropyLogStat                       (alg);
      else if (logsnames[i] == GAEDAConfig::Evals)                         stat = new EvalsLogStat                         (alg);
      else if (logsnames[i] == GAEDAConfig::Fitness)                       stat = new FitnessLogStat                       (alg);
      else if (logsnames[i] == GAEDAConfig::FitInc)                        stat = new FitIncLogStat                        (alg);
      else if (logsnames[i] == GAEDAConfig::FSSAvg)                        stat = new FSSLogStat                           (alg, FSSLogStat::AVG);
      else if (logsnames[i] == GAEDAConfig::FSSBest)                       stat = new FSSLogStat                           (alg, FSSLogStat::BEST);
      else if (logsnames[i] == GAEDAConfig::FSSWorst)                      stat = new FSSLogStat                           (alg, FSSLogStat::WORST);
      else if (logsnames[i] == GAEDAConfig::Improve)                       stat = new ImprovementsLogStat                  (alg);
      else if (logsnames[i] == GAEDAConfig::Genealogy)                     stat = new GenealogyLogStat                     (alg);
      else if (logsnames[i] == GAEDAConfig::GrefBias)                      stat = new GrefenstetteLogStat                  (alg);
      else if (logsnames[i] == GAEDAConfig::NativePrcnt)                   stat = new NativePrcntLogStat                   (alg, comm.getMyRank());
      else if (logsnames[i] == GAEDAConfig::NumNewChild)                   stat = new NumNewChildrenLogStat                (alg);
      else if (logsnames[i] == GAEDAConfig::OptNcomps)                     stat = new OptNCompsLogStat                     (alg, cfg.getOptimum (), cfg.getNComponentsOpt ());
      else if (logsnames[i] == GAEDAConfig::OptDist)                       stat = new OptDistLogStat                       (alg, cfg.getOptimum (), cfg.getDistToOptFunc  ());
      else if (logsnames[i] == GAEDAConfig::Participation)                 stat = new ParticipationLogStat                 (alg);
      else if (logsnames[i] == GAEDAConfig::Quality)                       stat = new QualityLogStat                       (alg);
      else if (logsnames[i] == GAEDAConfig::Score)                         stat = new ScoreLogStat                         (alg);
      else if (logsnames[i] == GAEDAConfig::RoutingAlgScore)               stat = new RoutingAlgScoreLogStat               (alg);
      else if (logsnames[i] == GAEDAConfig::VNSSuccessfulShaker)           stat = new VNSSuccessfulShakerLogStat           (alg);
      else if (logsnames[i] == GAEDAConfig::VNSShakerpos)                  stat = new VNSShakerPosLogStat                  (alg);
      else if (logsnames[i] == GAEDAConfig::RoutingAlgBestSolNonPenScore)  stat = new RoutingAlgBestSolNonPenScore         (alg);
      else if (logsnames[i] == GAEDAConfig::RoutingAlgBestSolPickupTime)   stat = new RoutingAlgBestSolPickupTimeLogStat   (alg);
      else if (logsnames[i] == GAEDAConfig::RoutingAlgBestSolDeliveryTime) stat = new RoutingAlgBestSolDeliveryTimeLogStat (alg);
      else if (logsnames[i] == GAEDAConfig::RoutingAlgBestSolTWV)          stat = new RoutingAlgBestSolTWV                 (alg);
      else if (logsnames[i] == GAEDAConfig::RoutingAlgBestSolRideV)        stat = new RoutingAlgBestSolRideV               (alg);
      else if (logsnames[i] == GAEDAConfig::RoutingAlgBestSolLoadV)        stat = new RoutingAlgBestSolLoadV               (alg);
      else if (logsnames[i] == GAEDAConfig::EvalRouteCalls)                stat = new EvalRouteCallsLogStat                (alg);
      else if (logsnames[i] == GAEDAConfig::Runtime)                       stat = new RuntimeLogStat                       (alg,comm,startruntime);
      else if (logsnames[i] == GAEDAConfig::Freemem)                       stat = new FreememLogStat                       (alg);

/*
      else if (logsnames[i] == GAEDAConfig::AvgDivIncDynCoevLogStat) stat = new AvgDivIncDynCoevLogStat (alg);
      else if (logsnames[i] == GAEDAConfig::CoevPhaseLogStat)        stat = new CoevPhaseLogStat        (alg);
      else if (logsnames[i] == GAEDAConfig::Epistasis)               stat = new EpistasisVarLogStat     (alg);
      else if (logsnames[i] == GAEDAConfig::GenDiv)                  stat = new GenDivLogStat           (alg);
      else if (logsnames[i] == GAEDAConfig::InfDivDynCoevLogStat)    stat = new InfDivDynCoevLogStat    (alg);
      else if (logsnames[i] == GAEDAConfig::MaxDivDynCoevLogStat)    stat = new MaxDivDynCoevLogStat    (alg);
      else if (logsnames[i] == GAEDAConfig::SupDivDynCoevLogStat)    stat = new SupDivDynCoevLogStat    (alg);
*/
      else                                                           throw runtime_error ("Error: Unrecognized statistic to log.");

      logstats->push_back (stat);

   }

   return logstats;

}

static EAIslandsModel* createAndInitializeIslandsModel (const GAEDAConfig&        cfg,
                                                              GAGeneticAlgorithm* ga,
                                                              GAGenome&           gen,
                                                              CommManager&        comm_manager) {

   EAIslandsTopology*    top                 = cfg.getTopology ();
   EAEmigrantsSelector*  emigrants_selector  = cfg.getEmigrantsSelector (*ga);
   EAImmigrantsSelector* immigrants_selector = new EABestImmigrantsSelector (*ga);
   EAIslandsModel*       im                  = NULL;

   switch (cfg.getModel ()) {

      case GAEDAConfig::SYNC:
         im = new EAIslandsModelSync (ga,
                                      top,
                                      cfg.getMigrationFrequency (),
                                      emigrants_selector,
                                      immigrants_selector,
                                      comm_manager);
         break;

      case GAEDAConfig::ASYNC:
         throw runtime_error("Asyn not yet implemented");
         break;

//      case GAEDAConfig::CLUST:
//         im = new EAIslandsModelClusteredPops (ga,
//                                               top,
//                                               cfg.getMigrationFrequency (),
//                                               emigrants_selector,
//                                               immigrants_selector,
//                                               cfg.getClusteringPeriod (),
//                                               EAIslandsModelClusteredPops::kmeansCluster,
//                                               EAIslandsModelClusteredPops::resizeWithBestIndMutation,
//                                               cfg.getCentroidsGenerator (gen),
//                                               comm_manager);
//         break;

      default:
         throw runtime_error ("Error: Communication model not recognized.");

   }

   return im;

}

static Algorithm* createRoutingDistModel(const GAEDAConfig&  cfg,
                                         GAGeneticAlgorithm& ga,
                                         GAGenome&           gen,
                                         CommManager&        comm_manager) {
  RoutingAlg& ralg = dynamic_cast<RoutingAlg&>(ga); assert(&ralg);

  RoutingGenome& rgen = dynamic_cast<RoutingGenome&>(gen); assert(&rgen);

/*
  if ( cfg.getRouteDistributorName() == "random") {
    rdist = new RandomRouteDistributor(rgen,comm_manager.getNumIslands());
  } else{
    throw runtime_error("Unrecognized Route Distributor");
  }
  */

  Algorithm* alg = new DistRoutingAlg(ralg,
                                      cfg.getMigrationFrequency (),
                                      cfg.getEmigrantRouteType(),
                                      gen,
                                      comm_manager);

  return alg;
}

static Algorithm* createAlgorithm (const GAEDAConfig& cfg,
                                         GAGenome*    gen,
                                         CommManager& comm_manager,
                                         GAPopulation* oldPop,
                                         int repetition) {

   if (cfg.getAlgorithm () != GAEDAConfig::MOS          &&
       cfg.getAlgorithm () != GAEDAConfig::MOS2         &&
       cfg.getAlgorithm () != GAEDAConfig::MOSMultiDeme &&
       cfg.getAlgorithm () != GAEDAConfig::MOSRL           ) // MOS no usa centroides (por ahora ...)
      gen->initializer (getInitializer (cfg, comm_manager, *gen));

   GAGeneticAlgorithm* ga = createAndConfigureGA (cfg, gen, comm_manager, oldPop, repetition);

   switch (cfg.getModel()) {
   case GAEDAConfig::SEQ:
     return ga;
   case GAEDAConfig::ROUTINGDISTSYNC:
     return createRoutingDistModel(cfg,*ga,*gen,comm_manager);
   default:    // island model
     return createAndInitializeIslandsModel (cfg, ga, *gen, comm_manager);
   }

}

static void runAndDisplayResults (Algorithm& alg, GALogger& logger, CommManager& comm_manager, long double startruntime, GAEDAConfig& cfg) {
   alg.run();
   long double runtime = comm_manager.getTime () - startruntime;
   auto_ptr< GAGenome > best (alg.best ());

   if (comm_manager.isIslandMaster ()) {

      logger.appendExecTime (runtime);

      long double score = best->score ();

      std::cout << setiosflags (ios::fixed | ios::showpoint) << setprecision (14);
      std::cout << "-> RES:   " << score   << std::endl;
      std::cout << "-> TIME:  " << runtime << std::endl;
      if (cfg.getPrintBestSol() ) {
        std::cout << "-> BEST:  " << *best   << std::endl;

        if (GAGenealogy::handle () != NULL) {
           if (GAGenealogy::handle()->isGenealogyMemory() &&
               GAEDAConfig::handle()->printGenealogy   ()    ) {
              stringstream mes;
              mes << "\n#Best Genome\n" << best->getId () << " " << best->getIsland () << STD_ENDL;
              mes << *best;
              GALogger::instance ()->appendLogMessage ("", mes.str (), GALogger::only_stats);
           }
        }
      }
   }

}

void printTrace() {
  const int maxSize = 50;
  void *array[maxSize];

  // get void*'s for all entries on the stack
  size_t read_size = backtrace(array,maxSize);

  // print out all the frames to stderr
  backtrace_symbols_fd(array, read_size, 2);
}

void handler(int sig) {
  cerr << "Error: captured signal " << sig << endl;
  printTrace();
  exit(-1);
}

int main (int argc, char **argv) {

  signal(SIGSEGV, handler);   // install our handler

  CommManager comm_manager (argc, argv);

  long double startruntime = comm_manager.getTime ();
  TimeConvergence::setInitTime(startruntime);

  try {

    GAEDAConfig* cfg = GAEDAConfig::handle (comm_manager);

    if (!cfg->parse (argc, argv)) {
       if (comm_manager.isIslandMaster ()) {
          std::cerr << "Try: " << argv[0] << " --help for more information." << std::endl;
       }
       GAEDAConfig::destroy();
       return -1;
     }

    if (cfg->getDebug ()) {
       bool cond = true;
       while (cond);
    }

    GARandomSeed (cfg->getSeed(), comm_manager.getMyRank());

    GAProblemStruct* prob   = cfg->getProblemStruct ();
    GAGenome*        gen    = NULL;
    Algorithm*       alg    = NULL;
    GAGenealogy*     geneal = NULL;

    RLImprovPerDimManager impmanager (cfg->getProblemSize());

    GAPopulation* oldPop = (GAPopulation*) 0;
    std::map<MOSTechnique*, unsigned> oldTechsMap;

    for (int repetition = 0; repetition < cfg->getRepetitions (); repetition++) {

       gen = prob->defineProblem ();

       if (comm_manager.isIslandMaster () && repetition == 0) {
          std::cout << "Trying to solve: " << prob->describeProblem () << std::endl;
       }
       std::cout << "Node: "<< comm_manager.getMyRank() << " Seed: " << GAGetRandomSeed() << std::endl;

       if (cfg->getAlgorithm () == GAEDAConfig::MOS          ||
           cfg->getAlgorithm () == GAEDAConfig::MOS2         ||
           cfg->getAlgorithm () == GAEDAConfig::MOSMultiDeme ||
           cfg->getAlgorithm () == GAEDAConfig::MOSRL           ) {
          if (!cfg->getTechsRepository ()) {
             std::cerr << "Error: configuration file not defined for techniques." << std::endl;
             return -1;
          }

          if (comm_manager.isIslandMaster () && repetition == 0)
             std::cout << cfg->getTechsRepository() << std::endl;

          ConfigFileParser techParser (cfg->getTechsRepository(), repetition==0);
          geneal = cfg->getGenealogy ();
       }

       // Re-assign the appropriate techniques (based on the oldTechsMap) to the old population
       if (oldPop) {
         for (unsigned k = 0; k < oldPop->size(); k++) {
           MOSGenome& ind = dynamic_cast<MOSGenome&>(oldPop->individual(k));
           unsigned techId = oldTechsMap[ind.getTechnique()];
           ind.setTechnique(MOSTechniqueSet::handle()->getTechnique(techId));
         }
       }

       // Set optimization criterion
       GAGenome::optCriterion (cfg->getOptCriterion());

       alg = createAlgorithm (*cfg, gen, comm_manager, oldPop, repetition);
       comm_manager.setGenome (&alg->population ().individual (0));

       auto_ptr<GALogger> logger (generateLogger (*cfg, repetition, comm_manager));

       const Algorithm* logstatsalg = alg;

       if (cfg->getModel() == GAEDAConfig::SYNC) {
         EAIslandsModel* im = DYN_CAST (EAIslandsModel*, alg);
         logstatsalg = im->getEA();
       } else if (cfg->getModel() == GAEDAConfig::ROUTINGDISTSYNC) {
         DistRoutingAlg* ralg = dynamic_cast<DistRoutingAlg*>(alg);
         logstatsalg = ralg->getRoutingAlg();
       }

       logger->setLogStats ( createLogStats (*cfg, *logstatsalg, comm_manager, startruntime) );

       // As it can be seen this configuration of the logger needs to be done after we have the algorithm created
       logger->setStats (const_cast<GAStatistics*> (&(alg->statistics ())));

       if (prob->configureAlg) prob->configureAlg(*alg);

       runAndDisplayResults (*alg, *logger, comm_manager, startruntime,*cfg);

       if (oldPop) delete oldPop;
       oldPop = alg->population().clone();

       // Store the conversion from technique pointer to technique ID in a map
       // in order to be able to re-assign the newly created techniques to the
       // old population in the following iteration of the main loop
       for (MOSTechniqueSet::MOSTechniqueSetIterator it = MOSTechniqueSet::handle()->begin(); it != MOSTechniqueSet::handle()->end(); it++)
         oldTechsMap[it->second] = it->second->getId();

       // Check if some postprocessing is needed and, if so, do it.
       if (prob->postprocess)
          prob->postprocess (oldPop, comm_manager.getMyRank (), repetition);

       if (alg   ) delete alg;
       if (geneal) delete geneal;
       if (gen   ) delete gen;

       MOSConversion::destroy();
       MOSTechniqueSet::destroy();
       MOSGenomeFactory::destroy();

    } // End FOR

    free (prob);

    if (oldPop) delete oldPop;

    // Destroy Singleton objects
    GAEDAConfig::destroy();

  } // End Try
  catch (runtime_error& e) {
    std::cerr << "The following error has occured: " << e.what () << std::endl;
    printTrace();
    exit(-1);
  }
  catch (exception& e) {
    std::cerr << "The following exception has occured: " << e.what () << std::endl;
    printTrace();
    exit(-1);
  }
  catch (...) {
    std::cerr << "Unknown type exception!" << std::endl;
    exit(-1);
  }

  return 0;
}

#include <assert.h>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <memory>

#include <GAEDAlib.h>
#include <GAEDAConfig.h>
#include <ConfigFileParser.h>
#include <MOSDefineProblem.h>
#include <MOSConfig.h>
#include <MOSConversion.h>


static GALogger* generateLogger (const GAEDAConfig&   cfg,
                                 const CommManager&   comm_manager) {

   char* baseName = cfg.getLogFile ();
   char* logFile;

   if (baseName == NULL) {
     return new GANullLogger (logFile, cfg.getStatsDisplayFreq (), cfg.getLogLevel ());
   }
   else {
     logFile = (char*) malloc (strlen (baseName) + 4);
     sprintf (logFile, "%s-%d", baseName, comm_manager.getMyRank ());
   }

   return new GAFileLogger (logFile, cfg.getStatsDisplayFreq (), cfg.getLogLevel ());

}

static void setConvergenceCriterion (const GAEDAConfig& cfg, GAGeneticAlgorithm& ga) {

   const GAEDAConfig::convergencecriterion& criterion = cfg.getConvCrit ();
   int                                      conviters = 0;

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

      case GAEDAConfig::SCOREORGEN:
         ga.setTerminator (GAGeneticAlgorithm::TerminateUponScoreOrGen);
         ga.pConvergence (cfg.getConv  ());
         ga.nGenerations (cfg.getIters ());

      case GAEDAConfig::EVALS:
         ga.setTerminator (GAGeneticAlgorithm::TerminateUponFitnessEvaluations);
         ga.nEvaluations  (cfg.getEvals ());
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
                                                       CommManager& comm_manager) {

   GAGeneticAlgorithm* ga  = NULL;

   switch (cfg.getAlgorithm ()) {

      case GAEDAConfig::MOS:

         ga = genericDefineMOSProblemCentral (cfg.getPopSize ());
         break;

      default:
         throw runtime_error("Error: Unexpected algorithm.");

   }

   GASelectionScheme* sel = cfg.getSelector();
   GAScalingScheme*  scal = cfg.getScaling ();

   ga->populationSize           ( cfg.getPopSize         ());
   ga->pMutation                ( cfg.getMutProb         ());
   ga->pCrossover               ( cfg.getCrossProb       ());
   ga->setRank                  ( comm_manager.getMyRank ());
   ga->selector                 (*sel );
   ga->scaling                  (*scal);

   setConvergenceCriterion (cfg,*ga);

   delete sel;
   delete scal;

   return ga;

}

vector<LogStat*>* createLogStats (const GAEDAConfig& cfg, const Algorithm& alg, unsigned rank) {

   vector<GAEDAConfig::logstat_t>logsnames = cfg.getLogStatsNames ();
   vector<LogStat*>*             logstats  = new vector<LogStat*> ();
   LogStat*                      stat;

   for (unsigned i=0; i<logsnames.size(); i++){

           if (logsnames[i] == GAEDAConfig::Evals)                   stat = new EvalsLogStat            (alg);
      else if (logsnames[i] == GAEDAConfig::Fitness)                 stat = new FitnessLogStat          (alg);
      else if (logsnames[i] == GAEDAConfig::FitInc)                  stat = new FitIncLogStat           (alg);
      else if (logsnames[i] == GAEDAConfig::Improve)                 stat = new ImprovementsLogStat     (alg);
      else if (logsnames[i] == GAEDAConfig::Participation)           stat = new ParticipationLogStat    (alg);
      else if (logsnames[i] == GAEDAConfig::Quality)                 stat = new QualityLogStat          (alg);
      else if (logsnames[i] == GAEDAConfig::QualityFuncActive)       stat = new QualityFuncActiveLogStat(alg);
      else if (logsnames[i] == GAEDAConfig::Score)                   stat = new ScoreLogStat            (alg);
      else                                                           throw runtime_error ("Error: Unrecognized statistic to log.");

      logstats->push_back (stat);

   }

   return logstats;

}


static void runAndDisplayResults (Algorithm& alg, GALogger& logger, CommManager& comm_manager) {

   long double runtime = comm_manager.getTime ();
   alg.run ();
   runtime = comm_manager.getTime () - runtime;
   auto_ptr< GAGenome > best (alg.best ());

   if (comm_manager.isIslandMaster ()) {

      logger.appendExecTime (runtime);

      long double score = best->score ();

      std::cout << setiosflags (ios::fixed | ios::showpoint) << setprecision (20);
      std::cout << "-> RES:   " << score   << std::endl;
      std::cout << "-> TIME:  " << runtime << std::endl;
      std::cout << "-> BEST:  " << *best   << std::endl;

   }

}

int main (int argc, char **argv) {

   CommManager comm_manager (argc, argv);

   try {

      GAEDAConfig* cfg = GAEDAConfig::handle (comm_manager);

      if (!cfg->parse (argc, argv)) {
         if (comm_manager.isIslandMaster ()) {
            std::cerr << "Try: " << argv[0] << " --help for more information." << std::endl;
            std::cout << "Seed: " << GAGetRandomSeed() << std::endl;
         }
         GAEDAConfig::destroy();
         return -1;
       }

      if (cfg->getDebug ()) {
         bool cond = true;
         while (cond);
      }

      GARandomSeed (cfg->getSeed ());

      GAProblemStruct* prob   = cfg->getProblemStruct ();
      GAGenome*        gen    = NULL;
      Algorithm*       alg    = NULL;

      gen = prob->defineProblem ();

      if (comm_manager.isIslandMaster ()) {
         std::cout << "Trying to solve: " << prob->describeProblem () << std::endl;
         std::cout << "Seed: " << GAGetRandomSeed() << std::endl;
      }

      if (cfg->getAlgorithm () == GAEDAConfig::MOS) {
         if (!cfg->getTechsRepository ()) {
            std::cerr << "Error: configuration file not defined for techniques." << std::endl;
            return -1;
         }

         if (comm_manager.isIslandMaster ())
            std::cout << cfg->getTechsRepository() << std::endl;

         ConfigFileParser techParser ( cfg->getTechsRepository() );
      }

      // Set optimization criterion
      GAGenome::optCriterion (cfg->getOptCriterion());

      alg = createAndConfigureGA (*cfg, gen, comm_manager);
      comm_manager.setGenome (&alg->population ().individual (0));

      auto_ptr<GALogger> logger (generateLogger (*cfg, comm_manager));

      logger->setLogStats (createLogStats (*cfg, *alg, comm_manager.getMyRank ()));

      // As it can be seen this configuration of the logger needs to be done after we have the algorithm created
      logger->setStats (const_cast<GAStatistics*> (&(alg->statistics ())));

      runAndDisplayResults (*alg, *logger, comm_manager);

      // Check if some postprocessing is needed and, if so, do it.
      if (prob->postprocess)
         prob->postprocess ((GAPopulation*) &(alg->population ()), comm_manager.getMyRank ());

      if (alg) delete alg;
      if (gen) delete gen;

      free (prob);

      // Destroy Singleton objects
      GAEDAConfig::destroy();
      MOSConversion::destroy();
      MOSTechniqueSet::destroy();
      MOSGenomeFactory::destroy();

   } // End Try
   catch (runtime_error& e) {
      std::cerr << "The following error has occured: " << e.what () << std::endl;
      return -1;
   }
   catch (exception& e) {
      std::cerr << "The following exception has occured: " << e.what () << std::endl;
      return -1;
   }
   catch (...) {
      std::cerr << "Unknown type exception!" << std::endl;
      return -1;
   }

   return 0;

}

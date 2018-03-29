#include "mos_algorithm.h"

#include <assert.h>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <memory>

#include <genomes/GA1DArrayGenome.h>
#include <GAEDAlib.h>
#include <GAEDAConfig.h>
#include <ConfigFileParser.h>
#include <MOSDefineProblem.h>
#include <MOSConfig.h>
#include <MOSConversion.h>
#include <RLImprovPerDimManager.h>

class GAGenome;
class GAPopulation;

long double MIN_ALLELE_VALUE = -5;
long double MAX_ALLELE_VALUE =  5;

long double (*fitnessFunc)(long double*)  = NULL;

extern "C" long double objective (GAGenome& g) {
  GA1DArrayAlleleGenome<long double>& genome = DYN_CAST (GA1DArrayAlleleGenome<long double>&, g);
  return fitnessFunc(genome);
}

GAGenome* defineProblem () {

   GAEDAConfig* cfg = GAEDAConfig::handle();

   GAAlleleSet<long double> alleles (MIN_ALLELE_VALUE, MAX_ALLELE_VALUE);
   GA1DArrayAlleleGenome<long double>* genome = new GA1DArrayAlleleGenome<long double> (cfg->getProblemSize (), alleles, objective);

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

bool postprocess (GAPopulation* pop, int rank, int repetition) {
  return true;
}

GAGenome::OptCriterion optCriterion() {
  return GAGenome::MINIMIZATION;
}

void individualInit (GAGenome& g) {
   GA1DArrayAlleleGenome<long double>& gen = DYN_CAST (GA1DArrayAlleleGenome<long double>&, g);

   gen.resize (GAGenome::ANY_SIZE);

   for (int i = 0; i < gen.length (); i++)
      gen.gene (i, GARandomDouble (-4, 4));

//  return RealUniformInitializer (g);
}

GAGenome::Initializer getInitializer (const GAEDAConfig& cfg, CommManager& comm_manager, GAGenome& gen) {

   GAGenome::Initializer initializer = NULL;
   GAGenome::Initializer ind_init    = individualInit;

   if (ind_init == NULL)
      throw runtime_error("Error: Problem module does not have an initializer for the genome (function individualInit)");

/*
   auto_ptr<CentroidsGenerator> cent_gen (cfg.getCentroidsGenerator (gen)); // only neccesary for Voronoi but must be placed
                                                                            // here (crosses initialization error)
*/
   switch (cfg.getIslandInit ()) {

/*
      case GAEDAConfig::VoronoiIslandInit:
         if (cfg.isSequential ())
            throw runtime_error("Error: voronoi initialization called with a sequential model");
         VoronoiIndInit::setVoronoiIndInitData (comm_manager, *cent_gen, ind_init);
         initializer = VoronoiIndInit::initialize;
         break;
*/

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
        logFile = (char*) malloc (strlen (baseName) + 6);
        sprintf (logFile, "%s-r%d-%d", baseName, repetition, comm_manager.getMyRank ());
     }
     else {
        logFile = (char*) malloc (strlen (baseName) + 4);
        sprintf (logFile, "%s-%d", baseName, comm_manager.getMyRank ());
     }
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
         ga = new DE (*gen, cfg.getFFactor (), cfg.getCRFactor ());
         break;

      case GAEDAConfig::ES:
        ga = new EvolutionStrategy (*gen);
        break;

      case GAEDAConfig::MOS:

        if (cfg.getEvolutiveApproach() == CentralEvolution)
          ga = genericDefineMOSProblemCentral   (cfg.getPopSize ());
        else
          ga = genericDefineMOSProblemAutonomic (cfg.getPopSize ());

         break;

      case GAEDAConfig::MOS2:

         ga = genericDefineMOS2ProblemCentral (cfg.getPopSize ());
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

      if      (logsnames[i] == GAEDAConfig::Age)                     stat = new AgeLogStat              (alg);
      else if (logsnames[i] == GAEDAConfig::DistAvgPoint)            stat = new DistAvgPoint            (alg);
      else if (logsnames[i] == GAEDAConfig::Entropy)                 stat = new EntropyLogStat          (alg);
      else if (logsnames[i] == GAEDAConfig::Evals)                   stat = new EvalsLogStat            (alg);
      else if (logsnames[i] == GAEDAConfig::Fitness)                 stat = new FitnessLogStat          (alg);
      else if (logsnames[i] == GAEDAConfig::FitInc)                  stat = new FitIncLogStat           (alg);
      else if (logsnames[i] == GAEDAConfig::FSSAvg)                  stat = new FSSLogStat              (alg, FSSLogStat::AVG);
      else if (logsnames[i] == GAEDAConfig::FSSBest)                 stat = new FSSLogStat              (alg, FSSLogStat::BEST);
      else if (logsnames[i] == GAEDAConfig::FSSWorst)                stat = new FSSLogStat              (alg, FSSLogStat::WORST);
      else if (logsnames[i] == GAEDAConfig::Improve)                 stat = new ImprovementsLogStat     (alg);
      else if (logsnames[i] == GAEDAConfig::Genealogy)               stat = new GenealogyLogStat        (alg);
      else if (logsnames[i] == GAEDAConfig::GrefBias)                stat = new GrefenstetteLogStat     (alg);
      else if (logsnames[i] == GAEDAConfig::NativePrcnt)             stat = new NativePrcntLogStat      (alg, rank);
      else if (logsnames[i] == GAEDAConfig::NumNewChild)             stat = new NumNewChildrenLogStat   (alg);
      else if (logsnames[i] == GAEDAConfig::OptNcomps)               stat = new OptNCompsLogStat        (alg, cfg.getOptimum (), cfg.getNComponentsOpt ());
      else if (logsnames[i] == GAEDAConfig::OptDist)                 stat = new OptDistLogStat          (alg, cfg.getOptimum (), cfg.getDistToOptFunc  ());
      else if (logsnames[i] == GAEDAConfig::Participation)           stat = new ParticipationLogStat    (alg);
      else if (logsnames[i] == GAEDAConfig::Quality)                 stat = new QualityLogStat          (alg);
      else if (logsnames[i] == GAEDAConfig::Score)                   stat = new ScoreLogStat            (alg);
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

static Algorithm* createAlgorithm (const GAEDAConfig& cfg,
                                         GAGenome*    gen,
                                         CommManager& comm_manager) {

   if (cfg.getAlgorithm () != GAEDAConfig::MOS          &&
       cfg.getAlgorithm () != GAEDAConfig::MOS2         &&
       cfg.getAlgorithm () != GAEDAConfig::MOSMultiDeme &&
       cfg.getAlgorithm () != GAEDAConfig::MOSRL           ) // MOS no usa centroides (por ahora ...)
      gen->initializer (getInitializer (cfg, comm_manager, *gen));

   GAGeneticAlgorithm* ga = createAndConfigureGA (cfg, gen, comm_manager);

   if (cfg.isSequential ())
      return ga;
   else
      return createAndInitializeIslandsModel (cfg, ga, *gen, comm_manager);

}

static void runAndDisplayResults (Algorithm& alg, GALogger& logger, CommManager& comm_manager) {

//   long double runtime = comm_manager.getTime ();
   alg.run ();
/*
   runtime = comm_manager.getTime () - runtime;
   auto_ptr< GAGenome > best (alg.best ());

   if (comm_manager.isIslandMaster ()) {

      logger.appendExecTime (runtime);

      long double score = best->score ();

      std::cout << setiosflags (ios::fixed | ios::showpoint) << setprecision (14);
      std::cout << "-> RES:   " << score   << std::endl;
      std::cout << "-> TIME:  " << runtime << std::endl;
      std::cout << "-> BEST:  " << *best   << std::endl;

      if (GAGenealogy::handle () != NULL)
         if (GAGenealogy::handle()->isGenealogyMemory() &&
             GAEDAConfig::handle()->printGenealogy   ()    ) {
            stringstream mes;
            mes << "\n#Best Genome\n" << best->getId () << " " << best->getIsland () << STD_ENDL;
            mes << *best;
            GALogger::instance ()->appendLogMessage ("", mes.str (), GALogger::only_stats);
         }

   }
*/
}

int mos_optimizer (long double(*fitnessfunction)(long double*), CommManager& comm_manager) {

   GAEDAConfig* cfg = GAEDAConfig::handle();

   GARandomSeed (cfg->getSeed());

   try {

      if (cfg->getDebug ()) {
         bool cond = true;
         while (cond);
      }

      GARandomSeed (GAGetRandomSeed() + 1);

      GAGenome*        gen    = NULL;
      Algorithm*       alg    = NULL;
      GAGenealogy*     geneal = NULL;

      RLImprovPerDimManager impmanager (cfg->getProblemSize());

      for (int repetition = 0; repetition < cfg->getRepetitions (); repetition++) {

         fitnessFunc = fitnessfunction;
         gen = defineProblem ();

         if (comm_manager.isIslandMaster () && !GAEDAConfig::handle()->quiet()) {
            std::cout << "MOS Seed: " << GAGetRandomSeed() << std::endl;
         }

         if (cfg->getAlgorithm () == GAEDAConfig::MOS          ||
             cfg->getAlgorithm () == GAEDAConfig::MOS2         ||
             cfg->getAlgorithm () == GAEDAConfig::MOSMultiDeme ||
             cfg->getAlgorithm () == GAEDAConfig::MOSRL           ) {
            if (!cfg->getTechsRepository ()) {
               std::cerr << "Error: configuration file not defined for techniques." << std::endl;
               return -1;
            }

            if (comm_manager.isIslandMaster () && !GAEDAConfig::handle()->quiet())
               std::cout << cfg->getTechsRepository() << std::endl;

            ConfigFileParser techParser (cfg->getTechsRepository(), repetition==0);
            geneal = cfg->getGenealogy ();
         }

         // Set optimization criterion
         GAGenome::optCriterion (optCriterion());

         alg = createAlgorithm (*cfg, gen, comm_manager);
         comm_manager.setGenome (&alg->population ().individual (0));

         auto_ptr<GALogger> logger (generateLogger (*cfg, repetition, comm_manager));

         if (cfg->isSequential ())
            logger->setLogStats (createLogStats (*cfg, *alg, comm_manager.getMyRank ()));
         else {
            EAIslandsModel* im = DYN_CAST (EAIslandsModel*, alg);
            logger->setLogStats (createLogStats (*cfg, *(im->getEA ()), comm_manager.getMyRank ()));
         }

         // As it can be seen this configuration of the logger needs to be done after we have the algorithm created
         logger->setStats (const_cast<GAStatistics*> (&(alg->statistics ())));

         runAndDisplayResults (*alg, *logger, comm_manager);

         // Check if some postprocessing is needed and, if so, do it.
         if (postprocess)
            postprocess ((GAPopulation*) &(alg->population ()), comm_manager.getMyRank (), repetition);

         if (alg   ) delete alg;
         if (geneal) delete geneal;
         if (gen   ) delete gen;

      } // End FOR

      // Destroy Singleton objects
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

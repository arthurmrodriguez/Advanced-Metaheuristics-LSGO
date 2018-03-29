#include <dlfcn.h>

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#include <libconfig.h++>

#include "gaid.h"
#include "GABayesianNetwork.h"
#include "GAGaussianNetwork.h"
#include "GARealOps.h"
#include "genomes/GAGenome.h"
#include "MOSTechnique.h"
#include "MOSTechniqueGA.h"
#include "MOSTechniqueMemeticGA.h"
#include "MOSTechniqueEDA.h"
#include "MOSTechniqueDE.h"
#include "MOSTechniqueDEPopSize.h"
#include "MOSTechniqueSTSDE.h"
#include "MOSTechniqueAdapDE.h"
#include "MOSTechniqueAdapSTSDE.h"
//#include "MOSTechniqueCTB2DE.h"
#include "MOSTechniqueCMAES.h"
#include "MOSTechniqueES.h"
#include "MOSTechniqueLS.h"
#include "MOSTechniqueLSObj.h"
#include "MOSTechniqueNMA.h"
#include "MOSTechniquePopLS.h"
#include "MOSTechniqueSet.h"

#include "islands/CommManager.h" //TODO quitar

using namespace std;

class ConfigFileParser {

   public:

   ConfigFileParser (const char* cfgFile, bool printFirstTime) : mCfgFile (cfgFile) {

      // ////////////////////////// //
      // Initialize class variables //
      // ////////////////////////// //
      mID = 0;
      buildMappingTable ();

      mGenFactory = MOSGenomeFactory::handle ();
      mTechSet = MOSTechniqueSet::handle ();

      mHandle = dlopen (NULL, RTLD_NOW | RTLD_GLOBAL);

      if (mHandle == NULL) {
        std::cerr << "[ConfigFileParser] Error: Shared library not found." << std::endl;
        std::cerr << dlerror () << std::endl;
        exit (1);
      }

      // ///////////////////// //
      // Parse techniques file //
      // ///////////////////// //
      libconfig::Config cfg;

      // To avoid execeptions when casting directly from integer to float
      cfg.setAutoConvert (true);

      try {
         cfg.readFile (cfgFile);
      }
      catch (libconfig::ParseException e) {
         std::cerr << "[ConfigFileParser] Error: parsing the configuration file '" << mCfgFile
                   << "' in line (" << e.getLine () << "): " << e.getError() << "." << std::endl;
         exit (-1);
      }
      catch (libconfig::FileIOException e) {
         std::cerr << "[ConfigFileParser] Error: configuration file not found or could not be open '"
                   << mCfgFile << "'." << std::endl;
         exit (-1);
      }

      if (GAEDAConfig::handle()->getHetMOSCfgFile() == "")
         parseHomogeneousMOS(cfg, printFirstTime);
      else
         parseHeterogeneousMOS(cfg, printFirstTime);

   }


   virtual ~ConfigFileParser () {};


   bool parseHomogeneousMOS(libconfig::Config& cfg, bool printFirstTime) {

      try {

         libconfig::Setting& techs = cfg.lookup ("techs");

         unsigned nTechs = techs.getLength ();
         const std::vector<std::string>* techsToUse = GAEDAConfig::handle()->getTechsToUse();
         unsigned nTechsUse = techsToUse->size();

         unsigned j = 0;

         if (CommManager::instance()->isIslandMaster() && !GAEDAConfig::handle()->quiet() && printFirstTime)
            std::cout << "[ConfigFileParser] Using MOS techniques: ";

         for (unsigned i = 0; i < nTechs || j < nTechsUse; i++) {

            std::string ids = (const char*) techs [i]["ids"];

            if( find (techsToUse->begin(), techsToUse->end(), ids) != techsToUse->end () ) {

              if (CommManager::instance()->isIslandMaster() && !GAEDAConfig::handle()->quiet() && printFirstTime)
                cout << ids << "  ";

               j++;

               std::string alg = (const char*) techs [i]["alg"];

               MOSTechnique* tech;

               if (alg == "ga")
                  tech = parseGATechnique(techs[i]);
               else if (alg == "memeticga")
                  tech = parseMemeticGATechnique(techs[i]);
               else if (alg == "eda")
                  tech = parseEDATechnique(techs[i]);
               else if (alg == "de" || alg == "ctb2de" || alg == "stsde" || alg == "depop" || alg == "adapde" || alg == "adapstsde")
                  tech = parseDETechnique(techs[i],alg);
               else if (alg == "cmaes")
                  tech = parseCMAESTechnique(techs[i]);
               else if (alg == "es")
                  tech = parseESTechnique(techs[i]);
               else if (alg == "mts")
                  tech = parseMTSTechnique(techs[i]);
               else if (alg == "ls")
                  tech = parseLSTechnique(techs[i]);
               else if (alg == "popls")
                  tech = parsePopLSTechnique(techs[i]);
               else if (alg == "nma")
                  tech = parseNMATechnique(techs[i]);
               else {

                 std::cerr << "[ConfigFileParser] Error: Unknown technique type. Aborting..." << std::endl;
                 exit (-1);

               }

               //tech->setInUse (true);

               mID++;
               mTechSet->registerTechnique (tech);

            }

         }
         if( j != nTechsUse ) {

           std::cerr << "[ConfigFileParser] Error: Some of the techniques couldn't be loaded "
                     << "probably because the technique id doesn't exists. Aborting..." << std::endl;

           exit (-1);

         }

      }
      catch (libconfig::SettingNotFoundException e) {

         std::cerr << "[ConfigFileParser] Error: one or more mandatory parameter(s) were not found." << std::endl;
         exit (-1);

      }

      if (CommManager::instance()->isIslandMaster() && !GAEDAConfig::handle()->quiet() && printFirstTime)
         std::cout << std::endl;

      return true;
   }


   bool parseHeterogeneousMOS(libconfig::Config& cfg, bool printFirstTime) {

      // Parse heterogeneous configuration file
      std::string hetMOScfg = GAEDAConfig::handle()->getHetMOSCfgFile();

      libconfig::Config hetCfg;

      try {
         hetCfg.readFile (hetMOScfg.c_str());
      }
      catch (libconfig::ParseException e) {
         std::cerr << "[ConfigFileParser] Error: parsing the configuration file '" << hetMOScfg
                   << "' in line (" << e.getLine () << "): " << e.getError() << "." << std::endl;
         exit (-1);
      }
      catch (libconfig::FileIOException e) {
         std::cerr << "[ConfigFileParser] Error: configuration file not found or could not be open '"
                   << hetMOScfg << "'." << std::endl;
         exit (-1);
      }

      // Parse techniques file and store only those defined in the configuration of this island
      try {

         // Retrieve the whole list of techniques to be used
         const std::vector<std::string>* techsToUse = GAEDAConfig::handle()->getTechsToUse();
         unsigned nTechsUse = techsToUse->size();

         // Retrieve the techniques to be used in this island
         int rank = CommManager::instance()->getMyRank();

         libconfig::Setting& islands = hetCfg.lookup ("islands");
         unsigned nIslands = islands.getLength ();

         if (rank < 0        ||
             rank > nIslands   ) {

            std::cerr << "[ConfigFileParser] Error: No config for island " << rank
                      << " in config file '" << hetMOScfg << "'" << std::endl;
            exit(-1);

         }

         unsigned nTechsIsland = islands[rank].getLength();

         std::vector<std::string> techsIsland;

         for (unsigned i = 0; i < nTechsIsland; i++)
            techsIsland.push_back (islands[rank][i]);

         // Retrieve the available techniques
         libconfig::Setting& techs = cfg.lookup ("techs");
         unsigned nTechs = techs.getLength ();

         // Add those techniques meaningful to this island to the techniques set
         if (CommManager::instance()->isIslandMaster() && !GAEDAConfig::handle()->quiet() && printFirstTime)
            std::cout << "[ConfigFileParser] Island " << rank << " using MOS techniques: ";

         unsigned j = 0;

         for (unsigned i = 0; i < nTechs && j < nTechsUse; i++) {

            std::string ids = (const char*) techs [i]["ids"];

            if (find (techsToUse->begin(), techsToUse->end(), ids) != techsToUse->end ()) {

              j++;

               std::string alg = (const char*) techs [i]["alg"];

               MOSTechnique* tech;

               if (alg == "ga")
                  tech = parseGATechnique(techs[i]);
               else if (alg == "memeticga")
                  tech = parseMemeticGATechnique(techs[i]);
               else if (alg == "eda")
                  tech = parseEDATechnique(techs[i]);
               else if (alg == "de" || alg == "ctb2de" || alg == "stsde" || alg == "depop" || alg == "adapde" || alg == "adapstsde")
                  tech = parseDETechnique(techs[i],alg);
               else if (alg == "cmaes")
                  tech = parseCMAESTechnique(techs[i]);
               else if (alg == "es")
                  tech = parseESTechnique(techs[i]);
               else if (alg == "mts")
                  tech = parseMTSTechnique(techs[i]);
               else if (alg == "ls")
                  tech = parseLSTechnique(techs[i]);
               else if (alg == "popls")
                  tech = parsePopLSTechnique(techs[i]);
               else if (alg == "nma")
                  tech = parseNMATechnique(techs[i]);
               else {

                 std::cerr << "[ConfigFileParser] Error: Unknown technique type. Aborting..." << std::endl;
                 exit (-1);

               }

               mID++;
               mTechSet->registerTechnique (tech);

//               if (find (techsIsland.begin(), techsIsland.end(), ids) != techsIsland.end ()) {
                    if (CommManager::instance()->isIslandMaster() && !GAEDAConfig::handle()->quiet() && printFirstTime)
                       cout << ids << "  ";
//                  tech->setInUse (true);
//               }
//               else
//                  tech->setInUse (false);

            }

         }

         if( j != nTechsUse ) {

           std::cerr << "[ConfigFileParser] Error: Some of the techniques couldn't be loaded "
                     << "probably because the technique id doesn't exists. Aborting..." << std::endl;

           exit (-1);

         }

      }
      catch (libconfig::SettingNotFoundException e) {

         std::cerr << "[ConfigFileParser] Error: one or more mandatory parameter(s) were not found." << std::endl;
         exit (-1);

      }

      if (CommManager::instance()->isIslandMaster() && !GAEDAConfig::handle()->quiet() && printFirstTime)
         std::cout << std::endl;

      return true;
   }


   MOSTechnique* parseGATechnique (libconfig::Setting& t) {

      // Parse crossover
      void* cx = dlsym (mHandle, (const char*) t["cx"]);

      // Parse mutator
      void* mut = dlsym (mHandle, (const char*) t["mut"]);

      // Parse comparator
      void* cmp = dlsym (mHandle, (const char*) t["cmp"]);

      // Parse initializer
      void* ini = dlsym (mHandle, (const char*) t["ini"]);

      // Parse objective function
      void* obj = dlsym (mHandle, (const char*) t["obj"]);

      // Parse crossover probability
      double pcx = (double) t["pcx"];

      // Parse crossover probability
      double pmut = (double) t["pmut"];

      // Parse encoding
      std::string coding = (const char*) t["coding"];

      // Parse technique name and description
      std::string ids  = (const char*) t["ids" ];
      std::string desc = (const char*) t["desc"];

      // Parse technique selector
      std::string strSelect = (const char*) t["select"];
      GASelectionScheme* select = stringToSelector( strSelect );

      // Check for missing operators or wrong values
      if (!cx || !mut || !cmp || !ini || !obj || pcx < 0.0 || pcx > 1.0 ||
          pmut < 0.0 || pmut > 1.0 || (stringToEnum (coding) == -1) || !select ) {

         std::cerr << "Error: Some of the parameters of the GA technique " << ids
                   << " is not properly defined in " << mCfgFile <<". Aborting..." << std::endl;
         exit (-1);

      }
      else {

         MOSTechnique* tech = new MOSTechniqueGA (mID,
                                                 desc,
                                                 (GAGenome::Mutator) mut,
                                                 (GAGenome::SexualCrossover) cx,
                                                 (GAGenome::Comparator) cmp,
                                                 (GAGenome::Initializer) ini,
                                                 (GAGenome::Evaluator) obj,
                                                 pcx,
                                                 pmut,
                                                 stringToEnum (coding),
                                                 mGenFactory->getGenome (stringToEnum (coding)),
                                                 select);

         return tech;

      }

   }


   MOSTechnique* parseMemeticGATechnique (libconfig::Setting& t) {

      // Parse crossover
      void* cx = dlsym (mHandle, (const char*) t["cx"]);

      // Parse mutator
      void* mut = dlsym (mHandle, (const char*) t["mut"]);

      // Parse comparator
      void* cmp = dlsym (mHandle, (const char*) t["cmp"]);

      // Parse initializer
      void* ini = dlsym (mHandle, (const char*) t["ini"]);

      // Parse objective function
      void* obj = dlsym (mHandle, (const char*) t["obj"]);

      // Parse crossover probability
      double pcx = (double) t["pcx"];

      // Parse crossover probability
      double pmut = (double) t["pmut"];

      // Parse encoding
      std::string coding = (const char*) t["coding"];

      // Parse technique name and description
      std::string ids  = (const char*) t["ids" ];
      std::string desc = (const char*) t["desc"];

      // Parse technique selector
      std::string strSelect = (const char*) t["select"];
      GASelectionScheme* select = stringToSelector( strSelect );

      // Parse LS methods and policies
      libconfig::Setting& LSSet = t["lss"];
      std::vector<LSType> lss;

      for (unsigned i = 0; i < LSSet.getLength(); i++) {
        LSType ls;
        ls.ls = (GAGenome::Mutator) dlsym (mHandle, (const char*) LSSet[i]["ls"]);
        ls.strategy = (int) LSSet[i]["strategy"];
        ls.freq = (int) LSSet[i]["freq"];
        lss.push_back(ls);
      }

      // Check for missing operators or wrong values
      if (!cx || !mut || !cmp || !ini || !obj || pcx < 0.0 || pcx > 1.0 ||
          pmut < 0.0 || pmut > 1.0 || (stringToEnum (coding) == -1) || !select || lss.size() == 0) {

         std::cerr << "Error: Some of the parameters of the GA technique " << ids
                   << " is not properly defined in " << mCfgFile <<". Aborting..." << std::endl;
         exit (-1);

      }
      else {

         MOSTechnique* tech = new MOSTechniqueMemeticGA (mID,
                                                 desc,
                                                 (GAGenome::Mutator) mut,
                                                 (GAGenome::SexualCrossover) cx,
                                                 (GAGenome::Comparator) cmp,
                                                 (GAGenome::Initializer) ini,
                                                 (GAGenome::Evaluator) obj,
                                                 pcx,
                                                 pmut,
                                                 stringToEnum (coding),
                                                 mGenFactory->getGenome (stringToEnum (coding)),
                                                 select,
                                                 lss);

         return tech;

      }

   }


  MOSTechnique* parseEDATechnique (libconfig::Setting& t) {

      // Parse initializer
      void* ini = dlsym (mHandle, (const char*) t["ini"]);

      // Parse objective function
      void* obj = dlsym (mHandle, (const char*) t["obj"]);

      // Parse encoding
      std::string coding = (const char*) t["coding"];

      // Parse technique name and description
      std::string ids  = (const char*) t["ids" ];
      std::string desc = (const char*) t["desc"];

      // Parse Select Per
      double sper = (double) t["sper"];

      // Parse network type
      std::string network((const char*)t["network"]);

      GAGraphModel::ModelType networkId = (GAGraphModel::ModelType)stringToEnum(network);

      std::string learnMethodString((const char*)t["learnMethod"]);
      GABayesianNetwork::LearningMethod learnMethodBayes;
      GAGaussianNetwork::LearningMethod learnMethodGauss;

      std::string localScoringString;
      GABayesianNetwork::EBNALocalScoring localScoring;
      std::string simMethodString;
      GABayesianNetwork::SimulationMethod simMethod;

      std::string scoreMethodString;
      GAGaussianNetwork::ScoreMethod scoreMethod;

      MOSTechnique* tech = NULL;

      switch (networkId) {

      case GAGraphModel::BAYESIAN_NETWORK:

         localScoringString = (const char*) t["localScoring"];
         localScoring = (GABayesianNetwork::EBNALocalScoring)stringToEnum(localScoringString);

         simMethodString = (const char*) t["simMethod"];
         simMethod = (GABayesianNetwork::SimulationMethod)stringToEnum(simMethodString);

         learnMethodBayes = (GABayesianNetwork::LearningMethod)stringToEnum(learnMethodString);

         if (!ini || !obj || (stringToEnum (coding) == -1) ){

            std::cerr << "[ConfigFileParser] Error: Some of the parameters of the EDA technique "
                      << ids <<" is not properly defined. Aborting..." << std::endl;
            exit (-1);

         }
         else {

            tech = new MOSTechniqueEDA (mID,
                                        desc,
                                        sper,
                                        networkId,
                                        learnMethodBayes,
                                        localScoring,
                                        simMethod,
                                        stringToEnum (coding),
                                        mGenFactory->getGenome (stringToEnum (coding)),
                                        (GAGenome::Initializer) ini,
                                        (GAGenome::Evaluator) obj);

         }

         break;

      case GAGraphModel::GAUSSIAN_NETWORK:

         scoreMethodString = (const char*) t["scoreMethod"];
         scoreMethod = (GAGaussianNetwork::ScoreMethod)stringToEnum(scoreMethodString);

         learnMethodGauss = (GAGaussianNetwork::LearningMethod)stringToEnum(learnMethodString);

         if (!ini || !obj || (stringToEnum (coding) == -1) ){

            std::cerr << "[ConfigFileParser] Error: Some of the parameters of the EDA technique "
                      << ids <<" is not properly defined. Aborting..." << std::endl;

            exit (-1);

         }
         else {

            tech = new MOSTechniqueEDA (mID,
                                       desc,
                                       sper,
                                       networkId,
                                       learnMethodGauss,
                                       scoreMethod,
                                       stringToEnum (coding),
                                       mGenFactory->getGenome (stringToEnum (coding)),
                                       (GAGenome::Initializer) ini,
                                       (GAGenome::Evaluator) obj);

         }

         break;

      default:

         std::cerr << "[ConfigFileParser] Error: EDA Graph model not recognized" << std::endl;

         exit (-1);
         break;

      }

      return tech;

   }


   MOSTechnique* parseDETechnique (libconfig::Setting& t, std::string de_alg) {

      // Parse crossover
      void* cx = dlsym (mHandle, (const char*) t["cx"]);

      // Parse initializer
      void* ini = dlsym (mHandle, (const char*) t["ini"]);

      // Parse objective function
      void* obj = dlsym (mHandle, (const char*) t["obj"]);

      // F
      double f = (double) t["f"];

      // CR
      double cr = (double) t["cr"];

      // Parse encoding
      std::string coding = (const char*) t["coding"];

      // Parse technique name and description
      std::string ids  = (const char*) t["ids" ];
      std::string desc = (const char*) t["desc"];

      // Parse technique selector
      std::string strSelect = (const char*) t["select"];
      GASelectionScheme* select = stringToSelector( strSelect );

      // Check for missing operators or wrong values
      if (!cx || !ini || !obj || !f || !cr || (stringToEnum (coding) == -1) || !select ) {

         std::cerr << "Error: Some of the parameters of the DE technique "<< ids
                   <<" is not properly defined in " << mCfgFile << ". Aborting..." << std::endl;

         exit (-1);

      }
      else {
        MOSTechnique* tech = NULL;

        if (de_alg == "de") {
          tech = new MOSTechniqueDE (mID,
                                     desc,
                                     (GAGenome::Initializer) ini,
                                     (GAGenome::Evaluator) obj,
                                     stringToEnum (coding),
                                     mGenFactory->getGenome (stringToEnum (coding)),
                                     (GAGenome::DECrossover) cx,
                                     f,
                                     cr,
                                     select);
        }
        else if (de_alg == "stsde") {
          tech = new MOSTechniqueSTSDE (mID,
                                        desc,
                                        (GAGenome::Initializer) ini,
                                        (GAGenome::Evaluator) obj,
                                        stringToEnum (coding),
                                        mGenFactory->getGenome (stringToEnum (coding)),
                                        (GAGenome::DECrossover) cx,
                                        f,
                                        cr,
                                        select);
        }
        else if (de_alg == "depop") {
          tech = new MOSTechniqueDEPopSize (mID,
                                            desc,
                                            (GAGenome::Initializer) ini,
                                            (GAGenome::Evaluator) obj,
                                            stringToEnum (coding),
                                            mGenFactory->getGenome (stringToEnum (coding)),
                                            (GAGenome::DECrossover) cx,
                                            f,
                                            cr,
                                            select);
        }
        else if (de_alg == "adapde") {
          tech = new MOSTechniqueAdapDE (mID,
                                         desc,
                                         (GAGenome::Initializer) ini,
                                         (GAGenome::Evaluator) obj,
                                         stringToEnum (coding),
                                         mGenFactory->getGenome (stringToEnum (coding)),
                                         (GAGenome::DECrossover) cx,
                                         f,
                                         cr,
                                         select);
        }
        else if (de_alg == "adapstsde") {
          tech = new MOSTechniqueAdapSTSDE (mID,
                                            desc,
                                            (GAGenome::Initializer) ini,
                                            (GAGenome::Evaluator) obj,
                                            stringToEnum (coding),
                                            mGenFactory->getGenome (stringToEnum (coding)),
                                            (GAGenome::DECrossover) cx,
                                            f,
                                            cr,
                                            select);
          }
/*         else if (de_alg == "ctb2de") { */
/*           tech = new MOSTechniqueCTB2DE (mID, */
/*                                          desc, */
/*                                          (GAGenome::Initializer) ini, */
/*                                          (GAGenome::Evaluator) obj, */
/*                                          stringToEnum (coding), */
/*                                          mGenFactory->getGenome (stringToEnum (coding)), */
/*                                          (GAGenome::DECrossover) cx, */
/*                                          f, */
/*                                          cr, */
/*                                          select); */
/*         } */
        else throw runtime_error("Unrecognaised DE algorithm in parseDETechnique");

        return tech;
      }
   }

   MOSTechnique* parseCMAESTechnique (libconfig::Setting& t) {

      // Parse initializer
      void* ini = dlsym (mHandle, (const char*) t["ini"]);

      // Parse objective function
      void* obj = dlsym (mHandle, (const char*) t["obj"]);

      // Parse encoding
      std::string coding = (const char*) t["coding"];

      // Parse technique name and description
      std::string ids  = (const char*) t["ids" ];
      std::string desc = (const char*) t["desc"];

      // Parse technique selector
      std::string strSelect = (const char*) t["select"];
      GASelectionScheme* select = stringToSelector( strSelect );

      // Check for missing operators or wrong values
      if (!ini || !obj || (stringToEnum (coding) == -1) || !select ) {

         std::cerr << "Error: Some of the parameters of the CMAES technique "<< ids
                   <<" is not properly defined in " << mCfgFile << ". Aborting..." << std::endl;

         exit (-1);

      }
      else {
        MOSTechnique* tech = NULL;

        tech = new MOSTechniqueCMAES (mID,
                                      desc,
                                      (GAGenome::Initializer) ini,
                                      (GAGenome::Evaluator) obj,
                                      stringToEnum (coding),
                                      mGenFactory->getGenome (stringToEnum (coding)),
                                      select);
        return tech;
      }
   }

   MOSTechnique* parseESTechnique (libconfig::Setting& t) {

      // Parse crossover
      void* cx = dlsym (mHandle, (const char*) t["cx"]);

      // Parse crossover update
      std::string cxu = t["cx"];
      cxu+="UpdateOnly";
      void* cx_upd = dlsym (mHandle, cxu.c_str());

      // Parse mutator
      void* mut = dlsym (mHandle, (const char*) t["mut"]);

      // Parse mutator update
      std::string mutu = t["mut"];
      mutu+="UpdateOnly";
      void* mut_upd = dlsym (mHandle, mutu.c_str());

      // Parse comparator
      void* cmp = dlsym (mHandle, (const char*) t["cmp"]);

      // Parse initializer
      void* ini = dlsym (mHandle, (const char*) t["ini"]);

      // Parse objective function
      void* obj = dlsym (mHandle, (const char*) t["obj"]);

      // Mu
      unsigned mu = (unsigned) t["mu"];

      // Ro
      unsigned ro = (unsigned) t["ro"];

      // Lambda
      unsigned lambda = (unsigned) t["lambda"];

      // Parse encoding
      std::string coding = (const char*) t["coding"];

      // Parse technique name and description
      std::string ids  = (const char*) t["ids" ];
      std::string desc = (const char*) t["desc"];

      // Parse technique selector
      std::string strSelect = (const char*) t["select"];
      GASelectionScheme* select = stringToSelector (strSelect);

      // Check for missing operators or wrong values
      if (!cx || !mut || !cmp || !ini || !obj || !mu || !ro || !lambda ||
          (stringToEnum (coding) == -1) || !select ) {

         std::cerr << "Error: Some of the parameters of the ES technique "<< ids
                   << " is not properly defined in " << mCfgFile << ". Aborting..." << std::endl;

         exit (-1);

      }
      else {

         MOSTechnique* tech = new MOSTechniqueES (mID,
                                                  desc,
                                                  (GAGenome::ESMutator) mut,
                                                  (GAGenome::ESMutator) mut_upd,
                                                  (GAGenome::ESCrossover) cx,
                                                  (GAGenome::ESCrossover) cx_upd,
                                                  (GAGenome::Comparator) cmp,
                                                  (GAGenome::Initializer) ini,
                                                  (GAGenome::Evaluator) obj,
                                                  mu,
                                                  ro,
                                                  lambda,
                                                  stringToEnum (coding),
                                                  mGenFactory->getGenome (stringToEnum (coding)),
                                                  select);

         return tech;

      }

   }

   MOSTechnique* parseMTSTechnique (libconfig::Setting& t) {

      // Parse LS method
      void* ls = dlsym (mHandle, (const char*) t["ls"]);

      // Parse comparator
      void* cmp = dlsym (mHandle, (const char*) t["cmp"]);

      // Parse initializer
      void* ini = dlsym (mHandle, (const char*) t["ini"]);

      // Parse objective function
      void* obj = dlsym (mHandle, (const char*) t["obj"]);

      // Parse encoding
      std::string coding = (const char*) t["coding"];

      // Parse technique name and description
      std::string ids  = (const char*) t["ids" ];
      std::string desc = (const char*) t["desc"];

      // Parse technique selector
      std::string strSelect = (const char*) t["select"];
      GASelectionScheme* select = stringToSelector( strSelect );

      // Check for missing operators or wrong values
      if (!ls || !cmp || !ini || !obj || (stringToEnum (coding) == -1) || !select ) {

         std::cerr << "Error: Some of the parameters of the MTS technique " << ids
                   << " is not properly defined in " << mCfgFile <<". Aborting..." << std::endl;
         exit (-1);

      }
      else {

         MOSTechnique* tech = new MOSTechniqueLS (mID,
                                                  desc,
                                                  (LocalSearch) ls,
                                                  (GAGenome::Comparator) cmp,
                                                  (GAGenome::Initializer) ini,
                                                  (GAGenome::Evaluator) obj,
                                                  stringToEnum (coding),
                                                  mGenFactory->getGenome (stringToEnum (coding)),
                                                  select);

         return tech;

      }

   }


   MOSTechnique* parseLSTechnique (libconfig::Setting& t) {

      // Parse LS method
      std::string ls = (const char*) t["ls"];

      // Parse comparator
      void* cmp = dlsym (mHandle, (const char*) t["cmp"]);

      // Parse initializer
      void* ini = dlsym (mHandle, (const char*) t["ini"]);

      // Parse objective function
      void* obj = dlsym (mHandle, (const char*) t["obj"]);

      // Parse encoding
      std::string coding = (const char*) t["coding"];

      // Parse technique name and description
      std::string ids  = (const char*) t["ids" ];
      std::string desc = (const char*) t["desc"];

      // Parse technique selector
      std::string strSelect = (const char*) t["select"];
      GASelectionScheme* select = stringToSelector( strSelect );

      // Check for missing operators or wrong values
      if (ls.length() == 0 || !cmp || !ini || !obj || (stringToEnum (coding) == -1) || !select ) {

         std::cerr << "Error: Some of the parameters of the MTS technique " << ids
                   << " is not properly defined in " << mCfgFile <<". Aborting..." << std::endl;
         exit (-1);

      }
      else {

         MOSTechnique* tech = new MOSTechniqueLSObj (mID,
                                                     desc,
                                                     ls,
                                                     (GAGenome::Comparator) cmp,
                                                     (GAGenome::Initializer) ini,
                                                     (GAGenome::Evaluator) obj,
                                                     stringToEnum (coding),
                                                     mGenFactory->getGenome (stringToEnum (coding)),
                                                     select);

         return tech;

      }

   }


   MOSTechnique* parsePopLSTechnique (libconfig::Setting& t) {

      // Parse LS method
      void* ls = dlsym (mHandle, (const char*) t["ls"]);

      // Parse comparator
      void* cmp = dlsym (mHandle, (const char*) t["cmp"]);

      // Parse initializer
      void* ini = dlsym (mHandle, (const char*) t["ini"]);

      // Parse objective function
      void* obj = dlsym (mHandle, (const char*) t["obj"]);

      // Parse encoding
      std::string coding = (const char*) t["coding"];

      // Parse technique name and description
      std::string ids  = (const char*) t["ids" ];
      std::string desc = (const char*) t["desc"];

      // Parse technique selector
      std::string strSelect = (const char*) t["select"];
      GASelectionScheme* select = stringToSelector( strSelect );

      // Check for missing operators or wrong values
      if (!ls || !cmp || !ini || !obj || (stringToEnum (coding) == -1) || !select ) {

         std::cerr << "Error: Some of the parameters of the PopLS technique " << ids
                   << " is not properly defined in " << mCfgFile <<". Aborting..." << std::endl;
         exit (-1);

      }
      else {

         MOSTechnique* tech = new MOSTechniquePopLS (mID,
                                                     desc,
                                                     (PopLocalSearch) ls,
                                                     (GAGenome::Comparator) cmp,
                                                     (GAGenome::Initializer) ini,
                                                     (GAGenome::Evaluator) obj,
                                                     stringToEnum (coding),
                                                     mGenFactory->getGenome (stringToEnum (coding)),
                                                     select);

         return tech;

      }

   }


   MOSTechnique* parseNMATechnique (libconfig::Setting& t) {

      // Parse comparator
      void* cmp = dlsym (mHandle, (const char*) t["cmp"]);

      // Parse initializer
      void* ini = dlsym (mHandle, (const char*) t["ini"]);

      // Parse objective function
      void* obj = dlsym (mHandle, (const char*) t["obj"]);

      // Parse encoding
      std::string coding = (const char*) t["coding"];

      // Parse technique name and description
      std::string ids  = (const char*) t["ids" ];
      std::string desc = (const char*) t["desc"];

      // Parse technique selector
      std::string strSelect = (const char*) t["select"];
      GASelectionScheme* select = stringToSelector( strSelect );

      // Check for missing operators or wrong values
      if (!cmp || !ini || !obj || (stringToEnum (coding) == -1) || !select ) {

         std::cerr << "Error: Some of the parameters of the NMA technique " << ids
                   << " is not properly defined in " << mCfgFile <<". Aborting..." << std::endl;
         exit (-1);

      }
      else {

         MOSTechnique* tech = new MOSTechniqueNMA (mID,
                                                   desc,
                                                   (GAGenome::Comparator) cmp,
                                                   (GAGenome::Initializer) ini,
                                                   (GAGenome::Evaluator) obj,
                                                   stringToEnum (coding),
                                                   mGenFactory->getGenome (stringToEnum (coding)),
                                                   select);

         return tech;

      }

   }


   void buildMappingTable () {

      mappingTable.insert (std::make_pair ("BAYESIAN_NETWORK",                    GAGraphModel::BAYESIAN_NETWORK));
      mappingTable.insert (std::make_pair ("GAUSSIAN_NETWORK",                    GAGraphModel::GAUSSIAN_NETWORK));

      mappingTable.insert (std::make_pair ("GAGaussianNetwork::UMDA",             GAGaussianNetwork::UMDA));
      mappingTable.insert (std::make_pair ("GAGaussianNetwork::EGNA_B",           GAGaussianNetwork::EGNA_B));
      mappingTable.insert (std::make_pair ("GAGaussianNetwork::EGNA_LOCAL",       GAGaussianNetwork::EGNA_LOCAL));
      mappingTable.insert (std::make_pair ("GAGaussianNetwork::MIMIC",            GAGaussianNetwork::MIMIC));
      mappingTable.insert (std::make_pair ("GAGaussianNetwork::EE",               GAGaussianNetwork::EE));
      mappingTable.insert (std::make_pair ("GAGaussianNetwork::EMNA",             GAGaussianNetwork::EMNA));

      mappingTable.insert (std::make_pair ("GAGaussianNetwork::BGe_SCORE",        GAGaussianNetwork::BGe_SCORE));
      mappingTable.insert (std::make_pair ("GAGaussianNetwork::BIC_SCORE",        GAGaussianNetwork::BIC_SCORE));

      mappingTable.insert (std::make_pair ("GABayesianNetwork::UMDA",             GABayesianNetwork::UMDA));
      mappingTable.insert (std::make_pair ("GABayesianNetwork::EBNA_B",           GABayesianNetwork::EBNA_B));
      mappingTable.insert (std::make_pair ("GABayesianNetwork::EBNA_LOCAL",       GABayesianNetwork::EBNA_LOCAL));
      mappingTable.insert (std::make_pair ("GABayesianNetwork::PBIL",             GABayesianNetwork::PBIL));
      mappingTable.insert (std::make_pair ("GABayesianNetwork::TREE",             GABayesianNetwork::TREE));
      mappingTable.insert (std::make_pair ("GABayesianNetwork::MIMIC",            GABayesianNetwork::MIMIC));
      mappingTable.insert (std::make_pair ("GABayesianNetwork::EBNA_K2",          GABayesianNetwork::EBNA_K2));

      mappingTable.insert (std::make_pair ("GABayesianNetwork::BIC_SCORE",        GABayesianNetwork::BIC_SCORE));
      mappingTable.insert (std::make_pair ("GABayesianNetwork::K2_SCORE",         GABayesianNetwork::K2_SCORE));

      mappingTable.insert (std::make_pair ("GABayesianNetwork::PLS",              GABayesianNetwork::PLS));
      mappingTable.insert (std::make_pair ("GABayesianNetwork::PLS_ALL_VALUES_1", GABayesianNetwork::PLS_ALL_VALUES_1));
      mappingTable.insert (std::make_pair ("GABayesianNetwork::PLS_ALL_VALUES_2", GABayesianNetwork::PLS_ALL_VALUES_2));
      mappingTable.insert (std::make_pair ("GABayesianNetwork::PLS_CORRECT",      GABayesianNetwork::PLS_CORRECT));
      mappingTable.insert (std::make_pair ("GABayesianNetwork::PENALIZATION",     GABayesianNetwork::PENALIZATION));

      mappingTable.insert (std::make_pair ("Integer", (int) GAID::IntegerEncoding));
      mappingTable.insert (std::make_pair ("Real",    (int) GAID::RealEncoding));

   }


   int stringToEnum (std::string& str) {

      std::map< std::string, int >::const_iterator it = mappingTable.find (str);

      if (it == mappingTable.end ())
         return -1;
      else
         return it->second;

   }

   GASelectionScheme* stringToSelector (std::string& str) {

     GASelectionScheme *selector = NULL, *sel_tmp = NULL;
     int selector_degree;

     if (sscanf(str.c_str(), "tournament:%d", &selector_degree)) {
       sel_tmp  = new GARouletteWheelSelector();
       selector = new GATournamentSelector(selector_degree, *sel_tmp);
       delete sel_tmp;
     }
     else if (sscanf(str.c_str(), "unif_tournament:%d", &selector_degree)) {
       sel_tmp  = new GAUniformSelector();
       selector = new GATournamentSelector(selector_degree, *sel_tmp);
       delete sel_tmp;
     }
     else if (str == "roulette")
       selector = new GARouletteWheelSelector();
     else if (str == "random")
       selector = new GAUniformSelector();

     return selector;

   }


   protected:

      void* mHandle;
      unsigned mID;

      std::string mCfgFile;

      MOSGenomeFactory* mGenFactory;
      MOSTechniqueSet*  mTechSet;

      std::map< std::string, int > mappingTable;

};

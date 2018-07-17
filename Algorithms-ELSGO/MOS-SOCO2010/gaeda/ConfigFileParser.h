#include <dlfcn.h>

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#include "gaid.h"
#include "GARealOps.h"
#include "config/libconfigcpp.h"
#include "genomes/GAGenome.h"
#include "MOSTechnique.h"
#include "MOSTechniqueDE.h"
#include "MOSTechniqueLS.h"
#include "MOSTechniqueSet.h"

#include "islands/CommManager.h" //TODO quitar

using namespace std;

class ConfigFileParser {

   public:

   ConfigFileParser (const char* cfgFile) : mCfgFile (cfgFile) {

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

      parseHomogeneousMOS(cfg);

   }


   virtual ~ConfigFileParser () {};


   bool parseHomogeneousMOS(libconfig::Config& cfg) {

      try {

         libconfig::Setting& techs = cfg.lookup ("techs");

         unsigned nTechs = techs.getLength ();
         const std::vector<std::string>* techsToUse = GAEDAConfig::handle()->getTechsToUse();
         unsigned nTechsUse = techsToUse->size();

         unsigned j = 0;

         if (CommManager::instance()->isIslandMaster() && !GAEDAConfig::handle()->quiet())
            std::cout << "[ConfigFileParser] Using MOS techniques: ";

         for (unsigned i = 0; i < nTechs || j < nTechsUse; i++) {

            std::string ids = (const char*) techs [i]["ids"];

            if( find (techsToUse->begin(), techsToUse->end(), ids) != techsToUse->end () ) {

              if (CommManager::instance()->isIslandMaster() && !GAEDAConfig::handle()->quiet())
                cout << ids << "  ";

               j++;

               std::string alg = (const char*) techs [i]["alg"];

               MOSTechnique* tech;

               if (alg == "de")
                  tech = parseDETechnique(techs[i],alg);
               else if (alg == "mts")
                  tech = parseMTSTechnique(techs[i]);
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

      if (CommManager::instance()->isIslandMaster() && !GAEDAConfig::handle()->quiet())
         std::cout << std::endl;

   }


   MOSTechnique* parseDETechnique (libconfig::Setting& t, std::string de_alg) {

      // Parse crossover
      void* cx = dlsym (mHandle, (const char*) t["cx"]);

      // Parse initializer
      void* ini = dlsym (mHandle, (const char*) t["ini"]);

      // Parse objective function
      void* obj = dlsym (mHandle, (const char*) t["obj"]);

      // F
      long double f = (double) t["f"];

      // CR
      long double cr = (double) t["cr"];

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


   void buildMappingTable () {

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

#include "DARPCostMatrix.h"
#include "OnlyDistCostMatrix.h"
#include "DistMatrix.h"
#include "VerticesList.h"
#include "VNSDARPGenome.h"
#include "TSDARPGenome.h"
#include "DARPEvaluator.h"
#include "darpOps.h"
#include "darpInit.h"
#include "RandomReqDistributor.h"
#include "ParReqDistMatrix.h"
#include "RandomRouteDistributor.h"
#include "KMedoidsReqDistributor.h"
#include "KMedoidsRouteDistributor.h"
#include "ReqDistMatrix.h"

#include <GAEDAlib.h>
#include <GAEDAConfig.h>
#include <libconfig.h++>

#include <algorithm>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <sstream>

#include <time.h>

extern "C" {
  #include "extras/cluster.h"
}

static DARPGenome* gen = 0;

static double maxRideValue;
static long vehicleCapacity;
static long nVehicles;
static int timeSumExp;
static long nStops;
static long maxDelay;
static int relativeDelayConstant = 0;
static DARPEvaluator::darpOptCriterionType optcrit;
static bool relativeDelays;

static int numRouteCallsStatsDispFreq = -1;

static double minEvalWeightsUpdateValue;
static double maxEvalWeightsUpdateValue;

static double loadVWeight;
static double TWVWeight;
static double rideVWeight;

static double initloadVWeight = -1;
static double initTWVWeight   = -1;
static double initrideVWeight = -1;

static std::string DistMatrixFile;
static std::string RequestsFile;

static CostMatrix*    costMatrix   = 0;
static DistMatrix*    distMatrix   = 0;
static VerticesList*  vertlist     = 0;
static DARPEvaluator* evaluator    = 0;

static string          genType;
static string          initname          = "";
static darpInitFuncT   initfunc          = 0;
static bool            postinitls        = true;
static bool            printBestSol      = true;
static string          objfuncinitMethod = "best";
static string          reqdistname       = "";
static int             initkmedoidsnpass = -1;
static kmedoidsMode    initkmedoidsmode  = KMEDNORMAL;
static string          localsearch       = "DARPLocalSearch";
static int             maxlsevals        = -1;
static bool            lsrandomroutes    = false;

static bool            kmedoidsuseconstrpen        = true;
static bool            kmedoidsusewaitingpen       = false;
static long            kmedoidswaitingpenthreshold = 0;
static double          kmedoidswaitingpenconstant  = 0;

static string            routedistname = "";
static RouteDistributor* routedist     = 0;
static int               kmedoidsroutedistnpass = -1;
static kmedoidsMode      kmedoidsroutedistmode  = KMEDNORMAL;


std::vector< VNSShaker* > shakers;

// Distributed algorithm options
static ReqDistributor* reqdist       = 0;
static ReqDistMatrix*  reqdistmatrix = 0;
static bool            isdistributed = false;

ReqDistributor* createReqDist(string name, DARPGenome& gen, CostMatrix& costmatrix, DARPEvaluator& eval) {
  ReqDistributor* dist = 0;
  if (name == "random") {
    dist = new RandomReqDistributor();
  }
  else if  (name == "kmedoids") {
    if (initkmedoidsnpass <0) initkmedoidsnpass = 10;

    dist = new KMedoidsReqDistributor(initkmedoidsnpass,initkmedoidsmode);
  }
  else {
    throw runtime_error("unrecognized req distributor name");
  }

  return dist;
}

kmedoidsMode getKMedoidsMode(string name) {
  kmedoidsMode mode;

  if      ( name == "normal"       ) mode = KMEDNORMAL;
  else if ( name == "forcebalanced") mode = KMEDFORCEBALANCE;
  else if ( name == "bestbalanced" ) mode = KMEDSELBESTBALANCED;
  else throw runtime_error("Error: unrecognized darp kmedoids mode criterion");

  return  mode;
}


bool parse_config (char* fich) {

   libconfig::Config cfg;
   cfg.setAutoConvert (true);

   try {
      cfg.readFile (fich);
   }
   catch (libconfig::ParseException e) {
     std::cerr << "Error: parsing the configuration file in line (" << e.getLine () << "): " << e.getError() << ". File name: " << fich << std::endl;
      exit (-1);
   }
   catch (libconfig::FileIOException e) {
      std::cerr << "Error: file not found or could not be open." << std::endl;
      exit (-1);
   }

   // First, we read mandatory parameters
   try {
      genType = (const char*) cfg.lookup("genomeType");
      if (genType.compare("VNS") == 0) {
        minEvalWeightsUpdateValue = (double) cfg.lookup("minEvalsWeightUpdateValue");
        maxEvalWeightsUpdateValue = (double) cfg.lookup("maxEvalsWeightUpdateValue");
      }
      else if (genType.compare("TS")  != 0) {throw runtime_error("Unrecognized Genome Type");}

      maxRideValue    = (double) cfg.lookup("maxRideValue");
      vehicleCapacity = (long) cfg.lookup("vehicleCapacity");
      nVehicles       = (long) cfg.lookup("nVehicles");
      timeSumExp      = (int) cfg.lookup("timeSumExp");
      nStops          = (long) cfg.lookup("nStops");
      maxDelay        = (long) cfg.lookup("maxDelay");

      relativeDelays        = (bool) cfg.lookup("relativeDelays");
      relativeDelayConstant = (long) cfg.lookup("relativeDelayConstant");

      loadVWeight = (double) cfg.lookup("loadVWeight");
      TWVWeight   = (double) cfg.lookup("TWVWeight");
      rideVWeight = (double) cfg.lookup("rideVWeight");

      DistMatrixFile = (const char*) cfg.lookup("distMatrixFile");
      RequestsFile   = (const char*) cfg.lookup("requestsFile");

      const char* darpOptCritRead = (const char*) cfg.lookup("darpOptCrit");
      if      (strcmp(darpOptCritRead,"cost")       == 0) optcrit = DARPEvaluator::TRAVELCOST;
      else if (strcmp(darpOptCritRead,"bothdelays") == 0) optcrit = DARPEvaluator::BOTHDELAYS;
      else if (strcmp(darpOptCritRead,"delivery")   == 0) optcrit = DARPEvaluator::DELIVERY;
      else throw runtime_error("Error: unrecognized darp opt criterion");

      initname = (const char*) cfg.lookup("initialization");
      initfunc = initMethod( initname );
      if (initfunc == DARPRegretInsertionInitC::DARPRegretInsertionInit) {
        double dencentr_alpha   = (double) cfg.lookup("decentrAlpha");
        double maxreqsperinitit = (int)    cfg.lookup("maxRequestsPetInitIt");
        DARPRegretInsertionInitC::setInitParameters(dencentr_alpha,maxreqsperinitit);
      }
      if (initfunc == DARPSlackInitC::DARPSlackInit) {
        double maxreqsperinitit = (int) cfg.lookup("maxRequestsPerInitIt");
        DARPSlackInitC::setInitParameters(maxreqsperinitit);
      }

      // reading the shakers lists
      libconfig::Setting& shakergroup = cfg.lookup("shakers");
      for (int i=0; i<shakergroup.getLength(); i++) {
        string name = shakergroup[i]["name"];
        int    size = shakergroup[i]["size"];

        shakers.push_back(DARPVNSShaker::createNewShaker(name,size));
      }
   }
   catch (libconfig::SettingNotFoundException e) {
      std::cerr << "[DARP] Error: one or more mandatory parameter(s) were not found." << std::endl;
      exit (-1);
   }

   // Optional stats

   if (cfg.exists("numRouteCallsStatsDispFreq")) numRouteCallsStatsDispFreq = (int)          cfg.lookup("numRouteCallsStatsDispFreq");
   if (cfg.exists("relativeDelayConstant"))      relativeDelayConstant      = (int)          cfg.lookup("relativeDelayConstant");
   if (cfg.exists("initloadVWeight") )           initloadVWeight            = (double)       cfg.lookup("initloadVWeight");
   if (cfg.exists("initTWVWeight") )             initTWVWeight              = (double)       cfg.lookup("initTWVWeight");
   if (cfg.exists("initrideVWeight") )           initrideVWeight            = (double)       cfg.lookup("initrideVWeight");
   if (cfg.exists("reqDist"))                    reqdistname                = (const char *) cfg.lookup("reqDist");
   if (cfg.exists("localsearch"))                localsearch                = (const char *) cfg.lookup("localsearch");
   if (cfg.exists("routedistributor"))           routedistname              = (const char *) cfg.lookup("routedistributor");
   if (cfg.exists("initkmedoidsnpass"))          initkmedoidsnpass          = (int)          cfg.lookup("initkmedoidsnpass");
   if (cfg.exists("initkmedoidsmode"))           initkmedoidsmode           = getKMedoidsMode( (const char *) cfg.lookup("initkmedoidsmode") );
   if (cfg.exists("kmedoidsroutedistnpass"))     kmedoidsroutedistnpass     = (int)          cfg.lookup("kmedoidsroutedistnpass");
   if (cfg.exists("kmedoidsroutedistmode"))      kmedoidsroutedistmode      = getKMedoidsMode( (const char *) cfg.lookup("kmedoidsroutedistmode") );
   if (cfg.exists("kmedoidsuseconstrpen"))       kmedoidsuseconstrpen       = (bool)         cfg.lookup("kmedoidsuseconstrpen");
   if (cfg.exists("kmedoidsusewaitingpen"))      kmedoidsusewaitingpen      = (bool)         cfg.lookup("kmedoidsusewaitingpen");
   if (cfg.exists("maxlsevals"))                 maxlsevals                 = (int)          cfg.lookup("maxlsevals");
   if (cfg.exists("lsrandomroutes"))             lsrandomroutes             = (bool)         cfg.lookup("lsrandomroutes");
   if (cfg.exists("objfuncinitMethod"))          objfuncinitMethod          = (const char *) cfg.lookup("objfuncinitMethod");
   if (cfg.exists("postinitls"))                 postinitls                 = (bool)         cfg.lookup("postinitls");
   if (cfg.exists("printbestsol"))               printBestSol               = (bool)         cfg.lookup("printbestsol");

   if (kmedoidsusewaitingpen) {
     kmedoidswaitingpenthreshold = (long)   cfg.lookup("waitingpenthreshold");
     kmedoidswaitingpenconstant  = (double) cfg.lookup("waitingpenconstant");
   }

   if (reqdistname != "") isdistributed = CommManager::instance()->isDistributed();

   return true;
}



VerticesList* createInitVertList(VerticesList& vertlist, ReqDistributor& reqdist, ReqDistMatrix& reqdistmatrix) {
  VerticesList* initvertlist = 0;
  if (isdistributed) {
    assert(& reqdist);

    CommManager& comm = * CommManager::instance();

    if (comm.isIslandMaster()) {
      vector<Request>           reqlist = vertlist.getReqsList();
      vector< vector<Request> > selreqs = reqdist.distributeReqs(reqlist, reqdistmatrix, comm.getNumIslands()  );

      //Serialize and send to each slave nodes
      for (int islandpos=1; islandpos<selreqs.size(); islandpos++) {
        vector<int>      serdata;
        vector<Request>& islandsreqs = selreqs[islandpos];
        for (int reqpos=0; reqpos<islandsreqs.size(); reqpos++) {
          serdata.push_back(islandsreqs[reqpos].pickup_vert->id_);
          serdata.push_back(islandsreqs[reqpos].delivery_vert->id_);
        }

        comm.sendIntVectorToIsland(islandpos,serdata);
      }

      initvertlist = new VerticesList( selreqs[0] );
    }
    else { // Slave
      // Receive the list of reqs ids and construct the initvertlist
      vector<int> vertids = comm.receiveIntVectorFromIsland(comm.MASTER_RANK,vertlist.size());
      assert(vertids.size() % 2 == 0);

      vector<Request> islandsreqs = vertlist.getReqsFromVertIds(vertids);

      initvertlist = new VerticesList( islandsreqs );
    }

    assert(initvertlist->size() > 0);
    GARandomSeed( GAGetRandomNoRankSeed(), comm.getMyRank() ); // Hack to avoid the bug that after calling MPI_Recv the random
                                                               // function gets called and all the random calls loose the
                                                               // synchronization

  }
  else{
    initvertlist = & vertlist;
  }



  return initvertlist;
}

extern "C" VNSOp *VNSLocalSearch() {
  VNSOp* ls = 0;
  if      (localsearch == "DARPLocalSearch"){
    if (maxlsevals > 0) ls = new DARPLocalSearch(maxlsevals,lsrandomroutes);
    else                ls = new DARPLocalSearch();
  }
  else if (localsearch == "DARPNOLocalSearch") ls = new DARPNOLocalSearch();
  else throw runtime_error("unrecognized local search");

  return ls;
}

extern "C" void individualInit (GAGenome& g) {
  assert(distMatrix!=0 && vertlist); // Using this global variables, no other way to avoid this
  assert(evaluator != 0);

  // Since they shouldnt have negative values we are using that knowledge to check if the values have been set
  if (initloadVWeight != -1 && initTWVWeight != -1 && initrideVWeight!= -1) {
    assert(initloadVWeight > 0 && initTWVWeight > 0 && initrideVWeight > 0);
    evaluator->setAllWeights(initloadVWeight,initTWVWeight,initrideVWeight);
  }

  DARPGenome& gen = dynamic_cast<DARPGenome&>(g); assert(&gen);

  VerticesList* initvertlist = createInitVertList(*vertlist,*reqdist,*reqdistmatrix);

  if (initfunc == DARPObjFuncBasedInitC::DARPObjFuncBasedInit) {
    if      (objfuncinitMethod == "best")  DARPObjFuncBasedInitC::setMethod(DARPObjFuncBasedInitC::BEST);
    else if (objfuncinitMethod == "first") DARPObjFuncBasedInitC::setMethod(DARPObjFuncBasedInitC::FIRST);
    else {
      throw runtime_error("unrecognized DARPObjFuncBasedInitC init method");
    }
  }

  initfunc(g,*distMatrix,*initvertlist);

  // If specified we apply the ls after the initialization
  if (postinitls) {
    DARPLocalSearch* ls = dynamic_cast<DARPLocalSearch*>( VNSLocalSearch());
    (*ls)(g);
    delete ls;
  }

  int gen_vertices=0;for (int i=0; i<gen.size();i++) gen_vertices += gen.routeLength(i);
  if (gen_vertices != initvertlist->size()) {
    stringstream msg; msg << "Error the constructed genome has different number of vertices, expected=" << initvertlist->size() << " gen=" << gen_vertices;
    throw runtime_error(msg.str());
  }

  if (initfunc == DARPObjFuncBasedInitC::DARPObjFuncBasedInit) {
    evaluator->setAllWeights(loadVWeight,TWVWeight,rideVWeight);
  }

  // TODO: delete here initvertlist but before we need to create a clone method in VerticesList
}

extern "C" void populationInit (GAPopulation& pop, double perc) {
  throw runtime_error("[populationInit] Error: initialization method not defined.");
}

extern "C" double objective(GAGenome& g) {
  DARPGenome& gen = dynamic_cast<DARPGenome&>(g);

  return evaluator->score(gen);
}

int getNVehicles4Node() {
  vector<int> nvehicles_assignment;

  CommManager& comm = * CommManager::instance();
  if (isdistributed) { // Distributed algorithm
    nvehicles_assignment.resize(comm.getNumIslands());
    int nvehiclespernode = nVehicles/comm.getNumIslands();

    for (int i=0; i<comm.getNumIslands(); i++) {
      nvehicles_assignment[i] = nvehiclespernode;
    }

    int nvehicles_assigned = nvehiclespernode*comm.getNumIslands(); assert(nvehicles_assigned<=nVehicles);

    // We assign the possible remaining vehicles to the first nodes. We cannot use a random criterion since
    // each node has a different seed. If we try to use a random approach and them send the values to each island
    // could be a little bit overkilling since later, the assignment is made from a random selection of vertices.
    for (int i=nvehicles_assigned; i<nVehicles; i++) {
      nvehicles_assignment[ nVehicles - i - 1 ] += 1;
    }
  }
  else {         // Sequential Algorithm
    for (int i=0; i<comm.getNumIslands(); i++) nvehicles_assignment.push_back(nVehicles);
  }

#ifdef DEBUG
  int sum_vehicles=0; for(int i=0;i<nvehicles_assignment.size(); i++) sum_vehicles += nvehicles_assignment[i];
  assert(sum_vehicles == nVehicles);
#endif

  return nvehicles_assignment[comm.getMyRank()];
}


void printOptions(int nvehicles4node) {
  if (CommManager::instance()->isIslandMaster()) {

    cout << endl << "Specific parameter values: " << endl;
    cout << "Genome:                   = " << genType << endl;
    if (genType.compare("VNS") == 0) {
      cout << "minEvalWeightsUpdateValue = " << minEvalWeightsUpdateValue << endl;
      cout << "maxEvalWeightsUpdateValue = " << maxEvalWeightsUpdateValue << endl;
    }
    cout << "maxRideValue               = " << maxRideValue << endl;
    cout << "vehicleCapacity            = " << vehicleCapacity << endl;
    cout << "Total nVehicles            = " << nVehicles << endl;
    cout << "nVehicles of Node          = " << nvehicles4node << endl;
    cout << "nStops                     = " << nStops << endl;
    cout << "maxDelay                   = " << maxDelay << endl;
    cout << "darpOptCrit                = " << optcrit << endl;
    cout << "relativeDelays             = " << relativeDelays << endl;
    cout << "relativeDelayConstant      = " << relativeDelayConstant << endl;
    cout << "initialization             = " << initname << endl;
    cout << "postinitls                 = " << postinitls << endl;
    cout << "printBestSol               = " << printBestSol << endl;
    cout << "objfuncinitMethod          = " << objfuncinitMethod << endl;
    cout << "loadVWeight                = " << loadVWeight << endl;
    cout << "TWVWeight                  = " << TWVWeight << endl;
    cout << "rideVWeight                = " << rideVWeight << endl;

    if (initloadVWeight != -1 && initTWVWeight != -1 && initrideVWeight != -1) {
      cout << "initloadVWeight            = " << initloadVWeight << endl;
      cout << "initTWVWeight              = " << initTWVWeight << endl;
      cout << "initrideVWeight            = " << initrideVWeight << endl;
    }
    if (initfunc == DARPRegretInsertionInitC::DARPRegretInsertionInit ) {
      cout << "decentralization Alpha     = " << DARPRegretInsertionInitC::decentrAlpha() << endl;
      cout << "maxReqsUsedInEachInitIt    = " << DARPRegretInsertionInitC::maxReqsForRegretIns() << endl;
    }
    if (initfunc == DARPSlackInitC::DARPSlackInit) {
      cout << "maxReqsUsedInEachInitIt    = " << DARPSlackInitC::maxReqsPerInitIt() << endl;
    }

    if (numRouteCallsStatsDispFreq>0) {
      cout << "numRouteCallsStatsDispFreq = " << numRouteCallsStatsDispFreq << endl;
    }

    if (reqdistname != "" ) {
      cout << "RequestDistributor         = " << reqdistname << endl;
      if (reqdistname == "kmedoids") {
        cout << " KMedoids Request distributor mode = "  << initkmedoidsmode  << endl;
        cout << " KMedoids Request distributor npass = " << initkmedoidsnpass << endl;
        cout << " KMedoids use constraint penalizations for matrix computation = " << kmedoidsuseconstrpen << endl;
      }

      if (kmedoidsusewaitingpen) {
        cout << "KMedoids waiting penalization threshold = " << kmedoidswaitingpenthreshold << endl;
        cout << "KMedoids waiting penalization constant = "  << kmedoidswaitingpenconstant << endl;
      }
    }

    if (routedistname != "") {
      cout << "RouteDistributor          =  " << routedistname << endl;
      if (routedistname == "kmedoids") {
        cout << "KMedoids Route distributor mode = "  << kmedoidsroutedistmode << endl;
        cout << "KMedoids Route distributor npass = " << kmedoidsroutedistnpass << endl;
      }
    }

    if (maxlsevals != -1) cout << "maxlsevals                = " << maxlsevals << endl;
    cout << "ls random routes          = " << lsrandomroutes << endl;

    cout << endl;

  }
}

extern "C" GAGenome* defineProblem () {
  // Parse the config file for the problem
  if (!parse_config (GAEDAConfig::handle()->getProblemData ())) {
    throw runtime_error ("Error: A valid configuration file must be provided for the TSP problem.");
  }

  Vertex::setTimeConstants(relativeDelays, relativeDelayConstant, maxRideValue, maxDelay);

  distMatrix = new DistMatrix  (DistMatrixFile.c_str(), nStops);
  vertlist   = new VerticesList(RequestsFile.c_str (), *distMatrix);
  costMatrix = new OnlyDistCostMatrix  (*distMatrix, *vertlist);
  evaluator  = new DARPEvaluator(vehicleCapacity,*vertlist,*costMatrix,timeSumExp,loadVWeight,TWVWeight,rideVWeight,optcrit);

  if (numRouteCallsStatsDispFreq > 0) evaluator->displayStatsEachNumRouteCalls(numRouteCallsStatsDispFreq);

  //costMatrix->pruneArcs(*evaluator,*vertlist);

  DARPGenome::setDARPEvaluator(evaluator);
  DARPVNSOp::setDARPEvaluator(evaluator);

  // Note the genome uses the evaluator for initializing its values so it should be place after the evaluator initialization

  int nvehicles4node = getNVehicles4Node();

  if      (genType.compare("VNS") == 0) gen = new VNSDARPGenome(nvehicles4node, minEvalWeightsUpdateValue, maxEvalWeightsUpdateValue, objective);
  else if (genType.compare("TS")  == 0) gen = new TSDARPGenome (nvehicles4node, objective);

  gen->initializer (individualInit);

  printOptions(nvehicles4node);

  return gen;
}

extern "C" void configureAlg (Algorithm& alg) {

  // This cannot be done in defineproblem because the optcriterion has not been read yet
  // nor in individualInit because it is called later than configureAlg and we need to have
  // the reqdistmatrix for the distrouting alg

  if (isdistributed) {
    reqdist                 = createReqDist(reqdistname, *gen, *costMatrix, *evaluator);
    vector<Request> reqlist = vertlist->getReqsList();
    CommManager*    comm    = CommManager::instance(); assert(comm);


    if (reqdistname == "kmedoids" || routedistname == "kmedoids") {

      reqdistmatrix = new ParReqDistMatrix(reqlist, *costMatrix, *evaluator,
                                           initTWVWeight, initrideVWeight, initloadVWeight, kmedoidsuseconstrpen,
                                           kmedoidsusewaitingpen, kmedoidswaitingpenthreshold, kmedoidswaitingpenconstant,*comm);
    }

    if (routedistname != "" ) {
      if      ( routedistname == "random" )   routedist = new RandomRouteDistributor(*gen,comm->getNumIslands());
      else if ( routedistname == "kmedoids" ) {
        assert(reqdistmatrix);
        routedist = new KMedoidsRouteDistributor(reqdistmatrix,kmedoidsroutedistnpass,kmedoidsroutedistmode,*gen,comm->getNumIslands());
      }
      else throw runtime_error("Unrecognized route distributor");

      DistRoutingAlg& distroutingalg = dynamic_cast<DistRoutingAlg&>(alg); assert(&distroutingalg);
      distroutingalg.setRouteDistributor(*routedist);
    }

  }
}


extern "C" const char *describeProblem (void) {
   return "DARP-TW problem.";
}


extern "C" GAGenome::OptCriterion optCriterion(){
  return GAGenome::MINIMIZATION;
}

void writeBIFile(DARPGenome& gen) {
  stringstream filename;
  filename << dynamic_cast<GAFileLogger*>(GALogger::instance())->path() << ".BI";

  ofstream os;
  os.open(filename.str().c_str(),ofstream::out | ofstream::trunc);

  for (int routepos=0; routepos<gen.length(); routepos++) {
    std::vector<int> route (gen.routeLength(routepos), 0);
    for (int vertpos=0; vertpos<gen.routeLength(routepos); vertpos++) route[vertpos] = gen.gene(routepos, vertpos);

    std::map<int,long> arrivalTimes;
    std::map<int,long> waitingTimes;
    std::map<int,long> beginningServiceTimes;
    std::map<int,long> departureTimes;
    std::map<int,long> rideTimes;
    std::map<int,long> forwardTimeSlacks;
    std::map<int,long> loadsWhenLeavingVertex;
    long cost;
    long loadViolation, TWViolation, rideViolation;
    double pickupDelay, deliveryDelay;

    evaluator->evalRoute(route,arrivalTimes,waitingTimes,beginningServiceTimes,departureTimes,rideTimes,forwardTimeSlacks,loadsWhenLeavingVertex,cost,loadViolation,TWViolation,rideViolation,pickupDelay,deliveryDelay);

    for (int vertpos=0; vertpos<gen.routeLength(routepos); vertpos++)  {
      int     vert_id = gen.gene(routepos,vertpos);
      Vertex& vert    = vertlist->getVertex(vert_id);

      long        readTime = vert.isPickUp() ? departureTimes[vert_id] : arrivalTimes[vert_id];
      struct tm * gtime    = localtime(&readTime);
      long        time     = gtime->tm_hour * 3600 + gtime->tm_min * 60 + gtime->tm_sec;

      int pos     = vert.pos_;
      string type = vert.isPickUp() ? "subida" : "bajada";

      os << routepos << "," << time << "," << pos << "," << abs(vert_id) << "," << type << endl;
    }
  }

  os.close();
}

extern "C" bool postprocess (GAPopulation* pop, int rank, int repetition) {
  if (printBestSol && rank == CommManager::MASTER_RANK) {

    DARPGenome& gen1 = dynamic_cast<DARPGenome&>(pop->individual(RoutingAlg::SBESTPOS));
    DARPGenome& gen2 = dynamic_cast<DARPGenome&>(pop->individual(RoutingAlg::SPOS));

    DARPGenome* sel_gen   = & gen1;
    DARPGenome* other_gen = & gen2;

    if (GAGenome::compareScores( gen2.score(), gen1.score() ) == GAGenome::BETTER )  {
      sel_gen   = & gen2;
      other_gen = & gen1;
    }

    if (!sel_gen->feasible() and other_gen->feasible() ) {
      sel_gen = other_gen;
    }

    cout << "best selected solution:" << endl;
    cout << *sel_gen << endl;

    writeBIFile(*sel_gen);

    for (unsigned i = 0; i < sel_gen->length(); i++) {
      std::vector<int> route (sel_gen->routeLength(i), 0);
      for (unsigned j = 0; j < sel_gen->routeLength(i); j++)
        route[j] = sel_gen->gene(i, j);

      std::cout << "***************************************************************************" << std::endl;
      std::cout << " ==> Route " << i << ": " << std::endl;
      std::cout << "***************************************************************************" << std::endl;

      {
        long cost;
        long loadViolation, TWViolation, rideViolation;
        double pickupDelay, deliveryDelay;
        evaluator->evalRoute(route, cost, loadViolation, TWViolation, rideViolation, pickupDelay, deliveryDelay, true,true);
      }
    }
  }

  // Free up memory
  delete distMatrix;
  delete vertlist;
  delete costMatrix;

  return true;
}

extern "C"  std::vector< VNSShaker* > getVNSShakersList () {
  return shakers;

//  shakers.push_back( new SwapNeighborhood(1) );
//  shakers.push_back( new ChainNeighborhood(1) );
//
//  shakers.push_back( new SwapNeighborhood(2) );
//  shakers.push_back( new ChainNeighborhood(2) );
//
//  shakers.push_back( new GreedyMoveNeighborhood(2) );
//  shakers.push_back( new GreedySwapNeighborhood(2) );
//
//  shakers.push_back( new SwapNeighborhood(3) );
//  shakers.push_back( new ChainNeighborhood(3) );
//
//  //shakers.push_back( new SwapNeighborhood(4) );
//  shakers.push_back( new ChainNeighborhood(4) );
//
//  shakers.push_back( new GreedyMoveNeighborhood(4) );
//  shakers.push_back( new GreedySwapNeighborhood(4) );
//
//  shakers.push_back( new SwapNeighborhood(5) );
//  shakers.push_back( new ChainNeighborhood(5) );
//
//  //shakers.push_back( new SwapNeighborhood(6) );
//  //shakers.push_back( new ChainNeighborhood(6) );
//  shakers.push_back( new GreedyMoveNeighborhood(6) );
//  shakers.push_back( new GreedySwapNeighborhood(6) );
//
//  shakers.push_back( new GreedyMoveNeighborhood(8) );
//  shakers.push_back( new GreedySwapNeighborhood(8) );
//
//  shakers.push_back( new CheckAllNaturalSeqsCombsNeighborhood(0) );
////  //shakers.push_back( new ZeroSplitNeighborhood(0) );

  return shakers;
}

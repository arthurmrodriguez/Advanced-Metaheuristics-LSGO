#include "DistRoutingAlg.h"
#include "../logger/GALogger.h"
#include "../../gaeda/garandom.h"
#include <sstream>
#include <memory>

using namespace std;

DistRoutingAlg::DistRoutingAlg(RoutingAlg&       alg,
                               int               mig_period,
                               emigrantRouteType emigrant_route,
                               GAGenome&         gen,
                               CommManager&      comm_manager) : alg_(alg.clone()),
                                                                 rdistributor_(0),
                                                                 mig_period_(mig_period),
                                                                 emigrant_route_(emigrant_route),
                                                                 comm_manager_(comm_manager) {

  if (comm_manager.isIslandMaster()) {
    cout << endl;
    cout << "Routing Distributed Algorithm parameter values "<< endl;
    cout << "Migration Period  : " << mig_period << endl;
    cout << "Migration         : "; if (emigrant_route == BEST) cout << "best solution"; else cout << "local solution"; cout << endl;
  }
}

DistRoutingAlg::DistRoutingAlg(const DistRoutingAlg& other) : alg_(0),
                                                              rdistributor_(0),
                                                              comm_manager_(other.comm_manager_) {
  copy(other);
}

DistRoutingAlg::~DistRoutingAlg() {
  delete alg_;
  if (rdistributor_) delete rdistributor_;
}

/*
 * TODO: Refactor from EAIslandsModel since it has many parts in common
 */
void DistRoutingAlg::evolve() {
  assert(&(alg_->bestSol()));
  assert(&(alg_->actualSol()));

//  /*LOG*/ GALogger::instance()->appendLogMessage("","before starting  test",GALogger::only_stats);
//
//  {stringstream msg; msg << "best ind is " << endl << alg_->bestSol() << endl;
//  /*LOG*/ GALogger::instance()->appendLogMessage("",msg.str(),GALogger::only_stats);
//  }

  while (true) {
    if ( !alg_->done() ) {
      do {
        ///*LOG*/ GALogger::instance()->appendLogMessage("","before step inner alg",GALogger::only_stats);
        alg_->step();
        ///*LOG*/ GALogger::instance()->appendLogMessage("","after step inner alg",GALogger::only_stats);

      } while ( alg_->statistics().generation() % mig_period_ != 0 && !alg_->done() );
    }
    else {
      ///*LOG*/ GALogger::instance()->appendLogMessage("","before update stats",GALogger::only_stats);
     alg_->updateNoStepsStats();
     ///*LOG*/ GALogger::instance()->appendLogMessage("","after update stats",GALogger::only_stats);
    }

    if (haveAllIslandsConverged()) break; // Not put it in the while condition to avoid having always a migration at the end 

    ///*LOG*/ GALogger::instance()->appendLogMessage("","before migration",GALogger::only_stats);
    doMigration();
    ///*LOG*/ GALogger::instance()->appendLogMessage("","after migration",GALogger::only_stats);
  }
  assert(&(alg_->bestSol()));
  assert(&(alg_->actualSol()));

  ///*LOG*/ GALogger::instance()->appendLogMessage("","before unite",GALogger::only_stats);
  uniteRoutes(); // Finally we group all the routes in the master node
  ///*LOG*/ GALogger::instance()->appendLogMessage("","after unite",GALogger::only_stats);
  if (comm_manager_.isIslandMaster() ) GALogger::instance()->appendStats("End of evolution stats");
}

void DistRoutingAlg::doMigration() {
  assert(rdistributor_);

  /*LOG*/ GALogger::instance()->appendLogMessage("","doing migration",GALogger::only_stats);

  if (comm_manager_.isIslandMaster() ) {

    auto_ptr< vector<GAGenome*> > received_inf ( comm_manager_.receiveOneIndFromAllSlaveIslands() );

    //{stringstream msg; msg << "received following best inds: " << endl;
    //for (int i=0; i<received_inf->size(); i++) msg << "score value is " << (*received_inf)[i]->score() << endl <<  *( (*received_inf)[i]) << endl << endl;
    ///*LOG*/ GALogger::instance()->appendLogMessage("DistRouting",msg.str(),GALogger::only_stats);}

    received_inf->push_back(alg_->bestSol().cloneWithoutEmptyRoutes()); // We receive all the routes from the slave islands and we add the master's route

    //{stringstream msg; msg << "including best ind: " << endl;
    //msg << alg_->bestSol() << endl;
    ///*LOG*/ GALogger::instance()->appendLogMessage("DistRouting",msg.str(),GALogger::only_stats);}

    vector<RoutingGenome*> inds2send = rdistributor_->assignRoutes2Nodes(*received_inf); // No need to delete inds2send genomes!!
    assert(inds2send.size() == comm_manager_.getNumIslands());

    //{stringstream msg; msg << "printing assignations: " << endl;
    //for (int i=0; i<inds2send.size(); i++) {
    //  int a = i;
    //  msg << "to node " << i << endl;
    //  msg << * inds2send[i] << endl; }
    ///*LOG*/ GALogger::instance()->appendLogMessage("DistRouting",msg.str(),GALogger::only_stats);}

    alg_->setNewRoutes( * inds2send[comm_manager_.MASTER_RANK] ); // Set routes for the master rank

//    {stringstream msg; msg << "new routes are: " << endl;
//    msg << "best sol: " << endl << alg_->bestSol() << endl;
//    msg << "actual sol: " << endl << alg_->actualSol() << endl;
//    /*LOG*/ GALogger::instance()->appendLogMessage("DistRouting",msg.str(),GALogger::only_stats);}


    for (int i=0; i<inds2send.size(); i++) {                      // Send the routes for the remaining nodes
      if (i==comm_manager_.MASTER_RANK) continue;
      comm_manager_.sendIndividual( * inds2send[i], i);
    }

    for (int i=0; i<received_inf->size(); i++) delete (*received_inf)[i];
  }
  else {   // Slave node
//    {stringstream msg; msg << "sending selected ind to master " << selectIndForSendingRoutes()  << endl;
//    /*LOG*/ GALogger::instance()->appendLogMessage("DistRouting",msg.str(),GALogger::only_stats);}

    RoutingGenome&          selectedInd = selectIndForSendingRoutes();
    auto_ptr<RoutingGenome> genNoEmptyRoutes ( selectedInd.cloneWithoutEmptyRoutes() );
    RoutingGenome&          gen2Send = *genNoEmptyRoutes;
    comm_manager_.sendIndToMaster( gen2Send );

    GAGenome* received_inf = comm_manager_.receiveIndFromMaster();

//    {stringstream msg; msg << "received ind from master: " << endl;
//    msg << *received_inf << endl;
//    /*LOG*/ GALogger::instance()->appendLogMessage("DistRouting",msg.str(),GALogger::only_stats);}

    alg_->setNewRoutes(*received_inf);

//    {stringstream msg; msg << "new routes are: " << endl;
//    msg << "best sol: " << endl << alg_->bestSol() << endl;
//    msg << "actual sol: " << endl << alg_->actualSol() << endl;
//    /*LOG*/ GALogger::instance()->appendLogMessage("DistRouting",msg.str(),GALogger::only_stats);}


    delete received_inf;
  }

  GARandomSeed( GAGetRandomNoRankSeed(), comm_manager_.getMyRank() ); // Hack to avoid the bug that after calling MPI_Recv the random function gets called
                                                                      // and all the random calls loose the synchronization
}

RoutingGenome& DistRoutingAlg::selectIndForSendingRoutes() {
  switch(emigrant_route_) {
    case BEST:   return alg_->bestSol();   break;
    case ACTUAL: return alg_->actualSol(); break;
    default: throw runtime_error("emigrant route not recognized");
  }
}

bool DistRoutingAlg::haveAllIslandsConverged(){
  bool haveIFinished = alg_->done();
  auto_ptr< list<int> > converged_islands_pos ( comm_manager_.sendAndReceiveConvergence( haveIFinished ) );

  GARandomSeed( GAGetRandomNoRankSeed(), comm_manager_.getMyRank() ); // Hack to avoid the bug that after calling MPI_Recv the random function gets called
                                                                      // and all the random calls loose the synchronization


//  /*LOG*/ stringstream message; message << "n# converged islands=" << converged_islands_pos->size();
//  GALogger::instance()->appendLogMessage("Convergence", message.str());

  // Note: we need to to do it with the haveIFinished variables because if we used the 
  // alg->done() method for both functions (sendAndReceiveConvergence and here below) we
  // could generate a running condition if we use the stop criterion based on time
  return ( (int) converged_islands_pos->size() == comm_manager_.getNumIslands() - 1  && haveIFinished);  // since we are not receiving a notification from ourselves
}

void DistRoutingAlg::copy(const DistRoutingAlg& other) {
  const DistRoutingAlg& otheralg = dynamic_cast<const DistRoutingAlg&>(other);
  assert(&otheralg);

  if (alg_) delete alg_;
  alg_ = otheralg.alg_->clone();

  if (rdistributor_) delete rdistributor_;
  rdistributor_ = otheralg.rdistributor_->clone();

  mig_period_     = otheralg.mig_period_;
  emigrant_route_ = otheralg.emigrant_route_;

  comm_manager_ = otheralg.comm_manager_;
}

DistRoutingAlg* DistRoutingAlg::clone() {
  return new DistRoutingAlg(*this);
}

void DistRoutingAlg::uniteRoutes() {
//  {stringstream msg; msg << "en distroutingalg uniteroutes : " << endl;
//  /*LOG*/ GALogger::instance()->appendLogMessage("DistRouting",msg.str(),GALogger::only_stats);}

//  {stringstream msg; msg << "en distroutingalg uniteroutes best solutions are: " << endl;
//  msg << "best sol: " << endl << alg_->bestSol() << endl;
//  msg << "actual sol: " << endl << alg_->actualSol() << endl;
//  /*LOG*/ GALogger::instance()->appendLogMessage("DistRouting",msg.str(),GALogger::only_stats);}


  if ( comm_manager_.isIslandSlave() ) {
//    {stringstream msg; msg << "sending the individual to the master: " << endl;
//    /*LOG*/ GALogger::instance()->appendLogMessage("DistRouting",msg.str(),GALogger::only_stats);}
    comm_manager_.sendIndToMaster(alg_->bestSol());
  }
  else {
//    {stringstream msg; msg << "master: receving the routes " << endl;
//    /*LOG*/ GALogger::instance()->appendLogMessage("DistRouting",msg.str(),GALogger::only_stats);}

    auto_ptr< vector<GAGenome*> > inds ( comm_manager_.receiveOneIndFromAllSlaveIslands() );
    inds->push_back(alg_->bestSol().clone()); // We receive all the routes from the slave islands and we add the master's route

//    {stringstream msg; msg << "master: received the routes " << endl;
//    /*LOG*/ GALogger::instance()->appendLogMessage("DistRouting",msg.str(),GALogger::only_stats);}

//    {stringstream msg; msg << "masteR: before setting routes: " << endl;
//    msg << "best sol: " << endl << alg_->bestSol() << endl;
//    msg << "actual sol: " << endl << alg_->actualSol() << endl;
//    /*LOG*/ GALogger::instance()->appendLogMessage("DistRouting",msg.str(),GALogger::only_stats);}

    alg_->setNewRoutes(*inds);

//    {stringstream msg; msg << "masteR: after setting routes: " << endl;
//    msg << "best sol: " << endl << alg_->bestSol() << endl;
//    /*LOG*/ GALogger::instance()->appendLogMessage("DistRouting",msg.str(),GALogger::only_stats);}


//    {stringstream msg; msg << "masteR: after setting routes: " << endl;
//    msg << "actual sol: " << endl << alg_->actualSol() << endl;
//    /*LOG*/ GALogger::instance()->appendLogMessage("DistRouting",msg.str(),GALogger::only_stats);}

    for (int i=0; i<inds->size(); i++) delete (*inds)[i];

//    {stringstream msg; msg << "master: after deleting individuals" << endl;
//    msg << "best sol: " << endl << alg_->bestSol() << endl;
//    msg << "actual sol: " << endl << alg_->actualSol() << endl;
//    /*LOG*/ GALogger::instance()->appendLogMessage("DistRouting",msg.str(),GALogger::only_stats);}
  }

}

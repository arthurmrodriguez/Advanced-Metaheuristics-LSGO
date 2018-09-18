#include "CommManager.h"
#include "../logger/GALogger.h"
#include <mpi.h>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <memory>
#include <algorithm>
#include <stdexcept>
#include "../GAGenealogy.h"
#include "../GAEDAConfig.h"

using namespace std;

CommManager* CommManager::comm_instance_ptr = NULL; //TODO QUITAR

const int MAX_INDS_SENT = 100000; // The limit of individuals that can be sent

const unsigned int TAG_SEND_SIZE          = 3;
const unsigned int TAG_SEND_DATA          = 4;
const unsigned int TAG_SEND_NOTIFICATION  = 5;
const unsigned int TAG_SEND_FIT           = 6;
const int          YES_NOTIFICATION_VALUE = 1;
const int          NO_NOTIFICATION_VALUE  = 0;

int argc_ = 0;
char** argv_ = (char**) 0;

CommManager::CommManager(int argc, char* argv[])
  :genome_(NULL){
  comm_instance_ptr = this;  //TODO QUITAR
  argc_ = argc;
  argv_ = argv;
  initialized_ = -1;
}

CommManager::~CommManager(){
  if (initialized_ == GAEDAConfig::SYNC)
    MPI_Finalize ();
  if (genome_)
    delete (genome_);
}

void CommManager::setGenome(GAGenome* genome){
  assert(initialized_ == GAEDAConfig::SYNC || initialized_ == GAEDAConfig::SEQ);
  if (genome_)
    delete genome_;
  genome_ = genome->clone();
}

double CommManager::getTime () const{
  return MPI_Wtime();
}

bool CommManager::isDistributed()  const {
  return getNumIslands() > 1;
}

bool CommManager::isIslandMaster() const{
  return ( my_rank_ == MASTER_RANK );
}

bool CommManager::isIslandSlave() const{
  assert(initialized_ == GAEDAConfig::SYNC || initialized_ == GAEDAConfig::SEQ);
  return ( my_rank_ != MASTER_RANK );
}

int CommManager::getMyRank() const{
  assert(initialized_ == GAEDAConfig::SYNC || initialized_ == GAEDAConfig::SEQ);
  return my_rank_;
}

int CommManager::getFirstRank() const {
  assert(initialized_ == GAEDAConfig::SYNC || initialized_ == GAEDAConfig::SEQ);
  return 0;
}

int CommManager::getLastRank() const {
  assert(initialized_ == GAEDAConfig::SYNC || initialized_ == GAEDAConfig::SEQ);
  return getNumIslands() - 1;
}

int CommManager::getNumIslands () const{
  assert(initialized_ == GAEDAConfig::SYNC || initialized_ == GAEDAConfig::SEQ);
  return num_islands_;
}

// This method is called from the master node
void CommManager::sendIndsToSlaveIslands(vector<GAGenome*>& individuals) const{
  assert(initialized_ == GAEDAConfig::SYNC);
  assert(individuals.size() > 0);
  assert( my_rank_ == MASTER_RANK );

  GAPopulation  pop;         // Note: the destructor of population calls the destructor of every individual of the pop
  unsigned int  medoids_size = individuals.size();
  for (unsigned int i=0; i<medoids_size; i++) {
    const GAGenome& medoid_i = *(individuals[i]);
    pop.add(medoid_i);      // This creates  a cloned population of medoids
  }

  // Send population of inds to every island
  for (int i = MASTER_RANK+1; i < num_islands_; i++) {                         // We don't send the pop to ourselves (master)
    sendPopulation(pop,i);
  }
}

void CommManager::sendIndToMaster(GAGenome& individual) const{
  assert(initialized_ == GAEDAConfig::SYNC);
  assert( my_rank_ != MASTER_RANK );
  sendIndividual(individual, MASTER_RANK);
}

void CommManager::sendIndividual(GAGenome& ind, int dest_island_rank) const{
  assert(initialized_ == GAEDAConfig::SYNC);
  // We have to do this little maneuver to make a 1-genome population, because
  // the constructor of GAPopulation will only copy the genome attributes
  GAPopulation pop;
  pop.add (ind);

  sendPopulation(pop,dest_island_rank);
  {/*LOG*/ stringstream message; message << "Sended ind " << ind << "to island " << dest_island_rank << endl;
           GALogger::instance()->appendLogMessage("CommManager:sendingIndividual", message.str() ); }
}

void CommManager::sendPopulation(GAPopulation& pop, int dest_island_rank) const {
  assert(initialized_ == GAEDAConfig::SYNC);

  // BEGIN: Genealogy
  // Compute the age of each genome
  pop.computeAge();
  // END: Genealogy

  /*LOG*/ stringstream message; message << "Sending population to island " << dest_island_rank << endl;
  /*LOG*/ GALogger::instance()->appendPopulation("CommManager:sendingPopulation", message.str(), pop );

  ostringstream os(ios::binary);

  // serialize the number of individuals
  int pop_size = pop.size();
  os.write((char *) &pop_size, sizeof( pop_size ) );

  for (int i=0; i<pop_size; i++) pop.individual(i).writeObject(os);

//   ostringstream p(ios::binary);
//   pop.individual(0).writeObject(p);

  // for (int i=0; i < 1000; i++) {
//     std::cout << "Antes de nuestro guiri9gay" << std::endl;
//     const string os_tmp = p.str();
//     istringstream is;
//     is.str(os_tmp);
//     GAGenome& ref = const_cast<GAGenome&> (pop.individual(0));
//     ref.readObject(is);
//     std::cout << "Lectura realizada correctamente" << std::endl;
//   }

  //std::cout << "sendPopultaion after WriteObject " << os.str() << std::endl;

  const string os_str = os.str();  // Note if we dont save in a string and we directly
                                   // call os.str().data() and then os.str().size(), the string
                                   // gets corrupted by the size method which writes the size value
                                   // in the same region as the first os.str().data() returned

  //std::cout << "sendPopulation after cast " << os.str() << std::endl;
  char* data = const_cast<char*> ( os_str.data() );
  int   size = os_str.size();

 //  std::cout << "sendPopulation " << size << std::endl;
//   for (int i=0; i < size; i++)
//     std::cout << (int)data[i];
//   std::cout << std::endl;

  MPI_Send (&size, 1, MPI_INT, dest_island_rank, TAG_SEND_SIZE, MPI_COMM_WORLD);

  if (size > 0 ) MPI_Send (data, size, MPI_BYTE, dest_island_rank, TAG_SEND_DATA, MPI_COMM_WORLD);
}

void CommManager::sendPopToMaster(GAPopulation& pop) const{
  assert(initialized_ == GAEDAConfig::SYNC);
  sendPopulation(pop,MASTER_RANK);
}

void CommManager::sendPopsToSlaveIslands( const vector<GAPopulation*>& pops) const{
  assert(initialized_ == GAEDAConfig::SYNC);
  // pops should contain ALL the islands populations although we are only sending to the slave islands. The reason for this is
  // to simpligy the code in the ea island cluster model
  assert( (int) pops.size() == num_islands_ );

  for (int i = MASTER_RANK+1; i < num_islands_; i++) {                         // We don't send the pop to ourselves (master)
    CommManager::sendPopulation(*pops[i],i);
  }
}

void CommManager::sendDoubleVectorToMaster( vector<double>& data ) const {
  // trick since the C++ spec now guarantees that vectors store their elements contiguously:
  MPI_Send ( &(data[0]), data.size(), MPI_DOUBLE, MASTER_RANK, TAG_SEND_NOTIFICATION, MPI_COMM_WORLD);
}

void CommManager::sendIntVectorToIsland( int rank, vector<int>& data ) const {
  MPI_Send ( &(data[0]), data.size(), MPI_INTEGER, rank, TAG_SEND_NOTIFICATION, MPI_COMM_WORLD);
}

GAPopulation* CommManager::receivePopFromMaster() const{
  assert(initialized_ == GAEDAConfig::SYNC);
  return receivePopulation(MASTER_RANK);
}

GAGenome* CommManager::receiveIndFromMaster() {
  return receiveIndividual(MASTER_RANK);
}

vector<double> CommManager::receiveDoubleVectorFromIsland(int rank, int maxsize) {
  vector<double> data(maxsize);

  MPI_Status status;
  // trick since the C++ spec now guarantees that vectors store their elements contiguously:
  MPI_Recv( &(data[0]), maxsize, MPI_DOUBLE, rank, TAG_SEND_NOTIFICATION, MPI_COMM_WORLD, &status);

  // We check how many were received and adapt the vector structure size
  int numReceived;
  MPI_Get_count(&status, MPI_DOUBLE, &numReceived);
  data.resize(numReceived);

  return data;
}

vector<int> CommManager::receiveIntVectorFromIsland(int rank, int maxsize) {
  vector<int> data(maxsize);

  MPI_Status status;
  // trick since the C++ spec now guarantees that vectors store their elements contiguously:
  MPI_Recv( &(data[0]), maxsize, MPI_INTEGER, rank, TAG_SEND_NOTIFICATION, MPI_COMM_WORLD, &status);

  // We check how many were received and adapt the vector structure size
  int numReceived;
  MPI_Get_count(&status, MPI_INTEGER, &numReceived);
  data.resize(numReceived);


  return data;
}

vector<GAGenome*>* CommManager::receiveOneIndFromAllSlaveIslands(){
  assert(initialized_ == GAEDAConfig::SYNC);
  assert( my_rank_ == MASTER_RANK );

  vector<GAGenome*>* inds = new vector<GAGenome*>(num_islands_-1);

  for (int i = 1; i < num_islands_; i++)  (*inds)[i-1] = receiveIndividual(i);

  /*LOG*/ GALogger::instance()->appendPopulation( "CommManager::receiveIndFromAllIslands", "the inds received are:", *inds );

  return inds;
}

GAPopulation* CommManager::receivePopsFromAllSlaveIslands (){
  assert(initialized_ == GAEDAConfig::SYNC);
  GAPopulation* all_pops = new GAPopulation();

  for (int rank = MASTER_RANK+1; rank < num_islands_; rank++)  {
    auto_ptr<GAPopulation> pop_received ( receivePopulation(rank) );
    for (unsigned i=0; i<pop_received->size(); i++) all_pops->add( pop_received->individual(i) ); // not very efficient lots of clones
  }

  return all_pops;
}

GAGenome* CommManager::receiveIndividual(int source){
  assert(initialized_ == GAEDAConfig::SYNC);
  auto_ptr<GAPopulation> pop ( receivePopulation(source) );         // Convert the serialized data to a GAPopulation object

  GAGenome* received_individual = pop->individual(0).clone();

  return received_individual;
}

GAPopulation* CommManager::receivePopulation(int source) const{
  assert(initialized_ == GAEDAConfig::SYNC);

  int   size = 0;

  MPI_Recv(&size, 1, MPI_INT,source,TAG_SEND_SIZE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  GAPopulation* pop = NULL;
  if (size > 0) {
    char* data = (char *) malloc( size );
    MPI_Recv(data, size, MPI_BYTE, source, TAG_SEND_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

//     std::cout << "receivePopultaion " << size << std::endl;
//     for (int i=0; i < size; i++)
//       std::cout << (int)data[i];
//     std::cout << std::endl;

    // Little hack. We cannot do str_data(data) since the string has escape characters that are interpreted in the
    // constructor of the string. Since the string has binary data and probably lots of escaped characters, we need to read
    // them by hand.
    string str_data(size,size); for (int i=0; i<size; i++) str_data[i] = data[i];

    istringstream is(str_data, ios::binary);

    int pop_size;
    is.read( (char*) &pop_size, sizeof(int) );

    assert(pop_size >= 0 && pop_size < MAX_INDS_SENT);

    pop = new GAPopulation();
    GAGenome* gen_tmp;

    for (int i=0; i<pop_size; i++){
      gen_tmp = genome_->clone();
      gen_tmp->readObject(is);
      pop->add(gen_tmp);
    }
    free (data);
  }
  else {                                              // We received an empty population
    pop = new GAPopulation();
  }

  /*LOG*/ stringstream message; message << "Population received from island " <<  source << " size is " << size <<  endl;
  /*LOG*/ GALogger::instance()->appendPopulation("CommManager:receivingPopulation", message.str(), *pop );

  // BEGIN: Genealogy
  // Add the genomes to the genealogy
  if (GAGenealogy::handle() != NULL)
    for (unsigned i = 0; i < pop->size(); i++)
      GAGenealogy::handle()->addImmigrant(pop->individual(i));
  // END: Genealogy

  return pop;
}

// This is a generic method used to avoid duplication of code. It sends the booleans values passed as argument
// to all islands and return the islands positions that have sent a YES_VALUE
list<int>* CommManager::sendAndReceiveNotifications ( vector<unsigned int>& values ) {

  assert(initialized_ == GAEDAConfig::SYNC);
  assert ((int) values.size() == num_islands_);

  list<int>* islands_have_notified_pos = new list<int>();

  // Notify the islands that have lower rank than us
  for (int i=MASTER_RANK; i<my_rank_; i++){
    MPI_Send ( &(values[i]), 1, MPI_INT, i, TAG_SEND_NOTIFICATION, MPI_COMM_WORLD);
  }

  // Receive notifications from islands that have higher rank than us
  for (int i=num_islands_-1; i>my_rank_; i--){
    int value_received;
    MPI_Recv (&value_received, 1, MPI_INT, i, TAG_SEND_NOTIFICATION, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if ( value_received == YES_NOTIFICATION_VALUE ) islands_have_notified_pos->push_back(i);
  }

  // Notify the islands that have higher rank than us
  for (int i=num_islands_-1; i>my_rank_; i--){
    MPI_Send ( &(values[i]), 1, MPI_INT, i, TAG_SEND_NOTIFICATION, MPI_COMM_WORLD);
  }

  // Receive notifications from islands that have lower rank than us
  for (int i=MASTER_RANK; i<my_rank_; i++){
    int value_received;
    MPI_Recv (&value_received, 1, MPI_INT, i, TAG_SEND_NOTIFICATION, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if ( value_received == YES_NOTIFICATION_VALUE ) islands_have_notified_pos->push_back(i);
  }
  return islands_have_notified_pos;
}

// This methos is used in the ReceiversOriented migrator_. It notifies the islands passed as an argument of this method
// its intention to receive individuals from them and returns the islands that want from the actual island to send them its individuals.
list<int>* CommManager::notifySendersAndGetReceivers ( list<int>& islands_to_notify ) {
  assert(initialized_ == GAEDAConfig::SYNC);
  assert (islands_to_notify.size() > 0);

  {/*LOG*/ stringstream message; message << "We are going to notify the following islands : ";
           for (list<int>::iterator it=islands_to_notify.begin(); it!=islands_to_notify.end(); it++) message << *it << " ," ;
           GALogger::instance()->appendLogMessage("CommManager::notifySelectedIslands",message.str());}

  vector<unsigned int> notifications(num_islands_);
  for (int i=0; i<num_islands_; i++) notifications[i] = NO_NOTIFICATION_VALUE;

  for (list<int>::iterator it=islands_to_notify.begin(); it!=islands_to_notify.end(); it++){
    notifications[*it] = YES_NOTIFICATION_VALUE;
  }

  return sendAndReceiveNotifications ( notifications );
}


GAPopulation* CommManager::sendAndReceiveInds(GAPopulation& emigrant_pop, list<int>& receivers_pos){
  assert(initialized_ == GAEDAConfig::SYNC);

  GAPopulation* inmigrants = new GAPopulation();
  GAPopulation  empty_pop;                                                                                       // Sent to no selected islands

  list<int>     islands_pos;                                                                                     // Just for logging the senders islands

  {/*LOG*/ stringstream message; message << "Before sending to the following islands: ";
           for (list<int>::iterator it=receivers_pos.begin(); it!=receivers_pos.end(); it++) message << *it << " ," ;
           GALogger::instance()->appendLogMessage("CommManager::sendAndReceiveIndividuals",message.str());}

  // Send to the islands that have lower id than us
  for (int i=MASTER_RANK; i<my_rank_; i++){
    if ( find(receivers_pos.begin(), receivers_pos.end(), i) != receivers_pos.end() ) sendPopulation(emigrant_pop,i);
    else                                                                              sendPopulation(empty_pop,i);
  }

  {/*LOG*/ stringstream message; message << "After sending to islands with lower id than my rank: " << my_rank_;
           GALogger::instance()->appendLogMessage("CommManager::sendAndReceiveIndividuals",message.str());}

  // Receive from islands that have higher id than us
  for (int i=num_islands_-1; i>my_rank_; i--){
    auto_ptr<GAPopulation> island_inmigrants ( receivePopulation(i) );                                        // First we receive the inds from each island

    if ( island_inmigrants->size() > 0 ) islands_pos.push_back(i);                                            // We log the island with inds

    for (unsigned i=0; i< island_inmigrants->size(); i++) inmigrants->add( island_inmigrants->individual(i) );     // Then we add them to the global inmigrant
                                                                                                              // The add function clones the genome
  }

  {/*LOG*/ stringstream message; message << "After receiving from islands with higher id than my rank: " << my_rank_;
           GALogger::instance()->appendLogMessage("CommManager::sendAndReceiveIndividuals",message.str());}

  // Send to the islands that have higher id than us
  for (int i=num_islands_-1; i>my_rank_; i--){
    if ( find(receivers_pos.begin(), receivers_pos.end(), i) != receivers_pos.end() ) sendPopulation(emigrant_pop,i);
    else                                                                              sendPopulation(empty_pop,i);
  }

  {/*LOG*/ stringstream message; message << "After sending to islands with higher id than my rank: " << my_rank_;
           GALogger::instance()->appendLogMessage("CommManager::sendAndReceiveIndividuals",message.str());}
  // Receive from islands that have lower id than us
  for (int i=MASTER_RANK; i<my_rank_; i++){
    auto_ptr<GAPopulation> island_inmigrants ( receivePopulation(i) );                                        // First we receive the inds from each island

    if ( island_inmigrants->size() > 0 ) islands_pos.push_back(i);                                            // We log the island with inds

    for (unsigned i=0; i< island_inmigrants->size(); i++) inmigrants->add( island_inmigrants->individual(i) );     // Then we add it to the global inmigrants
                                                                                                              // The add function clones the genome
  }
  {/*LOG*/ stringstream message; message << "After receiving from islands with lower id than my rank: " << my_rank_;
           GALogger::instance()->appendLogMessage("CommManager::sendAndReceiveIndividuals",message.str());}

  {/*LOG*/ stringstream message; message << "After receiving individuals from the following islands: ";
           for (list<int>::iterator it=islands_pos.begin(); it!=islands_pos.end(); it++) message << *it << " ," ;
           GALogger::instance()->appendPopulation("CommManager::sendAndReceiveIndividuals",message.str(), *inmigrants );}

  return inmigrants;
}


// We send our convergence value to all islands and we return the islands that have converged
list<int>*  CommManager::sendAndReceiveConvergence( bool value ){
  assert(initialized_ == GAEDAConfig::SYNC);

  vector<unsigned int> notifications(num_islands_);

  for (int i=0; i<num_islands_; i++)  notifications[i] = (value) ? YES_NOTIFICATION_VALUE : NO_NOTIFICATION_VALUE;

  return sendAndReceiveNotifications ( notifications );
}

// For a future refactoring, the problem is with the receive method since the signature for overloading does not include the return types
// template<typename T> list<T>* sendAndReceive(T value, list<int>& receivers_pos){
//  list<T>* values_received = new list<T>();
//
//  GAPopulation* inmigrants = new GAPopulation();
//  GAPopulation  empty_pop;                                                                                       // Sent to no selected islands
//
//  list<int>     islands_pos;                                                                                     // Just for logging the senders islands
//
//  {/*LOG*/ stringstream message; message << "Before sending to the following islands: ";
//           for (list<int>::iterator it=receivers_pos.begin(); it!=receivers_pos.end(); it++) message << *it << " ," ;
//           GAFileTracer::instance()->appendLogMessage("CommManager::sendAndReceive",message.str());}
//
//
//  MPI_Send ( NO_NOTIFICATION_VALUE, 1, MPI_INT, i, TAG_SEND_NOTIFICATION, MPI_COMM_WORLD);
//
//
//  MPI_Recv (&value_received, 1, MPI_INT, i, TAG_SEND_NOTIFICATION, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//
//
//
//
//
//  // Send to the islands that have lower id than us
//  for (int i=MASTER_RANK; i<my_rank; i++){
//    if ( find(receivers_pos.begin(), receivers_pos.end(), i) != receivers_pos.end() ){
//      MPI_Send ( YES_NOTIFICATION_VALUE, 1, MPI_INT, i, TAG_SEND_NOTIFICATION, MPI_COMM_WORLD);
//      send(emigrant_pop,i);
//    }
//    else{
//      MPI_Send ( NO_NOTIFICATION_VALUE, 1, MPI_INT, i, TAG_SEND_NOTIFICATION, MPI_COMM_WORLD);
//      send(empty_pop,i);
//    }
//  }
//
//  {/*LOG*/ stringstream message; message << "After sending to islands with lower id than my rank: " << my_rank;
//           GAFileTracer::instance()->appendLogMessage("CommManager::sendAndReceiveIndividuals",message.str());}
//
//  // Receive from islands that have higher id than us
//  for (int i=mpi_num_nodes-1; i>my_rank; i--){
//    int data_notification;
//    MPI_Recv (&data_notification, 1, MPI_INT, i, TAG_SEND_NOTIFICATION, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//
//    switch(data_notification){
//    case YES_NOTIFICATION_VALUE:
//      receive()
//      break;
//    case NO_NOTIFICATION_VALUE:
//      break;
//    default:
//      throw runtime_exception("Error in sendAndReceive notification value unrecognised");
//    }
//
//    auto_ptr<GAPopulation> island_inmigrants ( receivePopulation(i) );                                        // First we receive the inds from each island
//
//    if ( island_inmigrants->size() > 0 ) islands_pos.push_back(i);                                            // We log the island with inds
//
//    for (int i=0; i< island_inmigrants->size(); i++) inmigrants->add( island_inmigrants->individual(i) );     // Then we add them to the global inmigrant
//                                                                                                              // The add function clones the genome
//  }
//
//  {/*LOG*/ stringstream message; message << "After receiving from islands with higher id than my rank: " << my_rank;
//           GAFileTracer::instance()->appendLogMessage("CommManager::sendAndReceiveIndividuals",message.str());}
//
//  // Send to the islands that have higher id than us
//  for (int i=mpi_num_nodes-1; i>my_rank; i--){
//    if ( find(receivers_pos.begin(), receivers_pos.end(), i) != receivers_pos.end() ) sendPopulation(emigrant_pop,i);
//    else                                                                              sendPopulation(empty_pop,i);
//  }
//
//  {/*LOG*/ stringstream message; message << "After sending to islands with higher id than my rank: " << my_rank;
//           GAFileTracer::instance()->appendLogMessage("CommManager::sendAndReceiveIndividuals",message.str());}
//  // Receive from islands that have lower id than us
//  for (int i=MASTER_RANK; i<my_rank; i++){
//    auto_ptr<GAPopulation> island_inmigrants ( receivePopulation(i) );                                        // First we receive the inds from each island
//
//    if ( island_inmigrants->size() > 0 ) islands_pos.push_back(i);                                            // We log the island with inds
//
//    for (int i=0; i< island_inmigrants->size(); i++) inmigrants->add( island_inmigrants->individual(i) );     // Then we add it to the global inmigrants
//                                                                                                              // The add function clones the genome
//  }
//  {/*LOG*/ stringstream message; message << "After receiving from islands with lower id than my rank: " << my_rank;
//           GAFileTracer::instance()->appendLogMessage("CommManager::sendAndReceiveIndividuals",message.str());}
//
//  {/*LOG*/ stringstream message; message << "After receiving individuals from the following islands: ";
//           for (list<int>::iterator it=islands_pos.begin(); it!=islands_pos.end(); it++) message << *it << " ," ;
//           GAFileTracer::instance()->appendPopulation("CommManager::sendAndReceiveIndividuals",message.str(), *inmigrants );}
//
//  return inmigrants;
//
//}
//

//
//static vector<GAGenome*>* CommManager::receiveIndsFromMaster(){
//  vector<GAGenome*>*  inds = new vector<GAGenome *>;
//  auto_ptr< GAPopulation > pop_received ( CommManager::receivePopFromMaster() );
//  for (int i=0; i<pop_received->size(); i++) inds->push_back( &( pop_received->individual(i).clone() ) );
//  return inds;
//}

//static vector<GAGenome*>* CommManager::receiveSingleIndsFromSlaveIslands(){
//  vector<GAGenome*>*  inds = new vector<GAGenome *>;
//  auto_ptr< GAPopulation > pop_received ( CommManager::receivePopFromMaster() );
//  for (int i=0; i<pop_received->size(); i++) inds->push_back( &( pop_received->individual(i).clone() ) );
//  return inds;
//
//
//}

//vector<GAGenome*>* CommManager::receiveIndividualsFromMaster(){
//  int rank; MPI_Comm_rank (MPI_COMM_WORLD, &rank); assert( rank != MASTER_RANK );
//
//  auto_ptr<GAPopulation> pop ( receivePopulation(MASTER_RANK) );
//
//  int count; MPI_Comm_size (MPI_COMM_WORLD, &count); assert (pop->size() > 0); assert( pop->size() <= count);  // Some checks
//
//  vector<GAGenome*>* medoids = new vector<GAGenome*>();                                                        // Create the medoids vector and fill it
//
//  int pop_size = pop->size();
//
//  for (int i=0; i<pop_size; i++) {
//    GAGenome& gen_temp = *( pop->individual(i).clone() );
//    medoids->push_back(&gen_temp);
//  }
//
//  assert ( (int) medoids->size() >= 0 && (int) medoids->size() <= count);
//
//  { /*LOG*/ stringstream message; message << "After receiving medoids from master in island: " << rank;
//            GAFileTracer::instance()->appendPopulation( "CommManager::receiveAllMedoidsFromMaster",message.str(), *medoids );}
//
//  return medoids;
//}

//vector<GAGenome*>* CommManager::receiveSingleIndividualFromAllIslands(){
//  int rank; MPI_Comm_rank (MPI_COMM_WORLD, &rank); assert( rank == MASTER_RANK );
//
//  vector<GAGenome* >* medoids = new vector<GAGenome *>();
//
//  for (int i = MASTER_RANK+1; i < mpi_num_nodes; i++)  {
//    {/*LOG*/ stringstream log_messg; log_messg << "receiving medoid from island " << i << endl;
//             GAFileTracer::instance()->appendLogMessage("CommManager:receiveAllMedoidsFromALlIslands", log_messg.str()); }
//
//    medoids->push_back( receiveIndividual(i) );
//  }
//
//  int count; MPI_Comm_size (MPI_COMM_WORLD, &count);  assert ( (int) medoids->size() >= 0 && (int)medoids->size() <= count);
//
//  /*LOG*/ GAFileTracer::instance()->appendPopulation( "CommManager::receiveAllMedoidsFromAllIslands", "the medoids received are:", *medoids );
//
//  return medoids;
//}

//void CommManager::sendFitness(double fitness, int dest_island_rank){
//  MPI_Send (&fitness, 1, MPI_DOUBLE,d est_island_rank, TAG_SEND_DATA, MPI_COMM_WORLD);
//}

//double CommManager::receiveFitness(int source){
//  double result;
//  MPI_Recv(&result, 1, MPI_DOUBLE, source,TAG_SEND_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//  return result;
//}


int CommManager::initialize(int inittype) {
  switch(inittype) {
  case GAEDAConfig::SYNC:
  case GAEDAConfig::ASYNC:
  case GAEDAConfig::ASYNCFAST:
  case GAEDAConfig::CLUST:
  case GAEDAConfig::ROUTINGDISTSYNC:
    initialized_ = GAEDAConfig::SYNC; // We consider all the parallel approaches the same, for MPI Initia√lization purposes
    MPI_Init (&argc_, &argv_);
    MPI_Comm_size (MPI_COMM_WORLD, &num_islands_);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank_);
    break;
  case GAEDAConfig::SEQ:
    initialized_ = GAEDAConfig::SEQ;
    num_islands_ = 1;
    my_rank_ = 0;
    break;
  default:
    throw runtime_error("[CommManager] Error: call to initialize with an invalid communication model. Aborting...");
    break;
  }

  return initialized_;
}

int CommManager::isInitialized() const {
  return initialized_;
}

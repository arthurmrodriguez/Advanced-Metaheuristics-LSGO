#ifndef COMMMANAGER_H_
#define COMMMANAGER_H_

#include "../genomes/GAGenome.h"
#include "../GAPopulation.h"
#include <vector>
#include <list>

using namespace std;

// TODO change not onlu receiveOneIndFromAllSlaveIslads but all the receive methods in order to receive always a vector or
// array of genomes. This way, the methods would be more efficient

//GAGenome* genome_; //TODO Cambiar a atributos

class CommManager {
  GAGenome* genome_;
  int       num_islands_;
  int       my_rank_;

  void          sendIndividual              ( GAGenome& ind, int dest_island_rank           ) const;
  void          sendPopulation              ( GAPopulation& pop, int dest_island_rank       ) const;
  GAGenome*     receiveIndividual           ( int source                                    );
  GAPopulation* receivePopulation           ( int source                                    ) const;
  list<int>*    sendAndReceiveNotifications ( vector<unsigned int>& values                  );

  static CommManager* comm_instance_ptr;//TODO QUITAR
public:
  static CommManager* instance(){return comm_instance_ptr;} //TODO QUITAR

  GAGenome* getGenome (){return genome_;}  //TODO QUITAR

  static const int MASTER_RANK = 0;

  CommManager (int argc, char* argv[]);
  ~CommManager(                      );

  void               setGenome                        (GAGenome* genome                                            );
  long double             getTime                          (                                                            ) const;
  int                getMyRank                        (                                                            ) const;
  int                getNumIslands                    (                                                            ) const;
  bool               isIslandMaster                   (                                                            ) const;
  bool               isIslandSlave                    (                                                            ) const;

  void               sendIndsToSlaveIslands           ( vector<GAGenome*>& medoids                                 ) const;
  void               sendIndToMaster                  ( GAGenome& ind                                              ) const;
  void               sendPopToMaster                  ( GAPopulation& pop                                          ) const;
  void               sendPopsToSlaveIslands           ( const vector<GAPopulation*>& pops                          ) const;

  GAPopulation*      receivePopFromMaster             (                                                            ) const;
  vector<GAGenome*>* receiveOneIndFromAllSlaveIslands (                                                            );
  GAPopulation*      receivePopsFromAllSlaveIslands   (                                                            );
  GAPopulation*      sendAndReceiveInds               ( GAPopulation& emigrant_pop, list<int>& receivers_pos );

  list<int>*         notifySendersAndGetReceivers     ( list<int>& islands_to_notify                               );
  list<int>*         sendAndReceiveConvergence        ( bool value                                                 );
};


#endif /*COMMMANAGER_H_*/

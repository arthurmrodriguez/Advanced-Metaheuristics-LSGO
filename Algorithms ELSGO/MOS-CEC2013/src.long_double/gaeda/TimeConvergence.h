#ifndef TIMECONVERGENCE_H_
#define TIMECONVERGENCE_H_

#include "GAGeneticAlgorithm.h"
#include "islands/CommManager.h"

/*
 * Hack static class. Due to the original design of GAEDAlib, all the methods have been modeled as
 * belonging to a static class
 */

class TimeConvergence {
  static long double       start_time__;
  static long double       max_time__;
  static CommManager* comm_manager__;

public:
  static GABoolean TerminateUponTime (GAGeneticAlgorithm&) {
    assert(comm_manager__ != 0 && start_time__ > 0 && max_time__ > 0);

    return (comm_manager__->getTime() - start_time__ > max_time__) ? gaTrue : gaFalse;
  }

  static void setInitTime   (long double       init_time)    { start_time__    = init_time; assert(init_time >= 0); }
  static void setMaxTime    (long double       max_time)     { max_time__     = max_time;  assert(max_time  >= 0); }
  static void setCommManager(CommManager* comm_manager) { comm_manager__ = comm_manager; }

};


#endif /* TIMECONVERGENCE_H_ */

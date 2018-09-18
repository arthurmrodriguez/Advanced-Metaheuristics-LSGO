#ifndef GALOGGER_H_
#define GALOGGER_H_

#include <iostream>
#include <string>
#include <vector>

#include "LogStat.h"

class GAGenome;
class GAPopulation;

using namespace std;

class GALogger {
 public:
  enum LogLevel {debug, debug_lite, normal, only_stats, none};

protected:
  static GALogger* instance_ptr__;

  GALogger(unsigned int stats_display_freq, LogLevel log_level);

  unsigned int       stats_display_freq_;
  LogLevel           level_;
  GAStatistics*      stats_;
  vector<LogStat*>*  logstats_;                              // This class takes ownership of this vector so it deletes the logstats

public:
  static GALogger* instance();

  virtual ~GALogger();

  void setStats    (GAStatistics* stats);                    // Cannot be set on constructor since the logger is needed before
  void setLogStats (vector<LogStat*>* logstats);

  vector<LogStat*>* getLogStats (){return logstats_;} //TODO borrar

  virtual void appendLogMessage    (string title, string message, LogLevel level=debug)                                     = 0;
  virtual void appendPopulation    (string title, string message, const GAPopulation& pop, LogLevel level=debug)            = 0;
  virtual void appendPopulation    (string title, string message, const vector<GAGenome*>& pop, LogLevel level=debug)       = 0;
  virtual void appendMedoidsStats  (string title, const GAPopulation& pop, vector<GAGenome*>& medoids, GAGenome& my_medoid) = 0;
  virtual void appendStats         (string title)                                                                           = 0;
  virtual void appendStats         (string title, const GAPopulation& pop)                                                  = 0;
  virtual void appendOperatorResult(string title, GAGenome& father, GAGenome& mother, GAGenome& child1, GAGenome& child2)   = 0;
  virtual void appendOperatorResult(string title, GAGenome& father, GAGenome& mother, GAGenome& child1)                     = 0;
  virtual void appendOperatorResult(string title, GAGenome& gen1, GAGenome& gen2)                                           = 0;
  virtual void appendInd           (string title, GAGenome& gen1)                                                           = 0;
  virtual void appendExecTime      (double time)                                                                            = 0;
};

#endif /*GALOGGER_H_*/

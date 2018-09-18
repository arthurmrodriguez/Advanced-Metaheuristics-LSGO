#ifndef GAFILELOGGER_H
#define GAFILELOGGER_H

#include <string>
#include <fstream>
#include <vector>

#include "GALogger.h"

class GAStatistics;

using namespace std;

class GAFileLogger : public GALogger {
  ofstream           os_;
  string             path_;

  void   printStatsLine (                   );
  string getIndLogStr   (const GAGenome& gen);
 public:

  GAFileLogger (char* path, unsigned int stats_display_freq, GALogger::LogLevel log_level);
  ~GAFileLogger();

  void appendLogMessage    (string title, string message, LogLevel level=debug);
  void appendPopulation    (string title, string message, const GAPopulation& pop, LogLevel level=debug);
  void appendPopulation    (string title, string message, const vector<GAGenome*>& pop, LogLevel level=debug);
  void appendMedoidsStats  (string title, const GAPopulation& pop, vector<GAGenome*>& medoids, GAGenome& my_medoid);
  void appendStats         (string title);
  void appendStats         (string title, const GAPopulation& pop);
  void appendOperatorResult(string title, GAGenome& father, GAGenome& mother, GAGenome& child1, GAGenome& child2);
  void appendOperatorResult(string title, GAGenome& father, GAGenome& mother, GAGenome& child1);
  void appendOperatorResult(string title, GAGenome& gen1, GAGenome& gen2);
  void appendInd           (string title, GAGenome& gen1);
  void appendExecTime      (double time);
  void setStatsObj(GAStatistics* stats);
  void updateStats(const GAPopulation& pop);

  string path() const { return path_; }
};

#endif

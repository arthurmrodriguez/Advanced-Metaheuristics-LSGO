#include "GALogger.h"

GALogger* GALogger::instance_ptr__ = NULL;

GALogger::GALogger(unsigned stats_display_freq, LogLevel log_level):
                                                        stats_display_freq_(stats_display_freq),
                                                        level_(log_level),
                                                        logstats_(NULL) {
  instance_ptr__ = this;
}

GALogger::~GALogger() {

  if (logstats_) {

    vector<LogStat*>::iterator it;
    for ( it=logstats_->begin(); it<logstats_->end(); it++ ) {
      delete *it;
    }

    delete logstats_;
  }
}

void GALogger::setStats(GAStatistics* stats){
  stats_ = stats;
}

void GALogger::setLogStats(vector<LogStat*>* logstats){
  logstats_ = logstats;
}

GALogger* GALogger::instance(){
  assert (instance_ptr__ != NULL);
  return instance_ptr__;
}


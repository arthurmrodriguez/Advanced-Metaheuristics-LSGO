#include <iomanip>
#include <sstream>

#include "GAFileLogger.h"

#include "../GAPopulation.h"
#include "../GAStatistics.h"
#include "../genomes/GAGenome.h"

using namespace std;

static const string HEADER_PADDING = "============================ ";
static const string END_PADDING    = "========================================================";


GAFileLogger::GAFileLogger (char*    path,
                            unsigned stats_display_freq,
                            LogLevel log_level): GALogger(stats_display_freq,log_level) {
  os_.open (path, ofstream::out | ofstream::trunc);
  if (level_ < only_stats) os_ << "Starting logging with log Mode set to " << level_ << endl;
}

GAFileLogger::~GAFileLogger () {
    os_.close ();
}

void GAFileLogger::appendLogMessage(string title, string message, LogLevel level){
  if ( level_ <= level ) {
    if (level_ < only_stats ) os_ << HEADER_PADDING << title << HEADER_PADDING << endl << endl;
    os_ << message  << endl;
    if (level_ < only_stats ) os_ << END_PADDING << endl << endl;;
  }
}

// Necesitamos recuperar las estadísticas para optener el numero de generación pero puede que no queramos recuperar todas las estadísticas
void GAFileLogger::setStatsObj(GAStatistics* stats){
  stats_ = stats;
}


void GAFileLogger::appendPopulation(string title, string message,const GAPopulation& pop, LogLevel level){
  if ( level_ <= level ) {
    os_ << endl;
    os_ << HEADER_PADDING << title << HEADER_PADDING << endl << endl;
    os_ << message << endl << endl;
    stringstream  log_messg; log_messg << showbase << showpoint << showpos << fixed << setw(20) << setprecision(10) << right;
    log_messg << "The population is:" << endl << endl;
    for ( unsigned i=0; i<pop.size(); i++) log_messg << getIndLogStr(pop.individual(i,GAPopulation::RAW)) << endl;
    os_ << log_messg.str() << endl;
    for (unsigned  i=0; i<title.size(); i++) os_ << "=";
    os_ << END_PADDING << endl << endl;;
  }
}

void GAFileLogger::appendPopulation(string title, string message, const vector<GAGenome*>& pop, LogLevel level){
  if ( level_ <= level ) {
     os_ << endl;
     os_ << HEADER_PADDING << title << HEADER_PADDING << endl << endl;
     os_ << message << endl << endl;
     stringstream  log_messg; log_messg << showbase << showpoint << showpos << fixed << setw(20) << setprecision(10) << right;
     log_messg << "The population is:" << endl << endl;
     for (unsigned int i=0; i<pop.size(); i++) log_messg << getIndLogStr(*pop[i]) << endl;
     os_ << log_messg.str() << endl;
    for (unsigned int i=0; i<title.size(); i++) os_ << "=";
    os_ << END_PADDING << endl << endl;;
  }
}

void GAFileLogger::appendMedoidsStats(string title, const GAPopulation& pop, vector<GAGenome*>& medoids, GAGenome& my_medoid){
  if (level_ <= GALogger::debug ) {
    long double min_dist_medoid_to_pop = my_medoid.compare( pop.individual(0) );
    long double max_dist_medoid_to_pop = my_medoid.compare( pop.individual(0) );
    long double avg_dist_medoid_to_pop = my_medoid.compare( pop.individual(0) );

    for (unsigned i=1; i<pop.size(); i++){
      if ( my_medoid == pop.individual(i) ) continue;
      long double dist_medoid_to_ind = my_medoid.compare( pop.individual(i) );

      if (dist_medoid_to_ind < min_dist_medoid_to_pop) min_dist_medoid_to_pop = dist_medoid_to_ind;
      if (dist_medoid_to_ind > max_dist_medoid_to_pop) max_dist_medoid_to_pop = dist_medoid_to_ind;
      avg_dist_medoid_to_pop += dist_medoid_to_ind;
    }

    avg_dist_medoid_to_pop = avg_dist_medoid_to_pop / pop.size();

    long double min_dist_between_medoids = ( medoids[0] )->compare( *( medoids[1] ) );
    long double max_dist_between_medoids = ( medoids[0] )->compare( *( medoids[1] ) );
    long double avg_dist_between_medoids = 0.0;

    for (unsigned int i=0; i<medoids.size(); i++){
      for (unsigned int j=i+1; j<medoids.size(); j++){
        long double dist_between_medoids = medoids[i]->compare( *( medoids[j] ) );

        if (dist_between_medoids < min_dist_between_medoids) min_dist_between_medoids = dist_between_medoids;
        if (dist_between_medoids > max_dist_between_medoids) max_dist_between_medoids = dist_between_medoids;
        avg_dist_between_medoids += dist_between_medoids;
      }
    }
    avg_dist_between_medoids = avg_dist_between_medoids / medoids.size();

    appendPopulation( title, "Medoids stats: ", medoids, GALogger::debug );
    appendPopulation( title, "The Medoid  was computed from the following population: ",    pop, GALogger::debug );

    stringstream message;
    message << "Medoids stats ";
    message << " min_dist_medoid_to_pop=" << min_dist_medoid_to_pop;
    message << " max_dist_medoid_to_pop=" << max_dist_medoid_to_pop << " avg_dist_medoid_to_pop=" << avg_dist_medoid_to_pop;
    message << " max_dist_between_medoids=" << max_dist_between_medoids << " min_dist_between_medoids=" << min_dist_between_medoids;
    message << " avg_dist_between_medoids=" << avg_dist_between_medoids;

    appendLogMessage(title, message.str(),GALogger::normal);
  }
}

void GAFileLogger::appendStats(string title, const GAPopulation& pop){
  if (level_ <= GALogger::only_stats ) {
    if (stats_->generation() % stats_display_freq_ == 0) {
      stringstream title_pop; title_pop << "Stats on generation " << stats_->generation();
      appendPopulation( title, title_pop.str() , pop, GALogger::normal);

      stringstream message;
      message << stats_->generation() << " ";

      for (unsigned int i=0; i<logstats_->size(); i++){
        (*logstats_)[i]->update();
        message << (*logstats_)[i]->to_string();
      }

      appendLogMessage(title, message.str(), GALogger::only_stats);
    }
  }
}

void GAFileLogger::appendOperatorResult(string title, GAGenome& father, GAGenome& mother, GAGenome& child1, GAGenome& child2){
  if (level_ <= GALogger::debug_lite ) {
    stringstream inds_str;

    inds_str << "Parents are:" << endl;
    inds_str << getIndLogStr(father) << endl << getIndLogStr(mother) << endl << endl;
    inds_str << "children are:" << endl;
    inds_str << getIndLogStr(child1) << endl << getIndLogStr(child2) << endl;

    appendLogMessage(title,inds_str.str());
  }
}
void GAFileLogger::appendOperatorResult(string title, GAGenome& father, GAGenome& mother, GAGenome& child1){
  if (level_ <= GALogger::debug_lite ) {
    stringstream inds_str;

    inds_str << "Parents are:" << endl;
    inds_str << getIndLogStr(father) << endl << getIndLogStr(mother) << endl << endl;
    inds_str << "Son is:" << endl;
    inds_str << getIndLogStr(child1) << endl;

    appendLogMessage(title,inds_str.str(),debug_lite);
  }
}


void GAFileLogger::appendOperatorResult(string title, GAGenome& gen1, GAGenome& gen2){
  if (level_ <= GALogger::debug_lite ) {
    stringstream inds_str;  inds_str << getIndLogStr(gen1) << endl << getIndLogStr(gen2) << endl;
    appendLogMessage(title,inds_str.str(),debug_lite);
  }
}

void GAFileLogger::appendInd(string title, GAGenome& gen1){
  if (level_ <= GALogger::debug ) {
    appendLogMessage(title,getIndLogStr(gen1),debug_lite);
  }
}


string GAFileLogger::getIndLogStr(const GAGenome& gen){
  stringstream str;
  str << showbase << showpoint << showpos << scientific << setw(20) << setprecision(20) << right ;
  if (level_ <= GALogger::debug )
     str << gen;
  str << " score: " << gen.score() << " age: " << gen.age() << " origin: " << gen.origin();
  str << " id: " << gen.getId() << " island: " << gen.getIsland();
  return str.str();
}

void GAFileLogger::appendExecTime(long double time){
  os_ << showbase << showpoint << fixed << setw(6) << setprecision(2) << right ;
  os_ << "time: "<< time << endl;
}

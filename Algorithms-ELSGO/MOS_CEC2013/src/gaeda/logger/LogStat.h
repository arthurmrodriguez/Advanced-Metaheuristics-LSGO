#ifndef LOGSTAT_H_
#define LOGSTAT_H_

#include <assert.h>
#include <string>

#include "../Algorithm.h"

static const int MESSAGE_SIZE = 550; //TODO: parametrizar el tamanyo

/**
 * This class is used to represent a population stat. By using a class, the GAFileTracer or LogClass is completely
 * decoupled from the number of stats that the GAPopulation could return, making it
 * easier to increase or decrease this number in a future.
 */
class LogStat {

protected:

  std::string name_;
  const Algorithm& alg_;
  char message_[MESSAGE_SIZE];

public:

  LogStat(std::string name, const Algorithm& alg) :
    name_(name), alg_(alg) {
    assert(&(alg) != NULL);
  }

  virtual ~LogStat() {
  }

  std::string getName() {
    return name_;
  }

  virtual const char* to_string() {
    return message_;
  }

  virtual void update() = 0;

};

#endif /*LOGSTAT_H_*/

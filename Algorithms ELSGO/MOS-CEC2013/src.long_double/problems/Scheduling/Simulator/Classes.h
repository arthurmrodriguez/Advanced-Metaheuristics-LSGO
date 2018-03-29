#ifndef _CLASSES_H
#define _CLASSES_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// For file streams
#include <fstream>

// For setting precission in output streams
#include <iomanip>
#include <iostream>
#include <iterator>

#include "EventContainer.h"

struct PenaltyStruct {

   char* class_name;
   float penalty_limit;
   float penalty_per_unit;
   float probability_percent;

};


class Classes {

   public:

      Classes () {}
      Classes (const char* classesfile) {parse (classesfile);}

      virtual ~Classes () {}

      bool parse (const char* classesFileName);

      float calculate_penalty (const char* jobClass, float arrival_time,
                               float completion_time, float job_time) const;

      float get_factor (const char* jobClass) const;


   private:

      std::vector<PenaltyStruct*> m_classes;

};

#endif

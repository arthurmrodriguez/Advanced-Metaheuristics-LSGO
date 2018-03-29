#include "Classes.h"

#define MAX_PENALTY_LIMIT    10
#define MAX_PENALTY_PER_UNIT 10

//-----------------------------------------------------------------------------------------------------------------

bool Classes::parse (const char* classesfile) {

   std::ifstream f (classesfile);

   if (f.bad ()) {

      std::cerr << "ifstream: error opening file: " << classesfile << std::endl;
      return false;

   }
   else {

      char line [1024];

      // Read the header
      f.getline (line, 1024);

      if (line [0] != '#') {

         std::cerr << "Error: cabecera no encontrada en " << classesfile << std::endl;
         return false;

      }

      m_classes.clear();

      while (!f.eof ()) {

         PenaltyStruct* ps = new PenaltyStruct ();

         ps->class_name = new char [128];

         f >> ps->class_name;

         if (f.eof ())
            break;

         f >> ps->penalty_limit;
         f >> ps->penalty_per_unit;
         f >> ps->probability_percent;

         m_classes.push_back (ps);

      }

   }

   return true;

}

//-----------------------------------------------------------------------------------------------------------------

float Classes::calculate_penalty (const char* jobClass, float arrival_time,
                                  float completion_time, float job_time) const {

   long double time_exceeded;
   long double no_penalty_time;

   for (unsigned i = 0; i < m_classes.size (); i++) {

      if (strcmp (m_classes [i]->class_name, jobClass) == 0) {

         no_penalty_time = arrival_time + (job_time * (1.0 + m_classes [i]->penalty_limit));

         time_exceeded = completion_time - no_penalty_time;

         if (time_exceeded > 0.0)
            return (time_exceeded * m_classes [i]->penalty_per_unit);

      }

   }

   return 0.0;

}

//-----------------------------------------------------------------------------------------------------------------

float Classes::get_factor(const char* jobClass) const {

   for (unsigned i = 0; i < m_classes.size (); i++)
      if (strcmp (m_classes [i]->class_name, jobClass) == 0)
         return m_classes [i]->penalty_limit;

   return 0.0;

}

//-----------------------------------------------------------------------------------------------------------------

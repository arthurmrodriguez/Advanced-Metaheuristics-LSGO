#ifndef _EVENTVECTOR_H
#define _EVENTVECTOR_H

#include <map>
#include <iostream>

#include "EventContainer.h"

class EventVector : public EventContainer {

   public:

      EventVector() : m_job_id (0) {}

      virtual ~EventVector() {}

      std::ostream& print (std::ostream& os) const;

      bool     empty () const {return m_vector.empty ();}
      unsigned size  () const {return m_vector.size  ();}

      SimulatorEvent* first (         ) const {return m_vector.begin ()->second;}
      SimulatorEvent* item  (int index) const {return m_vector.find (index)->second;}

      unsigned add (SimulatorEvent* e);
      unsigned add (SimulatorEvent* e, unsigned job_id);

      void remove (        ) {m_vector.erase (m_vector.begin ());}
      void remove (int item) {m_vector.erase (item);}

      void reset ();

      void read_events_from_file (const char* filename, EventType type);
      void print_to_file         (const char* filename                ) const;
      void print_events_to_file  (const char* filename, int nodes     ) const;

      const std::vector<unsigned> get_job_ids ();

      int   get_job_cpus (unsigned index) const {return   m_vector.find (index)->second->m_job->required_cpus ();}
      float get_job_time (unsigned index) const {return   m_vector.find (index)->second->m_job->expected_time ();}
      const Job& get_job (unsigned index) const {return *(m_vector.find (index)->second->m_job);};


   private:

      std::map<int, SimulatorEvent*> m_vector;
      unsigned m_job_id;

};

#endif

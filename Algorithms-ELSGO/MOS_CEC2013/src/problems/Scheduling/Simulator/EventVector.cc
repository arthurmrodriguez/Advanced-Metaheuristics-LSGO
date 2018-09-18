#include <fstream>

#include "EventVector.h"

//-----------------------------------------------------------------------------------------------------------------

std::ostream& EventVector::print (std::ostream& os) const {

   os << size () << " jobs:" << std::endl;

   std::map<int, SimulatorEvent*>::const_iterator it;

   for (it = m_vector.begin (); it != m_vector.end (); it++)
      it->second->print (os);

   return os;

}

//-----------------------------------------------------------------------------------------------------------------

unsigned EventVector::add (SimulatorEvent* e) {

   m_vector [m_job_id] = e;

   m_job_id++;

   return (m_job_id - 1);

}

//-----------------------------------------------------------------------------------------------------------------

unsigned EventVector::add (SimulatorEvent* e, unsigned job_id) {

   m_vector [job_id] = e;
   m_job_id = job_id;

   return m_job_id;

}

//-----------------------------------------------------------------------------------------------------------------

void EventVector::reset () {

   std::map<int, SimulatorEvent*>::iterator iter;

   for (iter = m_vector.begin(); iter != m_vector.end();) {

      delete iter->second;

      iter++;

   }

   m_vector.erase (m_vector.begin (), m_vector.end ());

   m_vector.clear ();

   m_job_id = 0;

}

//-----------------------------------------------------------------------------------------------------------------

void EventVector::read_events_from_file (const char* filename, EventType type) {

   std::ifstream f (filename);

   if (f.bad ())
      std::cerr << "ifstream: error opening file: " << filename << std::endl;
   else {

      char line [1024], label [128], jobClass [128];
      int cpus, job_id;
      float timestamp, memory, expect_time, real_time, start_time, penalty;

      // Leemos la cabecera
      f.getline (line, 1024);

      if (line [0] != '#')
         std::cerr << "Error: cabecera no encontrada en " << filename << std::endl;

      if (m_vector.size () > 0)
         reset ();

      while (!f.eof ()) {

         f >> label;
	 f >> jobClass;

         if (f.eof ())
            break;

         f >> timestamp;
         f >> cpus;
         f >> memory;
         f >> expect_time;
         f >> real_time;
         f >> start_time;
         f >> penalty;
         f >> job_id;

         Job* j = new Job (label, cpus, memory, expect_time, real_time, timestamp, jobClass);
         SimulatorEvent* event = new SimulatorEvent (type, NULL, false, j, start_time, penalty);
         add (event, job_id);

      }

   }

}

//-----------------------------------------------------------------------------------------------------------------

void EventVector::print_to_file (const char* filename) const {

   FILE* f = fopen (filename, "wb");

   if (!f)
      perror ("fopen");
   else {

      char line [1024];
      std::map<int, SimulatorEvent*>::const_iterator it;

      sprintf (line, "# LABEL, TIMESTAMP, CPUS, MEMORY, EXPECT_TIME, REAL_TIME, CLASS\n");
      fputs (line, f);

      for (it = m_vector.begin (); it != m_vector.end (); it++) {

         sprintf (line, "%s, %f, %i, %f, %f, %f, %s\n", it->second->m_job->label          (),
                                                    	  it->second->m_job->timestamp      (),
                                                    	  it->second->m_job->required_cpus  (),
                                                    	  it->second->m_job->memory_per_cpu (),
                                                    	  it->second->m_job->expected_time  (),
                                                    	  it->second->m_job->real_time      (),
                                                        it->second->m_job->jobClass       ());

         fputs (line, f);

      }

      fclose (f);

   }

}

//-----------------------------------------------------------------------------------------------------------------

void EventVector::print_events_to_file (const char* filename, int nodes) const {

   std::ofstream f (filename);

   if (f.bad ())
      std::cerr << "ofstream: error opening file: " << filename << std::endl;
   else {

      std::map<int, SimulatorEvent*>::const_iterator it;

      f << "# LABEL, CLASS, TIMESTAMP, CPUS, MEMORY, EXPECT_TIME, REAL_TIME, START_TIME, PENALTY JOB_ID" << std::endl;

      for (it = m_vector.begin (); it != m_vector.end (); it++) {

         f << it->second->m_job->label          () << " "
           << it->second->m_job->jobClass       () << " "
           << it->second->m_job->timestamp      () << " "
           << it->second->m_job->required_cpus  () << " "
           << it->second->m_job->memory_per_cpu () << " "
           << it->second->m_job->expected_time  () << " "
           << it->second->m_job->real_time      () << " "
           << it->second->m_start_time             << " "
           << it->second->m_penalty                << " "
           << it->first
           << std::endl;

/*
         if (it->second->m_usage)
            for (int j = 0; j < nodes; j++)
               if (j < nodes - 1)
                  f << it->second->m_usage [j] << " ";
               else
                  f << it->second->m_usage [j] << std::endl;
*/

      }

   }

}

//-----------------------------------------------------------------------------------------------------------------

const std::vector<unsigned> EventVector::get_job_ids () {

   std::vector<unsigned> ids;
   std::map<int, SimulatorEvent*>::iterator it;

   for (it = m_vector.begin (); it != m_vector.end (); it++)
      ids.push_back (it->first);

   return ids;
}

//-----------------------------------------------------------------------------------------------------------------

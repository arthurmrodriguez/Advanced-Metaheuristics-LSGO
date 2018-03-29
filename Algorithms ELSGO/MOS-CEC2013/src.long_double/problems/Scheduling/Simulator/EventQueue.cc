#include <fstream>

#include "EventQueue.h"

//-----------------------------------------------------------------------------------------------------------------

std::ostream& EventQueue::print (std::ostream& os) const {

   os << size () << " jobs:" << std::endl;

   SimulHeap::const_iterator it;

   for (it = m_queue.begin (); it != m_queue.end (); it++)
      (*it)->print (os);

   return os;

}

//-----------------------------------------------------------------------------------------------------------------

unsigned EventQueue::add (SimulatorEvent* e)  {

   m_queue.push (e);

   return (m_queue.size () - 1);

}

//-----------------------------------------------------------------------------------------------------------------

unsigned EventQueue::add (SimulatorEvent* e, unsigned job_id) {

   m_queue.push (e);

   return (m_queue.size () - 1);

}

//-----------------------------------------------------------------------------------------------------------------

void EventQueue::reset () {

   SimulHeap::const_iterator it;

   for (it = m_queue.begin (); it != m_queue.end (); it++)
      delete *it;

   m_queue.clear ();

}

//-----------------------------------------------------------------------------------------------------------------

void EventQueue::read_events_from_file (const char* filename, EventType type, int nodes) {

   std::ifstream f (filename);

   if (f.bad ())
      std::cerr << "ifstream: error opening file: " << filename << std::endl;
   else {

      char line [1024], label [128], jobClass [128];
      int cpus, job_id;
      float timestamp, memory, expect_time, real_time, start_time, penalty;

      // Leemos las cabeceras
      f.getline (line, 1024);
      f.getline (line, 1024);

      if (line [0] != '#')
         std::cerr << "Error: cabecera no encontrada en " << filename << std::endl;

      if (m_queue.size () > 0)
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

         int* usage = (int*) malloc (sizeof (int) * nodes);

         for (int i = 0; i < nodes; i++)
            f >> usage [i];

         Job* j = new Job (label, cpus, memory, expect_time, real_time, timestamp, jobClass);
         SimulatorEvent* event = new SimulatorEvent (type, usage, true, j, start_time, penalty);
         add (event, job_id);

      }

   }

}

//-----------------------------------------------------------------------------------------------------------------

void EventQueue::print_to_file (const char* filename) const {

   FILE* f = fopen (filename, "wb");

   if (!f)
      perror ("fopen");
   else {

      char line [1024];

      sprintf (line, "# LABEL, TIMESTAMP, CPUS, MEMORY, EXPECT_TIME, REAL_TIME, CLASS\n");
      fputs (line, f);

      SimulHeap::const_iterator it;

      for (it = m_queue.begin (); it != m_queue.end (); it++) {

         SimulatorEvent* e = *it;

         sprintf (line, "%s, %f, %i, %f, %f, %f, %s\n", e->m_job->label          (),
                                                    	  e->m_job->timestamp      (),
                                                        e->m_job->required_cpus  (),
                                                    	  e->m_job->memory_per_cpu (),
                                                    	  e->m_job->expected_time  (),
                                                    	  e->m_job->real_time      (),
                                                        e->m_job->jobClass       ());

         fputs (line, f);

      }

      fclose (f);

   }

}

//-----------------------------------------------------------------------------------------------------------------

void EventQueue::print_events_to_file (const char* filename, int nodes) const {

   std::ofstream f (filename);

   if (f.bad ())
      std::cerr << "ofstream: error opening file: " << filename << std::endl;
   else {

      f << "# LABEL, CLASS, TIMESTAMP, CPUS, MEMORY, EXPECT_TIME, REAL_TIME, START_TIME, PENALTY, JOB_ID" << std::endl;
      f << "# USAGE" << std::endl;

      SimulHeap::const_iterator it;
      int i = 0;

      for (it = m_queue.begin (); it != m_queue.end (); it++) {

         SimulatorEvent* e = *it;

         f << e->m_job->label          () << " "
           << e->m_job->jobClass       () << " "
           << e->m_job->timestamp      () << " "
           << e->m_job->required_cpus  () << " "
           << e->m_job->memory_per_cpu () << " "
           << e->m_job->expected_time  () << " "
           << e->m_job->real_time      () << " "
           << e->m_start_time             << " "
           << e->m_penalty                << " "
           << i                                  // There's no need for this actually...
           << std::endl;


         if (e->m_usage)
            for (int j = 0; j < nodes; j++)
               if (j < nodes - 1)
                  f << e->m_usage [j] << " ";
               else
                  f << e->m_usage [j] << std::endl;

         i++;

      }

   }

}

//-----------------------------------------------------------------------------------------------------------------

#include "EventContainer.h"

#include <stdio.h>

//-----------------------------------------------------------------------------------------------------------------

std::ostream& SimulatorEvent::print (std::ostream& os) const {

   os << "   Job " << m_job->label ()
      << ": Requirements ("
      << m_job->required_cpus () << "CPUs, "
      << m_job->memory_per_cpu () << "MB/CPU) "
      << " Execution time ("
      << "Expected " << m_job->expected_time () << ", "
      << "Real " << m_job->real_time () << ", "
      << "Timestamp " << m_job->timestamp () << ")" << std::endl;

   return os;

}

//-----------------------------------------------------------------------------------------------------------------

void EventContainer::read_from_file (const char* filename, EventType type) {

   FILE* f = fopen (filename, "r");

   if (!f)
      perror ("fopen");
   else {

      char line [1024], label [128], jobClass[128];
      int cpus;
      float timestamp, memory, expect_time, real_time;

      if (m_job_order.size () > 0)
         reset ();

      while (fgets (line, sizeof (line), f)) {

         if (line [0] == '#')
            continue;

         if (sscanf (line, "%[^ ,],%f,%i,%f,%f,%f,%s", label, &timestamp, &cpus, &memory, &expect_time, &real_time, jobClass) == 7) {

            Job* j = new Job (label, cpus, memory, expect_time, real_time, timestamp, jobClass);
            SimulatorEvent* event = new SimulatorEvent (type, NULL, false, j, 0, 0);
            add (event);

         }

      }

      fclose (f);

   }

}

//-----------------------------------------------------------------------------------------------------------------

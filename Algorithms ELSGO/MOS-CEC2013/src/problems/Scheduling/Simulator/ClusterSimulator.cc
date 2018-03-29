#include "ClusterSimulator.h"

//-----------------------------------------------------------------------------------------------------------------

int assign_matrix [1000][10000];

//-----------------------------------------------------------------------------------------------------------------

ClusterSimulator::ClusterSimulator () : m_machine       ( ),
                                        m_start_time    (0),
                                        m_current_time  (0),
                                        m_total_penalty (0),
                                        m_new_time      (0)  {

   m_job_list = new EventVector ();
   m_queue    = new EventQueue  ();
   m_classes  = new Classes     ();

}

//-----------------------------------------------------------------------------------------------------------------

ClusterSimulator::ClusterSimulator (const char* machinefile, const char* classesfile)
                     : m_machine       (machinefile),
                       m_start_time    (0),
                       m_current_time  (0),
                       m_total_penalty (0),
                       m_new_time      (0)            {

   m_job_list = new EventVector ();
   m_queue    = new EventQueue  ();
   m_classes  = new Classes     (classesfile);

}

//-----------------------------------------------------------------------------------------------------------------

void ClusterSimulator::add_event (SimulatorEvent* event) {

   event->load (waiting_event, NULL, false, event->m_job, event->m_job->timestamp (), event->penalty ());

   unsigned job_id = m_job_list->add (event);
   m_order.push_back (job_id);

   m_current_time = event->m_job->timestamp ();

}

//-----------------------------------------------------------------------------------------------------------------

float ClusterSimulator::run_to_next_event () {

   SimulatorEvent* e;

   if (m_queue->empty ())
      return m_current_time;

   e = m_queue->first ();
   m_queue->remove    ();

   int   cpus   = 0;
   float memory = 0.0;

   for (int i = 0; i < m_machine.total_nodes (); i++)
      if (e->m_usage [i]) {

         m_machine.node (i).release_resources (e->m_usage [i],
                                               e->m_usage [i] * e->memory_per_cpu ());
         cpus   += e->m_usage [i];
         memory += e->m_usage [i] * e->memory_per_cpu ();

      }

   m_machine.release_machine_resources (cpus, memory);
   m_current_time = e->end_time ();

   float penalty = m_classes->calculate_penalty (e->m_job->jobClass (), e->m_job->timestamp (),
                                                 m_current_time, e->m_job->real_time ());

   m_total_penalty += penalty;

   m_new_time = m_current_time + m_total_penalty;

   // BEGIN: DEBUG
/*
   std::cerr << "     => Current_time: " << m_current_time << std::endl;
   std::cerr << "     => Accumulated penalty: " << m_total_penalty << std::endl;
   std::cerr << "     => Current time + Accumulated penalty: " << m_current_time + m_total_penalty << std::endl;
   std::cerr << "     => Finished job: " << std::endl;
   std::cerr << "     ";
   e->print (std::cerr);
   std::cerr << "     => Resources freed: " << cpus << " CPUS and " << memory << " MB of RAM." << std::endl;
   std::cerr << "     => Resources (CPUs): " << m_machine.total_cpus () << " CPUs:" << m_machine.total_cpus () - m_machine.used_cpus () << " free and " << m_machine.used_cpus () << " used." << std::endl;
   std::cerr << "     => Resources (RAM): " << m_machine.total_memory () << " MB of RAM: " << m_machine.total_memory () - m_machine.used_memory () << " free and " << m_machine.used_memory () << " MB used." << std::endl;
   std::cerr << "     => Jobs executing: " << std::endl;

   for (unsigned i = 0; i < m_queue->size (); i++) {
      std::cerr << "     ";
      m_queue->item(i)->print (std::cerr);
   }

   std::cerr << "     => Jobs waiting (ordered list): " << std::endl;

   std::vector<unsigned> order = m_job_list->get_job_ids ();

   for (unsigned i = 0; i < order.size (); i++) {
      std::cerr << "     ";
      m_job_list->item(order [i])->print (std::cerr);
   }

   std::cerr << std::endl;

   std::cerr << "     => Machine state: " << std::endl;
   m_machine.print (std::cerr);
*/
   // END: DEBUG

   // m_job is a shared object. We delete it here and ONLY here.
   delete e->m_job;
   delete e;

   return m_current_time;

}

//-----------------------------------------------------------------------------------------------------------------

EventVector* ClusterSimulator::step () {

   EventVector* event_vector = new EventVector ();

   while (!m_order.empty ()) {

      int* assign = (int*) malloc (sizeof (int) * m_machine.total_nodes ());
      SimulatorEvent* e = m_job_list->item (m_order [0]);

      // BEGIN: DEBUG
      /*
      std::cerr << "     => Checking Job to fit in the machine: " << std::endl;
      std::cerr << "     ";
      e->print (std::cerr);
      */
      // END: DEBUG

      if (!m_machine.first_available (e->m_job->required_cpus (), e->m_job->memory_per_cpu (), assign)) {

         e->load (executing_event, assign, true, e->m_job , m_current_time, 0.0);

         // BEGIN: DEBUG
         /*
         std::cerr << "     => Success! " << std::endl;
         std::cerr << std::endl;
         */
         // END: DEBUG

         m_queue->add (e);

         m_job_list->remove (m_order [0]);
         m_order.erase (m_order.begin ());

         // Memory of assign vector will be freed by run_to_next_event
         SimulatorEvent* out = new SimulatorEvent (out_event, assign, false, e->m_job, m_current_time, 0.0);
         event_vector->add (out);

      }
      else {

         // BEGIN: DEBUG
         /*
         std::cerr << "     => Fail! " << std::endl;
         std::cerr << std::endl;
         */
         // END: DEBUG

         break;

      }

   }

   return event_vector;

}

//-----------------------------------------------------------------------------------------------------------------

float ClusterSimulator::simulate (const std::vector<unsigned>& order) {

   if (!order.empty ()) {

      unsigned i = 0;

      while (i < order.size ()) {

         //int* assign = (int*) malloc (sizeof (int) * m_machine.total_nodes ());
         int* assign = assign_matrix [i];

         SimulatorEvent* e = m_job_list->item (order [i]);

         if (!m_machine.first_available (e->m_job->required_cpus (), e->m_job->memory_per_cpu (), assign)) {

            m_job_list->remove (order [i]);

            e->load (executing_event, assign, false, e->m_job, m_current_time, 0.0);
            m_queue->add (e);

            i++;

         }
         else
            run_to_next_event ();

      }

      while (keep_running ())
         run_to_next_event ();

   }

   return m_current_time;

}

//-----------------------------------------------------------------------------------------------------------------

float ClusterSimulator::simulate_backfilling (std::vector<unsigned>& order) {

   int i = 0;
   std::vector<unsigned> order2 = order;
   std::vector<unsigned>::iterator it;

   if (!order2.empty ()) {

      while (order2.size () > 0) {

         SimulatorEvent* e = NULL;
         //int* assign = (int*) malloc (sizeof (int) * m_machine.total_nodes ());
         int* assign = assign_matrix [i];

         it = order2.begin ();

         while (it != order2.end ()) {

            e = m_job_list->item (*it);

            if (!m_machine.first_available (e->m_job->required_cpus (), e->m_job->memory_per_cpu (), assign))
               break;

            it++;

         }

         if (it != order2.end ()) {

            order [i] = *it;
            i++;

            m_job_list->remove (*it);
            order2.erase (it);

            e->load (executing_event, assign, false, e->m_job, m_current_time, 0.0);
            m_queue->add (e);

         }
         else
            run_to_next_event ();

      }

      while (keep_running ())
         run_to_next_event ();

   }

   return m_current_time;

}

//-----------------------------------------------------------------------------------------------------------------

float ClusterSimulator::simulate_backfilling_res (std::vector<unsigned>& order) {

   int i = 0;
   std::vector<unsigned> order2 = order;
   std::vector<unsigned>::iterator it;

   if (!order2.empty ()) {

      while (order2.size () > 0) {

         SimulatorEvent* e     = NULL;
         SimulatorEvent* first = NULL;
         //int* assign = (int*) malloc (sizeof (int) * m_machine.total_nodes ());
         int* assign = assign_matrix [i];

         it = order2.begin ();
         e = first = m_job_list->item (*it);

         if (m_machine.first_available (first->m_job->required_cpus (), first->m_job->memory_per_cpu (), assign)) {

            it++;

            float res_time = getReservationTime (first->m_job->required_cpus (), first->m_job->memory_per_cpu ());

            while (it != order2.end ()) {

               e = m_job_list->item (*it);

               if (res_time > e->end_time ())
                  if (!m_machine.first_available (e->m_job->required_cpus (), e->m_job->memory_per_cpu (), assign))
                     break;

               it++;

            }

         }

         if (it != order2.end ()) {

            order [i] = *it;
            i++;

            m_job_list->remove (*it);
            order2.erase (it);

            e->load (executing_event, assign, false, e->m_job, m_current_time, 0.0);
            m_queue->add (e);

         }
         else
            run_to_next_event ();

      }

      while (keep_running ())
         run_to_next_event ();

   }

   return m_current_time;

}

//-----------------------------------------------------------------------------------------------------------------

int ClusterSimulator::get_job_cpus (unsigned index) const {

   return m_job_list->get_job_cpus (index);

}

//-----------------------------------------------------------------------------------------------------------------

float ClusterSimulator::get_job_time (unsigned index) const {

   return m_job_list->get_job_time (index);

}

//-----------------------------------------------------------------------------------------------------------------

float ClusterSimulator::get_job_effective_deadline (unsigned index) const {

   const Job& job = m_job_list->get_job (index);

   return job.timestamp () + (job.real_time () * m_classes->get_factor (job.jobClass ()));

}

//-----------------------------------------------------------------------------------------------------------------

float ClusterSimulator::get_job_deadline (unsigned index) const {

   const Job& job = m_job_list->get_job (index);

   return job.timestamp () + (job.real_time () * (1.0 + m_classes->get_factor (job.jobClass ())));

}

//-----------------------------------------------------------------------------------------------------------------

double ClusterSimulator::total_processor_time () const {

   double time = 0.0;

   std::vector<unsigned> job_ids = m_job_list->get_job_ids ();

   // We first compute the accumulated time of running jobs

   EventQueue::const_iterator it;

   for (it = m_queue->begin (); it != m_queue->end (); it++) {

      SimulatorEvent* e = *it;
      time += e->m_job->required_cpus () * (e->end_time () - m_current_time);

   }

   // And now the time of waiting events

   for (unsigned i = 0; i < job_ids.size (); i++) {

      Job* job = m_job_list->item (job_ids [i])->m_job;
      time += job->required_cpus () * job->real_time ();

   }

   return time;

}

//-----------------------------------------------------------------------------------------------------------------

void ClusterSimulator::read_simulator_state  (const char* filesimulator) {

   FILE* f;

   f = fopen (filesimulator, "r");

   if (!f)
      perror ("fopen");
   else {

      char line [1024];

      while (fgets (line, sizeof (line), f)) {

         if (line [0] == '#')
            continue;

         if (sscanf (line, "%f,%f,%f", &m_start_time, &m_current_time, &m_total_penalty) == 3)
            break;

      }

      fclose (f);

   }

}

//-----------------------------------------------------------------------------------------------------------------

void ClusterSimulator::print_simulator_state (const char* filesimulator) const {

   FILE* f = fopen (filesimulator, "wb");

   if (!f)
      perror ("fopen");
   else {

      char line [1024];

      sprintf (line, "#START_TIME, CURRENT_TIME, PENALTY_TIME\n");
      fputs (line, f);

      sprintf (line, "%f,%f,%f\n", m_start_time, m_current_time, m_total_penalty);
      fputs (line, f);

      fclose (f);

   }

}

//-----------------------------------------------------------------------------------------------------------------

float ClusterSimulator::getReservationTime (int cpus, float memory) const {

   int total_cpus = 0;
   int size = m_queue->size ();
   SimulatorEvent* e = NULL;

   for (int i = 0; total_cpus < cpus && i < size; i++) {

      e = m_queue->item (i);

      for (int j = 0; j < m_machine.total_nodes (); j++)
         if (e->m_usage [j]) {

            float free_memory = m_machine.node (j).total_memory () -
                                m_machine.node (j).used_memory  () +
                                e->memory_per_cpu ();

            if (free_memory >= memory)
               total_cpus   += e->m_usage [j];

         }

   }

   if (total_cpus >= cpus)
      return e->end_time ();
   else
      return -1.0;

}

//-----------------------------------------------------------------------------------------------------------------

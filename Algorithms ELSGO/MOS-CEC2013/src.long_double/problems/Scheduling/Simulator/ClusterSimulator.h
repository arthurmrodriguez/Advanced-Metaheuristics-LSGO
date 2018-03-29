#ifndef _CLUSTERSIMULATOR_H
#define _CLUSTERSIMULATOR_H

#include "ClusterMachine.h"
#include "Job.h"
#include "EventQueue.h"
#include "EventVector.h"
#include "Classes.h"

class ClusterSimulator {

   public:

      ClusterSimulator  ();
      ClusterSimulator  (const char* machinefile, const char* classesfile);

      virtual ~ClusterSimulator () {

         delete m_job_list;
         delete m_queue;
         delete m_classes;

      }

      const ClusterMachine& cluster_machine () {return  m_machine;};
      const EventVector&    job_list        () {return *m_job_list;};

      void add_event (SimulatorEvent* event);

      bool  keep_running () const {return !m_job_list->empty () || !m_queue->empty ();};
      float current_time () const {return m_current_time;};

      float run_to_next_event ();

      void resetSimulation () {

         m_current_time  = m_start_time;
         m_total_penalty = 0.0;
         m_job_list->reset ();
         m_queue->reset    ();

      }

      EventVector* step ();

      float simulate                 (const std::vector<unsigned>& order);
      float simulate_backfilling     (      std::vector<unsigned>& order);
      float simulate_backfilling_res (      std::vector<unsigned>& order);

      void set_order (const std::vector<unsigned>& order) {m_order = order;}

      std::ostream& print (std::ostream& os) const {

         m_job_list->print (os);
         m_queue->print    (os);

         return os;

      }

      void read_from_files (const char* filejobs, const char* filemachine) {

         if (filemachine)
            m_machine.read_from_file (filemachine);

         if (filejobs)
            m_job_list->read_from_file (filejobs,  waiting_event);

      }

      void read_from_files (const char* filejobs,    const char* filequeue,
                            const char* filemachine, const char* filesimulator) {

         if (filemachine)
            m_machine.read_from_file (filemachine);

         if (filejobs)
            m_job_list->read_events_from_file (filejobs,  waiting_event);

         if (filequeue)
            m_queue->read_events_from_file (filequeue, executing_event, m_machine.total_nodes ());

         if (filesimulator)
            read_simulator_state (filesimulator);

      }

      void print_to_files (const char* filejobs,    const char* filequeue,
                           const char* filemachine, const char* filesimulator) const {

         if (filejobs)
            m_job_list->print_events_to_file (filejobs,  m_machine.total_nodes ());

         if (filequeue)
            m_queue->print_events_to_file (filequeue, m_machine.total_nodes ());

         if (filemachine)
            m_machine.print_to_file (filemachine);

         if (filesimulator)
            print_simulator_state (filesimulator);

      }

      int num_events_waiting () const {return m_job_list->size ();}
      int num_events_running () const {return m_queue->size    ();}

      const std::vector<unsigned> get_job_ids () const {return m_job_list->get_job_ids ();}

      int   get_job_cpus               (unsigned index) const;
      float get_job_time               (unsigned index) const;
      float get_job_deadline           (unsigned index) const;
      float get_job_effective_deadline (unsigned index) const;

      long double total_processor_time () const;

      bool parse_classes (const char* classesfile) {return m_classes->parse (classesfile);}

      float get_total_penalty () const  {return m_total_penalty;}

   private:

      void read_simulator_state  (const char* filesimulator);
      void print_simulator_state (const char* filesimulator) const;

      float getReservationTime (int cpus, float memory) const;

      // Attributes
      EventVector*          m_job_list;
      EventQueue*           m_queue;
      ClusterMachine        m_machine;
      float                 m_start_time;
      float                 m_current_time;
      float                 m_total_penalty;
      std::vector<unsigned> m_order;
      Classes*              m_classes;

   public:

      float	      m_new_time;

};

#endif

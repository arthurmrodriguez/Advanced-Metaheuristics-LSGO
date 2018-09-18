#ifndef _EVENTQUEUEBASE_H
#define _EVENTQUEUEBASE_H

#include <iostream>
#include <stdlib.h>

#include "Job.h"

#define END_TIME_WEIGHT 1.00

enum EventType {in_event = 0, out_event = 1, executing_event = 2, waiting_event = 3};


struct SimulatorEvent {

   SimulatorEvent () {}
   SimulatorEvent (EventType type, int usage [], bool utype, Job* job, float start_time, float penalty)
                     : m_usage      (usage     ),
                       m_utype      (utype     ),
                       m_job        (job       ),
                       m_start_time (start_time),
							  m_penalty    (penalty   ),
                       m_steps      (1         ) {

      m_type = type;

   };

   virtual ~SimulatorEvent () {

      //delete m_job;

      if (m_utype)
         free (m_usage);

   }

   void load (EventType type, int* usage, bool utype, Job* job, float start_time, float penalty) {

      m_type       = type;
      m_utype      = utype;
      m_usage      = usage;
      m_job        = job;
      m_start_time = start_time;
      m_penalty    = penalty;

   }

   EventType type           () const {return m_type;}
   int       steps          () const {return m_steps;}
   float     memory_per_cpu () const {return m_job->memory_per_cpu ();}
   float     end_time       () const {return m_start_time + m_job->real_time ();}
   float     start_time     () const {return m_start_time;}
   float     penalty        () const {return m_penalty;}
   void      add_step       ()       {m_steps++;}

   std::ostream& print (std::ostream& os) const;

   EventType m_type;
   int*      m_usage;
   bool      m_utype;
   Job*      m_job;
   float     m_start_time;
   float     m_penalty;
   int       m_steps;

};


class EventContainer {

   public:

      virtual ~EventContainer () {}

      const std::vector<int>& job_order () const {return m_job_order;}

      virtual std::ostream& print (std::ostream& os) const = 0;

      void read_from_file (const char* filename, EventType type);

      virtual void print_to_file         (const char* filename           ) const = 0;
      virtual void print_events_to_file  (const char* filename, int nodes) const = 0;

      virtual bool      empty () const = 0;
      virtual unsigned  size  () const = 0;

      virtual SimulatorEvent* first (         ) const = 0;
      virtual SimulatorEvent* item  (int index) const = 0;

      virtual unsigned add (SimulatorEvent* e)                  = 0;
      virtual unsigned add (SimulatorEvent* e, unsigned job_id) = 0;

      virtual void remove () = 0;
      virtual void reset  () = 0;


   private:

      std::vector<int> m_job_order;

};


class EventCompare : public std::binary_function<SimulatorEvent*, SimulatorEvent*, bool> {

   public:

      bool operator () (SimulatorEvent*& e1, SimulatorEvent*& e2) {

         if (e1->type () == executing_event && e2->type () == executing_event) {

            float he1 = e1->end_time () * END_TIME_WEIGHT - e1->steps () * (1.00 - END_TIME_WEIGHT);
            float he2 = e2->end_time () * END_TIME_WEIGHT - e2->steps () * (1.00 - END_TIME_WEIGHT);

            return he1 > he2;

         }
         else if (e1->type () == in_event && e2->type () == in_event)
            return e1->m_job->timestamp () > e2->m_job->timestamp ();
         else if (e1->type () == out_event && e2->type () == out_event)
            return e1->end_time () > e2->end_time ();
         else if(e1->type () == in_event && e2->type () == out_event)
            return e1->m_job->timestamp () > e2->end_time ();
         else if(e1->type () == out_event && e2->type () == in_event)
            return e1->end_time () > e2->m_job->timestamp ();
         else
            return e1->end_time () > e2->end_time ();

      }

};

#endif

#ifndef _EVENTQUEUE_H
#define _EVENTQUEUE_H

#include <iostream>

#include "Heap.h"
#include "EventContainer.h"

class EventQueue : public EventContainer {

   protected:

      typedef Heap<SimulatorEvent*, EventCompare> SimulHeap;


   public:

      typedef SimulHeap::const_iterator const_iterator;

      EventQueue () {}
      virtual ~EventQueue () {}

      std::ostream& print (std::ostream& os) const;

      const_iterator begin () const {return m_queue.begin ();}
      const_iterator end   () const {return m_queue.end   ();}

      bool     empty () const {return m_queue.empty ();}
      unsigned size  () const {return m_queue.size  ();}

      SimulatorEvent* first (         ) const {return m_queue.top ();}
      SimulatorEvent* item  (int index) const {return *(m_queue.begin () + index);}

      unsigned add (SimulatorEvent* e);
      unsigned add (SimulatorEvent* e, unsigned job_id);

      void remove () {m_queue.pop ();}

      void reset ();

      void read_events_from_file (const char* filename, EventType type, int nodes);
      void print_to_file         (const char* filename                           ) const;
      void print_events_to_file  (const char* filename, int nodes                ) const;


   private:

      SimulHeap m_queue;

};

#endif

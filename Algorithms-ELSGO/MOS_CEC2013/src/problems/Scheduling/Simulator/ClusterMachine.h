#ifndef _CLUSTERMACHINE_H
#define _CLUSTERMACHINE_H

#include <iostream>
#include <vector>

#define RATIO_LOWER_THRESHOLD 0.75
#define RATIO_UPPER_THRESHOLD 1.25
#define MEM_ZERO 0.01


class ClusterNode {

   public:

      ClusterNode (int cpus, float memory) : m_total_cpus   (cpus  ),
                                             m_used_cpus    (0     ),
                                             m_total_memory (memory),
                                             m_used_memory  (0.0   ) {}

      int check_resources (int cpus, float memory) const {

         return ((m_used_cpus   + cpus   <= m_total_cpus  ) &&
                 (m_used_memory + memory <= m_total_memory)    );

      }

      void take_resources (int cpus, float memory) {

         m_used_cpus   += cpus;
         m_used_memory += memory;

      }

      void release_resources (int cpus, float memory) {

         m_used_cpus   -= cpus;
         m_used_memory -= memory;

         if (m_used_memory < MEM_ZERO)
            m_used_memory = 0.0;

      }

      float tentative_ratio (int cpus, float memory) const {

         return ((m_used_memory + memory) / m_total_memory) /
                 ((float) (m_used_cpus + cpus) / (float) m_total_cpus);

      }

      int total_cpus     () const {return m_total_cpus;}
      int used_cpus      () const {return m_used_cpus;}

      float total_memory () const {return m_total_memory;}
      float used_memory  () const {return m_used_memory;}

      float cpus_usage   () const {return (float) m_used_cpus / (float) m_total_cpus;}
      float memory_usage () const {return m_used_memory / m_total_memory;}


   private:

      int   m_total_cpus;
      int   m_used_cpus;
      float m_total_memory;
      float m_used_memory;

};

//-----------------------------------------------------------------------------------------------------------------

class ClusterMachine {

   public:

      ClusterMachine  () : m_total_nodes  (0),
                           m_total_cpus   (0),
                           m_used_cpus    (0),
                           m_total_memory (0),
                           m_used_memory  (0) {;}

      ClusterMachine  (const char* machinefilename);

      virtual ~ClusterMachine ();

      ClusterNode& node (int i) const {return *(m_nodes [i]);};

      bool first_available (int cpus, float memory_per_cpu, int assign []);
      bool best_available  (int cpus, float memory_per_cpu, int assign []);
      bool worst_available (int cpus, float memory_per_cpu, int assign []);

      int check_machine_resources (int cpus, float memory) const {

         return ((m_used_cpus   + cpus   <= m_total_cpus  ) &&
                 (m_used_memory + memory <= m_total_memory)    );

      }

      void take_machine_resources (int cpus, float memory) {

         m_used_cpus   += cpus;
         m_used_memory += memory;

      }

      void release_machine_resources (int cpus, float memory) {

         m_used_cpus   -= cpus;
         m_used_memory -= memory;

         if (m_used_memory < MEM_ZERO)
            m_used_memory = 0.0;

      }

      int total_nodes    () const {return m_total_nodes;}

      int total_cpus     () const {return m_total_cpus;}
      int used_cpus      () const {return m_used_cpus;}

      float total_memory () const {return m_total_memory;}
      float used_memory  () const {return m_used_memory;}

      float cpu_usage    () const {return (float) m_used_cpus / (float) m_total_cpus;}
      float memory_usage () const {return m_used_memory / m_total_memory;}

      std::ostream& print (std::ostream& os) const;

      void read_from_file (const char* filename);
      void print_to_file  (const char* filename) const;

      void destroy ();

      typedef std::vector<ClusterNode*>::iterator iterator;


   private:

      int                       m_total_nodes;
      int                       m_total_cpus;
      int                       m_used_cpus;
      float                     m_total_memory;
      float                     m_used_memory;
      std::vector<ClusterNode*> m_nodes;

};

#endif

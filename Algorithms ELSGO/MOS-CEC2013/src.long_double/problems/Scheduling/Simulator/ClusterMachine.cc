#include "ClusterMachine.h"

#include <string.h>
#include <stdio.h>

//-----------------------------------------------------------------------------------------------------------------

ClusterMachine::ClusterMachine (const char* filename) {

   read_from_file (filename);

}

//-----------------------------------------------------------------------------------------------------------------

ClusterMachine::~ClusterMachine () {

   this->destroy ();

}

//-----------------------------------------------------------------------------------------------------------------

bool ClusterMachine::first_available (int cpus, float memory_per_cpu, int use []) {

   int n = 0;
   int myCPUs = cpus;

   memset ((void*) use, 0, m_total_nodes * sizeof (int));

   if (check_machine_resources (cpus, cpus * memory_per_cpu)) {

      take_machine_resources (cpus, cpus * memory_per_cpu);

      for (ClusterMachine::iterator i = m_nodes.begin (); myCPUs && i < m_nodes.end (); n++, i++)
         while (myCPUs && (*i)->check_resources (1, memory_per_cpu)) {

            (*i)->take_resources (1, memory_per_cpu);
            use [n]++;
            myCPUs--;

         }

   }
   else
     return true;

   if (myCPUs) {

      for (n--; n >= 0; n--)
         if (use [n])
            m_nodes [n]->release_resources (use [n], use [n] * memory_per_cpu);

      release_machine_resources (cpus, cpus * memory_per_cpu);
      return true;

   }

   return false;

}

//-----------------------------------------------------------------------------------------------------------------

bool ClusterMachine::best_available (int cpus, float memory_per_cpu, int use []) {

   int n = 0;
   int myCPUs = cpus;

   memset ((void*) use, 0, m_total_nodes * sizeof (int));

   if (check_machine_resources (cpus, cpus * memory_per_cpu)) {

      take_machine_resources (cpus, cpus * memory_per_cpu);

      for (iterator i = m_nodes.begin (); myCPUs && i < m_nodes.end (); n++, i++)
         while (myCPUs && (*i)->check_resources (1, memory_per_cpu) &&
                (*i)->tentative_ratio (1, memory_per_cpu) > RATIO_UPPER_THRESHOLD) {

            (*i)->take_resources (1, memory_per_cpu);
            use [n]++;
            myCPUs--;

         }

      n = 0;

      for (iterator i = m_nodes.begin (); myCPUs && i < m_nodes.end (); n++, i++)
         while (myCPUs && (*i)->check_resources (1, memory_per_cpu)) {

            (*i)->take_resources (1, memory_per_cpu);
            use [n]++;
            myCPUs--;

         }

   }
   else
      return true;

   if (myCPUs) {

      for (n--; n >= 0; n--)
         if (use [n])
            m_nodes [n]->release_resources (use [n], use [n] * memory_per_cpu);

      release_machine_resources (cpus, cpus * memory_per_cpu);
      return true;

   }

   return false;

}

//-----------------------------------------------------------------------------------------------------------------

bool ClusterMachine::worst_available (int cpus, float memory_per_cpu, int use []) {

   int n = 0;
   int myCPUs = cpus;

   memset ((void*) use, 0, m_total_nodes * sizeof (int));

   if (check_machine_resources (cpus, cpus * memory_per_cpu)) {

      for (iterator i = m_nodes.begin (); myCPUs && i < m_nodes.end (); n++, i++)
         while (myCPUs && (*i)->check_resources (1, memory_per_cpu) &&
                (*i)->tentative_ratio (1, memory_per_cpu) < RATIO_LOWER_THRESHOLD) {

            (*i)->take_resources (1, memory_per_cpu);
            use [n]++;
            myCPUs--;

         }

      n = 0;

      for (iterator i = m_nodes.begin (); myCPUs && i < m_nodes.end (); n++, i++)
         while (myCPUs && (*i)->check_resources (1, memory_per_cpu)) {

            (*i)->take_resources (1, memory_per_cpu);
            use [n]++;
            myCPUs--;

         }

   }
   else
     return true;

   if(myCPUs) {

      for (n--; n >= 0; n--)
         if (use [n])
            m_nodes [n]->release_resources (use [n], use [n] * memory_per_cpu);

      release_machine_resources (cpus, cpus * memory_per_cpu);
      return true;

   }

   return false;

}

//-----------------------------------------------------------------------------------------------------------------

std::ostream& ClusterMachine::print (std::ostream& os) const {

   os << "Cluster Machine (" << total_nodes () << " nodos)" << std::endl;

   for (int i = 0; i < total_nodes (); i++)
      os << "   Node " << i << " : "
         << "Memory " << m_nodes [i]->used_memory ()
         << "/" << m_nodes [i]->total_memory () << "MB "
         << "Processors " << m_nodes [i]->used_cpus ()
         << "/" << m_nodes [i]->total_cpus () << "CPUs "
         << "Ratio (Mem:" << m_nodes [i]->memory_usage ()
         << ", CPU:" << m_nodes [i]->cpus_usage () << ")"
         << std::endl;

   return os;

}

//-----------------------------------------------------------------------------------------------------------------

void ClusterMachine::read_from_file (const char* filename) {

   FILE* f;

   f = fopen (filename, "r");

   if (!f)
      perror ("fopen");
   else {

      char line [1024];
      int nodes, cpus, used_cpus;
      float memory, used_memory;

      this->destroy ();

      m_total_nodes  = 0;
      m_total_cpus   = 0;
      m_used_cpus    = 0;
      m_total_memory = 0.0;
      m_used_memory  = 0.0;

      while (fgets (line, sizeof (line), f)) {

         if (line [0] == '#')
            continue;

         if (sscanf (line, "%d,%d,%f", &nodes, &cpus, &memory) == 3) {
            for (int n = 0; n < nodes; n++)
               m_nodes.push_back (new ClusterNode (cpus, memory));

            m_total_nodes  += nodes;
            m_total_cpus   += (nodes * cpus);
            m_total_memory += (nodes * memory);

         }
         else if (sscanf (line, "%d,%f,%d,%f", &cpus, &memory, &used_cpus, &used_memory) == 4) {

            ClusterNode* node = new ClusterNode (cpus, memory);
            node->take_resources (used_cpus, used_memory);
            m_nodes.push_back (node);

            m_total_nodes++;

            m_total_cpus   += cpus;
            m_used_cpus    += used_cpus;
            m_total_memory += memory;
            m_used_memory  += used_memory;

         }

      }

      fclose (f);

   }

}

//-----------------------------------------------------------------------------------------------------------------

void ClusterMachine::print_to_file (const char* filename) const {

   FILE* f = fopen (filename, "wb");

   if (!f)
      perror ("fopen");
   else {

      char line [1024];

      sprintf (line, "#CPUS, MEMORY, CPUS_USED, MEMORY_USED\n");
      fputs (line, f);

      for (int i = 0; i < total_nodes (); i++) {

         sprintf (line, "%d,%f,%d,%f\n", m_nodes [i]->total_cpus   (),
                                         m_nodes [i]->total_memory (),
                                         m_nodes [i]->used_cpus    (),
                                         m_nodes [i]->used_memory  ());

         fputs (line, f);

      }

      fclose (f);

   }

}

//-----------------------------------------------------------------------------------------------------------------

void ClusterMachine::destroy () {

   for (int i = 0; i < m_total_nodes; i++)
      if (m_nodes [i])
         delete m_nodes [i];

   m_nodes.clear ();

}

//-----------------------------------------------------------------------------------------------------------------

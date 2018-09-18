#include "Job.h"

#include <string.h>
#include <stdlib.h>

//-----------------------------------------------------------------------------------------------------------------

Job::Job (const char* label, int required_cpus, float memory_per_cpu, float expect_exec_time,
          float real_exec_time, float timestamp, const char* jobClass) : m_required_cpus  (required_cpus),
                                                                         m_memory_per_cpu (memory_per_cpu),
                                                                         m_expect_time    (expect_exec_time),
                                                                         m_real_time      (real_exec_time),
                                                                         m_timestamp      (timestamp)         {

   if (strlen (label))
      m_label = strdup (label);
   else
      m_label = NULL;

   if (strlen (jobClass))
      m_class = strdup (jobClass);
   else
      m_class = NULL;

}

//-----------------------------------------------------------------------------------------------------------------

Job::~Job () {

   if (m_label)
      free (m_label);

   if (m_class)
      free (m_class);

}

//-----------------------------------------------------------------------------------------------------------------

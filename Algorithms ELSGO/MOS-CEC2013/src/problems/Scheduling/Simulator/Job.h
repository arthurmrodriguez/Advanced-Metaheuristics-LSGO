#ifndef _JOB_H
#define _JOB_H

#include <vector>

class Job {

   public:

      Job (const char* label,
                 int   required_cpus,
                 float memory_per_cpu,
                 float expect_exec_time,
                 float real_exec_time,
                 float timestamp,
	        const char* jobClass);

      virtual ~Job();

      const char* label    () const {return m_label;}
      int   required_cpus  () const {return m_required_cpus;}
      float memory_per_cpu () const {return m_memory_per_cpu;}
      float expected_time  () const {return m_expect_time;}
      float real_time      () const {return m_real_time;}
      float timestamp      () const {return m_timestamp;}
      const char* jobClass	() const {return m_class;}


   private:

      char*   m_label;
      int     m_required_cpus;
      float   m_memory_per_cpu;
      float   m_expect_time;
      float   m_real_time;
      float   m_timestamp;
      char*   m_class;

};

#endif

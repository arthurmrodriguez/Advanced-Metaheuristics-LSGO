/**
 * @mainpage
 * @author Antonio Latorre, Alberto Molina.
 * @date February-2008
 *
 * This module generates valid files of jobs and 
 * classes to use as input data in the simulator
 */

/**
 * @file
 * @brief Input data generator to the simulator.
 *
 * This module generates valid files of jobs and 
 * classes to use as input data in the simulator
 */

// For file streams
#include <fstream>

#include <iostream>

// Guess why...
#include <vector>

// For setting precission in output streams
#include <iomanip>

// For sort
#include <algorithm>

// For rand, srand and gettimeofday
#include <sys/time.h>
#include <time.h>

// For assert (used in normal_random)
#include <assert.h>

#include <math.h>
#include <errno.h>
#include <limits.h>
#include <string.h>

//-----------------------------------------------------------------------------------------------------------------

/**
 * @struct Job
 * @brief Represents a job that execute into simulator.
 */

typedef struct Job {

   std::string name;
   std::string jobClass;
   double      timestamp;
   int         cpus;
   double      mem;
   double      time;

} TJob;

/**
 * @struct PenaltyStruct
 * @brief Represents a class penalty, that can be assigned to a job.
 */

struct PenaltyStruct {

      char* class_name;
      float penalty_limit;
      float penalty_per_unit;
      float probability_percent;

};

//-----------------------------------------------------------------------------------------------------------------

/**
 * @defgroup old_defines Old defines of datasets generation
 * @brief Defines for old datasets generation.
 */

//-------------------------------------------

/**
 * @def SMALL_CPU_AVG
 * @brief Small average of cpus.
 * @ingroup old_defines
 */
#define SMALL_CPU_AVG     10

/**
 * @def BIG_CPU_AVG
 * @brief Big average of cpus.
 * @ingroup old_defines
 */
#define BIG_CPU_AVG      100

/**
 * @def NORMAL_CPU_AVG
 * @brief Normal average of cpus.
 * @ingroup old_defines
 */
#define NORMAL_CPU_AVG    50

/**
 * @def WEIRD_CPU_AVG
 * @brief Weird average of cpus.
 * @ingroup old_defines
 */
#define WEIRD_CPU_AVG    120

//-------------------------------------------

/**
 * @def SMALL_MEM_AVG
 * @brief Small average of memory.
 * @ingroup old_defines
 */
#define SMALL_MEM_AVG   1024

/**
 * @def BIG_MEM_AVG
 * @brief Big average of memory.
 * @ingroup old_defines
 */
#define BIG_MEM_AVG     3072

/**
 * @def NORMAL_MEM_AVG
 * @brief Normal average of memory.
 * @ingroup old_defines
 */
#define NORMAL_MEM_AVG  2048

/**
 * @def WEIRD_MEM_AVG
 * @brief Weird average of memory.
 * @ingroup old_defines
 */
#define WEIRD_MEM_AVG   2048

//-------------------------------------------

/**
 * @def SMALL_TIME_AVG
 * @brief Small average of time.
 * @ingroup old_defines
 */
#define SMALL_TIME_AVG    10

/**
 * @def BIG_TIME_AVG
 * @brief Big average of time.
 * @ingroup old_defines
 */
#define BIG_TIME_AVG      30

/**
 * @def NORMAL_TIME_AVG
 * @brief Normal average of time.
 * @ingroup old_defines
 */
#define NORMAL_TIME_AVG   20

/**
 * @def WEIRD_TIME_AVG
 * @brief Weird average of time.
 * @ingroup old_defines
 */
#define WEIRD_TIME_AVG    50

//-------------------------------------------

/**
 * @def SMALL_VAR
 * @brief Small variable.
 * @ingroup old_defines
 */
#define SMALL_VAR      1.256

/**
 * @def BIG_VAR
 * @brief Big variable.
 * @ingroup old_defines
 */
#define BIG_VAR        3.024

/**
 * @def NORMAL_VAR
 * @brief Normal variable.
 * @ingroup old_defines
 */
#define NORMAL_VAR     2.576

/**
 * @def WEIRD_VAR
 * @brief Weird variable.
 * @ingroup old_defines
 */
#define WEIRD_VAR      6.600

//-----------------------------------------------------------------------------------------------------------------

/**
 * @defgroup new_defines New defines of datasets generation
 * @brief Defines for new datasets generation.
 */

//-------------------------------------------

/**
 * @def SMALL_PROC_MIN
 * @brief Minimum processor time small.
 * @ingroup new_defines
 */
#define SMALL_PROC_MIN           4

/**
 * @def SMALL_PROC_MAX
 * @brief Maximum processor time small.
 * @ingroup new_defines
 */
#define SMALL_PROC_MAX          12

/**
 * @def SMALL_PROC_STEP
 * @brief Step processor time small.
 * @ingroup new_defines
 */
#define SMALL_PROC_STEP          1

/**
 * @def NORMAL_PROC_MIN
 * @brief Minimum processor time normal.
 * @ingroup new_defines
 */
#define NORMAL_PROC_MIN         10

/**
 * @def NORMAL_PROC_MAX
 * @brief Maximum processor time normal.
 * @ingroup new_defines
 */
#define NORMAL_PROC_MAX         70

/**
 * @def NORMAL_PROC_STEP
 * @brief Step processor time normal.
 * @ingroup new_defines
 */
#define NORMAL_PROC_STEP         2

/**
 * @def BIG_PROC_MIN
 * @brief Minimum processor time big.
 * @ingroup new_defines
 */
#define BIG_PROC_MIN            64

/**
 * @def BIG_PROC_MAX
 * @brief Maximum processor time big.
 * @ingroup new_defines
 */
#define BIG_PROC_MAX           1024

/**
 * @def BIG_PROC_STEP
 * @brief Step processor time big.
 * @ingroup new_defines
 */
#define BIG_PROC_STEP           32

//-------------------------------------------

/**
 * @def SMALL_MEM_MIN
 * @brief Minimum memory small.
 * @ingroup new_defines
 */
#define SMALL_MEM_MIN            5

/**
 * @def SMALL_MEM_MAX
 * @brief Maximum memory small.
 * @ingroup new_defines
 */
#define SMALL_MEM_MAX           12

/**
 * @def NORMAL_MEM_MIN
 * @brief Minimum memory normal.
 * @ingroup new_defines
 */
#define NORMAL_MEM_MIN          10

/**
 * @def NORMAL_MEM_MAX
 * @brief Maximum memory normal.
 * @ingroup new_defines
 */
#define NORMAL_MEM_MAX        1000

/**
 * @def BIG_MEM_MIN
 * @brief Minimum memory big.
 * @ingroup new_defines
 */
#define BIG_MEM_MIN           1000

/**
 * @def BIG_MEM_MAX
 * @brief Maximum memory big.
 * @ingroup new_defines
 */
#define BIG_MEM_MAX           4000

//-------------------------------------------

/**
 * @def SMALL_LEN_MIN
 * @brief Minimum length small.
 * @ingroup new_defines
 */
#define SMALL_LEN_MIN         1

/**
 * @def SMALL_LEN_MAX
 * @brief Maximum length small.
 * @ingroup new_defines
 */
#define SMALL_LEN_MAX        7

/**
 * @def NORMAL_LEN_MIN
 * @brief Minimum length normal.
 * @ingroup new_defines
 */
#define NORMAL_LEN_MIN       7

/**
 * @def NORMAL_LEN_MAX
 * @brief Maximum length normal.
 * @ingroup new_defines
 */
#define NORMAL_LEN_MAX      48

/**
 * @def BIG_LEN_MIN
 * @brief Maximum length big.
 * @ingroup new_defines
 */
#define BIG_LEN_MIN         40

/**
 * @def BIG_LEN_MAX
 * @brief Minimum length big.
 * @ingroup new_defines
 */
#define BIG_LEN_MAX        168

//-----------------------------------------------------------------------------------------------------------------

/**
 * @defgroup prob_proc_class Defines of probabilities
 * @brief Probabilities of the class processor.
 */

//-------------------------------------------

/**
 * @def PROB_PS
 * @brief Small probability.
 * @ingroup prob_proc_class
 */
#define PROB_PS      0.20

/**
 * @def PROB_PN
 * @brief Normal probability.
 * @ingroup prob_proc_class
 */
#define PROB_PN      0.55

/**
 * @def PROB_PB
 * @brief Big probability.
 * @ingroup prob_proc_class
 */
#define PROB_PB      0.24

//-----------------------------------------------------------------------------------------------------------------

/**
 * @defgroup prob_mem_class_small Defines of probabilities of memory for classes with small processor time 
 * @brief Probabilities of memory for classes conditioned with processor_class = SMALL.
 */

//-------------------------------------------

/**
 * @def PROB_MS_PS
 * @brief Small processor time, small probability of memory.
 * @ingroup prob_mem_class_small
 */
#define PROB_MS_PS   0.05

/**
 * @def PROB_MN_PS
 * @brief Small processor time, normal probability of memory.
 * @ingroup prob_mem_class_small
 */
#define PROB_MN_PS   0.35

/**
 * @def PROB_MB_PS
 * @brief Small processor time, big probability of memory.
 * @ingroup prob_mem_class_small
 */
#define PROB_MB_PS   0.60

//-----------------------------------------------------------------------------------------------------------------

/**
 * @defgroup prob_mem_class_normal Defines of probabilities of memory for classes with normal processor time 
 * @brief Probabilities of memory for classes conditioned with processor_class = NORMAL.
 */

//-------------------------------------------

/**
 * @def PROB_MS_PN
 * @brief Normal processor time, small probability of memory.
 * @ingroup prob_mem_class_normal
 */
#define PROB_MS_PN   0.05

/**
 * @def PROB_MN_PN
 * @brief Normal processor time, normal probability of memory.
 * @ingroup prob_mem_class_normal
 */
#define PROB_MN_PN   0.55

/**
 * @def PROB_MB_PN
 * @brief Normal processor time, big probability of memory.
 * @ingroup prob_mem_class_normal
 */
#define PROB_MB_PN   0.40

//-----------------------------------------------------------------------------------------------------------------

/**
 * @defgroup prob_mem_class_big Defines of probabilities of memory for classes with big processor time 
 * @brief Probabilities of memory for classes conditioned with processor_class = BIG.
 */

//-------------------------------------------

/**
 * @def PROB_MS_PB
 * @brief Big processor time, small probability of memory.
 * @ingroup prob_mem_class_big
 */
#define PROB_MS_PB   0.35

/**
 * @def PROB_MN_PB
 * @brief Big processor time, normal probability of memory.
 * @ingroup prob_mem_class_big
 */
#define PROB_MN_PB   0.30

/**
 * @def PROB_MB_PB
 * @brief Big processor time, big probability of memory.
 * @ingroup prob_mem_class_big
 */
#define PROB_MB_PB   0.35

//-----------------------------------------------------------------------------------------------------------------

/**
 * @defgroup prob_len_class_small Defines of probabilities of length for classes with small memory 
 * @brief Probabilities of length for classes conditioned with memory = SMALL.
 */

//-------------------------------------------

/**
 * @def PROB_LS_MS
 * @brief Small memory, small probability of length.
 * @ingroup prob_len_class_small
 */
#define PROB_LS_MS   0.65

/**
 * @def PROB_LN_MS
 * @brief Small memory, normal probability of length.
 * @ingroup prob_len_class_small
 */
#define PROB_LN_MS   0.25

/**
 * @def PROB_LB_MS
 * @brief Small memory, big probability of length.
 * @ingroup prob_len_class_small
 */
#define PROB_LB_MS   0.10

//-----------------------------------------------------------------------------------------------------------------

/**
 * @defgroup prob_len_class_normal Defines of probabilities of length for classes with normal memory 
 * @brief Probabilities of length for classes conditioned with memory = NORMAL.
 */

//-------------------------------------------

/**
 * @def PROB_LS_MN
 * @brief Normal memory, small probability of length.
 * @ingroup prob_len_class_normal
 */
#define PROB_LS_MN   0.35

/**
 * @def PROB_LN_MN
 * @brief Normal memory, normal probability of length.
 * @ingroup prob_len_class_normal
 */
#define PROB_LN_MN   0.50

/**
 * @def PROB_LB_MN
 * @brief Normal memory, big probability of length.
 * @ingroup prob_len_class_normal
 */
#define PROB_LB_MN   0.15

//-----------------------------------------------------------------------------------------------------------------

/**
 * @defgroup prob_len_class_big Defines of probabilities of length for classes with big memory 
 * @brief Probabilities of length for classes conditioned with memory = BIG.
 */

//-------------------------------------------

/**
 * @def PROB_LS_MB
 * @brief Big memory, small probability of length.
 * @ingroup prob_len_class_big
 */
#define PROB_LS_MB   0.05

/**
 * @def PROB_LN_MB
 * @brief Big memory, normal probability of length.
 * @ingroup prob_len_class_big
 */
#define PROB_LN_MB   0.65

/**
 * @def PROB_LB_MB
 * @brief Big memory, big probability of length.
 * @ingroup prob_len_class_big
 */
#define PROB_LB_MB   0.30

//-----------------------------------------------------------------------------------------------------------------

/**
 * @defgroup penalty_lim_man New defines for manage penalty limits.
 * @brief Probabilities of penalty limits.
 */

//-------------------------------------------

/**
 * @def MAX_PENALTY_LIMIT
 * @brief Maximum penalty limit.
 * @ingroup penalty_lim_man
 */
#define MAX_PENALTY_LIMIT      10

/**
 * @def MIN_PERCENT_LOW_PENALTY_LIMIT
 * @brief Minimum percent of low penalty limit.
 * @ingroup penalty_lim_man
 */
#define MIN_PERCENT_LOW_PENALTY_LIMIT      0.00

/**
 * @def MAX_PERCENT_LOW_PENALTY_LIMIT
 * @brief Maximum percent of low penalty limit.
 * @ingroup penalty_lim_man
 */
#define MAX_PERCENT_LOW_PENALTY_LIMIT      0.30

/**
 * @def MIN_PERCENT_MIDDLE_PENALTY_LIMIT
 * @brief Minimum percent of middle penalty limit.
 * @ingroup penalty_lim_man
 */
#define MIN_PERCENT_MIDDLE_PENALTY_LIMIT   0.31

/**
 * @def MAX_PERCENT_MIDDLE_PENALTY_LIMIT
 * @brief Maximum percent of middle penalty limit.
 * @ingroup penalty_lim_man
 */
#define MAX_PERCENT_MIDDLE_PENALTY_LIMIT   0.60

/**
 * @def MIN_PERCENT_HIGH_PENALTY_LIMIT
 * @brief Minimum percent of high penalty limit.
 * @ingroup penalty_lim_man
 */
#define MIN_PERCENT_HIGH_PENALTY_LIMIT     0.61

/**
 * @def MAX_PERCENT_HIGH_PENALTY_LIMIT
 * @brief Maximum percent of high penalty limit.
 * @ingroup penalty_lim_man
 */
#define MAX_PERCENT_HIGH_PENALTY_LIMIT     1.00

//-----------------------------------------------------------------------------------------------------------------

/**
 * @defgroup penalty_per_unit New defines for magane penalties per unit
 * @brief Probabilities of penalties per unit.
 */

//-------------------------------------------

/**
 * @def MAX_PENALTY_PER_UNIT
 * @brief Maximum penalty per unit.
 * @ingroup penalty_per_unit
 */
#define MAX_PENALTY_PER_UNIT   10

/**
 * @def MIN_LOW_PENALTY_PER_UNIT_LIMIT
 * @brief Minimum penalty per unit of low penalty limit.
 * @ingroup penalty_per_unit
 */
#define MIN_LOW_PENALTY_PER_UNIT_LIMIT      0.05

/**
 * @def MAX_LOW_PENALTY_PER_UNIT_LIMIT
 * @brief Maximum penalty per unit of low penalty limit.
 * @ingroup penalty_per_unit
 */
#define MAX_LOW_PENALTY_PER_UNIT_LIMIT      0.25

/**
 * @def MIN_MIDDLE_PENALTY_PER_UNIT_LIMIT
 * @brief Minimum penalty per unit of middle penalty limit.
 * @ingroup penalty_per_unit
 */
#define MIN_MIDDLE_PENALTY_PER_UNIT_LIMIT   0.26

/**
 * @def MAX_MIDDLE_PENALTY_PER_UNIT_LIMIT
 * @brief Maximum penalty per unit of middle penalty limit.
 * @ingroup penalty_per_unit
 */
#define MAX_MIDDLE_PENALTY_PER_UNIT_LIMIT   1.00

/**
 * @def MIN_HIGH_PENALTY_PER_UNIT_LIMIT
 * @brief Minimum penalty per unit of high penalty limit.
 * @ingroup penalty_per_unit
 */
#define MIN_HIGH_PENALTY_PER_UNIT_LIMIT     1.01

/**
 * @def MAX_HIGH_PENALTY_PER_UNIT_LIMIT
 * @brief Maximum penalty per unit of high penalty limit.
 * @ingroup penalty_per_unit
 */
#define MAX_HIGH_PENALTY_PER_UNIT_LIMIT     5.00

//-----------------------------------------------------------------------------------------------------------------

/**
 * @def INTERVAL_LEN
 * @brief Average of number of jobs for each time of the day, and length interval.
 */
// FIXME: This should be automatically calculated.
#define INTERVAL_LEN 4

/**
 * @def ALPHA
 * @brief Average of number of jobs for the adjustment factor.
 */
#define ALPHA        2.8

/**
 * @var weekday
 * @brief Frequency of jobs worked on week days divided in intervals.
 */
double weekday [] = {0.25, 0.125, 2.0, 1.0, 4.0, 0.5};

/**
 * @var weekend
 * @brief Frequency of jobs worked on weekend days divided in intervals.
 */
double weekend [] = {0.0625, 0.03125, 0.25, 0.25, 0.25, 0.125};

unsigned special_days [] = {0, 4};
unsigned num_special_days = sizeof (special_days) / sizeof (unsigned);

//-------------------------------------------

/**
 * @var DaysOfWeek
 * @brief Order of the week days
 */
const char* DaysOfWeek[] = {"Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday", "Monday"};

/**
 * @var WeekendDays
 * @brief Weekend days (for unsorted DaysOfWeek arrays)
 */
int WeekendDays[] = {4, 5};

/**
 * @var MonthAndYear
 * @brief Month and year string, for printing purposes only
 */
std::string MonthAndYear = "jan. 2008";

//-------------------------------------------

/**
 * @var available_procs
 * @brief Available processors in the system
 */
// FIXME: This should be read from a machine configuration file.
unsigned available_procs = 1024;

/**
 * @var max_time
 * @brief Maximum simulation time
 */
double max_time;

/**
 * @var act_time
 * @brief Current simulation time
 */
double act_time = 0.001;

/**
 * @var hugeJob
 * @brief Has a Huge job been generated?
 */
bool hugeJob = false;

/**
 * @var classes
 * @brief Array of classes of jobs
 */
std::vector<PenaltyStruct*> classes;

//-----------------------------------------------------------------------------------------------------------------

/**
 * @brief Calculate probability with poisson distribution.
 * @param avg Average
 * @param n N param for poisson probability
 * @return Poisson probability
 */
double poisson_prob (double avg, double n) {

   unsigned fact = 1;

   for (int i = 1; i <= n; i++)
      fact *= i;

   return (((exp (avg * (-1))) * (pow (avg, n))) / fact);

}

//-----------------------------------------------------------------------------------------------------------------

/**
 * @brief Generate a random number using normal distribution
 * @param mean Mean for normal distribution
 * @param variance Variance for normal distribution
 * @return Return random number generated
 */
double normal_random (double mean, double variance) {

   int i;
   double x = 0;
   double normal_01;
   double normal;

   // Sanity check
   assert (variance >= 0);

   // First, we generate a N(0,1) random number based on the Central

   // Average 12 U(0,1) variables.
   for (i = 0; i < 12; i++)
      x += rand ();

   x /= RAND_MAX;

   // Normalize.
   normal_01 = x - 6;
   normal = mean + normal_01 * sqrt (variance * mean);

   if (normal < 0)
      normal = mean / 10;

   return normal;

}

//-----------------------------------------------------------------------------------------------------------------

// strtol function with error checking
long strtol2 (char* str, int base) {

   char* endptr;
   long val;

   errno = 0;  // To distinguish success/failure after call
   val = strtol (str, &endptr, base);

   // Check for various possible errors
   if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) ||
       (errno != 0 && val == 0)) {

      std::cerr << "strtol2: Error while converting string." << std::endl;
      return 0;

   }

   if (endptr == str) {

      std::cerr << "strtol2: No digits were found." << std::endl;
      return 0;

   }

   if (*endptr != '\0') // Not necessarily an error...
      std::cerr << "strtol2: Further characters after number: " << endptr << std::endl;

   return val;

}

//-----------------------------------------------------------------------------------------------------------------

// Random number between an interval
double rand2 (double min, double max) {

   double diff = max - min;
   double n = (rand() / (RAND_MAX + 1.0));

   return (min + n * diff);

}

//-----------------------------------------------------------------------------------------------------------------

/**
 * @brief Calculate actual day with a timestamp.
 * @param act_time Total ejecution time
 * @return Return calculate day
 */
unsigned day (double timestamp) {

   return ((unsigned) timestamp) / 24;

}

//-----------------------------------------------------------------------------------------------------------------

/**
 * @brief Calculate the number of jobs with a total ejecution time.
 * @param act_time Total ejecution time
 * @return Return calculate number of jobs
 */
unsigned number_of_jobs (double act_time) {

   // Determine current time slot
   int index = ((int) (act_time - (double) ((((unsigned) act_time) / 24) * 24.0))) / INTERVAL_LEN;

   unsigned mday = day (act_time);

   double factor = 1.0;

   for (unsigned i = 0; i < num_special_days; i++)
      if (mday == special_days [i])
         factor = 6.0;

   // Obtain average number of jobs for the current time slot
   int wday = mday % 7;

   // Average number oj jobs for a time slot
   double avg = 0.0;

   // Obtain average number of jobs for current time slot
   if (wday == WeekendDays [0] || wday == WeekendDays [1])  // Weekend
      avg = weekend [index]  * ALPHA * factor;
   else  // Week day
      avg = weekday [index] * ALPHA * factor;

   // Calculate number of jobs
   int    jobs = -1; // First posibility is to have 0 jobs
   double acum = 0.0;
   double j = rand () / (RAND_MAX + 1.0);

   while (j > acum) {

      jobs++;
      acum += poisson_prob (avg, jobs);

   }

   return jobs;

}

//-----------------------------------------------------------------------------------------------------------------

/**
 * @brief Check order between two jobs
 * @param j1 First job
 * @param j2 Second job
 * @return Return if order it's ok
 */
bool jobs_order (TJob j1, TJob j2) {

   return j1.timestamp < j2.timestamp;

}

//-----------------------------------------------------------------------------------------------------------------

/**
 * @brief Generate a class with random values
 * @param random Random param used to generate random values of the class.
 * @return Return name of class generated
 */
char* random_class(float random)
{
   float last_percent = 0.0;

   for (unsigned int i = 0; i < classes.size (); i++) {

      if ((random >= last_percent) && (random < last_percent + classes [i]->probability_percent))
         return classes [i]->class_name;
      else
         last_percent += classes [i]->probability_percent;

   }

   return random_class (rand2 (0.0, 1.0));

}

//-----------------------------------------------------------------------------------------------------------------

/**
 * @brief Generate a job with random values.
 * @return Return the job generated.
 */
TJob generate_job (void) {

   TJob job;

   double j;
   // double avg, var;

   double prob_ms, prob_mn, prob_mb;
   double prob_ls, prob_ln, prob_lb;


   // Random number of processors
   j = rand () / (RAND_MAX + 1.0);

   if (j >= 0 && j <= PROB_PS) {

      prob_ms = PROB_MS_PS;
      prob_mn = PROB_MN_PS;
      prob_mb = PROB_MB_PS;

      job.name += 'S';

      job.cpus = (int) rand2 (0, (SMALL_PROC_MAX - SMALL_PROC_MIN) / SMALL_PROC_STEP);
      job.cpus = SMALL_PROC_MIN + ((job.cpus + 1) * SMALL_PROC_STEP);

   }
   else if (j > PROB_PS && j <= PROB_PS + PROB_PN) {

      prob_ms = PROB_MS_PN;
      prob_mn = PROB_MN_PN;
      prob_mb = PROB_MB_PN;

      job.name += 'N';

      job.cpus = (int) rand2 (0, (NORMAL_PROC_MAX - NORMAL_PROC_MIN) / NORMAL_PROC_STEP);
      job.cpus = NORMAL_PROC_MIN + ((job.cpus + 1) * NORMAL_PROC_STEP);

   }
   else if ((j > PROB_PS + PROB_PN && j <= PROB_PS + PROB_PN + PROB_PB) || hugeJob) {

      prob_ms = PROB_MS_PB;
      prob_mn = PROB_MN_PB;
      prob_mb = PROB_MB_PB;

      job.name += 'B';

      job.cpus = (int) rand2 (0, (BIG_PROC_MAX - BIG_PROC_MIN) / BIG_PROC_STEP);
      job.cpus = BIG_PROC_MIN + ((job.cpus + 1) * BIG_PROC_STEP);

   }
   else {

      job.name += "BBH";

      job.cpus = (int) rand2 (0, (BIG_PROC_MAX - BIG_PROC_MIN) / BIG_PROC_STEP);
      job.cpus = BIG_PROC_MIN + ((job.cpus + 1) * BIG_PROC_STEP);

      job.mem = rand2 (BIG_MEM_MIN, BIG_MEM_MAX);
      job.time = 1.5 * (rand2 (act_time, max_time));

      hugeJob = true;

      job.jobClass = random_class(j);

      return job;

   }

   // Random amount of memory
   j = rand () / (RAND_MAX + 1.0);

   if (j >= 0 && j <= prob_ms) {

      prob_ls = PROB_LS_MS;
      prob_ln = PROB_LN_MS;
      prob_lb = PROB_LB_MS;

      job.name += 'S';

      job.mem = rand2 (SMALL_MEM_MIN, SMALL_MEM_MAX);

   }
   else if (j > prob_ms && j <= prob_ms + prob_mn) {

      prob_ls = PROB_LS_MN;
      prob_ln = PROB_LN_MN;
      prob_lb = PROB_LB_MN;

      job.name += 'N';

      job.mem = rand2 (NORMAL_MEM_MIN, NORMAL_MEM_MAX);

   }
   else {

      prob_ls = PROB_LS_MB;
      prob_ln = PROB_LN_MB;
      prob_lb = PROB_LB_MB;

      job.name += 'B';

      job.mem = rand2 (BIG_MEM_MIN, BIG_MEM_MAX);

   }

   // Random execution time
   j = rand () / (RAND_MAX + 1.0);

   if (j >= 0 && j <= prob_ls) {

      job.name += 'S';

      job.time = rand2 (SMALL_LEN_MIN, SMALL_LEN_MAX);

   }
   else if (j > prob_ls && j <= prob_ls + prob_ln) {

      job.name += 'N';

      job.time = rand2 (NORMAL_LEN_MIN, NORMAL_LEN_MAX);

   }
   else {

      job.name += 'B';

      job.time = rand2 (BIG_LEN_MIN, BIG_LEN_MAX);

   }

   job.jobClass = random_class(j);

   return job;

}

//-----------------------------------------------------------------------------------------------------------------

/**
 * @brief Parse classes file to memory
 * @param fileclasses Name of file that contents the classes
 */
void parse (char* classesfile) {

   std::ifstream f (classesfile);

   if (f.bad ())
      std::cerr << "ifstream: error opening file: " << classesfile << std::endl;
   else {

      char line [1024];

      // Leemos la cabecera
      f.getline (line, 1024);

      if (line [0] != '#')
         std::cerr << "Error: cabecera no encontrada en " << classesfile << std::endl;

      classes.clear();

      while (!f.eof ()) {

         PenaltyStruct* ps = new PenaltyStruct ();
         ps->class_name = new char [128];

         f >> ps->class_name;

         if (f.eof ())
            break;

         f >> ps->penalty_limit;
         f >> ps->penalty_per_unit;
         f >> ps->probability_percent;

         classes.push_back (ps);

      }

   }

}

//-----------------------------------------------------------------------------------------------------------------

/**
 * @brief Generate a classes file to use in datasets generations.
 * @param fileclasses Name of file where will saved the new classes
 * @param classes_numer Number of classes generated
 * @see generate_classes_from_file
 */
void generate_classes_file (char* fileclasses, int classes_number) {

   int type = 0;
   float free_percent = 1.0 / classes_number;

   std::ofstream of (fileclasses);

   of << std::setiosflags (std::ios::fixed | std::ios::showpoint)
      << std::setprecision (5);

   of << "#CLASS, LIMIT, PENALTY_PER_UNIT, PROBABILITY_PERCENT" << std::endl;

   classes.clear ();

   for (int j = 0; j < (classes_number / 26) + 1; j++) {

      for (int i = 65; (i <= 90) && (26 * j + i - 65 < classes_number); i++) {

         PenaltyStruct* ps = new PenaltyStruct ();
         float limit, penalty_per_unit;

         ps->class_name = new char [128];

	      switch (type % 9) {

            case 0 :

               limit = rand2 (MIN_PERCENT_HIGH_PENALTY_LIMIT, MAX_PERCENT_HIGH_PENALTY_LIMIT);
               penalty_per_unit = rand2 (MIN_LOW_PENALTY_PER_UNIT_LIMIT, MAX_LOW_PENALTY_PER_UNIT_LIMIT);

               break;

            case 1 :

               limit = rand2 (MIN_PERCENT_HIGH_PENALTY_LIMIT, MAX_PERCENT_HIGH_PENALTY_LIMIT);
               penalty_per_unit = rand2 (MIN_MIDDLE_PENALTY_PER_UNIT_LIMIT, MAX_MIDDLE_PENALTY_PER_UNIT_LIMIT);

               break;

            case 2 :

               limit = rand2 (MIN_PERCENT_HIGH_PENALTY_LIMIT, MAX_PERCENT_HIGH_PENALTY_LIMIT);
               penalty_per_unit = rand2 (MIN_HIGH_PENALTY_PER_UNIT_LIMIT, MAX_HIGH_PENALTY_PER_UNIT_LIMIT);

               break;

            case 3 :

               limit = rand2 (MIN_PERCENT_MIDDLE_PENALTY_LIMIT, MAX_PERCENT_MIDDLE_PENALTY_LIMIT);
               penalty_per_unit = rand2 (MIN_LOW_PENALTY_PER_UNIT_LIMIT, MAX_LOW_PENALTY_PER_UNIT_LIMIT);

               break;

            case 4 :

               limit = rand2 (MIN_PERCENT_MIDDLE_PENALTY_LIMIT, MAX_PERCENT_MIDDLE_PENALTY_LIMIT);
               penalty_per_unit = rand2 (MIN_MIDDLE_PENALTY_PER_UNIT_LIMIT, MAX_MIDDLE_PENALTY_PER_UNIT_LIMIT);

               break;

            case 5 :

               limit = rand2 (MIN_PERCENT_MIDDLE_PENALTY_LIMIT, MAX_PERCENT_MIDDLE_PENALTY_LIMIT);
               penalty_per_unit = rand2 (MIN_HIGH_PENALTY_PER_UNIT_LIMIT, MAX_HIGH_PENALTY_PER_UNIT_LIMIT);

               break;

            case 6 :

               limit = rand2 (MIN_PERCENT_LOW_PENALTY_LIMIT, MAX_PERCENT_LOW_PENALTY_LIMIT);
               penalty_per_unit = rand2 (MIN_LOW_PENALTY_PER_UNIT_LIMIT, MAX_LOW_PENALTY_PER_UNIT_LIMIT);

               break;

            case 7 :

               limit = rand2 (MIN_PERCENT_LOW_PENALTY_LIMIT, MAX_PERCENT_LOW_PENALTY_LIMIT);
               penalty_per_unit = rand2 (MIN_MIDDLE_PENALTY_PER_UNIT_LIMIT, MAX_MIDDLE_PENALTY_PER_UNIT_LIMIT);

               break;

            default:
            case 8 :

               limit = rand2 (MIN_PERCENT_LOW_PENALTY_LIMIT, MAX_PERCENT_LOW_PENALTY_LIMIT);
               penalty_per_unit = rand2 (MIN_HIGH_PENALTY_PER_UNIT_LIMIT, MAX_HIGH_PENALTY_PER_UNIT_LIMIT);

               break;
         }

         for (int z = 0; z < j + 1; z++) {

            of << (char) i;
            ps->class_name [z] = (char) i;

         }

         of <<  " " << limit << " " << penalty_per_unit << " " << free_percent << std::endl;// probability_percent << std::endl;

         ps->penalty_limit = limit;
         ps->penalty_per_unit = penalty_per_unit;
         ps->probability_percent = free_percent; //probability_percent;

         classes.push_back(ps);

         type++;

      }

   }

   of.close ();

}

//-----------------------------------------------------------------------------------------------------------------

/**
 * @brief Generate a classes file, to use in datasets generations, with a names file
 * @param fileclasses Name of file where will saved the new classes
 * @param classes_numer Name of file that contents the name
 * @see generate_classes_file
 */
void generate_classes_from_file (char* fileclasses, char* filenames) {

   FILE* f = fopen (filenames, "r");
   std::vector<char*> names;
   float free_percent = 1.0;
   int classes_number;

   std::ofstream of (fileclasses);

   if (!f) {

      perror ("fopen");

   }
   else {

      char line [1024];

      while (fgets (line, sizeof (line), f)) {

         char* class_name = new char [128];

         if (line [0] == '#')
            continue;

         if (sscanf (line, "%[^ ,]", class_name) == 1) {

            for (unsigned i = 0; i < 128; i++)
               if (class_name [i] == '\n' || class_name [i] == '\r') {

                  class_name [i] = '\0';
                  break;

               }

            names.push_back (class_name);

         }
         else {

            perror ("sscanf");
            fclose (f);

         }

      }

      fclose (f);

   }

   of << std::setiosflags (std::ios::fixed | std::ios::showpoint)
      << std::setprecision (5);

   of << "#CLASS, LIMIT, PENALTY_PER_UNIT, PROBABILITY_PERCENT" << std::endl;

   classes_number = names.size ();
   classes.clear ();

   for (int j = 0; j < classes_number; j++) {

      PenaltyStruct* ps = new PenaltyStruct ();
      float limit, penalty_per_unit, probability_percent;
      ps->class_name = new char [128];

      strcpy(ps->class_name, names [j]);
      limit = rand2 (0.0, MAX_PENALTY_LIMIT);
      penalty_per_unit = rand2 (0.0, MAX_PENALTY_PER_UNIT);

      probability_percent = free_percent;

      if (j != classes_number) {

         while (probability_percent >= free_percent)
            probability_percent = normal_random (free_percent / classes_number, 1.0 / classes_number);

         free_percent -= probability_percent;

      }

      of <<  ps->class_name << " " << limit << " " << penalty_per_unit << " " << probability_percent << std::endl;

      ps->penalty_limit = limit;
      ps->penalty_per_unit = penalty_per_unit;
      ps->probability_percent = probability_percent;
      classes.push_back (ps);

   }

   of.close ();

}

//-----------------------------------------------------------------------------------------------------------------

/**
 * @brief Old version of generate datasets function.
 * @param argv Array of of necessary parameters.
 * @see generate_dataset_new
 */
void generate_dataset_old (char** argv) {

   float cpu, mem, t_time;
   int jobs = atoi (argv [1]);

   int n_jobs = 0;
   float total_time = 0.0;

   if (jobs <= 0)
      std::cerr << "Usage: " << argv [0] << " num_jobs file_dest" << std::endl;


   std::ofstream of (argv [2]);


   of << "#LABEL, CPUS, MEMORY, EXPECT_TIME, REAL_TIME" << std::endl;


   for (int i = 1; i <= jobs; i++) {

      int j;
      char name [4];

      name [3] = '\0';

      // Random number of CPU's
      j = 1 + (int) (100.0 * (rand() / (RAND_MAX + 1.0)));

      if (j >= 0 and j <= 15) {

         cpu = normal_random (SMALL_CPU_AVG,  SMALL_VAR);
         name [0] = 'S';

      }
      else if (j > 15 and j <= 45) {

         cpu = normal_random (BIG_CPU_AVG,  BIG_VAR);
         name [0] = 'B';

      }
      else if (j > 45 and j <= 95) {

         cpu = normal_random (NORMAL_CPU_AVG,  NORMAL_VAR);
         name [0] = 'N';

      }
      else {

         cpu = normal_random (WEIRD_CPU_AVG,  WEIRD_VAR);
         name [0] = 'W';

      }



      // Random number of CPU's

      j = 1 + (int) (100.0 * (rand() / (RAND_MAX + 1.0)));

      if (j >= 0 and j <= 15) {

         mem = normal_random (SMALL_MEM_AVG,  SMALL_VAR);
         name [1] = 'S';

      }
      else if (j > 15 and j <= 45) {

         mem = normal_random (BIG_MEM_AVG,  BIG_VAR);
         name [1] = 'B';

      }
      else if (j > 45 and j <= 95) {

         mem = normal_random (NORMAL_MEM_AVG,  NORMAL_VAR);
         name [1] = 'N';

      }
      else {

         mem = normal_random (WEIRD_MEM_AVG,  WEIRD_VAR);
         name [1] = 'W';

      }


      // Random number of CPU's

      j = 1 + (int) (100.0 * (rand() / (RAND_MAX + 1.0)));

      if (j >= 0 and j <= 15) {

         t_time = normal_random (SMALL_TIME_AVG,  SMALL_VAR);
         name [2] = 'S';

      }
      else if (j > 15 and j <= 45) {

         t_time = normal_random (BIG_TIME_AVG,  BIG_VAR);
         name [2] = 'B';

      }
      else if (j > 45 and j <= 95) {

         t_time = normal_random (NORMAL_TIME_AVG,  NORMAL_VAR);
         name [2] = 'N';

      }
      else {

         t_time = normal_random (WEIRD_TIME_AVG,  WEIRD_VAR);
         name [2] = 'W';

      }

      n_jobs++;

      of << name << n_jobs << ", " << (int) cpu << ", " << mem << ", " << t_time << ", " << t_time << std::endl;

      total_time += t_time;

   }

   of << "#Total sequential time: " << total_time << std::endl;

   of.close ();

}

//-----------------------------------------------------------------------------------------------------------------

/**
 * @brief New version of generate datasets function, generate a simulator jobs file.
 * @param filename Name of file where will saved the dataset
 * @param fileclasses Name of classes file used to generate the dataset
 * @see generate_classes_file, generate_classes_from_file
 */
void generate_dataset_new (char* filename, char* fileclasses) {

   std::vector<TJob> jobs;

   double total_time = 0.0;
   double total_proc_time = 0.0;

   parse (fileclasses);

   std::ofstream of (filename);

   of << std::setiosflags (std::ios::fixed | std::ios::showpoint)
      << std::setprecision (5);

   of << "#LABEL, TIMESTAMP, CPUS, MEMORY, EXPECT_TIME, REAL_TIME, CLASS" << std::endl;


   do {

      unsigned n = number_of_jobs (act_time);

      for (unsigned i = 0; i < n; i++) {

         TJob job = generate_job ();

         job.timestamp = rand2 ((act_time - 0.001), (act_time - 0.001) + INTERVAL_LEN);
         total_time += job.time;
         total_proc_time += (job.time * job.cpus);

         jobs.push_back (job);

      }

      act_time += INTERVAL_LEN;

   } while (act_time <= max_time);


   // Order jobs vector according to the job timestamp
   std::sort (jobs.begin (), jobs.end (), jobs_order);

   int day1 = 0, day2 = 0;
   double last_day = 0.0;

   for (unsigned i = 0; i< jobs.size (); i++) {

      day2 = day (jobs [i].timestamp);

      if (day1 != day2) {

         of << "# End " << DaysOfWeek [(day1 % 7)] << " " << day1 + 1 << " "
            << MonthAndYear << " (" << (last_day / 24) << " procs/day)" << std::endl;

         of << "#######################################################################################" << std::endl;

         day1 = day2;
         last_day = 0.0;

      }

      of << jobs [i].name << i << ", " << jobs [i].timestamp << ", " << jobs [i].cpus << ", " << jobs [i].mem       << ", " << jobs [i].time      << ", " << jobs [i].time << ", " << jobs [i].jobClass << std::endl;
      of << "# " << jobs [i].cpus * jobs [i].time << std::endl;

      last_day += (jobs [i].time * jobs [i].cpus);

   }

   of << "# End " << DaysOfWeek [(day1 % 7)] << " " << day1 + 1 << " "
      << MonthAndYear << " (" << (last_day / 24) << " procs/day)" << std::endl;

   of << "# Total sequential time: "     << total_time                        << std::endl;
   of << "# Total processor time: "      << total_proc_time                   << std::endl;
   of << "# Total processors required: " << total_proc_time / max_time << std::endl;
   of << "# Total machine time: "        << total_proc_time / available_procs << std::endl;

   of.close ();

}

//-----------------------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------
// Main program
//-------------------------------------------------------------------------------

int main (int argc, char** argv) {

   // Check for correct usage of the program
   if (argc < 4) {

      std::cerr << "Usage: " << argv [0]  << " -j" << " max_time output_file classes_file" << std::endl;
      std::cerr << "                      " << " -c" << " classes_file num_classes" << std::endl;
      std::cerr << "                      " << " -f" << " classes_file names_file" << std::endl;
      exit (EXIT_FAILURE);

   }

   if (strcmp("-c",argv[1]) == 0) {

      // Generate classes file with n classes
      generate_classes_file(argv[2], strtol2 (argv [3], 10));

   }
   else if (strcmp("-f",argv[1]) == 0) {

      // Generate classes file from names file
      generate_classes_from_file(argv[2],argv[3]);

   }
   else if (strcmp("-j",argv[1]) == 0) {

      if ((max_time = strtol2 (argv [2], 10)) <= 0) {

         std::cerr << "Error: max_time must be a positive number." << std::endl;
         exit (EXIT_FAILURE);

      }

      // Random seed initialization
      struct timeval tv;
      gettimeofday (&tv, NULL);
      srand ((unsigned) ((tv.tv_sec * 1000000) + tv.tv_usec));

      // Generate jobs file
      generate_dataset_new (argv [3], argv [4]);

   }
   else {

      std::cerr << "Error: posible options -j -c -f." << std::endl;
      exit (EXIT_FAILURE);

   }

   return 0;

}

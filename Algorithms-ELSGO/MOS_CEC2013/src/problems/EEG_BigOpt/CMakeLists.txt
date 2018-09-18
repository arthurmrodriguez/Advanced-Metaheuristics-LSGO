#Added commands
cmake_minimum_required (VERSION 2.8)
cmake_policy(SET CMP0042 NEW)


#EEG_problem only
set (PROBLEMS_SRCS eeg_problem)		
include_directories (${GAEDA_INCLUDE})

foreach (PROBLEM ${PROBLEMS_SRCS})
   add_library (${PROBLEM} SHARED ${PROBLEM}.cc)
   target_link_libraries (${PROBLEM} ${GAEDA_LIBRARY})
   set_target_properties (${PROBLEM} PROPERTIES SUFFIX ".so")
endforeach (PROBLEM)
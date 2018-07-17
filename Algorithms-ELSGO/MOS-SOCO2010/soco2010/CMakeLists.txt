project (soco2010)

include_directories (../gaeda)
link_directories (../gaeda)

set (PROBLEMS_SRCS eeg_problem) #schwefel_problem shifted_ackley shifted_griewank shifted_rastrigin shifted_rosembrock shifted_sphere f7 f8 f9 f10 f11)
#set (PROBLEMS_HYBRID_SRCS f12 f13 f14 f15 f16 f17 f18 f19)
set(CMAKE_MACOS_RPATH OFF)

foreach (PROBLEM ${PROBLEMS_SRCS})
  add_library (${PROBLEM} SHARED ${PROBLEM}.cc)
  target_link_libraries (${PROBLEM} 'gaeda')
  set_target_properties (${PROBLEM} PROPERTIES SUFFIX ".so")
endforeach (PROBLEM)

#foreach (PROBLEM ${PROBLEMS_HYBRID_SRCS})
#  add_library (${PROBLEM} SHARED ${PROBLEM}.cc aux.cc)
#  target_link_libraries (${PROBLEM} 'gaeda')
#  set_target_properties (${PROBLEM} PROPERTIES SUFFIX ".so")
#endforeach (PROBLEM)

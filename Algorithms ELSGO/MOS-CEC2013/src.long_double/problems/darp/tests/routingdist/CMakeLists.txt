find_library(GOOGLE_TEST_LIB gtest $ENV{CMAKE_PREFIX_PATH})
find_path   (GOOGLE_TEST_INCLUDE gtest/gtest.h $ENV{CMAKE_PREFIX_PATH})

include_directories (${GAEDA_INCLUDE} ${GOOGLE_TEST_INCLUDE})

link_directories (${CMAKE_CURRENT_SOURCE_DIR})

set (TESTS_SRCS routedistributortests routinggenomenewmethods randomroutedistributortests reqdistmatrixtests routedistmatrixtests )

foreach(TEST ${TESTS_SRCS})
  add_executable (${TEST} ${TEST} ../../aux.cc)
  target_link_libraries (${TEST} ${GAEDA_LIBRARY} ${GOOGLE_TEST_LIB})
endforeach (TEST)

set (TESTS_CLUST_SRCS kmedoidsreqdistributortests kmedoidsroutedisttests)
# These files must be compiled with the cluster file
foreach(TEST ${TESTS_CLUST_SRCS})
  add_executable (${TEST} ${TEST} ../../aux.cc ../../extras/cluster.cc)
  target_link_libraries (${TEST} ${GAEDA_LIBRARY} ${GOOGLE_TEST_LIB})
endforeach (TEST)

find_library(GOOGLE_TEST_LIB gtest $ENV{CMAKE_PREFIX_PATH})
find_path   (GOOGLE_TEST_INCLUDE gtest/gtest.h $ENV{CMAKE_PREFIX_PATH})

include_directories (${GAEDA_INCLUDE} ${GOOGLE_TEST_INCLUDE})

link_directories (${CMAKE_CURRENT_SOURCE_DIR})

set (TESTS_SRCS findposofvalues_unittests insertions_tests swap_unittests chain_unittests zero_unittests ls_unittests init_unittests regretinit_unittests slackinit_unittests checknatseqs_unittests greedymove_unittests greedyswap_unittests darpactiontests tsCordeauTests darpgenserializationtests)

foreach(TEST ${TESTS_SRCS})
  add_executable (${TEST} ${TEST} ../aux.cc)
  target_link_libraries (${TEST} ${GAEDA_LIBRARY} ${GOOGLE_TEST_LIB})
  #target_link_libraries(${TEST} 'gtest_main')
endforeach (TEST)

add_subdirectory(routingdist)

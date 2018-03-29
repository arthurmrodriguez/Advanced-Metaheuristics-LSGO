find_path(SPIRIT_INCLUDE "boost/spirit.hpp")
if ("${SPIRIT_INCLUDE}" STREQUAL "SPIRIT_INCLUDE-NOTFOUND")
  MESSAGE(FATAL_ERROR "Couldn't find Boost SPIRIT header files.")
endif()

# Used in the genealogy
include_directories (${SPIRIT_INCLUDE})

set (UTILS_SRCS genealogy)

foreach (UTIL ${UTILS_SRCS})
  add_executable (${UTIL} ${UTIL}.cc)
endforeach (UTIL)

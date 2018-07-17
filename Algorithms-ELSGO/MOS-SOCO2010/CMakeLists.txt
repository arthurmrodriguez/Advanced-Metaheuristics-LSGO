project (GAEDALib)

# We require at least CMake V2.6
cmake_minimum_required (VERSION 2.6)

#set (CMAKE_VERBOSE_MAKEFILE true)
cmake_policy(SET CMP0042 NEW)
SET(CMAKE_MACOSX_RPATH OFF)
#execute_process (COMMAND ctags-exuberant --recurse ./*)

add_subdirectory (gaeda)        # generate libgaeda
add_subdirectory (gaedaexec)    # generate gaedaexec
#add_subdirectory (cec2012)      # generate problems
#add_subdirectory (cec2013_lsgo) # generate problems
add_subdirectory (soco2010)     # generate problems

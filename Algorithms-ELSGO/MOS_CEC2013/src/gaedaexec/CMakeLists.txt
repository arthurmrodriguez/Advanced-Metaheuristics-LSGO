include_directories(${LIBCONFIG_INCLUDE} ${GAEDA_INCLUDE})

set (GAEDAEXEC_SRCS main.cc)
add_executable (gaedaexec ${GAEDAEXEC_SRCS})

target_link_libraries (gaedaexec ${GAEDA_LIBRARY} ${LIBCONFIG} ${Boost_PROGRAM_OPTIONS_LIBRARY})

# Needed in order to be able to dlopen 'libvesicles.so', as libgomp.so can
# not be dlopened (at least in GCC 4.2)
# See: http://gcc.gnu.org/ml/gcc-help/2007-09/msg00050.html
# target_link_libraries (gaedaexec 'gomp')

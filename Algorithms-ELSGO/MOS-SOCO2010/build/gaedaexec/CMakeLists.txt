project (gaedaexec)

find_package (MPI REQUIRED)
find_package (Boost REQUIRED COMPONENTS program_options)

include (CMakeForceCompiler)
cmake_force_cxx_compiler (mpicxx "MPI C++ Compiler")
cmake_force_c_compiler (mpicc "MPI C Compiler")

include_directories (../gaeda)
link_directories (../gaeda)

set (GAEDAEXEC_SRCS main.cc)
add_executable (gaedaexec ${GAEDAEXEC_SRCS})

target_link_libraries (gaedaexec 'gaeda')
target_link_libraries (gaedaexec ${Boost_LIBRARIES})

# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/Arthur18/Documents/UGR-GII/4toCurso/TFG-AdvancedMetaheuristicsForBigOpt/UGR-GII-TFG/Algorithms-ELSGO/MOS-SOCO2010

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/Arthur18/Documents/UGR-GII/4toCurso/TFG-AdvancedMetaheuristicsForBigOpt/UGR-GII-TFG/Algorithms-ELSGO/MOS-SOCO2010

# Include any dependencies generated for this target.
include soco2010/CMakeFiles/eeg_problem.dir/depend.make

# Include the progress variables for this target.
include soco2010/CMakeFiles/eeg_problem.dir/progress.make

# Include the compile flags for this target's objects.
include soco2010/CMakeFiles/eeg_problem.dir/flags.make

soco2010/CMakeFiles/eeg_problem.dir/eeg_problem.cc.o: soco2010/CMakeFiles/eeg_problem.dir/flags.make
soco2010/CMakeFiles/eeg_problem.dir/eeg_problem.cc.o: soco2010/eeg_problem.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/Arthur18/Documents/UGR-GII/4toCurso/TFG-AdvancedMetaheuristicsForBigOpt/UGR-GII-TFG/Algorithms-ELSGO/MOS-SOCO2010/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object soco2010/CMakeFiles/eeg_problem.dir/eeg_problem.cc.o"
	cd /Users/Arthur18/Documents/UGR-GII/4toCurso/TFG-AdvancedMetaheuristicsForBigOpt/UGR-GII-TFG/Algorithms-ELSGO/MOS-SOCO2010/soco2010 && /Library/Developer/CommandLineTools/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/eeg_problem.dir/eeg_problem.cc.o -c /Users/Arthur18/Documents/UGR-GII/4toCurso/TFG-AdvancedMetaheuristicsForBigOpt/UGR-GII-TFG/Algorithms-ELSGO/MOS-SOCO2010/soco2010/eeg_problem.cc

soco2010/CMakeFiles/eeg_problem.dir/eeg_problem.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/eeg_problem.dir/eeg_problem.cc.i"
	cd /Users/Arthur18/Documents/UGR-GII/4toCurso/TFG-AdvancedMetaheuristicsForBigOpt/UGR-GII-TFG/Algorithms-ELSGO/MOS-SOCO2010/soco2010 && /Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/Arthur18/Documents/UGR-GII/4toCurso/TFG-AdvancedMetaheuristicsForBigOpt/UGR-GII-TFG/Algorithms-ELSGO/MOS-SOCO2010/soco2010/eeg_problem.cc > CMakeFiles/eeg_problem.dir/eeg_problem.cc.i

soco2010/CMakeFiles/eeg_problem.dir/eeg_problem.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/eeg_problem.dir/eeg_problem.cc.s"
	cd /Users/Arthur18/Documents/UGR-GII/4toCurso/TFG-AdvancedMetaheuristicsForBigOpt/UGR-GII-TFG/Algorithms-ELSGO/MOS-SOCO2010/soco2010 && /Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/Arthur18/Documents/UGR-GII/4toCurso/TFG-AdvancedMetaheuristicsForBigOpt/UGR-GII-TFG/Algorithms-ELSGO/MOS-SOCO2010/soco2010/eeg_problem.cc -o CMakeFiles/eeg_problem.dir/eeg_problem.cc.s

soco2010/CMakeFiles/eeg_problem.dir/eeg_problem.cc.o.requires:

.PHONY : soco2010/CMakeFiles/eeg_problem.dir/eeg_problem.cc.o.requires

soco2010/CMakeFiles/eeg_problem.dir/eeg_problem.cc.o.provides: soco2010/CMakeFiles/eeg_problem.dir/eeg_problem.cc.o.requires
	$(MAKE) -f soco2010/CMakeFiles/eeg_problem.dir/build.make soco2010/CMakeFiles/eeg_problem.dir/eeg_problem.cc.o.provides.build
.PHONY : soco2010/CMakeFiles/eeg_problem.dir/eeg_problem.cc.o.provides

soco2010/CMakeFiles/eeg_problem.dir/eeg_problem.cc.o.provides.build: soco2010/CMakeFiles/eeg_problem.dir/eeg_problem.cc.o


# Object files for target eeg_problem
eeg_problem_OBJECTS = \
"CMakeFiles/eeg_problem.dir/eeg_problem.cc.o"

# External object files for target eeg_problem
eeg_problem_EXTERNAL_OBJECTS =

soco2010/libeeg_problem.so: soco2010/CMakeFiles/eeg_problem.dir/eeg_problem.cc.o
soco2010/libeeg_problem.so: soco2010/CMakeFiles/eeg_problem.dir/build.make
soco2010/libeeg_problem.so: soco2010/CMakeFiles/eeg_problem.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/Arthur18/Documents/UGR-GII/4toCurso/TFG-AdvancedMetaheuristicsForBigOpt/UGR-GII-TFG/Algorithms-ELSGO/MOS-SOCO2010/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library libeeg_problem.so"
	cd /Users/Arthur18/Documents/UGR-GII/4toCurso/TFG-AdvancedMetaheuristicsForBigOpt/UGR-GII-TFG/Algorithms-ELSGO/MOS-SOCO2010/soco2010 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/eeg_problem.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
soco2010/CMakeFiles/eeg_problem.dir/build: soco2010/libeeg_problem.so

.PHONY : soco2010/CMakeFiles/eeg_problem.dir/build

soco2010/CMakeFiles/eeg_problem.dir/requires: soco2010/CMakeFiles/eeg_problem.dir/eeg_problem.cc.o.requires

.PHONY : soco2010/CMakeFiles/eeg_problem.dir/requires

soco2010/CMakeFiles/eeg_problem.dir/clean:
	cd /Users/Arthur18/Documents/UGR-GII/4toCurso/TFG-AdvancedMetaheuristicsForBigOpt/UGR-GII-TFG/Algorithms-ELSGO/MOS-SOCO2010/soco2010 && $(CMAKE_COMMAND) -P CMakeFiles/eeg_problem.dir/cmake_clean.cmake
.PHONY : soco2010/CMakeFiles/eeg_problem.dir/clean

soco2010/CMakeFiles/eeg_problem.dir/depend:
	cd /Users/Arthur18/Documents/UGR-GII/4toCurso/TFG-AdvancedMetaheuristicsForBigOpt/UGR-GII-TFG/Algorithms-ELSGO/MOS-SOCO2010 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/Arthur18/Documents/UGR-GII/4toCurso/TFG-AdvancedMetaheuristicsForBigOpt/UGR-GII-TFG/Algorithms-ELSGO/MOS-SOCO2010 /Users/Arthur18/Documents/UGR-GII/4toCurso/TFG-AdvancedMetaheuristicsForBigOpt/UGR-GII-TFG/Algorithms-ELSGO/MOS-SOCO2010/soco2010 /Users/Arthur18/Documents/UGR-GII/4toCurso/TFG-AdvancedMetaheuristicsForBigOpt/UGR-GII-TFG/Algorithms-ELSGO/MOS-SOCO2010 /Users/Arthur18/Documents/UGR-GII/4toCurso/TFG-AdvancedMetaheuristicsForBigOpt/UGR-GII-TFG/Algorithms-ELSGO/MOS-SOCO2010/soco2010 /Users/Arthur18/Documents/UGR-GII/4toCurso/TFG-AdvancedMetaheuristicsForBigOpt/UGR-GII-TFG/Algorithms-ELSGO/MOS-SOCO2010/soco2010/CMakeFiles/eeg_problem.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : soco2010/CMakeFiles/eeg_problem.dir/depend


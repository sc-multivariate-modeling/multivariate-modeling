# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.11

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
CMAKE_COMMAND = /opt/ecrc/cmake/3.11.1/ub16/bin/cmake

# The command to remove a file.
RM = /opt/ecrc/cmake/3.11.1/ub16/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/abdullsm/develop/stars-h

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/abdullsm/develop/stars-h/build

# Include any dependencies generated for this target.
include examples/problem/CMakeFiles/example_problem_dense.dir/depend.make

# Include the progress variables for this target.
include examples/problem/CMakeFiles/example_problem_dense.dir/progress.make

# Include the compile flags for this target's objects.
include examples/problem/CMakeFiles/example_problem_dense.dir/flags.make

examples/problem/CMakeFiles/example_problem_dense.dir/dense.c.o: examples/problem/CMakeFiles/example_problem_dense.dir/flags.make
examples/problem/CMakeFiles/example_problem_dense.dir/dense.c.o: ../examples/problem/dense.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object examples/problem/CMakeFiles/example_problem_dense.dir/dense.c.o"
	cd /home/abdullsm/develop/stars-h/build/examples/problem && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/example_problem_dense.dir/dense.c.o   -c /home/abdullsm/develop/stars-h/examples/problem/dense.c

examples/problem/CMakeFiles/example_problem_dense.dir/dense.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/example_problem_dense.dir/dense.c.i"
	cd /home/abdullsm/develop/stars-h/build/examples/problem && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/abdullsm/develop/stars-h/examples/problem/dense.c > CMakeFiles/example_problem_dense.dir/dense.c.i

examples/problem/CMakeFiles/example_problem_dense.dir/dense.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/example_problem_dense.dir/dense.c.s"
	cd /home/abdullsm/develop/stars-h/build/examples/problem && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/abdullsm/develop/stars-h/examples/problem/dense.c -o CMakeFiles/example_problem_dense.dir/dense.c.s

# Object files for target example_problem_dense
example_problem_dense_OBJECTS = \
"CMakeFiles/example_problem_dense.dir/dense.c.o"

# External object files for target example_problem_dense
example_problem_dense_EXTERNAL_OBJECTS =

examples/problem/dense: examples/problem/CMakeFiles/example_problem_dense.dir/dense.c.o
examples/problem/dense: examples/problem/CMakeFiles/example_problem_dense.dir/build.make
examples/problem/dense: src/libstarsh.a
examples/problem/dense: /opt/ecrc/mkl/2018-update-1/mkl/lib/intel64/libmkl_intel_lp64.so
examples/problem/dense: /opt/ecrc/mkl/2018-update-1/mkl/lib/intel64/libmkl_sequential.so
examples/problem/dense: /opt/ecrc/mkl/2018-update-1/mkl/lib/intel64/libmkl_core.so
examples/problem/dense: /opt/ecrc/mkl/2018-update-1/mkl/lib/intel64/libmkl_intel_lp64.so
examples/problem/dense: /opt/ecrc/mkl/2018-update-1/mkl/lib/intel64/libmkl_sequential.so
examples/problem/dense: /opt/ecrc/mkl/2018-update-1/mkl/lib/intel64/libmkl_core.so
examples/problem/dense: /opt/ecrc/gsl/2.4-gcc-7.2.0/ub16/lib/libgsl.so
examples/problem/dense: /opt/ecrc/gsl/2.4-gcc-7.2.0/ub16/lib/libgslcblas.so
examples/problem/dense: examples/problem/CMakeFiles/example_problem_dense.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable dense"
	cd /home/abdullsm/develop/stars-h/build/examples/problem && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/example_problem_dense.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/problem/CMakeFiles/example_problem_dense.dir/build: examples/problem/dense

.PHONY : examples/problem/CMakeFiles/example_problem_dense.dir/build

examples/problem/CMakeFiles/example_problem_dense.dir/clean:
	cd /home/abdullsm/develop/stars-h/build/examples/problem && $(CMAKE_COMMAND) -P CMakeFiles/example_problem_dense.dir/cmake_clean.cmake
.PHONY : examples/problem/CMakeFiles/example_problem_dense.dir/clean

examples/problem/CMakeFiles/example_problem_dense.dir/depend:
	cd /home/abdullsm/develop/stars-h/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/abdullsm/develop/stars-h /home/abdullsm/develop/stars-h/examples/problem /home/abdullsm/develop/stars-h/build /home/abdullsm/develop/stars-h/build/examples/problem /home/abdullsm/develop/stars-h/build/examples/problem/CMakeFiles/example_problem_dense.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/problem/CMakeFiles/example_problem_dense.dir/depend

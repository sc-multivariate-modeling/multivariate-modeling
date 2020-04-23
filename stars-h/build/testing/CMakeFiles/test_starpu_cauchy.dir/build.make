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
include testing/CMakeFiles/test_starpu_cauchy.dir/depend.make

# Include the progress variables for this target.
include testing/CMakeFiles/test_starpu_cauchy.dir/progress.make

# Include the compile flags for this target's objects.
include testing/CMakeFiles/test_starpu_cauchy.dir/flags.make

testing/CMakeFiles/test_starpu_cauchy.dir/starpu_cauchy.c.o: testing/CMakeFiles/test_starpu_cauchy.dir/flags.make
testing/CMakeFiles/test_starpu_cauchy.dir/starpu_cauchy.c.o: ../testing/starpu_cauchy.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object testing/CMakeFiles/test_starpu_cauchy.dir/starpu_cauchy.c.o"
	cd /home/abdullsm/develop/stars-h/build/testing && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/test_starpu_cauchy.dir/starpu_cauchy.c.o   -c /home/abdullsm/develop/stars-h/testing/starpu_cauchy.c

testing/CMakeFiles/test_starpu_cauchy.dir/starpu_cauchy.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/test_starpu_cauchy.dir/starpu_cauchy.c.i"
	cd /home/abdullsm/develop/stars-h/build/testing && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/abdullsm/develop/stars-h/testing/starpu_cauchy.c > CMakeFiles/test_starpu_cauchy.dir/starpu_cauchy.c.i

testing/CMakeFiles/test_starpu_cauchy.dir/starpu_cauchy.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/test_starpu_cauchy.dir/starpu_cauchy.c.s"
	cd /home/abdullsm/develop/stars-h/build/testing && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/abdullsm/develop/stars-h/testing/starpu_cauchy.c -o CMakeFiles/test_starpu_cauchy.dir/starpu_cauchy.c.s

# Object files for target test_starpu_cauchy
test_starpu_cauchy_OBJECTS = \
"CMakeFiles/test_starpu_cauchy.dir/starpu_cauchy.c.o"

# External object files for target test_starpu_cauchy
test_starpu_cauchy_EXTERNAL_OBJECTS =

testing/starpu_cauchy: testing/CMakeFiles/test_starpu_cauchy.dir/starpu_cauchy.c.o
testing/starpu_cauchy: testing/CMakeFiles/test_starpu_cauchy.dir/build.make
testing/starpu_cauchy: src/libstarsh.a
testing/starpu_cauchy: /opt/ecrc/mkl/2018-update-1/mkl/lib/intel64/libmkl_intel_lp64.so
testing/starpu_cauchy: /opt/ecrc/mkl/2018-update-1/mkl/lib/intel64/libmkl_sequential.so
testing/starpu_cauchy: /opt/ecrc/mkl/2018-update-1/mkl/lib/intel64/libmkl_core.so
testing/starpu_cauchy: /opt/ecrc/mkl/2018-update-1/mkl/lib/intel64/libmkl_intel_lp64.so
testing/starpu_cauchy: /opt/ecrc/mkl/2018-update-1/mkl/lib/intel64/libmkl_sequential.so
testing/starpu_cauchy: /opt/ecrc/mkl/2018-update-1/mkl/lib/intel64/libmkl_core.so
testing/starpu_cauchy: /opt/ecrc/gsl/2.4-gcc-7.2.0/ub16/lib/libgsl.so
testing/starpu_cauchy: /opt/ecrc/gsl/2.4-gcc-7.2.0/ub16/lib/libgslcblas.so
testing/starpu_cauchy: testing/CMakeFiles/test_starpu_cauchy.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable starpu_cauchy"
	cd /home/abdullsm/develop/stars-h/build/testing && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_starpu_cauchy.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
testing/CMakeFiles/test_starpu_cauchy.dir/build: testing/starpu_cauchy

.PHONY : testing/CMakeFiles/test_starpu_cauchy.dir/build

testing/CMakeFiles/test_starpu_cauchy.dir/clean:
	cd /home/abdullsm/develop/stars-h/build/testing && $(CMAKE_COMMAND) -P CMakeFiles/test_starpu_cauchy.dir/cmake_clean.cmake
.PHONY : testing/CMakeFiles/test_starpu_cauchy.dir/clean

testing/CMakeFiles/test_starpu_cauchy.dir/depend:
	cd /home/abdullsm/develop/stars-h/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/abdullsm/develop/stars-h /home/abdullsm/develop/stars-h/testing /home/abdullsm/develop/stars-h/build /home/abdullsm/develop/stars-h/build/testing /home/abdullsm/develop/stars-h/build/testing/CMakeFiles/test_starpu_cauchy.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : testing/CMakeFiles/test_starpu_cauchy.dir/depend


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
include examples/approximation/CMakeFiles/example_approximation_minimal.dir/depend.make

# Include the progress variables for this target.
include examples/approximation/CMakeFiles/example_approximation_minimal.dir/progress.make

# Include the compile flags for this target's objects.
include examples/approximation/CMakeFiles/example_approximation_minimal.dir/flags.make

examples/approximation/CMakeFiles/example_approximation_minimal.dir/minimal.c.o: examples/approximation/CMakeFiles/example_approximation_minimal.dir/flags.make
examples/approximation/CMakeFiles/example_approximation_minimal.dir/minimal.c.o: ../examples/approximation/minimal.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object examples/approximation/CMakeFiles/example_approximation_minimal.dir/minimal.c.o"
	cd /home/abdullsm/develop/stars-h/build/examples/approximation && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/example_approximation_minimal.dir/minimal.c.o   -c /home/abdullsm/develop/stars-h/examples/approximation/minimal.c

examples/approximation/CMakeFiles/example_approximation_minimal.dir/minimal.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/example_approximation_minimal.dir/minimal.c.i"
	cd /home/abdullsm/develop/stars-h/build/examples/approximation && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/abdullsm/develop/stars-h/examples/approximation/minimal.c > CMakeFiles/example_approximation_minimal.dir/minimal.c.i

examples/approximation/CMakeFiles/example_approximation_minimal.dir/minimal.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/example_approximation_minimal.dir/minimal.c.s"
	cd /home/abdullsm/develop/stars-h/build/examples/approximation && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/abdullsm/develop/stars-h/examples/approximation/minimal.c -o CMakeFiles/example_approximation_minimal.dir/minimal.c.s

# Object files for target example_approximation_minimal
example_approximation_minimal_OBJECTS = \
"CMakeFiles/example_approximation_minimal.dir/minimal.c.o"

# External object files for target example_approximation_minimal
example_approximation_minimal_EXTERNAL_OBJECTS =

examples/approximation/minimal: examples/approximation/CMakeFiles/example_approximation_minimal.dir/minimal.c.o
examples/approximation/minimal: examples/approximation/CMakeFiles/example_approximation_minimal.dir/build.make
examples/approximation/minimal: src/libstarsh.a
examples/approximation/minimal: /opt/ecrc/mkl/2018-update-1/mkl/lib/intel64/libmkl_intel_lp64.so
examples/approximation/minimal: /opt/ecrc/mkl/2018-update-1/mkl/lib/intel64/libmkl_sequential.so
examples/approximation/minimal: /opt/ecrc/mkl/2018-update-1/mkl/lib/intel64/libmkl_core.so
examples/approximation/minimal: /opt/ecrc/mkl/2018-update-1/mkl/lib/intel64/libmkl_intel_lp64.so
examples/approximation/minimal: /opt/ecrc/mkl/2018-update-1/mkl/lib/intel64/libmkl_sequential.so
examples/approximation/minimal: /opt/ecrc/mkl/2018-update-1/mkl/lib/intel64/libmkl_core.so
examples/approximation/minimal: /opt/ecrc/gsl/2.4-gcc-7.2.0/ub16/lib/libgsl.so
examples/approximation/minimal: /opt/ecrc/gsl/2.4-gcc-7.2.0/ub16/lib/libgslcblas.so
examples/approximation/minimal: examples/approximation/CMakeFiles/example_approximation_minimal.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable minimal"
	cd /home/abdullsm/develop/stars-h/build/examples/approximation && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/example_approximation_minimal.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/approximation/CMakeFiles/example_approximation_minimal.dir/build: examples/approximation/minimal

.PHONY : examples/approximation/CMakeFiles/example_approximation_minimal.dir/build

examples/approximation/CMakeFiles/example_approximation_minimal.dir/clean:
	cd /home/abdullsm/develop/stars-h/build/examples/approximation && $(CMAKE_COMMAND) -P CMakeFiles/example_approximation_minimal.dir/cmake_clean.cmake
.PHONY : examples/approximation/CMakeFiles/example_approximation_minimal.dir/clean

examples/approximation/CMakeFiles/example_approximation_minimal.dir/depend:
	cd /home/abdullsm/develop/stars-h/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/abdullsm/develop/stars-h /home/abdullsm/develop/stars-h/examples/approximation /home/abdullsm/develop/stars-h/build /home/abdullsm/develop/stars-h/build/examples/approximation /home/abdullsm/develop/stars-h/build/examples/approximation/CMakeFiles/example_approximation_minimal.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/approximation/CMakeFiles/example_approximation_minimal.dir/depend


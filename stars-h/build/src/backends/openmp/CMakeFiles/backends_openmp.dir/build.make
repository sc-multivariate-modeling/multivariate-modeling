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
include src/backends/openmp/CMakeFiles/backends_openmp.dir/depend.make

# Include the progress variables for this target.
include src/backends/openmp/CMakeFiles/backends_openmp.dir/progress.make

# Include the compile flags for this target's objects.
include src/backends/openmp/CMakeFiles/backends_openmp.dir/flags.make

src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dqp3.c.o: src/backends/openmp/CMakeFiles/backends_openmp.dir/flags.make
src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dqp3.c.o: ../src/backends/openmp/blrm/dqp3.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dqp3.c.o"
	cd /home/abdullsm/develop/stars-h/build/src/backends/openmp && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/backends_openmp.dir/blrm/dqp3.c.o   -c /home/abdullsm/develop/stars-h/src/backends/openmp/blrm/dqp3.c

src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dqp3.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/backends_openmp.dir/blrm/dqp3.c.i"
	cd /home/abdullsm/develop/stars-h/build/src/backends/openmp && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/abdullsm/develop/stars-h/src/backends/openmp/blrm/dqp3.c > CMakeFiles/backends_openmp.dir/blrm/dqp3.c.i

src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dqp3.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/backends_openmp.dir/blrm/dqp3.c.s"
	cd /home/abdullsm/develop/stars-h/build/src/backends/openmp && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/abdullsm/develop/stars-h/src/backends/openmp/blrm/dqp3.c -o CMakeFiles/backends_openmp.dir/blrm/dqp3.c.s

src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/drsdd.c.o: src/backends/openmp/CMakeFiles/backends_openmp.dir/flags.make
src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/drsdd.c.o: ../src/backends/openmp/blrm/drsdd.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/drsdd.c.o"
	cd /home/abdullsm/develop/stars-h/build/src/backends/openmp && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/backends_openmp.dir/blrm/drsdd.c.o   -c /home/abdullsm/develop/stars-h/src/backends/openmp/blrm/drsdd.c

src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/drsdd.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/backends_openmp.dir/blrm/drsdd.c.i"
	cd /home/abdullsm/develop/stars-h/build/src/backends/openmp && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/abdullsm/develop/stars-h/src/backends/openmp/blrm/drsdd.c > CMakeFiles/backends_openmp.dir/blrm/drsdd.c.i

src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/drsdd.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/backends_openmp.dir/blrm/drsdd.c.s"
	cd /home/abdullsm/develop/stars-h/build/src/backends/openmp && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/abdullsm/develop/stars-h/src/backends/openmp/blrm/drsdd.c -o CMakeFiles/backends_openmp.dir/blrm/drsdd.c.s

src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dsdd.c.o: src/backends/openmp/CMakeFiles/backends_openmp.dir/flags.make
src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dsdd.c.o: ../src/backends/openmp/blrm/dsdd.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dsdd.c.o"
	cd /home/abdullsm/develop/stars-h/build/src/backends/openmp && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/backends_openmp.dir/blrm/dsdd.c.o   -c /home/abdullsm/develop/stars-h/src/backends/openmp/blrm/dsdd.c

src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dsdd.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/backends_openmp.dir/blrm/dsdd.c.i"
	cd /home/abdullsm/develop/stars-h/build/src/backends/openmp && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/abdullsm/develop/stars-h/src/backends/openmp/blrm/dsdd.c > CMakeFiles/backends_openmp.dir/blrm/dsdd.c.i

src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dsdd.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/backends_openmp.dir/blrm/dsdd.c.s"
	cd /home/abdullsm/develop/stars-h/build/src/backends/openmp && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/abdullsm/develop/stars-h/src/backends/openmp/blrm/dsdd.c -o CMakeFiles/backends_openmp.dir/blrm/dsdd.c.s

src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dmml.c.o: src/backends/openmp/CMakeFiles/backends_openmp.dir/flags.make
src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dmml.c.o: ../src/backends/openmp/blrm/dmml.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dmml.c.o"
	cd /home/abdullsm/develop/stars-h/build/src/backends/openmp && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/backends_openmp.dir/blrm/dmml.c.o   -c /home/abdullsm/develop/stars-h/src/backends/openmp/blrm/dmml.c

src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dmml.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/backends_openmp.dir/blrm/dmml.c.i"
	cd /home/abdullsm/develop/stars-h/build/src/backends/openmp && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/abdullsm/develop/stars-h/src/backends/openmp/blrm/dmml.c > CMakeFiles/backends_openmp.dir/blrm/dmml.c.i

src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dmml.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/backends_openmp.dir/blrm/dmml.c.s"
	cd /home/abdullsm/develop/stars-h/build/src/backends/openmp && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/abdullsm/develop/stars-h/src/backends/openmp/blrm/dmml.c -o CMakeFiles/backends_openmp.dir/blrm/dmml.c.s

src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dfe.c.o: src/backends/openmp/CMakeFiles/backends_openmp.dir/flags.make
src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dfe.c.o: ../src/backends/openmp/blrm/dfe.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dfe.c.o"
	cd /home/abdullsm/develop/stars-h/build/src/backends/openmp && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/backends_openmp.dir/blrm/dfe.c.o   -c /home/abdullsm/develop/stars-h/src/backends/openmp/blrm/dfe.c

src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dfe.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/backends_openmp.dir/blrm/dfe.c.i"
	cd /home/abdullsm/develop/stars-h/build/src/backends/openmp && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/abdullsm/develop/stars-h/src/backends/openmp/blrm/dfe.c > CMakeFiles/backends_openmp.dir/blrm/dfe.c.i

src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dfe.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/backends_openmp.dir/blrm/dfe.c.s"
	cd /home/abdullsm/develop/stars-h/build/src/backends/openmp && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/abdullsm/develop/stars-h/src/backends/openmp/blrm/dfe.c -o CMakeFiles/backends_openmp.dir/blrm/dfe.c.s

backends_openmp: src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dqp3.c.o
backends_openmp: src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/drsdd.c.o
backends_openmp: src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dsdd.c.o
backends_openmp: src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dmml.c.o
backends_openmp: src/backends/openmp/CMakeFiles/backends_openmp.dir/blrm/dfe.c.o
backends_openmp: src/backends/openmp/CMakeFiles/backends_openmp.dir/build.make

.PHONY : backends_openmp

# Rule to build all files generated by this target.
src/backends/openmp/CMakeFiles/backends_openmp.dir/build: backends_openmp

.PHONY : src/backends/openmp/CMakeFiles/backends_openmp.dir/build

src/backends/openmp/CMakeFiles/backends_openmp.dir/clean:
	cd /home/abdullsm/develop/stars-h/build/src/backends/openmp && $(CMAKE_COMMAND) -P CMakeFiles/backends_openmp.dir/cmake_clean.cmake
.PHONY : src/backends/openmp/CMakeFiles/backends_openmp.dir/clean

src/backends/openmp/CMakeFiles/backends_openmp.dir/depend:
	cd /home/abdullsm/develop/stars-h/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/abdullsm/develop/stars-h /home/abdullsm/develop/stars-h/src/backends/openmp /home/abdullsm/develop/stars-h/build /home/abdullsm/develop/stars-h/build/src/backends/openmp /home/abdullsm/develop/stars-h/build/src/backends/openmp/CMakeFiles/backends_openmp.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/backends/openmp/CMakeFiles/backends_openmp.dir/depend

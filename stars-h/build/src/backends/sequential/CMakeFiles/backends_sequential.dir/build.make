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
include src/backends/sequential/CMakeFiles/backends_sequential.dir/depend.make

# Include the progress variables for this target.
include src/backends/sequential/CMakeFiles/backends_sequential.dir/progress.make

# Include the compile flags for this target's objects.
include src/backends/sequential/CMakeFiles/backends_sequential.dir/flags.make

src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dqp3.c.o: src/backends/sequential/CMakeFiles/backends_sequential.dir/flags.make
src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dqp3.c.o: ../src/backends/sequential/dense/dqp3.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dqp3.c.o"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/backends_sequential.dir/dense/dqp3.c.o   -c /home/abdullsm/develop/stars-h/src/backends/sequential/dense/dqp3.c

src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dqp3.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/backends_sequential.dir/dense/dqp3.c.i"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/abdullsm/develop/stars-h/src/backends/sequential/dense/dqp3.c > CMakeFiles/backends_sequential.dir/dense/dqp3.c.i

src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dqp3.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/backends_sequential.dir/dense/dqp3.c.s"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/abdullsm/develop/stars-h/src/backends/sequential/dense/dqp3.c -o CMakeFiles/backends_sequential.dir/dense/dqp3.c.s

src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/drsdd.c.o: src/backends/sequential/CMakeFiles/backends_sequential.dir/flags.make
src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/drsdd.c.o: ../src/backends/sequential/dense/drsdd.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/drsdd.c.o"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/backends_sequential.dir/dense/drsdd.c.o   -c /home/abdullsm/develop/stars-h/src/backends/sequential/dense/drsdd.c

src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/drsdd.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/backends_sequential.dir/dense/drsdd.c.i"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/abdullsm/develop/stars-h/src/backends/sequential/dense/drsdd.c > CMakeFiles/backends_sequential.dir/dense/drsdd.c.i

src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/drsdd.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/backends_sequential.dir/dense/drsdd.c.s"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/abdullsm/develop/stars-h/src/backends/sequential/dense/drsdd.c -o CMakeFiles/backends_sequential.dir/dense/drsdd.c.s

src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dsdd.c.o: src/backends/sequential/CMakeFiles/backends_sequential.dir/flags.make
src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dsdd.c.o: ../src/backends/sequential/dense/dsdd.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dsdd.c.o"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/backends_sequential.dir/dense/dsdd.c.o   -c /home/abdullsm/develop/stars-h/src/backends/sequential/dense/dsdd.c

src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dsdd.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/backends_sequential.dir/dense/dsdd.c.i"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/abdullsm/develop/stars-h/src/backends/sequential/dense/dsdd.c > CMakeFiles/backends_sequential.dir/dense/dsdd.c.i

src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dsdd.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/backends_sequential.dir/dense/dsdd.c.s"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/abdullsm/develop/stars-h/src/backends/sequential/dense/dsdd.c -o CMakeFiles/backends_sequential.dir/dense/dsdd.c.s

src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dsvfr.c.o: src/backends/sequential/CMakeFiles/backends_sequential.dir/flags.make
src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dsvfr.c.o: ../src/backends/sequential/dense/dsvfr.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dsvfr.c.o"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/backends_sequential.dir/dense/dsvfr.c.o   -c /home/abdullsm/develop/stars-h/src/backends/sequential/dense/dsvfr.c

src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dsvfr.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/backends_sequential.dir/dense/dsvfr.c.i"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/abdullsm/develop/stars-h/src/backends/sequential/dense/dsvfr.c > CMakeFiles/backends_sequential.dir/dense/dsvfr.c.i

src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dsvfr.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/backends_sequential.dir/dense/dsvfr.c.s"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/abdullsm/develop/stars-h/src/backends/sequential/dense/dsvfr.c -o CMakeFiles/backends_sequential.dir/dense/dsvfr.c.s

src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dna.c.o: src/backends/sequential/CMakeFiles/backends_sequential.dir/flags.make
src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dna.c.o: ../src/backends/sequential/dense/dna.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dna.c.o"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/backends_sequential.dir/dense/dna.c.o   -c /home/abdullsm/develop/stars-h/src/backends/sequential/dense/dna.c

src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dna.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/backends_sequential.dir/dense/dna.c.i"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/abdullsm/develop/stars-h/src/backends/sequential/dense/dna.c > CMakeFiles/backends_sequential.dir/dense/dna.c.i

src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dna.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/backends_sequential.dir/dense/dna.c.s"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/abdullsm/develop/stars-h/src/backends/sequential/dense/dna.c -o CMakeFiles/backends_sequential.dir/dense/dna.c.s

src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dca.c.o: src/backends/sequential/CMakeFiles/backends_sequential.dir/flags.make
src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dca.c.o: ../src/backends/sequential/blrm/dca.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dca.c.o"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/backends_sequential.dir/blrm/dca.c.o   -c /home/abdullsm/develop/stars-h/src/backends/sequential/blrm/dca.c

src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dca.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/backends_sequential.dir/blrm/dca.c.i"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/abdullsm/develop/stars-h/src/backends/sequential/blrm/dca.c > CMakeFiles/backends_sequential.dir/blrm/dca.c.i

src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dca.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/backends_sequential.dir/blrm/dca.c.s"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/abdullsm/develop/stars-h/src/backends/sequential/blrm/dca.c -o CMakeFiles/backends_sequential.dir/blrm/dca.c.s

src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dfe.c.o: src/backends/sequential/CMakeFiles/backends_sequential.dir/flags.make
src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dfe.c.o: ../src/backends/sequential/blrm/dfe.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building C object src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dfe.c.o"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/backends_sequential.dir/blrm/dfe.c.o   -c /home/abdullsm/develop/stars-h/src/backends/sequential/blrm/dfe.c

src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dfe.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/backends_sequential.dir/blrm/dfe.c.i"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/abdullsm/develop/stars-h/src/backends/sequential/blrm/dfe.c > CMakeFiles/backends_sequential.dir/blrm/dfe.c.i

src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dfe.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/backends_sequential.dir/blrm/dfe.c.s"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/abdullsm/develop/stars-h/src/backends/sequential/blrm/dfe.c -o CMakeFiles/backends_sequential.dir/blrm/dfe.c.s

src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dmml.c.o: src/backends/sequential/CMakeFiles/backends_sequential.dir/flags.make
src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dmml.c.o: ../src/backends/sequential/blrm/dmml.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building C object src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dmml.c.o"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/backends_sequential.dir/blrm/dmml.c.o   -c /home/abdullsm/develop/stars-h/src/backends/sequential/blrm/dmml.c

src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dmml.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/backends_sequential.dir/blrm/dmml.c.i"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/abdullsm/develop/stars-h/src/backends/sequential/blrm/dmml.c > CMakeFiles/backends_sequential.dir/blrm/dmml.c.i

src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dmml.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/backends_sequential.dir/blrm/dmml.c.s"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/abdullsm/develop/stars-h/src/backends/sequential/blrm/dmml.c -o CMakeFiles/backends_sequential.dir/blrm/dmml.c.s

src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dqp3.c.o: src/backends/sequential/CMakeFiles/backends_sequential.dir/flags.make
src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dqp3.c.o: ../src/backends/sequential/blrm/dqp3.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building C object src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dqp3.c.o"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/backends_sequential.dir/blrm/dqp3.c.o   -c /home/abdullsm/develop/stars-h/src/backends/sequential/blrm/dqp3.c

src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dqp3.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/backends_sequential.dir/blrm/dqp3.c.i"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/abdullsm/develop/stars-h/src/backends/sequential/blrm/dqp3.c > CMakeFiles/backends_sequential.dir/blrm/dqp3.c.i

src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dqp3.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/backends_sequential.dir/blrm/dqp3.c.s"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/abdullsm/develop/stars-h/src/backends/sequential/blrm/dqp3.c -o CMakeFiles/backends_sequential.dir/blrm/dqp3.c.s

src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/drsdd.c.o: src/backends/sequential/CMakeFiles/backends_sequential.dir/flags.make
src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/drsdd.c.o: ../src/backends/sequential/blrm/drsdd.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building C object src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/drsdd.c.o"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/backends_sequential.dir/blrm/drsdd.c.o   -c /home/abdullsm/develop/stars-h/src/backends/sequential/blrm/drsdd.c

src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/drsdd.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/backends_sequential.dir/blrm/drsdd.c.i"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/abdullsm/develop/stars-h/src/backends/sequential/blrm/drsdd.c > CMakeFiles/backends_sequential.dir/blrm/drsdd.c.i

src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/drsdd.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/backends_sequential.dir/blrm/drsdd.c.s"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/abdullsm/develop/stars-h/src/backends/sequential/blrm/drsdd.c -o CMakeFiles/backends_sequential.dir/blrm/drsdd.c.s

src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dsdd.c.o: src/backends/sequential/CMakeFiles/backends_sequential.dir/flags.make
src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dsdd.c.o: ../src/backends/sequential/blrm/dsdd.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abdullsm/develop/stars-h/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building C object src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dsdd.c.o"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/backends_sequential.dir/blrm/dsdd.c.o   -c /home/abdullsm/develop/stars-h/src/backends/sequential/blrm/dsdd.c

src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dsdd.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/backends_sequential.dir/blrm/dsdd.c.i"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/abdullsm/develop/stars-h/src/backends/sequential/blrm/dsdd.c > CMakeFiles/backends_sequential.dir/blrm/dsdd.c.i

src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dsdd.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/backends_sequential.dir/blrm/dsdd.c.s"
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && /opt/ecrc/gcc/7.2.0/ub16/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/abdullsm/develop/stars-h/src/backends/sequential/blrm/dsdd.c -o CMakeFiles/backends_sequential.dir/blrm/dsdd.c.s

backends_sequential: src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dqp3.c.o
backends_sequential: src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/drsdd.c.o
backends_sequential: src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dsdd.c.o
backends_sequential: src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dsvfr.c.o
backends_sequential: src/backends/sequential/CMakeFiles/backends_sequential.dir/dense/dna.c.o
backends_sequential: src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dca.c.o
backends_sequential: src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dfe.c.o
backends_sequential: src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dmml.c.o
backends_sequential: src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dqp3.c.o
backends_sequential: src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/drsdd.c.o
backends_sequential: src/backends/sequential/CMakeFiles/backends_sequential.dir/blrm/dsdd.c.o
backends_sequential: src/backends/sequential/CMakeFiles/backends_sequential.dir/build.make

.PHONY : backends_sequential

# Rule to build all files generated by this target.
src/backends/sequential/CMakeFiles/backends_sequential.dir/build: backends_sequential

.PHONY : src/backends/sequential/CMakeFiles/backends_sequential.dir/build

src/backends/sequential/CMakeFiles/backends_sequential.dir/clean:
	cd /home/abdullsm/develop/stars-h/build/src/backends/sequential && $(CMAKE_COMMAND) -P CMakeFiles/backends_sequential.dir/cmake_clean.cmake
.PHONY : src/backends/sequential/CMakeFiles/backends_sequential.dir/clean

src/backends/sequential/CMakeFiles/backends_sequential.dir/depend:
	cd /home/abdullsm/develop/stars-h/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/abdullsm/develop/stars-h /home/abdullsm/develop/stars-h/src/backends/sequential /home/abdullsm/develop/stars-h/build /home/abdullsm/develop/stars-h/build/src/backends/sequential /home/abdullsm/develop/stars-h/build/src/backends/sequential/CMakeFiles/backends_sequential.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/backends/sequential/CMakeFiles/backends_sequential.dir/depend


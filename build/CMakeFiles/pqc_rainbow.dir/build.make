# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/flocau01/.local/lib/python3.8/site-packages/cmake/data/bin/cmake

# The command to remove a file.
RM = /home/flocau01/.local/lib/python3.8/site-packages/cmake/data/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/flocau01/projects/pqc_rainbow_signature

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/flocau01/projects/pqc_rainbow_signature/build

# Include any dependencies generated for this target.
include CMakeFiles/pqc_rainbow.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/pqc_rainbow.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/pqc_rainbow.dir/flags.make

CMakeFiles/pqc_rainbow.dir/src/main/signature_system.c.o: CMakeFiles/pqc_rainbow.dir/flags.make
CMakeFiles/pqc_rainbow.dir/src/main/signature_system.c.o: ../src/main/signature_system.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/flocau01/projects/pqc_rainbow_signature/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/pqc_rainbow.dir/src/main/signature_system.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/pqc_rainbow.dir/src/main/signature_system.c.o -c /home/flocau01/projects/pqc_rainbow_signature/src/main/signature_system.c

CMakeFiles/pqc_rainbow.dir/src/main/signature_system.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/pqc_rainbow.dir/src/main/signature_system.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/flocau01/projects/pqc_rainbow_signature/src/main/signature_system.c > CMakeFiles/pqc_rainbow.dir/src/main/signature_system.c.i

CMakeFiles/pqc_rainbow.dir/src/main/signature_system.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/pqc_rainbow.dir/src/main/signature_system.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/flocau01/projects/pqc_rainbow_signature/src/main/signature_system.c -o CMakeFiles/pqc_rainbow.dir/src/main/signature_system.c.s

CMakeFiles/pqc_rainbow.dir/src/main/test_functions.c.o: CMakeFiles/pqc_rainbow.dir/flags.make
CMakeFiles/pqc_rainbow.dir/src/main/test_functions.c.o: ../src/main/test_functions.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/flocau01/projects/pqc_rainbow_signature/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/pqc_rainbow.dir/src/main/test_functions.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/pqc_rainbow.dir/src/main/test_functions.c.o -c /home/flocau01/projects/pqc_rainbow_signature/src/main/test_functions.c

CMakeFiles/pqc_rainbow.dir/src/main/test_functions.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/pqc_rainbow.dir/src/main/test_functions.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/flocau01/projects/pqc_rainbow_signature/src/main/test_functions.c > CMakeFiles/pqc_rainbow.dir/src/main/test_functions.c.i

CMakeFiles/pqc_rainbow.dir/src/main/test_functions.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/pqc_rainbow.dir/src/main/test_functions.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/flocau01/projects/pqc_rainbow_signature/src/main/test_functions.c -o CMakeFiles/pqc_rainbow.dir/src/main/test_functions.c.s

# Object files for target pqc_rainbow
pqc_rainbow_OBJECTS = \
"CMakeFiles/pqc_rainbow.dir/src/main/signature_system.c.o" \
"CMakeFiles/pqc_rainbow.dir/src/main/test_functions.c.o"

# External object files for target pqc_rainbow
pqc_rainbow_EXTERNAL_OBJECTS =

pqc_rainbow: CMakeFiles/pqc_rainbow.dir/src/main/signature_system.c.o
pqc_rainbow: CMakeFiles/pqc_rainbow.dir/src/main/test_functions.c.o
pqc_rainbow: CMakeFiles/pqc_rainbow.dir/build.make
pqc_rainbow: librainbow.a
pqc_rainbow: /usr/lib/x86_64-linux-gnu/libssl.so
pqc_rainbow: /usr/lib/x86_64-linux-gnu/libcrypto.so
pqc_rainbow: CMakeFiles/pqc_rainbow.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/flocau01/projects/pqc_rainbow_signature/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C executable pqc_rainbow"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pqc_rainbow.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/pqc_rainbow.dir/build: pqc_rainbow

.PHONY : CMakeFiles/pqc_rainbow.dir/build

CMakeFiles/pqc_rainbow.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/pqc_rainbow.dir/cmake_clean.cmake
.PHONY : CMakeFiles/pqc_rainbow.dir/clean

CMakeFiles/pqc_rainbow.dir/depend:
	cd /home/flocau01/projects/pqc_rainbow_signature/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/flocau01/projects/pqc_rainbow_signature /home/flocau01/projects/pqc_rainbow_signature /home/flocau01/projects/pqc_rainbow_signature/build /home/flocau01/projects/pqc_rainbow_signature/build /home/flocau01/projects/pqc_rainbow_signature/build/CMakeFiles/pqc_rainbow.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/pqc_rainbow.dir/depend

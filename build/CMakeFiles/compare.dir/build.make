# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/mahenina/Documents/GitHub/bezier

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mahenina/Documents/GitHub/bezier/build

# Include any dependencies generated for this target.
include CMakeFiles/compare.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/compare.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/compare.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/compare.dir/flags.make

CMakeFiles/compare.dir/src/HornBez.cpp.o: CMakeFiles/compare.dir/flags.make
CMakeFiles/compare.dir/src/HornBez.cpp.o: ../src/HornBez.cpp
CMakeFiles/compare.dir/src/HornBez.cpp.o: CMakeFiles/compare.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mahenina/Documents/GitHub/bezier/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/compare.dir/src/HornBez.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/compare.dir/src/HornBez.cpp.o -MF CMakeFiles/compare.dir/src/HornBez.cpp.o.d -o CMakeFiles/compare.dir/src/HornBez.cpp.o -c /home/mahenina/Documents/GitHub/bezier/src/HornBez.cpp

CMakeFiles/compare.dir/src/HornBez.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/compare.dir/src/HornBez.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mahenina/Documents/GitHub/bezier/src/HornBez.cpp > CMakeFiles/compare.dir/src/HornBez.cpp.i

CMakeFiles/compare.dir/src/HornBez.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/compare.dir/src/HornBez.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mahenina/Documents/GitHub/bezier/src/HornBez.cpp -o CMakeFiles/compare.dir/src/HornBez.cpp.s

CMakeFiles/compare.dir/src/VS.cpp.o: CMakeFiles/compare.dir/flags.make
CMakeFiles/compare.dir/src/VS.cpp.o: ../src/VS.cpp
CMakeFiles/compare.dir/src/VS.cpp.o: CMakeFiles/compare.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mahenina/Documents/GitHub/bezier/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/compare.dir/src/VS.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/compare.dir/src/VS.cpp.o -MF CMakeFiles/compare.dir/src/VS.cpp.o.d -o CMakeFiles/compare.dir/src/VS.cpp.o -c /home/mahenina/Documents/GitHub/bezier/src/VS.cpp

CMakeFiles/compare.dir/src/VS.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/compare.dir/src/VS.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mahenina/Documents/GitHub/bezier/src/VS.cpp > CMakeFiles/compare.dir/src/VS.cpp.i

CMakeFiles/compare.dir/src/VS.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/compare.dir/src/VS.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mahenina/Documents/GitHub/bezier/src/VS.cpp -o CMakeFiles/compare.dir/src/VS.cpp.s

CMakeFiles/compare.dir/src/deCasteljau.cpp.o: CMakeFiles/compare.dir/flags.make
CMakeFiles/compare.dir/src/deCasteljau.cpp.o: ../src/deCasteljau.cpp
CMakeFiles/compare.dir/src/deCasteljau.cpp.o: CMakeFiles/compare.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mahenina/Documents/GitHub/bezier/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/compare.dir/src/deCasteljau.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/compare.dir/src/deCasteljau.cpp.o -MF CMakeFiles/compare.dir/src/deCasteljau.cpp.o.d -o CMakeFiles/compare.dir/src/deCasteljau.cpp.o -c /home/mahenina/Documents/GitHub/bezier/src/deCasteljau.cpp

CMakeFiles/compare.dir/src/deCasteljau.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/compare.dir/src/deCasteljau.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mahenina/Documents/GitHub/bezier/src/deCasteljau.cpp > CMakeFiles/compare.dir/src/deCasteljau.cpp.i

CMakeFiles/compare.dir/src/deCasteljau.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/compare.dir/src/deCasteljau.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mahenina/Documents/GitHub/bezier/src/deCasteljau.cpp -o CMakeFiles/compare.dir/src/deCasteljau.cpp.s

CMakeFiles/compare.dir/src/main.cpp.o: CMakeFiles/compare.dir/flags.make
CMakeFiles/compare.dir/src/main.cpp.o: ../src/main.cpp
CMakeFiles/compare.dir/src/main.cpp.o: CMakeFiles/compare.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mahenina/Documents/GitHub/bezier/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/compare.dir/src/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/compare.dir/src/main.cpp.o -MF CMakeFiles/compare.dir/src/main.cpp.o.d -o CMakeFiles/compare.dir/src/main.cpp.o -c /home/mahenina/Documents/GitHub/bezier/src/main.cpp

CMakeFiles/compare.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/compare.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mahenina/Documents/GitHub/bezier/src/main.cpp > CMakeFiles/compare.dir/src/main.cpp.i

CMakeFiles/compare.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/compare.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mahenina/Documents/GitHub/bezier/src/main.cpp -o CMakeFiles/compare.dir/src/main.cpp.s

# Object files for target compare
compare_OBJECTS = \
"CMakeFiles/compare.dir/src/HornBez.cpp.o" \
"CMakeFiles/compare.dir/src/VS.cpp.o" \
"CMakeFiles/compare.dir/src/deCasteljau.cpp.o" \
"CMakeFiles/compare.dir/src/main.cpp.o"

# External object files for target compare
compare_EXTERNAL_OBJECTS =

compare: CMakeFiles/compare.dir/src/HornBez.cpp.o
compare: CMakeFiles/compare.dir/src/VS.cpp.o
compare: CMakeFiles/compare.dir/src/deCasteljau.cpp.o
compare: CMakeFiles/compare.dir/src/main.cpp.o
compare: CMakeFiles/compare.dir/build.make
compare: CMakeFiles/compare.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/mahenina/Documents/GitHub/bezier/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX executable compare"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/compare.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/compare.dir/build: compare
.PHONY : CMakeFiles/compare.dir/build

CMakeFiles/compare.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/compare.dir/cmake_clean.cmake
.PHONY : CMakeFiles/compare.dir/clean

CMakeFiles/compare.dir/depend:
	cd /home/mahenina/Documents/GitHub/bezier/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mahenina/Documents/GitHub/bezier /home/mahenina/Documents/GitHub/bezier /home/mahenina/Documents/GitHub/bezier/build /home/mahenina/Documents/GitHub/bezier/build /home/mahenina/Documents/GitHub/bezier/build/CMakeFiles/compare.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/compare.dir/depend


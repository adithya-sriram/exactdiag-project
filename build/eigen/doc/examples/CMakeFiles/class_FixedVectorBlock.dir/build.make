# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

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
CMAKE_COMMAND = /opt/homebrew/lib/python3.11/site-packages/cmake/data/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/lib/python3.11/site-packages/cmake/data/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build"

# Include any dependencies generated for this target.
include Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/compiler_depend.make

# Include the progress variables for this target.
include Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/progress.make

# Include the compile flags for this target's objects.
include Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/flags.make

Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/class_FixedVectorBlock.cpp.o: Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/flags.make
Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/class_FixedVectorBlock.cpp.o: /Users/adithyasriram/Research\ Files/Khemani/repos/Exact\ Diagonalization\ Project/Eigen/doc/examples/class_FixedVectorBlock.cpp
Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/class_FixedVectorBlock.cpp.o: Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/class_FixedVectorBlock.cpp.o"
	cd "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/doc/examples" && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/class_FixedVectorBlock.cpp.o -MF CMakeFiles/class_FixedVectorBlock.dir/class_FixedVectorBlock.cpp.o.d -o CMakeFiles/class_FixedVectorBlock.dir/class_FixedVectorBlock.cpp.o -c "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/Eigen/doc/examples/class_FixedVectorBlock.cpp"

Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/class_FixedVectorBlock.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/class_FixedVectorBlock.dir/class_FixedVectorBlock.cpp.i"
	cd "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/doc/examples" && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/Eigen/doc/examples/class_FixedVectorBlock.cpp" > CMakeFiles/class_FixedVectorBlock.dir/class_FixedVectorBlock.cpp.i

Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/class_FixedVectorBlock.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/class_FixedVectorBlock.dir/class_FixedVectorBlock.cpp.s"
	cd "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/doc/examples" && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/Eigen/doc/examples/class_FixedVectorBlock.cpp" -o CMakeFiles/class_FixedVectorBlock.dir/class_FixedVectorBlock.cpp.s

# Object files for target class_FixedVectorBlock
class_FixedVectorBlock_OBJECTS = \
"CMakeFiles/class_FixedVectorBlock.dir/class_FixedVectorBlock.cpp.o"

# External object files for target class_FixedVectorBlock
class_FixedVectorBlock_EXTERNAL_OBJECTS =

Eigen/doc/examples/class_FixedVectorBlock: Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/class_FixedVectorBlock.cpp.o
Eigen/doc/examples/class_FixedVectorBlock: Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/build.make
Eigen/doc/examples/class_FixedVectorBlock: Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir="/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable class_FixedVectorBlock"
	cd "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/doc/examples" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/class_FixedVectorBlock.dir/link.txt --verbose=$(VERBOSE)
	cd "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/doc/examples" && ./class_FixedVectorBlock >/Users/adithyasriram/Research\ Files/Khemani/repos/Exact\ Diagonalization\ Project/build/Eigen/doc/examples/class_FixedVectorBlock.out

# Rule to build all files generated by this target.
Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/build: Eigen/doc/examples/class_FixedVectorBlock
.PHONY : Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/build

Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/clean:
	cd "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/doc/examples" && $(CMAKE_COMMAND) -P CMakeFiles/class_FixedVectorBlock.dir/cmake_clean.cmake
.PHONY : Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/clean

Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/depend:
	cd "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project" "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/Eigen/doc/examples" "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build" "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/doc/examples" "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/DependInfo.cmake" "--color=$(COLOR)"
.PHONY : Eigen/doc/examples/CMakeFiles/class_FixedVectorBlock.dir/depend

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

# Utility rule file for unsupported_examples.

# Include any custom commands dependencies for this target.
include Eigen/unsupported/doc/examples/CMakeFiles/unsupported_examples.dir/compiler_depend.make

# Include the progress variables for this target.
include Eigen/unsupported/doc/examples/CMakeFiles/unsupported_examples.dir/progress.make

unsupported_examples: Eigen/unsupported/doc/examples/CMakeFiles/unsupported_examples.dir/build.make
.PHONY : unsupported_examples

# Rule to build all files generated by this target.
Eigen/unsupported/doc/examples/CMakeFiles/unsupported_examples.dir/build: unsupported_examples
.PHONY : Eigen/unsupported/doc/examples/CMakeFiles/unsupported_examples.dir/build

Eigen/unsupported/doc/examples/CMakeFiles/unsupported_examples.dir/clean:
	cd "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/unsupported/doc/examples" && $(CMAKE_COMMAND) -P CMakeFiles/unsupported_examples.dir/cmake_clean.cmake
.PHONY : Eigen/unsupported/doc/examples/CMakeFiles/unsupported_examples.dir/clean

Eigen/unsupported/doc/examples/CMakeFiles/unsupported_examples.dir/depend:
	cd "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project" "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/Eigen/unsupported/doc/examples" "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build" "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/unsupported/doc/examples" "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/unsupported/doc/examples/CMakeFiles/unsupported_examples.dir/DependInfo.cmake" "--color=$(COLOR)"
.PHONY : Eigen/unsupported/doc/examples/CMakeFiles/unsupported_examples.dir/depend


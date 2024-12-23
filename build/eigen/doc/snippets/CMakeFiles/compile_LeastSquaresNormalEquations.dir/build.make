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
include Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/compiler_depend.make

# Include the progress variables for this target.
include Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/progress.make

# Include the compile flags for this target's objects.
include Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/flags.make

Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/compile_LeastSquaresNormalEquations.cpp.o: Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/flags.make
Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/compile_LeastSquaresNormalEquations.cpp.o: Eigen/doc/snippets/compile_LeastSquaresNormalEquations.cpp
Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/compile_LeastSquaresNormalEquations.cpp.o: /Users/adithyasriram/Research\ Files/Khemani/repos/Exact\ Diagonalization\ Project/Eigen/doc/snippets/LeastSquaresNormalEquations.cpp
Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/compile_LeastSquaresNormalEquations.cpp.o: Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/compile_LeastSquaresNormalEquations.cpp.o"
	cd "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/doc/snippets" && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/compile_LeastSquaresNormalEquations.cpp.o -MF CMakeFiles/compile_LeastSquaresNormalEquations.dir/compile_LeastSquaresNormalEquations.cpp.o.d -o CMakeFiles/compile_LeastSquaresNormalEquations.dir/compile_LeastSquaresNormalEquations.cpp.o -c "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/doc/snippets/compile_LeastSquaresNormalEquations.cpp"

Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/compile_LeastSquaresNormalEquations.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/compile_LeastSquaresNormalEquations.dir/compile_LeastSquaresNormalEquations.cpp.i"
	cd "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/doc/snippets" && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/doc/snippets/compile_LeastSquaresNormalEquations.cpp" > CMakeFiles/compile_LeastSquaresNormalEquations.dir/compile_LeastSquaresNormalEquations.cpp.i

Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/compile_LeastSquaresNormalEquations.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/compile_LeastSquaresNormalEquations.dir/compile_LeastSquaresNormalEquations.cpp.s"
	cd "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/doc/snippets" && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/doc/snippets/compile_LeastSquaresNormalEquations.cpp" -o CMakeFiles/compile_LeastSquaresNormalEquations.dir/compile_LeastSquaresNormalEquations.cpp.s

# Object files for target compile_LeastSquaresNormalEquations
compile_LeastSquaresNormalEquations_OBJECTS = \
"CMakeFiles/compile_LeastSquaresNormalEquations.dir/compile_LeastSquaresNormalEquations.cpp.o"

# External object files for target compile_LeastSquaresNormalEquations
compile_LeastSquaresNormalEquations_EXTERNAL_OBJECTS =

Eigen/doc/snippets/compile_LeastSquaresNormalEquations: Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/compile_LeastSquaresNormalEquations.cpp.o
Eigen/doc/snippets/compile_LeastSquaresNormalEquations: Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/build.make
Eigen/doc/snippets/compile_LeastSquaresNormalEquations: Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir="/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable compile_LeastSquaresNormalEquations"
	cd "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/doc/snippets" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/compile_LeastSquaresNormalEquations.dir/link.txt --verbose=$(VERBOSE)
	cd "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/doc/snippets" && ./compile_LeastSquaresNormalEquations >/Users/adithyasriram/Research\ Files/Khemani/repos/Exact\ Diagonalization\ Project/build/Eigen/doc/snippets/LeastSquaresNormalEquations.out

# Rule to build all files generated by this target.
Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/build: Eigen/doc/snippets/compile_LeastSquaresNormalEquations
.PHONY : Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/build

Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/clean:
	cd "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/doc/snippets" && $(CMAKE_COMMAND) -P CMakeFiles/compile_LeastSquaresNormalEquations.dir/cmake_clean.cmake
.PHONY : Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/clean

Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/depend:
	cd "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project" "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/Eigen/doc/snippets" "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build" "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/doc/snippets" "/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build/Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/DependInfo.cmake" "--color=$(COLOR)"
.PHONY : Eigen/doc/snippets/CMakeFiles/compile_LeastSquaresNormalEquations.dir/depend


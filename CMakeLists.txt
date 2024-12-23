cmake_minimum_required(VERSION 3.4...3.18)
set(CMAKE_CXX_STANDARD 17)
project(pyED)

# Add pybind and eigen subdirectories
add_subdirectory(pybind11)
add_subdirectory(Eigen)
INCLUDE_DIRECTORIES(pybind11/include)
INCLUDE_DIRECTORIES(Eigen)

# Add module
pybind11_add_module(pyED 
pyED.cpp 
lattice.cpp 
model.cpp 
op.cpp
construct_matrix.cpp
util.cpp)

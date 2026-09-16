cmake_minimum_required(VERSION 3.16)
project(bnct_clustering LANGUAGES CXX)

find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
execute_process(
  COMMAND "${Python3_EXECUTABLE}" -m pybind11 --cmakedir
  OUTPUT_VARIABLE pybind11_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  COMMAND_ERROR_IS_FATAL ANY
)
find_package(pybind11 REQUIRED CONFIG)

pybind11_add_module(clustering
  pyClustering.cc
  src/ClusteringAlgorithm.cc
  src/ClusterSBPoints_.cc
  src/SBPoint_.cc
)
target_include_directories(clustering PRIVATE include)
target_compile_features(clustering PRIVATE cxx_std_17)

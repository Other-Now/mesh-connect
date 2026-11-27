cmake_minimum_required(VERSION 3.20)
project(mesh_connect LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(MC_OPENMP "Parallel connector search with OpenMP" ON)
find_package(OpenMP)

if(MSVC)
  add_compile_options(/W4 /permissive- /utf-8)
  add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
else()
  add_compile_options(-Wall -Wextra -Wpedantic -Wshadow)
endif()

# Core library (static; linked into the CLI, the tests and the DLL)
add_library(mc_core STATIC
  src/geom.cpp src/bvh.cpp src/connect.cpp src/topology.cpp
  src/nastran.cpp src/gen.cpp src/ccx.cpp)
target_include_directories(mc_core PUBLIC include)
set_target_properties(mc_core PROPERTIES CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
if(MC_OPENMP AND OpenMP_CXX_FOUND)
  target_link_libraries(mc_core PUBLIC OpenMP::OpenMP_CXX)
endif()

# Shared library with the C ABI only (libmeshconnect.so / meshconnect.dll)
add_library(meshconnect SHARED src/c_api.cpp)
target_link_libraries(meshconnect PRIVATE mc_core)
target_compile_definitions(meshconnect PRIVATE MC_BUILDING_DLL)
set_target_properties(meshconnect PROPERTIES CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)

add_executable(mesh-connect src/main.cpp)
target_link_libraries(mesh-connect PRIVATE mc_core)

enable_testing()
add_executable(mc_tests tests/test_main.cpp)
target_link_libraries(mc_tests PRIVATE mc_core meshconnect)
add_test(NAME unit COMMAND mc_tests)

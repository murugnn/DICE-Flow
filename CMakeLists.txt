cmake_minimum_required(VERSION 3.10)
project(DiCE_Duan VERSION 1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Optimization flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")

# Set default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# OpenMP support (required for DNF pairwise parallelization)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
else()
    message(WARNING "OpenMP not found. DNF pairwise analysis will run single-threaded.")
endif()

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/cpp)

# OpenMP (engine.hpp uses omp pragmas / omp.h)
find_package(OpenMP)

# Main executable
add_executable(dice_analyzer
    src/cpp/main.cpp
)
if(OpenMP_CXX_FOUND)
    target_link_libraries(dice_analyzer PUBLIC OpenMP::OpenMP_CXX)
endif()

# Link OpenMP if available
if(OpenMP_CXX_FOUND)
    target_link_libraries(dice_analyzer PUBLIC OpenMP::OpenMP_CXX)
endif()

# Installation
install(TARGETS dice_analyzer DESTINATION bin)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ flags: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")

# ============================================================================
# Tests
# ============================================================================
enable_testing()
add_executable(test_engine src/cpp/tests/test_engine.cpp)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_engine PUBLIC OpenMP::OpenMP_CXX)
endif()
add_test(NAME engine_tests COMMAND test_engine)

cmake_minimum_required(VERSION 2.8)
project(EMPHF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

configure_file(
  ${EMPHF_SOURCE_DIR}/emphf_config.hpp.in
  ${EMPHF_SOURCE_DIR}/emphf_config.hpp)

option(EMPHF_USE_POPCOUNT
  "Use hardware popcount in hash computation (requires SSE4.2)"
  OFF)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif ()

if (UNIX)
   # C++11
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

   # Extensive warnings
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-missing-braces")

   if (EMPHF_USE_POPCOUNT)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2")
   endif ()
else ()
  message (FATAL_ERROR "Unsupported platform")
endif ()

OPTION (BUILD_EXECS "Build EMPHF Executabless" OFF)

IF (BUILD_EXECS)

# Sequential version
add_executable(compute_mphf_seq compute_mphf_seq.cpp)

# Sort-and-scan version (internal memory and mmap)
add_executable(compute_mphf_scan compute_mphf_scan.cpp)
add_executable(compute_mphf_scan_mmap compute_mphf_scan_mmap.cpp)

# HEM
add_executable(compute_mphf_hem compute_mphf_hem.cpp)

# Utilities
add_executable(test_mphf test_mphf.cpp)
add_executable(test_mphf_hem test_mphf_hem.cpp)
add_executable(gen_synthetic_data gen_synthetic_data.cpp)

ENDIF (BUILD_EXECS)

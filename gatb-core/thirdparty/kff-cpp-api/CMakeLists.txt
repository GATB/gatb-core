cmake_minimum_required(VERSION 3.16)

# set the project name
project(KFF)
project(KFF VERSION 1.0)
configure_file(${PROJECT_SOURCE_DIR}/kff_io.hpp.in ${CMAKE_BINARY_DIR}/kff_io.hpp)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# set(CMAKE_BUILD_TYPE Debug)

include_directories(${CMAKE_BINARY_DIR})

if (NOT CMAKE_CXX_FLAGS)
    set(CMAKE_CXX_FLAGS "-Wall -Ofast -g")
endif()

add_subdirectory("benchmarks")

add_library(kff kff_io.cpp)

target_include_directories(kff
    PUBLIC
    .
    ${PROJECT_SOURCE_DIR}/
)


# add the executable
add_executable(kff_example
    main.cpp
)

target_link_libraries(kff_example
    kff
)


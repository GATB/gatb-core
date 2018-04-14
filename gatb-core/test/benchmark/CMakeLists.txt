################################################################################
#  EXAMPLES 
################################################################################

include_directories (${gatb-core-includes})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") # needed for bench_mphf


list (APPEND PROGRAMS bench1 bench_bloom bench_mphf bench_minim bench_graph bench_bagfile) 

FOREACH (program ${PROGRAMS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} ${gatb-core-libraries})
ENDFOREACH (program)

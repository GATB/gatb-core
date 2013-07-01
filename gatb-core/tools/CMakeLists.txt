################################################################################
#  EXAMPLES 
################################################################################

set (LIBRARY_OUTPUT_PATH    lib/${CMAKE_BUILD_TYPE})

list (APPEND PROGRAMS ClearCache KmerChecksum)

FOREACH (program ${PROGRAMS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} ${gatb-core-libraries})
ENDFOREACH (program)

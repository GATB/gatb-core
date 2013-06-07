################################################################################
#  EXAMPLES 
################################################################################

set (LIBRARY_OUTPUT_PATH    lib/${CMAKE_BUILD_TYPE})

list (APPEND PROGRAMS BankConverter ClearCache DSK KmerChecksum Bloocoo)

FOREACH (program ${PROGRAMS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} GatbTools-static   dl  pthread  z  ${OMP})
ENDFOREACH (program)

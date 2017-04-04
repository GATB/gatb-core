################################################################################
message ("--  ---------- GATB TOOLS ----------")
################################################################################

################################################################################
#  TOOLS
################################################################################

# We add the compilation options for the library
add_definitions (${gatb-core-flags})

include_directories (${gatb-core-includes})

# We add the path for extra libraries
link_directories (${gatb-core-extra-libraries-path})

list (APPEND PROGRAMS dbgh5 dbginfo leon)

FOREACH (program ${PROGRAMS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} ${gatb-core-libraries})
  
  # (rayan) attempt at static linking, not sure if it works properly
  # but we need this to distribute binaries that are compatible on old system 
  # (see lh3's post on Biostar during the 2.0.1 dsk/minia release, he could not
  # run dsk even though i compiled it on gcc 4.4.7 on genocluster2)
  if (static)
     set_target_properties(${program} PROPERTIES LINK_SEARCH_START_STATIC w1)
      set_target_properties(${program} PROPERTIES LINK_SEARCH_END_STATIC 1)
      set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
      if (APPLE)
          set(CMAKE_EXE_LINKER_FLAGS "-static-libstdc++") #removed -static because doesn't work on mac's
      else()
          set(CMAKE_EXE_LINKER_FLAGS "-static-libstdc++ -static")
      endif()
      message ("STATIC BINARIES for ${program}")
  else()
      message ("DYNAMIC BINARIES for ${program}")
  endif()


ENDFOREACH (program)

################################################################################
#  INSTALLATION 
################################################################################
install (TARGETS dbgh5 dbginfo leon DESTINATION bin)

################################################################################
message ("--  ---------- GATB EXAMPLES ----------")
################################################################################

################################################################################
#  EXAMPLES 
################################################################################

# We add the compilation options for the library
add_definitions (${gatb-core-flags})

include_directories (${gatb-core-includes})

# We add the path for extra libraries
link_directories (${gatb-core-extra-libraries-path})

file (GLOB_RECURSE  PROGRAMS_FILES  *.cpp)

FOREACH (program ${PROGRAMS_FILES})

    # we get the base name of the current example
    get_filename_component(basename ${program} NAME_WE)
   
    # we add an executable for this example
    add_executable (${basename} ${program})
  
    # we use the gatb core libraries for the link
    target_link_libraries (${basename} ${gatb-core-libraries})

    list (APPEND PROGRAMS_NAMES ${basename})
      
ENDFOREACH (program)

# We add a custom target for making all examples
add_custom_target (examples ALL DEPENDS ${PROGRAMS_NAMES})

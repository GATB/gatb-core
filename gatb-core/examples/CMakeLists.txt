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

file (GLOB_RECURSE  PROGRAMS_FILES_EXTRA  
    bank13.cpp 
    debruijn7.cpp debruijn18.cpp debruijn19.cpp debruijn20.cpp debruijn21.cpp debruijn22.cpp debruijn23.cpp debruijn25.cpp
    kmer5.cpp kmer8.cpp kmer11.cpp
    storage4.cpp 
    multithreading2.cpp  multithreading4.cpp  multithreading5.cpp  multithreading6.cpp
    ToyTool.cpp  
)

# Some examples may require lambda expressions. If this option is not activated, we
# don't compile such examples.
if (NOT WITH_LAMBDA_EXPRESSION)
    message ("-- WARNING !!! Some examples won't be compiled (use of lambda expressions)")
    LIST (REMOVE_ITEM PROGRAMS_FILES ${PROGRAMS_FILES_EXTRA})
endif()

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

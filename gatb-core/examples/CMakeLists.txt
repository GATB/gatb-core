################################################################################
#  EXAMPLES 
################################################################################

include_directories (${gatb-core-includes})

list (APPEND PROGRAMS iterators1 iterators2 iterators3 iterators4 iterators5)
list (APPEND PROGRAMS bank1 bank2 bank3 bank4 bank5 bank6 bank7 bank8)
list (APPEND PROGRAMS kmer1 kmer2 kmer3 kmer4 kmer5 kmer6)

list (APPEND PROGRAMS 
    debruijn1  debruijn2  debruijn3  debruijn4  debruijn5  debruijn6  debruijn7  debruijn8  debruijn9
    debruijn10 debruijn11 debruijn12 debruijn13
)

FOREACH (program ${PROGRAMS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} ${gatb-core-libraries})
ENDFOREACH (program)

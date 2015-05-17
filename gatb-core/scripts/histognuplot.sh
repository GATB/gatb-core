#!/bin/bash

echo "-----------------------------------------------------------"
echo "DISPLAY kmers distribution for $1 with gnuplot"
echo "-----------------------------------------------------------"

h5dump -y -d histogram/histogram $1  | grep "[0-9]" | grep -v "[A-Z].*" | paste - - | gnuplot -p -e 'plot [][0:100] "-" with lines'

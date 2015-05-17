#!/bin/bash

h5dump -y -d histogram/histogram $1  | grep "[0-9]" | grep -v "[A-Z].*" | paste - - 

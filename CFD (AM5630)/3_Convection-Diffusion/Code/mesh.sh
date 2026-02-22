#!/bin/bash

# Utility to Plot the Mesh using gnuplot
# Author: Kanak Agarwal (kanakaero.github.io)
# Last Modified: 25th September 2025

gnuplot -persist -e "
    plot \
    'vlines.dat' with lines lc 'black' notitle, \
    'hlines.dat' with lines lc 'black' notitle, \
    'ccent.dat' with points pt 7 ps 0.4 lc 'red' notitle
"



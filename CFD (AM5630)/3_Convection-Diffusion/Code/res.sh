#!/bin/bash

gnuplot -persist <<-EOF
    set title "Residuals"
    set xlabel "Iteration"
    set ylabel "Residual" 
    set logscale y
    set grid

    # Format y-axis ticks as 10^n
    set format y "10^{%L}"
    
    plot "residuals_$1.dat" using 1:2 with lines linecolor rgb "red" linewidth 2 title "Residual"
EOF
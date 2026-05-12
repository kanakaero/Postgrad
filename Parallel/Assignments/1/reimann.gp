set terminal pngcairo size 900,600
set output "riemann_scaling.png"

set title "OpenMP Riemann Integration Performance"
set xlabel "Number of Threads"
set ylabel "Execution Time (seconds)"

set grid
set key top left

plot "riemann_timing.dat" using 1:2 with linespoints lw 2 pt 7 title "Reimann Integration"